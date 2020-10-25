#include "base.hxx"
#include "parameter_matcher.hxx"
#include "util/util.hxx"
#include <algorithm>
#include <cstring>
#include <sstream>


using namespace desres::viparr;

namespace {

    /* Match pattern string p to expression string e, where p may
     * contain wild characters */
    bool matchString(const std::string& e, const std::string& p, char wild) {

        if (p == "")
            return (e == "");

        /* Tokenize p by wild characters */
        std::vector<std::string> tokens;
        size_t i = 0;
        while (i < p.size()) {
            if (p[i] == wild) {
                ++i;
            }
            else {
                size_t pos = p.find(wild, i);
                if (pos == std::string::npos)
                    tokens.push_back(p.substr(i));
                else
                    tokens.push_back(p.substr(i,pos-i));
                i = pos;
            }
        }

        int n = tokens.size();
        if (n == 0) {
            /* p has all wild characters */
            return true;
        }

        /* Find first possible matching positions for tokens in e */
        std::vector<size_t> matches(n);
        matches[0] = e.find(tokens[0]);
        if (matches[0] == std::string::npos)
            return false;
        for (int i = 1; i < n; ++i) {
            matches[i] = e.find(tokens[i], matches[i-1] + tokens[i-1].size());
            if (matches[i] == std::string::npos)
                return false;
        }

        /* Check current match; if not good, find next match */
        while (true) {
            if (p.find(wild) != 0 && matches[0] != 0) {
                /* p starts with a non-wild, but there are no more matches
                 * where token[0] matches the start of e */
                return false;
            }
            if (p.find_last_of(wild) == p.size() - 1
                    || matches[n-1] + tokens[n-1].size() == e.size()) {
                /* p ends with a wild, or the current match extends to the end
                 * of e: match is good! */
                return true;
            }
            /* This match does not extend to the end of e; find next match */
            int i = n - 1;
            for (; i >= 0; --i) {
                /* Find next possible matching position for token i */
                size_t new_match = e.find(tokens[i], matches[i] + 1);
                if (new_match == std::string::npos)
                    continue;
                matches[i] = new_match;
                /* Find first possible matching positions for tokens i+1
                 * onward */
                bool found_next = true;
                for (int j = i + 1; j < n; ++j) {
                    new_match = e.find(tokens[j],
                            matches[j-1] + tokens[j-1].size());
                    if (new_match == std::string::npos) {
                        found_next = false;
                        break;
                    }
                    matches[j] = new_match;
                }
                if (found_next) {
                    /* Next possible matching positions were found */
                    break;
                }
                /* Else: Need to find next possible match for token i-1 */
            }
            if (i == -1) {
                /* No more possible matches */
                return false;
            }
        }
    }

    bool matchStringList(const std::vector<std::string>& expr,
            const std::vector<std::string>& pattern, char wild) {
        for (unsigned i = 0, n = expr.size(); i < n; ++i) {
            if (!matchString(expr[i], pattern[i], wild))
                return false;
        }
        return true;
    }

    bool isWild(const Pattern& pattern) {
        for (unsigned i = 0; i < pattern.atoms.size(); ++i)
            if (pattern.atoms[i].find(ParameterMatcher::ATOM_WILD)
                    != std::string::npos)
                return true;
        return false;
    }


}

namespace desres { namespace viparr {
    ParameterMatcher::ParameterMatcher(ForcefieldPtr ff,
            const std::string& table, SystemToPatternPtr sys_to_pattern,
            TypeToPatternPtr type_to_pattern,
            const std::vector<PermutationPtr>& perms) :
        _row_ids(ff->rowIDs(table).begin(), ff->rowIDs(table).end()),
        _sys_to_pattern(sys_to_pattern), _type_to_pattern(type_to_pattern),
        _perms(perms) {

        if (_row_ids.size() == 0)
            VIPARR_FAIL("Forcefield does not have param table " + table);
        _row_to_row_id = msys::IdList(
                Forcefield::ParamTable(table)->paramCount(), msys::BadId);
        for (unsigned i = 0; i < _row_ids.size(); ++i)
            if (_row_ids[i] < _row_to_row_id.size())
                _row_to_row_id[_row_ids[i]] = i;
        _param_table = Forcefield::ParamTable(table);
        init(type_to_pattern);
    }

    ParameterMatcher::ParameterMatcher(msys::ParamTablePtr param_table,
            const std::list<msys::Id>& rowIDs, 
            SystemToPatternPtr sys_to_pattern,
            TypeToPatternPtr type_to_pattern,
            const std::vector<PermutationPtr>& perms) :
        _param_table(param_table), _row_ids(rowIDs.begin(), rowIDs.end()),
        _sys_to_pattern(sys_to_pattern), _type_to_pattern(type_to_pattern),
        _perms(perms) {
      
        _row_to_row_id = msys::IdList(_param_table->paramCount(), msys::BadId);
        for (unsigned i = 0; i < _row_ids.size(); ++i)
            if (_row_ids[i] < _row_to_row_id.size())
                _row_to_row_id[_row_ids[i]] = i;
        init(type_to_pattern);
    }

    void ParameterMatcher::init(TypeToPatternPtr type_to_pattern) {

        unsigned nrows = _row_ids.size();
        _pattern_list.resize(nrows);
        if (nrows > 0)
            match_bonds = ((*type_to_pattern)(_param_table->value(_row_ids[0],
                            "type").asString()).bonds.size() > 0);
        for (unsigned i = 0; i < nrows; ++i) {
            _pattern_list[i] = (*type_to_pattern)(_param_table->value(
                        _row_ids[i], "type").asString());
            if (match_bonds != (_pattern_list[i].bonds.size() > 0)) {
                VIPARR_FAIL("Either all rows or no rows from a parameter "
                        "file can match bond types");
            }
            isWild(_pattern_list[i]) ?
                _wilds.push_back(i) : _exacts.push_back(i);
        }
    }

    ParameterMatcherPtr ParameterMatcher::create(ForcefieldPtr ff,
            const std::string& table, SystemToPatternPtr sys_to_pattern,
            TypeToPatternPtr type_to_pattern,
            const std::vector<PermutationPtr>& perms) {
        return ParameterMatcherPtr(new ParameterMatcher(ff, table,
                    sys_to_pattern, type_to_pattern, perms));
    }

    ParameterMatcherPtr ParameterMatcher::create(
            msys::ParamTablePtr param_table, const std::list<msys::Id>& rowIDs,
            SystemToPatternPtr sys_to_pattern, TypeToPatternPtr type_to_pattern,
            const std::vector<PermutationPtr>& perms) {
        return ParameterMatcherPtr(new ParameterMatcher(param_table, rowIDs,
                    sys_to_pattern, type_to_pattern, perms));
    }

    std::string determine_table_name(msys::SystemPtr sys,
                                     msys::ParamTablePtr pt) {
      std::vector<std::string> tableNames = sys->tableNames();
      for(unsigned i = 0; i < tableNames.size(); i++) {
        msys::TermTablePtr tt = sys->table(tableNames.at(i));
        if(tt->params().get() == pt.get()) {
          return tt->name();
        }
      }
      return "undetermined_functional_form";
    }

    void ParameterMatcher::print_fingerprint(TemplatedSystemPtr tsys, std::string table_name) {
      /* Stash the id's of virtuals and drudes in std::sets so we can print _ND_ and _A_. */
      auto pseudo_types = tsys->pseudoTypes();
      std::set<msys::Id> virtuals, drudes;
      for(const auto & pseudo_type : pseudo_types) {
        std::set<msys::Id> *destination = &virtuals;
        
        if(pseudo_type.name == "virtual_shift")
          destination = &virtuals;
        else if(pseudo_type.name == "drude_anharm")
          destination = &drudes;
        else
          VIPARR_ERR << "Encountered non-virtual, non-drude pseudo-type: " + pseudo_type.name << "; counting it as virtual.\n";

        for(const auto & ids : pseudo_type.sites_list)
          destination->insert(ids[0]);
      }
      
      for(const auto & kv : _fingerprint_counts) {
        msys::Id paramid = _row_ids[kv.first];

        std::string param_type;
        if(_param_table->propIndex("type") == msys::BadId)    
          param_type = _pattern_list[kv.first].print();
        else {
          param_type = _param_table->value(paramid, "type").asString();
          ViparrReplaceAll(param_type, " -", "-");
          ViparrReplaceAll(param_type, "- ", "-");
          ViparrReplaceAll(param_type, " =", "=");
          ViparrReplaceAll(param_type, "= ", "=");
        }

        std::cout << table_name << ":: " << param_type << " " << kv.second;

        /* Now, decide whether to apply one of _A_ or _ND_ types. For
         * now, this designation will be applied based on the
         * pseudo-typeness of the LAST atom in the tuple. */
        msys::Id designated_particle = _fingerprint_example_tuples[kv.first].back();
        if(table_name.substr(0,7) == "virtual")
          designated_particle = _fingerprint_example_tuples[kv.first][0];

        if(drudes.find(designated_particle) == drudes.end()) {
          if(virtuals.find(designated_particle) == virtuals.end())
            std::cout << " _A_,";
          std::cout << " _ND_";
        }        

        /*
          Now, print the values in the param corresponding to _row_ids[_matched_index].
        */
        std::vector<std::pair<std::string, double>> param_kv_pairs;
        for(unsigned ii = 0; ii < _param_table->propCount(); ii++)
          if(_param_table->propType(ii) != msys::ValueType::StringType)
            param_kv_pairs.push_back(std::make_pair(_param_table->propName(ii),  _param_table->value(paramid, ii).asFloat()));

        std::cout << " { ";
        for(unsigned ii = 0; ii < param_kv_pairs.size(); ++ii) {
          if(ii > 0)
            std::cout << " , ";
          std::cout << "\'" << param_kv_pairs[ii].first << "\' : " << param_kv_pairs[ii].second;
        }
        std::cout << " }\n";
      }
    }

    msys::Id ParameterMatcher::match(TemplatedSystemPtr sys,
            const msys::IdList& tuple, PermutationPtr* perm,
            bool allow_repeat) {
      /*
      for(const auto id : tuple)
        std::cout << id << " ";
      std::cout << "\n";
      */

        /* Generate Pattern from tuple */
        Pattern sys_pattern = (*_sys_to_pattern)(sys, tuple);
        StringList base_type = sys_pattern.atoms;

        /* If pattern/hierarchies is cached, return cached result */
        std::stringstream ss;
        ss << sys_pattern;
        Cache::const_iterator iter = _cache.find(ss.str());
//    std::cout << ss.str() << "\n";
      
        if (iter != _cache.end()) {
            if (perm != NULL)
                *perm = iter->second.second;
            //          std::cout << "Found it in the cache\n";
            msys::Id answer_row_id = iter->second.first;
            unsigned matched_index = std::find(_row_ids.begin(), _row_ids.end(), answer_row_id) - _row_ids.begin(); //index of answer_row_id in row_ids
            
            if(matched_index < _row_ids.size()) {
              _fingerprint_counts[matched_index]++;
              _fingerprint_example_tuples[matched_index] = tuple;
            }

            //std::cout << "Match found: " << _pattern_list[matched_index] << "\n";
            
            
            return answer_row_id;
        }

        /* Generate priority map from atom types for hierarchical matching */
        PriorityMap priority_map = getPriorityMap(sys, tuple);

        /*
        std::cout << "Priority levels for this tuple:";
        for(const auto & kv : priority_map) {
          auto priority_level = kv.second;
          std::cout << "\npriority_level ";

          for(auto & aa : kv.first)
            std::cout << aa << ",";               
          
          std::cout << " has " << priority_level.size() << " items: ";
          for(auto & aa : priority_level) {
            for(auto & bb : aa) {
              std::cout << bb << " ";
            }
          }
        }
        std::cout << "\n\n";
        */

        /* Match exact atom type patterns in priority order; make sure
         * there is only one match at the matched priority level. */
        for (auto iter = priority_map.rbegin(); 
                iter != priority_map.rend(); ++iter) {
            /* Loop through priority levels */
            std::set<StringList> priority_level = iter->second;
            int matched_index = -1;
            PermutationPtr matched_perm = Permutation::Null;
            for (unsigned i = 0; i < _exacts.size(); ++i) {
                /* Loop through exact atom type patterns */
                Pattern type_pattern = _pattern_list[_exacts[i]];
                std::set<StringList>::iterator level_iter;
                for (level_iter = priority_level.begin();
                        level_iter != priority_level.end(); ++level_iter) {
                    /* Loop through current priority level */
                    sys_pattern.atoms = *level_iter;
                    PermutationPtr tmp_perm = matchPattern(sys_pattern,
                            type_pattern);
                    if (tmp_perm != Permutation::Null) { /* Found a match */
                        if (matched_perm == Permutation::Null) {
                            /* First match found */
                            matched_index = _exacts[i];
                            matched_perm = tmp_perm;
                        } else if (matched_index != _exacts[i] &&
                                !allow_repeat) {
                            /* Matched a different parameter */
                            std::stringstream msg;
                            std::string tn = determine_table_name(sys->system()
                                                                  ,_param_table); 
                            sys_pattern.atoms = base_type;
                            msg << "For " << tn << ": "
                                << "Found two matches for " << sys_pattern 
                                << " (at the same priority): " << type_pattern
                                << " and " << _pattern_list[matched_index];
                            VIPARR_FAIL(msg.str());
                        }
                    }
                }
            }
            if (matched_perm != Permutation::Null) { /* Match was found */
                if (perm != NULL)
                    *perm = matched_perm;
                _cache.insert(std::make_pair(ss.str(),
                            CacheEntry(_row_ids[matched_index], matched_perm)));
                _fingerprint_counts[matched_index]++;
                _fingerprint_example_tuples[matched_index] = tuple;

                //std::cout << "Match found: " << _pattern_list[matched_index] << "\n";                  
                
                return _row_ids[matched_index];
            }
        }

        /* Match wildcard atom type patterns in the order in which they
         * appear, and return the first match */
        for (unsigned i = 0; i < _wilds.size(); ++i) {
            /* Loop through wildcard atom type patterns */
            Pattern type_pattern = _pattern_list[_wilds[i]];
            for (PriorityMap::reverse_iterator iter = priority_map.rbegin(); 
                    iter != priority_map.rend(); ++iter) {
                std::set<StringList>::iterator level_iter;
                for (level_iter = iter->second.begin();
                        level_iter != iter->second.end(); ++level_iter) {
                    /* Loop through priority levels and tuples in level */
                    sys_pattern.atoms = *level_iter;
                    PermutationPtr matched_perm = matchPattern(sys_pattern,
                            type_pattern);
                    if (matched_perm != Permutation::Null) { /* Found match */
                        if (perm != NULL)
                            *perm = matched_perm;
                        _cache.insert(std::make_pair(ss.str(),
                                    CacheEntry(_row_ids[_wilds[i]],
                                        matched_perm)));
//                        std::cout << "Wildcard whatever\n";
                        _fingerprint_counts[_wilds[i]]++;
                        _fingerprint_example_tuples[_wilds[i]] = tuple;
                        //std::cout << "Match found: " << _pattern_list[_wilds[i]] << "\n";                  

                        /* really, I should rename this function
                         * _match and write a wrapping function match
                         * that does the cache and fingerprint
                         * bookkeeping. also I should do fingerprint
                         * counts by row-ids since that's what this
                         * function returns, then in print_fingerprint
                         * I can compute the reverse map once and for
                         * all. this way I also don't burden callers
                         * that don't need the bookkeeping. */
                        
                        return _row_ids[_wilds[i]];
                    }
                }
            }
        }

        /* No match found */
        if (perm != NULL)
            *perm = Permutation::Null;
        _cache.insert(std::make_pair(ss.str(),
                    CacheEntry(msys::BadId, Permutation::Null)));
        //std::cout << "No match\n";
      
        return msys::BadId;
    }

    msys::IdList ParameterMatcher::matchMultiple(TemplatedSystemPtr sys,
            const msys::IdList& tuple) {

        /* Generate Pattern from tuple */
        Pattern sys_pattern = (*_sys_to_pattern)(sys, tuple);
        StringList base_type = sys_pattern.atoms;

        /* If pattern/hierarchies is cached, return cached result */
        std::stringstream ss;
        ss << sys_pattern << std::endl;
        MultipleCache::const_iterator iter = _multi_cache.find(ss.str());
        if (iter != _multi_cache.end())
            return iter->second;

        /* Generate priority map from atom types */
        PriorityMap priority_map = getPriorityMap(sys, tuple);

        /* Match type patterns; save all matches */
        std::set<msys::Id> matches;
        msys::IdList matches_list;
        for (PriorityMap::reverse_iterator iter = priority_map.rbegin(); 
                iter != priority_map.rend(); ++iter) {
            for (unsigned i = 0; i < _exacts.size(); ++i) {
                Pattern type_pattern = _pattern_list[_exacts[i]];
                for (std::set<StringList>::iterator level_iter
                        = iter->second.begin(); level_iter != iter->second.end();
                        ++level_iter) {
                    sys_pattern.atoms = *level_iter;
                    PermutationPtr tmp_perm = matchPattern(sys_pattern,
                            type_pattern);
                    if (tmp_perm != Permutation::Null) {
                        if (matches.find(_row_ids[_exacts[i]]) == matches.end()) {
                            matches_list.push_back(_row_ids[_exacts[i]]);
                            matches.insert(_row_ids[_exacts[i]]);
                        }
                    }
                }
            }
        }
        for (unsigned i = 0; i < _wilds.size(); ++i) {
            Pattern type_pattern = _pattern_list[_wilds[i]];
            for (PriorityMap::reverse_iterator iter = priority_map.rbegin(); 
                    iter != priority_map.rend(); ++iter) {
                for (std::set<StringList>::iterator level_iter
                        = iter->second.begin(); level_iter != iter->second.end();
                        ++level_iter) {
                    sys_pattern.atoms = *level_iter;
                    PermutationPtr tmp_perm = matchPattern(sys_pattern,
                            type_pattern);
                    if (tmp_perm != Permutation::Null) {
                        if (matches.find(_row_ids[_wilds[i]]) == matches.end()) {
                            matches_list.push_back(_row_ids[_wilds[i]]);
                            matches.insert(_row_ids[_wilds[i]]);
                        }
                    }
                }
            }
        }
        _multi_cache.insert(std::make_pair(ss.str(), matches_list));
        return matches_list;
    }

    ParameterMatcher::PriorityMap ParameterMatcher::getPriorityMap(
            TemplatedSystemPtr sys, const msys::IdList& tuple) const {

        PriorityMap priority_map;
        unsigned natoms = tuple.size();
        std::vector<int> priority(natoms, 0);

            std::set<StringList> patterns;
            patterns.insert((*_sys_to_pattern)(sys, tuple).atoms);
            priority_map.insert(std::make_pair(priority, patterns));
            return priority_map;
    }
    
    PermutationPtr ParameterMatcher::matchPattern(const Pattern& sys_pattern,
            const Pattern& type_pattern) const {
        for (unsigned i = 0; i < _perms.size(); ++i) {
            Pattern copy_pattern = (*_perms[i])(sys_pattern);
            if (copy_pattern.atoms.size() != type_pattern.atoms.size())
                continue;
            if (!matchStringList(copy_pattern.atoms, type_pattern.atoms,
                        ATOM_WILD))
                continue;
            if (type_pattern.bonds.size() != 0) {
                if (copy_pattern.bonds.size() != type_pattern.bonds.size())
                    continue;
                if (!matchStringList(copy_pattern.bonds, type_pattern.bonds,
                            BOND_WILD))
                    continue;
            }
            if (copy_pattern.flags.size() != type_pattern.flags.size())
                continue;
            if (!matchStringList(copy_pattern.flags, type_pattern.flags, '\0'))
                continue;
            return _perms[i];
        }
        return Permutation::Null;
    }

    void ParameterMatcher::writeMultiple(msys::Id row,
            const msys::IdList& term, msys::TermTablePtr table) const {
        if (row >= _row_to_row_id.size() ||
                _row_to_row_id[row] == msys::BadId ||
                _row_ids[_row_to_row_id[row]] != row)
            VIPARR_FAIL("Invalid row for writeMultiple");
        for (unsigned i = _row_to_row_id[row];
                i < _row_ids.size(); ++i) {
            if (_param_table->value(_row_ids[i], "type")
                        == _param_table->value(row, "type"))
                table->addTerm(term, _row_ids[i]);
            else
                break;
        }
    }
}}
