#include "add_system_tables.hxx"
#include "base.hxx"
#include "execute_viparr.hxx"
#include "postprocess/build_constraints.hxx"
#include "postprocess/compile_plugins.hxx"
#include "postprocess/fix_masses.hxx"

#include <viparr/version.hxx> /* Auto-generated in SConscript */

#include <msys/analyze.hxx>
#include <msys/clone.hxx>
#include <msys/elements.hxx>
#include <msys/override.hxx>
#include <sstream>
#include <map>
#include <algorithm>

#include <cstdlib>
#include <cmath>
#include <csignal>

namespace desres { namespace viparr {

    constexpr double charge_difference_threshold = 0.05;

    void ExecuteViparr(const msys::SystemPtr sys,
                       const std::vector<ForcefieldPtr>& fflist,
                       const msys::IdList& atoms, bool rename_atoms, bool rename_residues,
                       bool with_constraints, bool fix_masses, bool fatal,
                       bool compile_plugins, bool verbose, bool verbose_matching) {

      if (atoms.size() == 0)
        VIPARR_FAIL("No atoms selected for VIPARR parametrization");

      double original_charge = 0.0;
      for(const auto & atom_id : sys->atoms()) {
        original_charge += sys->atom(atom_id).charge;
      }

      /* If system was not loaded structure only, add system's existing
       * tables to global parameter tables */
      if (sys->tableNames().size() > 0) {
        if (verbose)
          VIPARR_OUT << "Adding existing tables in input system to "
            "viparr's global param tables" << std::endl;
        AddSystemTables(sys);
      }

      /* Set vdw_funct and vdw_rule */
      std::string& vdw_funct = sys->nonbonded_info.vdw_funct;
      std::string& vdw_rule = sys->nonbonded_info.vdw_rule;
      if (vdw_funct != "") {
        bool found = false;
        for (std::map<std::string, Rules::VDWFunc>::const_iterator
               vdw_func_iter = Rules::VDWFuncRegistry().begin();
             vdw_func_iter != Rules::VDWFuncRegistry().end();
             ++vdw_func_iter) {
          /* Assumes there is a 1-1 reverse mapping of VDW table name */
          if (vdw_func_iter->second.vdw_table_name == vdw_funct) {
            vdw_funct = vdw_func_iter->first;
            found = true;
            break;
          }
        }
        if (!found)
          VIPARR_FAIL("Unrecognized vdw_funct '" + vdw_funct
                      + "' in input chemical system");
      }
      for (ForcefieldPtr ff : fflist) {
        vdw_funct = Rules::MergeVDW(vdw_funct, ff->rules()->vdw_func);
        vdw_rule = Rules::MergeVDW(vdw_rule, ff->rules()->vdw_comb_rule);
        if (!fatal)
          ff->rules()->fatal = false;
      }

      /* Provide default values for vdw_funct and vdw_rule in system, check
       * compatibility */
      if (vdw_funct == "") {
        VIPARR_ERR << "WARNING: No vdw_func provided; assuming default "
          "of lj12_6_sig_epsilon" << std::endl;
        vdw_funct = "lj12_6_sig_epsilon";
      }
      std::map<std::string, Rules::VDWFunc>::const_iterator vdw_func_iter
        = Rules::VDWFuncRegistry().find(vdw_funct);
      if (vdw_func_iter == Rules::VDWFuncRegistry().end())
        VIPARR_FAIL("Unsupported VDW function: " + vdw_funct);
      const std::vector<std::string>& supported_rules
        = vdw_func_iter->second.supported_rules;
      if (vdw_rule == "") {
        VIPARR_ERR << "WARNING: No vdw_comb_rule provided; assuming "
          "default of " << supported_rules[0] << std::endl;
        vdw_rule = supported_rules[0];
      }
      bool supported = false;
      for (unsigned i = 0; i < supported_rules.size(); ++i) {
        if (vdw_rule == supported_rules[i]) {
          supported = true;
          break;
        }
      }
      if (!supported) {
        VIPARR_FAIL("Unsupported combine rule " + vdw_rule
                    + " for VDW function: " + vdw_funct);
      }

      /* Sort atoms into fragments */
      std::vector<msys::IdList> fragments;
      unsigned nfrags = sys->updateFragids(&fragments);
      unsigned nfrags_all = nfrags;

      /* Take only fragments with selected atoms */
      msys::IdList atoms_no_pseudos;
      msys::IdList pseudos;
      std::vector<bool> in_selection(sys->maxAtomId(), false);
      for (unsigned i = 0; i < atoms.size(); ++i)
        in_selection[atoms[i]] = true;
      std::vector<msys::IdList>::iterator frag_iter = fragments.begin();
      while (frag_iter != fragments.end()) {
        msys::IdList& frag = *frag_iter;
        unsigned j = 0;
        while (j < frag.size() && sys->atom(frag[j]).atomic_number == 0)
          ++j;
        if (j < frag.size() && in_selection[frag[j]]) {
          /* First non-pseudo atom in fragment is in selection */
          msys::IdList::iterator atom_iter = frag.begin();
          while (atom_iter != frag.end()) {
            if (sys->atom(*atom_iter).atomic_number == 0) {
              pseudos.push_back(*atom_iter);
              atom_iter = frag.erase(atom_iter);
            } else if (!in_selection[*atom_iter]) {
              VIPARR_FAIL("Atom selection contains a partial "
                          "fragment");
            } else {
              atoms_no_pseudos.push_back(*atom_iter);
              ++atom_iter;
            }
          }
          ++frag_iter;
        } else {
          /* Fragment is all pseudos, or first non-pseudo atom is not in
           * selection */
          for (unsigned i = 0; i < frag.size(); ++i) {
            if (sys->atom(frag[i]).atomic_number > 0
                && in_selection[frag[i]])
              VIPARR_FAIL("Atom selection contains a partial "
                          "fragment");
          }
          frag_iter = fragments.erase(frag_iter);
          --nfrags;
        }
      }

      if (verbose)
        VIPARR_OUT << "Parametrizing " << nfrags << " of " << nfrags_all
                   << " total fragments" << std::endl;

      /* Get atomic formulas for parametrized fragments */
      std::vector<std::string> formulas(nfrags);
      for (unsigned i = 0; i < nfrags; ++i) {
        int counts[128];
        for (unsigned j = 0; j < 128; ++j)
          counts[j] = 0;
        for (unsigned j = 0; j < fragments[i].size(); ++j)
          ++counts[sys->atom(fragments[i][j]).atomic_number];
        std::stringstream formula;
        for (unsigned j = 1; j < 128; ++j) {
          if (counts[j] > 1)
            formula << msys::AbbreviationForElement(j) << counts[j];
          else if (counts[j] == 1)
            formula << msys::AbbreviationForElement(j);
        }
        formulas[i] = formula.str();
      }

      /* Remove terms with selected atoms from all tables */
      std::vector<std::string> tables = sys->tableNames();
      for (unsigned i = 0; i < tables.size(); ++i) {
        msys::TermTablePtr table = sys->table(tables[i]);
        msys::IdList terms = table->findWithAny(atoms);
        for (unsigned j = 0; j < terms.size(); ++j)
          table->delTerm(terms[j]);
        /* Remove overrides for unused params */
        if (table->overrides()->count() != 0) {
          std::vector<bool> used(table->params()->paramCount(), false);
          terms = table->terms();
          for (msys::Id term : terms) {
            if (table->param(term) != msys::BadId)
              used[table->param(term)] = true;
          }
          std::vector<msys::IdPair> overrides
            = table->overrides()->list();
          for (const msys::IdPair& pair : overrides) {
            if (!used[pair.first] && !used[pair.second])
              table->overrides()->del(pair);
          }
        }
      }

      /* Remove pseudos in selected atoms from system */
      for (unsigned i = 0; i < pseudos.size(); ++i)
        sys->delAtom(pseudos[i]);

      /* Remove any extraneous cmap tables */
      std::set<std::string> cmapids;
      msys::TermTablePtr torsion = sys->table("torsiontorsion_cmap");
      if (torsion != msys::TermTablePtr()) {
        for (msys::Id term : torsion->terms())
          cmapids.insert(torsion->propValue(term, "cmapid"));
      }
      std::vector<std::string> aux_tables = sys->auxTableNames();
      for (const std::string& aux_table : aux_tables) {
        if (aux_table.substr(0,4) == "cmap"
            && cmapids.find(aux_table) == cmapids.end())
          sys->delAuxTable(aux_table);
      }

      for (msys::Id atom : atoms) {
          for (msys::Id bond : sys->bondsForAtom(atom)) {
            if(sys->bond(bond).order == 0)
              sys->bond(bond).order = 1;
          }
        }
      /* Stores the name of every plugin matched, for use during the
         compilation stage of plugin application. */
      std::set<std::string> all_plugins;

      /* Apply forcefields */
      std::vector<bool> assigned(nfrags, false);
      std::vector<std::vector<std::string> > whynot(nfrags);
      for (ForcefieldPtr ff : fflist) {
        std::set<std::string> matched_formulas;
        std::set<std::string> warned_formulas;
        int matched_frags = 0;
        if (verbose)
          VIPARR_OUT << "Applying forcefield " << ff->name << std::endl;

        /* Assign atomtypes */
        if (verbose)
          VIPARR_OUT << "  Matching fragments and assigning atom types"
                     << std::endl;
        TemplatedSystemPtr tsys = TemplatedSystem::create(sys);
        for (unsigned frag = 0; frag < nfrags; ++frag) {
          std::stringstream ss;
          ss << "Forcefield " << ff->name << " ";
          std::vector<std::pair<TemplatedSystemPtr, msys::IdList> > matches;
          bool matched = ff->typer()->matchFragment(tsys, fragments[frag],
                                                    matches, ss);
          if (!matched)
            whynot[frag].push_back(ss.str());
          else {
            if (assigned[frag]) {
              if (warned_formulas.find(formulas[frag])
                  == warned_formulas.end()) {
                VIPARR_ERR << "WARNING: Fragment " << frag << ": "
                           << formulas[frag]
                           << " was matched by multiple forcefields; first"
                           << " match takes precedence" << std::endl;
                warned_formulas.insert(formulas[frag]);
              }
            } else {
              if ((verbose || verbose_matching) && matched_formulas.find(formulas[frag])
                  == matched_formulas.end()) {
                std::string tpl_name
                  = matches[0].first->system()->residue(0).name;
                if (matches.size() > 1)
                  tpl_name += " ...";
                VIPARR_OUT << "    fragment " << frag << ": "
                           << formulas[frag] << " (matched by " << tpl_name
                           << ")" << std::endl;
                matched_formulas.insert(formulas[frag]);

                if(verbose_matching) {
                  for(unsigned i = 0; i < matches.size(); ++i) {
                    msys::Id res = sys->atom(matches[i].second[0]).residue;
                    VIPARR_OUT << "      Matched residue " << res <<
                      " (" << sys->residue(res).name << ") to " <<
                      matches[i].first->system()->residue(0).name << "\n";
                  }
                }
              }
              ff->typer()->assignMatch(tsys, matches, rename_atoms,
                                       rename_residues);
              ++matched_frags;
              assigned[frag] = true;
            }
          }
        }
        if (verbose)
          VIPARR_OUT << "  Matched " << matched_frags
                     << " total fragments" << std::endl;

        /* Apply plugins for parameter matching */
        std::vector<std::string> plugins = ff->rules()->plugins;
        std::vector<std::string> plugins_done;
        if (verbose)
          VIPARR_OUT << "  Applying plugins" << std::endl;
        for (unsigned i = 0; i < plugins.size(); ++i) {
          std::string name = plugins[i];
          all_plugins.insert(name);
          std::map<std::string, Forcefield::PluginPtr>::const_iterator
            iter = Forcefield::PluginRegistry().find(name);
          if (iter == Forcefield::PluginRegistry().end())
            VIPARR_FAIL("Unsupported plugin: " + name);
          const std::vector<std::string>& prerequisites
            = iter->second->prerequisites;
          for (unsigned j = 0; j < prerequisites.size(); ++j) {
            std::vector<std::string>::iterator done_iter = std::find(
              plugins_done.begin(), plugins_done.end(),
              prerequisites[j]);
            std::vector<std::string>::iterator all_iter = std::find(
              plugins.begin(), plugins.end(), prerequisites[j]);
            /* Fail if a prerequisite plugin is in the forcefield but
             * after this plugin */
            if (done_iter == plugins_done.end()
                && all_iter != plugins.end())
              VIPARR_FAIL("Plugin " + prerequisites[j] +
                          " must come before plugin " + name);
          }
          if (verbose)
            VIPARR_OUT << "    " << name << std::endl;
          iter->second->match(tsys, ff);
          plugins_done.push_back(name);
        }
      }
      for (unsigned frag = 0; frag < nfrags; ++frag) {
        if (!assigned[frag]) {
          std::stringstream msg;
          msg << "Forcefield assignment failed: No forcefield could "
              << "parametrize fragment " << frag << ": "
              << formulas[frag] << ". Reason for each "
              << "forcefield: " << std::endl;
          for (unsigned i = 0; i < whynot[frag].size(); ++i)
            msg << whynot[frag][i];
          VIPARR_FAIL(msg.str());
        }
      }

      if (compile_plugins) {
        if (verbose)
          VIPARR_OUT << "Compiling system for DMS output" << std::endl;

        CompilePlugins(sys,all_plugins);
        CleanupSystem(sys);
      }

      if (with_constraints) {
        if (verbose)
          VIPARR_OUT << "Building constraints" << std::endl;
        BuildConstraints(sys, atoms, false, std::set<std::string>(), verbose);
      }

      if (fix_masses) {
        if (verbose)
          VIPARR_OUT << "Fixing masses" << std::endl;
        FixMasses(sys, atoms_no_pseudos, verbose);
      }

      double end_charge = 0.0;
      for(const auto & atom_id : sys->atoms()) {
        end_charge += sys->atom(atom_id).charge;
      }

      if(std::fabs(end_charge-original_charge) > charge_difference_threshold) {
        VIPARR_ERR << "WARNING: Total charge changed substantially during parameterization, from ";
        VIPARR_ERR << original_charge << " to " << end_charge << ".\n";
      }


      /* Add ff meta info */
      msys::ParamTablePtr ffMetaTable;
      if (sys->auxTable("forcefield") == msys::ParamTablePtr()) {
        ffMetaTable = msys::ParamTable::create();
        ffMetaTable->addProp("path", msys::StringType);
        ffMetaTable->addProp("info", msys::StringType);
      } else
        ffMetaTable = sys->auxTable("forcefield");
      for (ForcefieldPtr ff : fflist) {
        msys::Id param = ffMetaTable->addParam();
        ffMetaTable->value(param, "path") = ff->name;
        std::stringstream info;
        for (unsigned i = 0; i < ff->rules()->info.size(); ++i)
          info << ff->rules()->info[i] << "\n";
        ffMetaTable->value(param, "info") = info.str();
      }
      sys->addAuxTable("forcefield", ffMetaTable);
    }

    msys::SystemPtr ReorderIDs(msys::SystemPtr sys) {
      /* Reorder so that pseudos are next to parents */
      sys = msys::Clone(sys, sys->orderedIds());

      /* Reorder pairs in pairs and exclusions tables, if necessary */
      std::vector<std::string> tables = sys->tableNames();
      for (unsigned i = 0; i < tables.size(); ++i) {
        if (tables[i] != "exclusion" && tables[i].find("pair") != 0)
          continue;
        msys::TermTablePtr table = sys->table(tables[i]);
        msys::IdList terms = table->terms();
        for (unsigned j = 0; j < terms.size(); ++j) {
          msys::IdList atoms = table->atoms(terms[j]);
          if (atoms[1] < atoms[0]) {
            std::swap(atoms[0], atoms[1]);
            table->addTerm(atoms, table->param(terms[j]));
            table->delTerm(terms[j]);
          }
        }
      }

      /* msys::Clone creates new ParamTables; merge them back into the
       * global shared tables */
      AddSystemTables(sys);

      return sys;
    }

  }}
