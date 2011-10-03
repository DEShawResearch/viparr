#ifndef desres_viparr_parameter_matcher_hxx
#define desres_viparr_parameter_matcher_hxx

#include "ff.hxx"
#include "pattern.hxx"
#include <set>
#include <map>
#include <msys/term_table.hxx>

namespace desres { namespace viparr {
    /* An instance of ParameterMatcher stores a Forcefield's param table and
     * corresponding rowIDs; it can match atom tuples from a system against
     * the param table's 'type' column, using hierarchical atom-types.
     * Typical use:
     *
     * ParameterMatcherPtr matcher
     *   = ParameterMatcher::create(ff, table_name, ...);
     * Id row = matcher->match(system, tuple);
     * term_table->addTerm(tuple, row);
     */
    class ParameterMatcher {

    public:
      ParameterMatcher(ForcefieldPtr ff, const std::string& table,
                       SystemToPatternPtr sys_to_pattern,
                       TypeToPatternPtr type_to_pattern,
                       const std::vector<PermutationPtr>& perms);
      ParameterMatcher(msys::ParamTablePtr param_table,
                       const std::list<msys::Id>& rowIDs,
                       SystemToPatternPtr sys_to_pattern,
                       TypeToPatternPtr type_to_pattern,
                       const std::vector<PermutationPtr>& perms);


      /* These characters in the pattern table are treated as wildcards
       * during string matching */
      static const char ATOM_WILD = '*';
      static const char BOND_WILD = '~';

      /* The param table, subset of row IDs, and atom type hierarchy
       * can be passed to create() directly or by specifying a
       * forcefield and table name. The matcher operates by constructing
       * a system Pattern object from the atom tuple and a type Pattern
       * object from the param table and matching these Pattern objects;
       * Pattern objects are constructed using instances of
       * SystemToPattern and TypeToPattern. The matcher can match various
       * permutations of the system Pattern, specified by perms. See
       * "pattern.hxx" for more documentation. */
      static std::shared_ptr<ParameterMatcher> create(ForcefieldPtr ff,
                                                        const std::string& table,
                                                        SystemToPatternPtr sys_to_pattern,
                                                        TypeToPatternPtr type_to_pattern,
                                                        const std::vector<PermutationPtr>& perms);
      static std::shared_ptr<ParameterMatcher> create(
        msys::ParamTablePtr param_table,
        const std::list<msys::Id>& rowIDs,
        SystemToPatternPtr sys_to_pattern,
        TypeToPatternPtr type_to_pattern,
        const std::vector<PermutationPtr>& perms);

      /* Matches a single atom tuple in a system to the contained param
       * table, by constructing the system Pattern from the tuple and
       * the type Pattern for each row.
       *
       * Rows of the param table are matched in the following order:
       * 1. Rows with no atom types containing ATOM_WILD
       *  Note: If the atom type hierarchy is non-empty, the returned
       *  match is the row whose atom types have the highest
       *  hierarchical-priority for the given atom tuple; see
       *  getPriorityMap() below for definition of priority. If
       *  there is more than one row matching at the same priority,
       *  an exception is thrown.
       * 2. Rows with wild-card atom types, in the order in which they
       *  appear in the pattern table
       * 
       * If the type Pattern specifies bond types, then bond types of the
       * system Pattern and type Pattern must match; otherwise bond types
       * of the system Pattern are ignored. Flags in the system Pattern
       * and type Pattern must always match.
       *
       * Returns the ID of the first matched row in the param table,
       * and optionally also the Permutation that generated the match
       * (currently needed for bci charge balancing). If no match is
       * found, returns msys::BadId and Permutation::Null. */
      msys::Id match(TemplatedSystemPtr sys, const msys::IdList& tuple,
                     PermutationPtr* perm = NULL, bool allow_repeat=false);

      /* Match a single atom tuple in a system to the contained param
       * table; return IDs of all rows that match. */
      msys::IdList matchMultiple(TemplatedSystemPtr sys, const
                                 msys::IdList& tuple);

      /* For a given tuple of atoms and matched row, finds all subsequent
       * consecutive rows in the param table with identical "type" and
       * writes one entry in the term table for each such row. (Used for
       * dihedral_trig tables with multiple rows for a single type
       * pattern.) Note that the behavior is different from
       * matchMultiple. */
      void writeMultiple(msys::Id row, const msys::IdList& tuple,
                         msys::TermTablePtr table) const;

      /* Access the contained param table */
      msys::ParamTablePtr paramTable() const { return _param_table; }

      /* Access the row IDs being matched by this matcher */
      std::list<msys::Id> rowIDs() const { return std::list<msys::Id>(
          _row_ids.begin(), _row_ids.end()); }

      /* Access the SystemToPattern and TypeToPattern objects for this
       * matcher */
      const SystemToPatternPtr sysToPattern() const {
        return _sys_to_pattern; }
      const TypeToPatternPtr typeToPattern() const {
        return _type_to_pattern; }

      //void set_match_callback(std::function<void(Pattern p, msys::Id row)> new_callback);
      void print_fingerprint(TemplatedSystemPtr tsys, std::string table_name);

      bool match_bonds;

    private:
      typedef std::vector<std::string> StringList;

      msys::ParamTablePtr _param_table;
      msys::IdList _row_ids;
      msys::IdList _row_to_row_id;
      SystemToPatternPtr _sys_to_pattern;
      TypeToPatternPtr _type_to_pattern;
      std::vector<PermutationPtr> _perms; 

      /* A map from _matched_index to a count of the number of times
         that index was the match found. DESRESCode#1849 */
      std::map<int, int> _fingerprint_counts;

      /* A map from _matched_index to a SINGLE example of a tuple (of
       * atom id's)that was matched to that parameter. This is
       * collected to assign _A_ and _ND_ types during
       * fingerprinting. */
      std::map<int, msys::IdList> _fingerprint_example_tuples;
      
      //std::function<void(Pattern p, msys::Id row)> _match_callback;

      /* We preprocess the param table using TypeToPattern during
       * construction, saving the constructed Pattern objects in
       * _pattern_list */
      std::vector<Pattern> _pattern_list;

      /* Keep track of which rows in the pattern table correspond to exact
       * atom types and which have wild atom type characters */
      std::vector<int> _exacts;
      std::vector<int> _wilds;

      /* Helper function for constructors */
      void init(TypeToPatternPtr type_to_pattern);

      /* Keep a cache of matched Patterns */
      typedef std::pair<msys::Id, PermutationPtr> CacheEntry;
      typedef std::map<std::string, CacheEntry> Cache;
      Cache _cache;
      typedef std::map<std::string, msys::IdList> MultipleCache;
      MultipleCache _multi_cache;

      /* Helper functions to match two Patterns; returns matching
       * Permutation if matched, or Permutation::Null if not matched */
      PermutationPtr matchPattern(const Pattern& sys_pattern,
                                  const Pattern& type_pattern) const;

      /* Used for hierarchical matching. Given a tuple of atoms,
       * generates a priority map of atom type tuples to match. The
       * priority map keys are tuples of ints, with each key associated
       * to a list of all atom type tuples at that priority. */
      typedef std::map<std::vector<int>,
                       std::set<StringList> > PriorityMap;
      PriorityMap getPriorityMap(TemplatedSystemPtr sys,
                                 const msys::IdList& tuple) const;
    };
    typedef std::shared_ptr<ParameterMatcher> ParameterMatcherPtr;

  }}

#endif
