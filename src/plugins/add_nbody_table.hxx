#ifndef desres_viparr_add_nbody_table_hxx
#define desres_viparr_add_nbody_table_hxx

#include "../parameter_matcher.hxx"

namespace desres { namespace viparr {

    /* Matches all atom tuples in the list nbodies to a Forcefield ff and
     * pattern/param table table_name, using the given SystemToPattern,
     * TypeToPattern, and Permutations, and writes the matches to a term table
     * table_name of given category. Throws an exception if any tuple cannot
     * be matched. */
    void AddNbodyTable(TemplatedSystemPtr sys, ForcefieldPtr ff, 
            const std::string& table_name, const std::string& plugin_name,
            int natoms, const std::vector<msys::IdList>& nbodies,
            SystemToPatternPtr sys_to_pattern, TypeToPatternPtr type_to_pattern,
            const std::vector<PermutationPtr>& perms, bool required,
            msys::Category category);

}}

#endif

