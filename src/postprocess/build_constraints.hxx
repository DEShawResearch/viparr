#ifndef desres_viparr_build_constraints_hxx
#define desres_viparr_build_constraints_hxx

#include <msys/system.hxx>
#include <string>
#include <set>

namespace desres { namespace viparr {

    /* Adds constraint parameter tables to the system and updates
     * "constrained" fields of stretch_harm and angle_harm tables. Supports
     * HOH and AHn constraints. If keep is true, sets "constrained" to 0
     * for all stretch_harm and angle_harm params; otherwise, sets
     * "constrained" to 1 for constrained params. Optionally takes a set
     * of constraint types to ignore; types can be "hoh", "ah1", "ah2",
     * etc. */
    void BuildConstraints(msys::SystemPtr sys, const msys::IdList& atoms,
            bool keep=false, const std::set<std::string>& 
            exclude=std::set<std::string>(), bool verbose=true);

}}

#endif
