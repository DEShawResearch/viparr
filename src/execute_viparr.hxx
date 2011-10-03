#ifndef desres_viparr_execute_viparr_hxx
#define desres_viparr_execute_viparr_hxx

#include "ff.hxx"
#include <msys/system.hxx>
#include <vector>

namespace desres { namespace viparr {

    /* The main viparr executable to parametrize a system with a list of
     * forcefields.
     * FIXME: this is a horror show.
     */
    void ExecuteViparr(const msys::SystemPtr input_sys,
            const std::vector<ForcefieldPtr>& ffs,
            const msys::IdList& atoms, bool rename_atoms=false,
            bool rename_residues=false, bool with_constraints=true,
            bool optimize_vsite_defs=true, bool fix_masses=true,
            bool fatal=true,
            bool compile_plugins=true, bool verbose=true,
            bool verbose_matching=false);

    /* Create a copy of the system in which IDs of pseudo atoms are adjacent
     * to their parent atoms */
    msys::SystemPtr ReorderIDs(msys::SystemPtr sys);

}}

#endif
