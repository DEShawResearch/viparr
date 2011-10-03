#ifndef desres_viparr_fix_masses_hxx
#define desres_viparr_fix_masses_hxx

#include <msys/system.hxx>

namespace desres { namespace viparr {

    void FixMasses(msys::SystemPtr sys, const msys::IdList& atoms,
            bool verbose=true);

}}

#endif
