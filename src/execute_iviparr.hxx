#ifndef desres_viparr_iviparr_hxx
#define desres_viparr_iviparr_hxx

#include "ff.hxx"
#include <msys/system.hxx>
#include <vector>

namespace desres { namespace viparr {

    /* WARNING -- Call this function only from the iviparr command-line
     * program. This function currently does not leave the Forcefield class 
     * shared param tables in a proper state for further use. */
    ForcefieldPtr ExecuteIviparr(msys::SystemPtr sys,
            const msys::IdList& atoms,
            RulesPtr rules, bool templateonly);

}}

#endif

