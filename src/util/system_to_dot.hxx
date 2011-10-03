#ifndef desres_viparr_system_to_dot_hxx
#define desres_viparr_system_to_dot_hxx

#include <msys/system.hxx>
#include <iostream>

namespace desres { namespace viparr {

    void SystemToDot(msys::SystemPtr sys, std::ostream& dotfile, msys::Id residue_id=msys::BadId);

}}

#endif
