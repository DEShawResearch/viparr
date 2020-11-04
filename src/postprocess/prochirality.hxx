#ifndef desres_viparr_prochirality_hxx
#define desres_viparr_prochirality_hxx

#include <msys/system.hxx>

namespace desres { namespace viparr {
    msys::IdList FixProchiralProteinAtomNames(msys::SystemPtr s, bool checkOnly=false);
}}

#endif
