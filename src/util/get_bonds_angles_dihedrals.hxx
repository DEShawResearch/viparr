#ifndef desres_viparr_get_bonds_angles_dihedrals_hxx
#define desres_viparr_get_bonds_angles_dihedrals_hxx

#include <msys/system.hxx>
#include <vector>

namespace desres { namespace viparr {

    /* Returns lists of non-pseudo bonds, pseudo bonds, angles, and dihedrals
     * in a given fragment or set of fragments */
    void GetBondsAnglesDihedrals(msys::SystemPtr sys, const msys::IdList& atoms,
            std::vector<msys::IdList>& non_pseudo_bonds,
            std::vector<msys::IdList>& pseudo_bonds,
            std::vector<msys::IdList>& angles,
            std::vector<msys::IdList>& dihedrals);
}}

#endif
