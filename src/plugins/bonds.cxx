#include "add_nbody_table.hxx"

using namespace desres;
using namespace desres::viparr;

static void apply_stretch_harm(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    /* Match forward and reverse permutations with bonds */
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    AddNbodyTable(sys, ff, "stretch_harm", "bonds", 2, sys->nonPseudoBonds(),
            SystemToPattern::Bonded, TypeToPattern::Default, perms, true,
            msys::BOND);
    /* Add "constrained" property */
    sys->system()->table("stretch_harm")->addTermProp("constrained",
            msys::IntType);
}

static Forcefield::RegisterPlugin _("bonds", apply_stretch_harm);
