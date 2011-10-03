#include "add_nbody_table.hxx"

using namespace desres;
using namespace desres::viparr;

static void apply_angle_harm(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    /* Match forward and reverse permutations with bonds */
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    AddNbodyTable(sys, ff, "angle_harm", "angles", 3, sys->angles(),
            SystemToPattern::Bonded, TypeToPattern::Default, perms, true,
            msys::BOND);
    sys->system()->table("angle_harm")->addTermProp("constrained",
            msys::IntType);
}

static Forcefield::RegisterPlugin _("angles", apply_angle_harm);
