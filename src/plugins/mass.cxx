#include "add_nbody_table.hxx"

using namespace desres;
using namespace desres::viparr;

/* If nonbonded, match mass based on nbtype; otherwise use btype */
template <bool nonbonded> 
static void match_mass(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    SystemToPatternPtr sys_to_pattern = SystemToPattern::BType;
    std::string plugin_name = "mass";
    if (nonbonded) {
        sys_to_pattern = SystemToPattern::NBType;
        plugin_name = "mass2";
    }
    std::vector<PermutationPtr> perms(1, Permutation::Identity);
    AddNbodyTable(sys, ff, "mass", plugin_name, 1, sys->typedAtoms(),
            sys_to_pattern, TypeToPattern::Default, perms, true,
            msys::NO_CATEGORY);
}

static void compile_mass(msys::SystemPtr sys) {
    msys::TermTablePtr table = sys->table("mass");
    msys::IdList terms = table->terms();
    for (msys::Id term : terms) {
        if (table->param(term) == msys::BadId)
            VIPARR_ERR << "WARNING: Missing mass term for atom "
                       << table->atoms(term)[0] << std::endl;
        else
            sys->atom(table->atoms(term)[0]).mass = table->propValue(term,
                    "amu").asFloat();
    }
}

static Forcefield::RegisterPlugin _("mass", match_mass<false>, std::vector<std::string>(),
                                    compile_mass);
static Forcefield::RegisterPlugin __("mass2", match_mass<true>, std::vector<std::string>(),
                                    compile_mass);
