#include "add_nbody_table.hxx"

using namespace desres;
using namespace desres::viparr;

static void match_charges_formal(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    std::vector<PermutationPtr> perms(1, Permutation::Identity);
    AddNbodyTable(sys, ff, "charges_formal", "charges_formal", 1,
            sys->typedAtoms(), SystemToPattern::Bonded, TypeToPattern::Default,
            perms, false, msys::NO_CATEGORY);
    msys::ParamTablePtr params = Forcefield::ParamTable("charges_formal");
    msys::Id p = params->addParam();
    params->value(p, "charge") = 0;
    msys::TermTablePtr charges = sys->system()->table("charges_formal");
    for (msys::IdList atom : sys->typedAtoms())
        if (charges->findExact(atom).size() == 0)
            charges->addTerm(atom, p);
}

static void compile_charges_formal(msys::SystemPtr sys) {
    msys::TermTablePtr table = sys->table("charges_formal");
    if (table == msys::TermTablePtr()) return;
    msys::IdList terms = table->terms();
    for (msys::Id term : terms)
        sys->atom(table->atoms(term)[0]).charge = table->propValue(term,
                "charge").asFloat();
}

/* Must recompute both charges_formal and charges_bci if either is updated */
static Forcefield::RegisterPlugin _("charges_formal", match_charges_formal,
                                    std::vector<std::string>(),
                                    compile_charges_formal);
