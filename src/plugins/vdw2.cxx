#include "../ff.hxx"

using namespace desres;
using namespace desres::viparr;

static void apply_vdw2(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    if (ff->rowIDs("vdw2").size() == 0)
        VIPARR_FAIL("Must have 'vdw2' table for 'vdw2' plugin");
    if (sys->system()->atoms().size() == 0) return;
    /* Create a term table with dummy terms to keep track of which vdw2 params
     * should be applied to this system */
    msys::TermTablePtr table = sys->system()->addTable("vdw2", 1,
            Forcefield::ParamTable("vdw2"));
    table->category = msys::NO_CATEGORY;
    for (msys::Id row : ff->rowIDs("vdw2")) {
        table->addTerm(msys::IdList(1, sys->system()->atoms()[0]), row);
    }
}

static Forcefield::RegisterPlugin _("vdw2", apply_vdw2);
