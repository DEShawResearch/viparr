#include "../parameter_matcher.hxx"

using namespace desres;
using namespace desres::viparr;

static void apply_vdw(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    if (ff->rowIDs("vdw1").size() == 0)
        VIPARR_FAIL("Must have 'vdw1' table for 'vdw1' plugin");
    msys::TermTablePtr table = sys->system()->addTable("nonbonded", 1,
            Forcefield::ParamTable("vdw1"));
    table->category = msys::NONBONDED;
    ParameterMatcherPtr matcher = ParameterMatcher::create(ff, "vdw1",
            SystemToPattern::NBType, TypeToPattern::Default,
            std::vector<PermutationPtr>(1, Permutation::Identity));
    unsigned old_size = table->termCount();
    const std::vector<msys::IdList>& atoms = sys->typedAtoms();
    for (unsigned i = 0, n = atoms.size(); i < n; ++i) {
        const msys::IdList& term = atoms[i];
        msys::Id row = matcher->match(sys, term);
        if (row == msys::BadId) {
            Pattern patt = (*SystemToPattern::NBType)(sys, term);
            if (!matcher->match_bonds)
                patt.bonds.clear();
            std::stringstream msg;
            msg << "No match found for table 'vdw1', pattern "
                << patt.print() << ", atom " << term[0];
            if (!ff->rules()->fatal) {
                VIPARR_ERR << "WARNING: " << msg.str() << std::endl;
                continue;
            } else
                VIPARR_FAIL(msg.str());
        }
        table->addTerm(term, row);
    }
    /* Double-check that all typed atoms were matched */
    if (table->termCount() - old_size != atoms.size()
            && ff->rules()->fatal)
        VIPARR_FAIL("VIPARR BUG--incorrect number of terms added");
}

static Forcefield::RegisterPlugin _("vdw1", apply_vdw);
