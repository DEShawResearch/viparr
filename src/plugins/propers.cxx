#include "../parameter_matcher.hxx"

using namespace desres;
using namespace desres::viparr;

/* If fatal, exits with an error if a dihedral term cannot be matched;
 * otherwise prints a warning. */
template <bool fatal> 
static void apply_dihedral_trig(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    if (ff->rowIDs("dihedral_trig").size() == 0)
        VIPARR_FAIL("Must have 'dihedral_trig' table for 'propers' plugin");
    msys::TermTablePtr table = sys->system()->addTable("dihedral_trig", 4,
            Forcefield::ParamTable("dihedral_trig"));
    table->category = msys::BOND;

    /* Match forward and reverse permutations with bonds */
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    ParameterMatcherPtr matcher = ParameterMatcher::create(ff, "dihedral_trig",
            SystemToPattern::Bonded, TypeToPattern::Default, perms);

    unsigned old_size = table->termCount();
    const std::vector<msys::IdList>& dihedrals = sys->dihedrals();
    for (unsigned i = 0, n = dihedrals.size(); i < n; ++i) {
        const msys::IdList& term = dihedrals[i];
        msys::Id row = matcher->match(sys, term, NULL, true);
        if (row == msys::BadId) {
            Pattern patt = (*SystemToPattern::Bonded)(sys, term);
            if (!matcher->match_bonds)
                patt.bonds.clear();
            std::stringstream msg;
            msg << "No match found for table 'dihedral_trig', pattern "
                << patt.print() << ", atoms (" << term[0] << "," << term[1]
                << "," << term[2] << "," << term[3] << ")";
            if (!ff->rules()->fatal || !fatal) {
                VIPARR_ERR << "WARNING: " << msg.str() << std::endl;
                continue;
            }
            else
                VIPARR_FAIL(msg.str());
        }
        /* In some forcefields, there may be multiple matches for a single
         * term */
        matcher->writeMultiple(row, term, table);
    }
    if (fatal && ff->rules()->fatal
            && (table->termCount() - old_size < dihedrals.size()))
        VIPARR_FAIL("VIPARR BUG--incorrect number of terms added");
}

static Forcefield::RegisterPlugin _("propers", apply_dihedral_trig<true>);
static Forcefield::RegisterPlugin __("propers_allowmissing", apply_dihedral_trig<false>);
