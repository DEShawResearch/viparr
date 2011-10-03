#include "../parameter_matcher.hxx"
#include <sstream>

using namespace desres;
using namespace desres::viparr;
using desres::msys::Id;
using desres::msys::IdList;

static void apply_virtuals_regular(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    const std::vector<TemplatedSystem::PseudoType>& pseudo_types
        = sys->pseudoTypes();
    for (unsigned i = 0; i < pseudo_types.size(); ++i) {
        std::string name = pseudo_types[i].name;
        if (name.compare(0, 8, "virtual_") != 0) continue;
        if (name == "virtual_shift") continue;
        const std::vector<IdList>& virtuals = pseudo_types[i].sites_list;
        if (virtuals.size() == 0) continue;
        int nsites = virtuals[0].size();

        /* Forcefield templates use "virtual_..." while param files use
         * "virtuals_..." (this is ridiculous...) */
        std::string vname = "virtuals_";
        vname += name.substr(8);
        msys::TermTablePtr vtable = sys->system()->addTable(name,
                nsites, Forcefield::ParamTable(vname));
        vtable->category = msys::VIRTUAL;

        /* Match pset of pseudo, btypes of site atoms, and bonds to
         * parent atom (or second site atom for virtual_fdat3) */
        SystemToPatternPtr sys_to_pattern = SystemToPattern::PseudoBType;
        if (name == "virtual_fdat3")
            sys_to_pattern = SystemToPattern::PseudoBondToSecond;
        std::vector<PermutationPtr> perms(1, Permutation::Identity);
        ParameterMatcherPtr matcher = ParameterMatcher::create(ff, vname,
                sys_to_pattern, TypeToPattern::Pseudo, perms);

        unsigned old_size = vtable->termCount();
        for (unsigned j = 0; j < virtuals.size(); ++j) {
            const IdList& term = virtuals[j];
            Id row = matcher->match(sys, term);
            if (row == msys::BadId) {
                Pattern patt = (*sys_to_pattern)(sys, term);
                if (!matcher->match_bonds)
                    patt.bonds.clear();
                std::stringstream msg;
                msg << "No match found for table '" << vname << "', pattern "
                    << patt.print() << ", atoms (" << term[0];
                for (unsigned k = 1; k < term.size(); ++k)
                    msg << "," << term[k];
                msg << ")";
                if (!ff->rules()->fatal) {
                    VIPARR_ERR << "WARNING: " << msg.str() << std::endl;
                    continue;
                } else
                    VIPARR_FAIL(msg.str());
            }
            vtable->addTerm(term, row);
        }
        if (ff->rules()->fatal &&
                vtable->termCount() - old_size != virtuals.size())
            VIPARR_FAIL("VIPARR BUG--incorrect number of terms added");
    }
}

static Forcefield::RegisterPlugin _("virtuals_regular", apply_virtuals_regular);
