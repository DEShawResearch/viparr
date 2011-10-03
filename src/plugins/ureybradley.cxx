#include "add_nbody_table.hxx"

using namespace desres;
using namespace desres::viparr;

/* urey_bradley terms are paired with stretch_harm params */
static void apply_ureybradley(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    if (ff->rowIDs("ureybradley_harm").size() == 0)
        VIPARR_FAIL("Must have 'ureybradley_harm' table for "
                "'ureybradley' plugin");
    msys::TermTablePtr table = sys->system()->addTable("ureybradley_harm",
            2, Forcefield::ParamTable("ureybradley_harm"));
    table->category = msys::NO_CATEGORY;
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    ParameterMatcherPtr matcher = ParameterMatcher::create(ff,
            "ureybradley_harm", SystemToPattern::Bonded, TypeToPattern::Default,
            perms);
    const std::vector<msys::IdList>& angles = sys->angles();
    for (unsigned i = 0, n = angles.size(); i < n; ++i) {
        const msys::IdList& term = angles[i];
        msys::Id row = matcher->match(sys, term);
        if (row != msys::BadId) {
            msys::IdList pair(2);
            pair[0] = term[0];
            pair[1] = term[2];
            table->addTerm(pair, row);
        }
    }
}

static void compile_ureybradley(msys::SystemPtr sys) {
    msys::TermTablePtr ureybradley = sys->table("ureybradley_harm");
    if (ureybradley == msys::TermTablePtr()) return;
    msys::ParamTablePtr stretch_params;
    if (!Forcefield::HasParamTable("stretch_harm")) {
        stretch_params = msys::ParamTable::create();
        Forcefield::AddParamTable("stretch_harm", stretch_params);
    } else
        stretch_params = Forcefield::ParamTable("stretch_harm");
    for (unsigned i = 0; i < ureybradley->params()->propCount(); ++i)
        stretch_params->addProp(ureybradley->params()->propName(i),
                ureybradley->params()->propType(i));
    msys::TermTablePtr stretch = sys->addTable("stretch_harm", 2,
            stretch_params);
    stretch->addTermProp("constrained", msys::IntType);
    msys::IdList terms = ureybradley->terms();
    for (msys::Id term : terms) {
        msys::IdList match = stretch->findWithAll(ureybradley->atoms(term));
        for (msys::IdList::iterator iter = match.begin();
                iter != match.end(); ) {
            if (stretch->propValue(*iter, "type").asString().substr(0,11)
                    != "ureybradley")
                iter = match.erase(iter);
            else
                ++iter;
        }
        if (match.size() > 1)
            VIPARR_FAIL("Multiple identical ureybradley terms found in "
                    "stretch_harm table");
        if (match.size() == 0)
            match.push_back(stretch->addTerm(ureybradley->atoms(term),
                        stretch_params->addParam()));
        for (unsigned i = 0; i < ureybradley->params()->propCount(); ++i) {
            if (ureybradley->params()->propName(i) == "type") {
                std::string value = "ureybradley ";
                value += ureybradley->propValue(term, i).asString();
                stretch->propValue(match[0], ureybradley->params()->propName(i))
                    = value;
            } else
                stretch->propValue(match[0], ureybradley->params()->propName(i))
                    = ureybradley->propValue(term, i);
        }
    }
}


static Forcefield::RegisterPlugin _("ureybradley", apply_ureybradley,
                                    std::vector<std::string>(),
                                    compile_ureybradley);
