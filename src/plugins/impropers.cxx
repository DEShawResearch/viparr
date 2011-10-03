#include "add_nbody_table.hxx"

using namespace desres;
using namespace desres::viparr;
using desres::msys::Id;
using desres::msys::IdList;

static void apply_improper_harm(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    /* Match forward and reverse permutations with no bonds. (We do not specify
     * the bond pattern for improper_harm because different forcefields have
     * different conventions for where to put the center atom.) */
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    AddNbodyTable(sys, ff, "improper_harm", "impropers", 4, sys->impropers(),
            SystemToPattern::BType, TypeToPattern::Default, perms, true,
            msys::BOND);
}

static void apply_improper_anharm(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    msys::TermTablePtr table = sys->system()->addTable("improper_anharm", 4,
            Forcefield::ParamTable("improper_anharm"));
    table->category = msys::BOND;
    /* If the center atom is last, the improper comes from an old-style
     * template, and we match only the identity permutation with no bonds for
     * backwards compatibility. */
    std::vector<PermutationPtr> perms_old(1, Permutation::Identity);
    ParameterMatcherPtr matcher_old = ParameterMatcher::create(ff,
            "improper_anharm", SystemToPattern::BType, TypeToPattern::Default,
            perms_old);
    /* If the center atom is first, the improper comes from SMARTS matching,
     * and we match all six permutations of the non-center atoms with bonds to
     * the center atom. */
    std::vector<PermutationPtr> perms_new(Permutation::Improper,
            Permutation::Improper + 6);
    ParameterMatcherPtr matcher_new = ParameterMatcher::create(ff,
            "improper_anharm", SystemToPattern::BondToFirst,
            TypeToPattern::Default, perms_new);
    const std::vector<IdList>& impropers = sys->impropers();
    for (unsigned i = 0, n = impropers.size(); i < n; ++i) {
        const IdList& term = impropers[i];
        if (sys->system()->findBond(term[3], term[0]) != msys::BadId
                && sys->system()->findBond(term[3], term[1]) != msys::BadId
                && sys->system()->findBond(term[3], term[2]) != msys::BadId) {
            /* Old style, center atom last */
            Id row = matcher_old->match(sys, term);
            if (row == msys::BadId) {
                Pattern patt = (*SystemToPattern::BType)(sys, term);
                std::stringstream msg;
                msg << "No match found for table improper_anharm, pattern "
                    << patt.print() << ", atoms (" << term[0] << "," << term[1]
                    << "," << term[2] << "," << term[3] << ")";
                if (!ff->rules()->fatal) {
                    VIPARR_ERR << "WARNING: " << msg.str() << std::endl;
                    continue;
                }
                else
                    VIPARR_FAIL(msg.str());
            }
            table->addTerm(term, row);
        } else if (sys->system()->findBond(term[0], term[1]) != msys::BadId
                && sys->system()->findBond(term[0], term[2]) != msys::BadId
                && sys->system()->findBond(term[0], term[3]) != msys::BadId) {
            /* New style, center atom first */
            Id row = matcher_new->match(sys, term);
            if (row == msys::BadId) {
                Pattern patt = (*SystemToPattern::BondToFirst)(sys, term);
                if (!matcher_new->match_bonds)
                    patt.bonds.clear();
                std::stringstream msg;
                msg << "No match found for table 'improper_anharm', pattern "
                    << patt.print() << ", atoms (" << term[0] << "," << term[1]
                    << "," << term[2] << "," << term[3] << ")";
                if (!ff->rules()->fatal) {
                    VIPARR_ERR << "WARNING: " << msg.str() << std::endl;
                    continue;
                }
                else
                    VIPARR_FAIL(msg.str());
            }
            table->addTerm(term, row);
        } else {
          /* DESRESCode#3487 Ignore this improper term. */
        }
    }
}

/* improper_trig terms are paired with dihedral_trig params */
static void apply_improper_trig(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    msys::TermTablePtr table = sys->system()->addTable("improper_trig", 4,
            Forcefield::ParamTable("improper_trig"));
    table->category = msys::NO_CATEGORY;
    /* Match forward and reverse permutations with no bonds. (We do not specify
     * the bond pattern for improper_trig because different forcefields have
     * different conventions for where to put the center atom.) */
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    ParameterMatcherPtr matcher = ParameterMatcher::create(ff, "improper_trig",
            SystemToPattern::BType, TypeToPattern::Default, perms);
    const std::vector<IdList>& impropers = sys->impropers();
    for (unsigned i = 0, n = impropers.size(); i < n; ++i) {
        const IdList& term = impropers[i];
        Id row = matcher->match(sys, term);
        if (row == msys::BadId) {
            Pattern patt = (*SystemToPattern::BType)(sys, term);
            std::stringstream msg;
            msg << "No match found for table 'improper_trig', pattern "
                << patt.print() <<", atoms (" << term[0] << "," << term[1]
                << "," << term[2] << "," << term[3] << ")";

                VIPARR_ERR << "WARNING: " << msg.str() << std::endl;
                continue;
        }
        table->addTerm(term, row);
    }
}

static void apply_impropers(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    if (ff->rowIDs("improper_harm").size() > 0)
        apply_improper_harm(sys, ff);
    if (ff->rowIDs("improper_anharm").size() > 0)
        apply_improper_anharm(sys, ff);
    if (ff->rowIDs("improper_trig").size() > 0)
        apply_improper_trig(sys, ff);
}

static void compile_improper_trig(msys::SystemPtr sys) {
    msys::TermTablePtr impropers = sys->table("improper_trig");
    if (impropers == msys::TermTablePtr()) return;
    msys::ParamTablePtr dihedral_params;
    if (!Forcefield::HasParamTable("dihedral_trig")) {
        dihedral_params = msys::ParamTable::create();
        Forcefield::AddParamTable("dihedral_trig", dihedral_params);
    } else
        dihedral_params = Forcefield::ParamTable("dihedral_trig");
    for (unsigned i = 0; i < impropers->params()->propCount(); ++i)
        dihedral_params->addProp(impropers->params()->propName(i),
                impropers->params()->propType(i));
    msys::TermTablePtr dihedrals = sys->addTable("dihedral_trig", 4,
            dihedral_params);
    dihedrals->category = msys::BOND;
    msys::IdList terms = impropers->terms();
    // sentinel value added to type field to indicate the added improper terms
    std::string sentinel="improper ";
    for (msys::Id term : terms) {
        msys::IdList match = dihedrals->findWithAll(impropers->atoms(term));
        for (msys::IdList::iterator iter = match.begin();
                iter != match.end(); ) {
            if (dihedrals->propValue(*iter, "type").asString().substr(0,sentinel.length())
                    != sentinel)
                iter = match.erase(iter);
            else
                ++iter;
        }
        if (match.size() > 1)
            VIPARR_FAIL("Multiple identical improper trig terms found in "
                    "dihedral_trig table");
        if (match.size() == 0)
            match.push_back(dihedrals->addTerm(impropers->atoms(term),
                        dihedral_params->addParam()));
        for (unsigned i = 0; i < impropers->params()->propCount(); ++i) {
            if (impropers->params()->propName(i) == "type") {
                std::string value = sentinel + impropers->propValue(term, i).asString();
                dihedrals->propValue(match[0], impropers->params()->propName(i))
                    = value;
            } else
                dihedrals->propValue(match[0], impropers->params()->propName(i))
                    = impropers->propValue(term, i);
        }
    }
    // Now remove the sentinel value
    Id typeIdx=dihedral_params->propIndex("type");
    for (unsigned i=0; i < dihedral_params->paramCount(); ++i){
        std::string value=dihedral_params->value(i,typeIdx).asString();
        if(value.substr(0,sentinel.length()) == sentinel)
            dihedral_params->value(i,typeIdx)=value.substr(sentinel.length());
    }

}

static Forcefield::RegisterPlugin _("impropers", apply_impropers, std::vector<std::string>(),
                                    compile_improper_trig);
