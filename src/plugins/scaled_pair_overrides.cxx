#include "../parameter_matcher.hxx"
#include "pairs_helper.hxx"

using namespace desres;
using namespace desres::viparr;
using msys::Id;
using msys::IdList;

/* Adds directly specified scaled pairs terms to the pairs table. If a pair that
 * is to be added already exists, the parameters of that pair term are
 * overridden with the new parameters.
 *
 * This plugin is used by iviparr forcefield outputs to adjust pairs terms that
 * may have been generated using vdw1_14 or charges_screening parameters. */

void match_scaled_pair_overrides(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    if (ff->rowIDs("scaled_pair_overrides").size() == 0)
        VIPARR_FAIL("Must have scaled_pair_overrides table for "
                "scaled_pair_overrides plugin");
    msys::TermTablePtr overrides = sys->system()->addTable("scaled_pair_overrides", 2,
            Forcefield::ParamTable("scaled_pair_overrides"));
    overrides->category = msys::NO_CATEGORY;
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    ParameterMatcherPtr matcher = ParameterMatcher::create(ff,
            "scaled_pair_overrides", SystemToPattern::NBType,
            TypeToPattern::Default, perms);
    int rule = ff->rules()->exclusions();
    if (rule > 4)
        VIPARR_FAIL("Maximum supported exclusion rule is 4");
    const std::vector<IdList>* alltuples[3];
    alltuples[0] = &sys->nonPseudoBonds();
    alltuples[1] = &sys->angles();
    alltuples[2] = &sys->dihedrals();
    for (int i = 2; i <= rule; ++i) {
        for (unsigned j = 0; j < alltuples[i-2]->size(); ++j) {
            IdList pair(2);
            pair[0] = (*alltuples[i-2])[j][0];
            pair[1] = (*alltuples[i-2])[j][i-1];
            if (overrides->findWithAll(pair).size() > 0) continue;
            /* Add pairs between ai/its pseudos and aj/its pseudos*/
            IdList atomsi;
            atomsi.push_back(pair[0]);
            IdList bonded = sys->system()->bondedAtoms(pair[0]);
            for (unsigned ki = 0; ki < bonded.size(); ++ki)
                if (sys->system()->atom(bonded[ki]).atomic_number == 0)
                    atomsi.push_back(bonded[ki]);
            IdList atomsj;
            atomsj.push_back(pair[1]);
            bonded = sys->system()->bondedAtoms(pair[1]);
            for (unsigned kj = 0; kj < bonded.size(); ++kj)
                if (sys->system()->atom(bonded[kj]).atomic_number == 0)
                    atomsj.push_back(bonded[kj]);
            for (unsigned ki = 0; ki < atomsi.size(); ++ki) {
                for (unsigned kj = 0; kj < atomsj.size(); ++kj) {
                    /* Find scaled pair override params */
                    IdList term(2);
                    term[0] = atomsi[ki];
                    term[1] = atomsj[kj];
                    IdList rows = matcher->matchMultiple(sys, term);
                    Id row = msys::BadId;
                    for (unsigned r = 0; r < rows.size(); ++r) {
                        if (matcher->paramTable()->value(rows[r], "separation")
                                == i) {
                            if (row != msys::BadId) {
                                std::stringstream msg;
                                msg << "Rows " << row << " and " << rows[r]
                                    << " of scaled_pair_overrides table "
                                    << " have the same types and separation";
                                VIPARR_FAIL(msg.str());
                            }
                            row = rows[r];
                        }
                    }
                    if (row != msys::BadId)
                        overrides->addTerm(term, row);
                }
            }
        }
    }
}

void compile_scaled_pair_overrides(msys::SystemPtr sys) {
    msys::TermTablePtr pairs = get_pairs_table(sys);
    if(pairs == msys::TermTablePtr())
        VIPARR_FAIL("Unable to find pairs table.\n");
        
    msys::TermTablePtr overrides = sys->table("scaled_pair_overrides");
    if (overrides == msys::TermTablePtr()) return;
    std::map<std::string, Rules::VDWFunc>::const_iterator iter =
        Rules::VDWFuncRegistry().find(sys->nonbonded_info.vdw_funct);
    if (iter == Rules::VDWFuncRegistry().end()) {
        VIPARR_FAIL("Unsupported VDW function '" + sys->nonbonded_info.vdw_funct
                + "'; make sure the VDW function was copied from the "
                "forcefield to the system and is in the VDWFuncRegistry");
    }
    const std::vector<std::string>& pair_props = iter->second.pair_param_names;
    msys::IdList terms = overrides->terms();
    for (msys::Id term : terms) {
        msys::IdList atoms = overrides->atoms(term);
        std::map<std::string, double> params;
        params.insert(std::make_pair("qij", overrides->propValue(term,
                        "qij").asFloat()));
        for (std::string prop : pair_props) {
            params.insert(std::make_pair(prop, overrides->propValue(term,
                            prop).asFloat()));
        }
        update_pair_params(pairs, sys, atoms[0], atoms[1], params);
    }
}

static Forcefield::RegisterPlugin _("scaled_pair_overrides",
                                    match_scaled_pair_overrides,
                                    std::vector<std::string>(),
                                    compile_scaled_pair_overrides);
