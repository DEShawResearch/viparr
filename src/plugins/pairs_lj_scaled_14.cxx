#include "add_nbody_table.hxx"
#include "pairs_helper.hxx"

using namespace desres;
using namespace desres::viparr;

/* Overwrite scaled lennard-jones pair interactions in the nonbonded pairs
 * table for 1-4 interactions using the parameters in table 'vdw1_14'. If
 * parameters are not found, parameters from the default nonbonded table are
 * used. */
void match_pairs_lj_scaled_14(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    std::vector<PermutationPtr> perms(1, Permutation::Identity);
    AddNbodyTable(sys, ff, "vdw1_14", "pairs_lj_scaled_14", 1,
            sys->typedAtoms(), SystemToPattern::Bonded, TypeToPattern::Default,
            perms, false, msys::NO_CATEGORY);
}
/* This compiles both pairs_lj_scaled and pairs_lj_scaled_14 */
void compile_pairs_lj_scaled(msys::SystemPtr sys) {
    msys::TermTablePtr pairs = get_pairs_table(sys);
    if(pairs == msys::TermTablePtr())
        VIPARR_FAIL("Unable to find pairs table.\n");

    msys::TermTablePtr exclusions = sys->table("exclusion");
    if (exclusions == msys::TermTablePtr()) {
        VIPARR_ERR << "WARNING: Missing exclusion table" << std::endl;
        return;
    }
    if (exclusions->params() == msys::ParamTablePtr()
            || exclusions->params()->propIndex("lj_scale") == msys::BadId)
        return;
    msys::TermTablePtr vdw = sys->table("nonbonded");
    if (vdw == msys::TermTablePtr()) {
        VIPARR_ERR << "WARNING: Missing nonbonded table" << std::endl;
        return;
    }
    msys::TermTablePtr vdw14 = sys->table("vdw1_14");
    std::map<std::string, Rules::VDWFunc>::const_iterator iter =
        Rules::VDWFuncRegistry().find(sys->nonbonded_info.vdw_funct);
    if (iter == Rules::VDWFuncRegistry().end()) {
        VIPARR_FAIL("Unsupported VDW function '" + sys->nonbonded_info.vdw_funct
                + "'; make sure the VDW function was copied from the "
                "forcefield to the system and is in the VDWFuncRegistry");
    }
    const Rules::VDWFunc& vdw_func = iter->second;
    const std::vector<std::string>& vdw_props = vdw_func.param_names;
    const std::vector<std::string>& pair_props = vdw_func.pair_param_names;
    std::string vdw_rule_name = sys->nonbonded_info.vdw_rule;
    std::map<std::string, Rules::VDWCombRulePtr>::const_iterator rules_iter
        = Rules::VDWCombRuleRegistry().find(vdw_rule_name);
    Rules::VDWCombRulePtr vdw_rule;
    std::stringstream missing_rule;
    if (rules_iter == Rules::VDWCombRuleRegistry().end()) {
        missing_rule << "Unsupported VDW combine rule '"
            << sys->nonbonded_info.vdw_rule
            << "'; make sure the VDW combine rule was copied from the "
            << "forcefield to the system and is in the VDWCombRuleRegistry"
            << std::endl;
    } else
        vdw_rule = rules_iter->second;

    msys::IdList terms = exclusions->terms();
    for (msys::Id term : terms) {
        if (exclusions->param(term) == msys::BadId) continue;
        double lj_scale = exclusions->propValue(term, "lj_scale").asFloat();
        if (lj_scale == 0) continue;
        msys::IdList atoms = exclusions->atoms(term);
        msys::Id row0 = msys::BadId;
        msys::TermTablePtr vdw_row0;
        msys::Id row1 = msys::BadId;
        msys::TermTablePtr vdw_row1;
        if (vdw14 != msys::TermTablePtr() &&
                separation_lookup(exclusions, atoms[0], atoms[1]) == 4) {
            vdw_row0 = vdw14;
            vdw_row1 = vdw14;
            row0 = vdw_lookup(vdw14, atoms[0], false);
            row1 = vdw_lookup(vdw14, atoms[1], false);
        }
        if (row0 == msys::BadId) {
            vdw_row0 = vdw;
            row0 = vdw_lookup(vdw, atoms[0], true);
        }
        if (row1 == msys::BadId) {
            vdw_row1 = vdw;
            row1 = vdw_lookup(vdw, atoms[1], true);
        }
        if (row0 == msys::BadId || row1 == msys::BadId)
            continue;
        std::vector<double> vi(vdw_props.size());
        for (unsigned p = 0; p < vdw_props.size(); ++p)
            vi[p] = vdw_row0->propValue(row0, vdw_props[p]).asFloat();
        std::vector<double> vj(vdw_props.size());
        for (unsigned p = 0; p < vdw_props.size(); ++p)
            vj[p] = vdw_row1->propValue(row1, vdw_props[p]).asFloat();
        /* Do not throw unrecognized comb rule error unless we actually need to
         * combine params */
        if (vdw_rule == Rules::VDWCombRulePtr())
            VIPARR_FAIL(missing_rule.str());
        std::vector<double> vcomb = (*vdw_rule)(vi, vj, lj_scale);
        if (vcomb.size() != pair_props.size()) {
            std::stringstream msg;
            msg << "VDW combine rule must return " << pair_props.size()
                << " properties; currently returns " << vcomb.size()
                << " properties" << std::endl;
            VIPARR_FAIL(msg.str());
        }
        std::map<std::string, double> params;
        for (unsigned p = 0; p < vcomb.size(); ++p)
            params.insert(std::make_pair(pair_props[p], vcomb[p]));
        update_pair_params(pairs, sys, atoms[0], atoms[1], params);
    }
}

static Forcefield::RegisterPlugin _("pairs_lj_scaled_14",
                                    match_pairs_lj_scaled_14,
                                    std::vector<std::string>(),
                                    compile_pairs_lj_scaled);
