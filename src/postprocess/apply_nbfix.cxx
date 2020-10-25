#include "compile_plugins.hxx"
#include "../add_system_tables.hxx"
#include "../ff.hxx"
#include "../util/util.hxx"
#include <msys/override.hxx>

void desres::viparr::ApplyNBFix(msys::SystemPtr sys) {
    if (sys->table("nonbonded") == msys::TermTablePtr()
            || sys->table("vdw2") == msys::TermTablePtr())
        return;
    msys::ParamTablePtr vdw2 = sys->table("vdw2")->params();
    msys::ParamTablePtr vdw1 = sys->table("nonbonded")->params();
    if (vdw1 == msys::ParamTablePtr())
        VIPARR_FAIL("nonbonded table is missing param table");
    if (vdw2 == msys::ParamTablePtr())
        VIPARR_FAIL("vdw2 term table is missing param table");
    if (sys->table("nonbonded")->overrides()->count() == 0)
        sys->table("nonbonded")->overrides()->resetParams(vdw2);
        /* If override count > 0, then the nonbonded override table should
         * already have been merged with the vdw2 table during application of
         * AddSystemTables */
    else {
        /* Reapply overrides that affect at least one previously parametrized
         * atom in the system, even if they are not specified in the new
         * forcefields */
        std::vector<msys::IdPair> override_pairs
            = sys->table("nonbonded")->overrides()->list();
        for (const msys::IdPair& pair : override_pairs) {
            sys->table("vdw2")->addTerm(msys::IdList(1, sys->atoms()[0]),
                    sys->table("nonbonded")->overrides()->get(pair));
        }
    }
    if (vdw2 != sys->table("nonbonded")->overrides()->params())
        VIPARR_FAIL("vdw2 term table does not point to the nonbonded overrides"
                " param table");
    msys::Id vdw1_nbfix_col = vdw1->propIndex("nbfix_identifier");
    if (vdw1_nbfix_col == msys::BadId)
        VIPARR_FAIL("Cannot apply NBFix: nonbonded table does not have "
                "'nbfix_identifier' column");
    if (vdw2->propIndex("nbfix_identifier") == msys::BadId)
        VIPARR_FAIL("Cannot apply NBFix: vdw2 table does not have "
                "'nbfix_identifier' column");
    if (vdw1->propIndex("type") == msys::BadId)
        VIPARR_FAIL("Cannot apply NBFix: nonbonded table does not have "
                "'type' column");
    if (vdw2->propIndex("type") == msys::BadId)
        VIPARR_FAIL("Cannot apply NBFix: vdw2 table does not have "
                "'type' column");
    msys::IdList vdw1_terms = sys->table("nonbonded")->terms();
    msys::IdList vdw2_terms = sys->table("vdw2")->terms();
    typedef std::map<std::string, msys::Id> Cache;
    Cache done;
    std::vector<bool> used(vdw1->paramCount(), false);
    for (msys::Id vdw1_term : vdw1_terms) {
        if (sys->table("nonbonded")->param(vdw1_term) != msys::BadId)
            used[sys->table("nonbonded")->param(vdw1_term)] = true;
    }
    for (msys::Id vdw2_term : vdw2_terms) {
        msys::Id vdw2_id = sys->table("vdw2")->param(vdw2_term);
        std::string nbfix_identifier = vdw2->value(vdw2_id,
                "nbfix_identifier").asString();
        std::string type = vdw2->value(vdw2_id, "type").asString();
        auto tokens = ViparrSplitString(type);
        if (tokens.size() != 2)
            VIPARR_FAIL("vdw2 parameters must have exactly two types");
        std::string type1 = tokens[0];
        std::string type2 = tokens[1];
        bool found = false;
        Cache::iterator iter = done.find(nbfix_identifier
                + " " + type1 + " " + type2);
        if (iter != done.end()) {
            found = true;
            for (unsigned j = 0; j < vdw2->propCount(); ++j) {
                if (vdw2->propName(j) == "type" || vdw2->propName(j) == "memo")
                    continue;
                if (vdw2->value(vdw2_id, j) != vdw2->value(iter->second, j)) {
                    std::stringstream msg;
                    msg << "Rows " << iter->second << " and " << vdw2_id
                        << " of global vdw2 table override the same types ("
                        << type1 << ", " << type2
                        << ") but have different parameters";
                    VIPARR_FAIL(msg.str());
                }
            }
        }
        iter = done.find(nbfix_identifier + " " + type2 + " " + type1);
        if (iter != done.end()) {
            found = true;
            for (unsigned j = 0; j < vdw2->propCount(); ++j) {
                if (vdw2->propName(j) == "type" || vdw2->propName(j) == "memo")
                    continue;
                if (vdw2->value(vdw2_id, j) != vdw2->value(iter->second, j)) {
                    std::stringstream msg;
                    msg << "Rows " << iter->second << " and " << vdw2_id
                        << " of global vdw2 table override the same types ("
                        << type1 << ", " << type2
                        << ") but have different parameters";
                    VIPARR_FAIL(msg.str());
                }
            }
        }
        if (!found) {
            done.insert(std::make_pair(nbfix_identifier + " " + type1
                        + " " + type2, vdw2_id));
            done.insert(std::make_pair(nbfix_identifier + " " + type2
                        + " " + type1, vdw2_id));
            msys::IdList all_vdw1_params = vdw1->findString(vdw1_nbfix_col,
                    nbfix_identifier);
            msys::IdList vdw1_params_1;
            msys::IdList vdw1_params_2;
            for (msys::Id p : all_vdw1_params) {
                /* Only override this parameter if it is used by the nonbonded
                 * term table of this system */
                if (used[p] && vdw1->value(p, "type").asString() == type1)
                    vdw1_params_1.push_back(p);
                if (used[p] && vdw1->value(p, "type").asString() == type2)
                    vdw1_params_2.push_back(p);
            }
            for (unsigned j = 0; j < vdw1_params_1.size(); ++j) {
                for (unsigned k = 0; k < vdw1_params_2.size(); ++k) {
                    sys->table("nonbonded")->overrides()->set(
                            std::make_pair(vdw1_params_1[j], vdw1_params_2[k]),
                            vdw2_id);
                }
            }
        }
    }
}
