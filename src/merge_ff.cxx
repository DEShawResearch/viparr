#include "base.hxx"
#include "merge_ff.hxx"
#include <vector>

namespace desres { namespace viparr {

    void MergeForcefields(ForcefieldPtr src, ForcefieldPtr patch,
            bool append_only, bool verbose) {
        if (src->rules() == RulesPtr()
                || patch->rules() == RulesPtr()
                || src->typer() == TemplateTyperPtr()
                || patch->typer() == TemplateTyperPtr())
            VIPARR_FAIL("Forcefields being merged must have valid Rules and "
                    "Typer objects");
        MergeRules(src->rules(), patch->rules(), verbose);
        MergeTemplates(src->typer(), patch->typer(), append_only, verbose);
        std::vector<std::string> tables = patch->paramTables();
        for (unsigned i = 0; i < tables.size(); ++i) {
            if (verbose)
                VIPARR_OUT << "Merging param table " << tables[i] << std::endl;
            msys::ParamTablePtr ptable = Forcefield::ParamTable(tables[i]);
            std::list<msys::Id> merged_rows = MergeParams(
                    src->rowIDs(tables[i]), patch->rowIDs(tables[i]), ptable,
                    append_only, verbose);
            src->clearParams(tables[i]);
            src->appendParams(tables[i], merged_rows);
        }
        if (patch->cmapTables().size() > 0) {
            bool same = (patch->cmapTables().size()
                    == src->cmapTables().size());
            if (same) {
                for (unsigned i = 1; i <= patch->cmapTables().size(); ++i) {
                    if (patch->cmapTable(i) != src->cmapTable(i)) {
                        same = false;
                        break;
                    }
                }
            }
            if (!same) {
                if (append_only && src->cmapTables().size() > 0)
                    VIPARR_FAIL("Merge error: Cannot overwrite SMARTS "
                            "atomtype tree in append-only mode");
                if (verbose)
                    VIPARR_OUT << "Replacing cmap tables" << std::endl;
                src->delCmapTables();
                for (unsigned i = 1; i <= patch->cmapTables().size(); ++i)
                    src->addCmapTable(patch->cmapTable(i));
            }
        }
        if (append_only)
            src->name += "__-a_" + patch->name;
        else
            src->name += "__-m_" + patch->name;
    }

    void MergeRules(RulesPtr src_rules, RulesPtr patch_rules, bool verbose) {
        if (patch_rules->info.size() == 0
                && patch_rules->vdw_func == ""
                && patch_rules->vdw_comb_rule == ""
                && patch_rules->plugins.size() == 0
                && patch_rules->exclusions() == 1
                && patch_rules->fatal)
            return;
        if (verbose)
            VIPARR_OUT << "Merging rules" << std::endl;
        src_rules->vdw_func = Rules::MergeVDW(src_rules->vdw_func,
                patch_rules->vdw_func);
        src_rules->vdw_comb_rule = Rules::MergeVDW(
                src_rules->vdw_comb_rule, patch_rules->vdw_comb_rule);
        if (patch_rules->exclusions() != 1) {
            if (src_rules->exclusions() != patch_rules->exclusions())
                VIPARR_FAIL("Cannot merge rules: Exclusion rules do not match");
            for (unsigned i = 2; i <= src_rules->exclusions(); ++i) {
                if (src_rules->es_scale(i) != patch_rules->es_scale(i))
                    VIPARR_FAIL("Cannot merge rules: ES scales do not match");
                if (src_rules->lj_scale(i) != patch_rules->lj_scale(i))
                    VIPARR_FAIL("Cannot merge rules: LJ scales do not match");
            }
        }
        if (src_rules->fatal && !patch_rules->fatal) {
            if (verbose)
                VIPARR_OUT << "...Setting fatal flag to false" << std::endl;
            src_rules->fatal = false;
        }
        if (src_rules->nbfix_identifier == ""
                && patch_rules->nbfix_identifier != "") {
            if (verbose)
                VIPARR_OUT << "...Setting NBFix identifier" << std::endl;
            src_rules->nbfix_identifier = patch_rules->nbfix_identifier;
        } else if (src_rules->nbfix_identifier != patch_rules->nbfix_identifier
                && patch_rules->nbfix_identifier != "") {
            VIPARR_FAIL("Cannot merge rules: NBFix identifiers do not match");
        }
        /* Take union of info strings */
        for (unsigned i = 0; i < patch_rules->info.size(); ++i) {
            bool found = false;
            for (unsigned j = 0; j < src_rules->info.size(); ++j) {
                if (src_rules->info[j] == patch_rules->info[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                if (verbose)
                    VIPARR_OUT << "...Adding info '" << patch_rules->info[i]
                        << "'" << std::endl;
                src_rules->info.push_back(patch_rules->info[i]);
            }
        }
        /* Take union of plugins */
        for (unsigned i = 0; i < patch_rules->plugins.size(); ++i) {
            bool found = false;
            for (unsigned j = 0; j < src_rules->plugins.size(); ++j) {
                if (src_rules->plugins[j] == patch_rules->plugins[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                if (verbose) {
                    VIPARR_OUT << "...Adding plugin '"
                        << patch_rules->plugins[i] << "'" << std::endl;
                }
                src_rules->plugins.push_back(patch_rules->plugins[i]);
            }
        }
    }

    void MergeTemplates(TemplateTyperPtr src_typer, TemplateTyperPtr
            patch_typer, bool append_only, bool verbose) {
        if (patch_typer->templates().size() == 0)
            return;
        if (verbose)
            VIPARR_OUT << "Merging templates" << std::endl;
        typedef std::map<std::string, std::vector<TemplatedSystemPtr> > TplMap;
        TplMap src_tpls;
        std::vector<TemplatedSystemPtr> templates = src_typer->templates();
        for (unsigned i = 0; i < templates.size(); ++i) {
            std::string name = templates[i]->system()->residue(0).name;
            if (src_tpls.find(name) == src_tpls.end())
                src_tpls[name] = std::vector<TemplatedSystemPtr>();
            src_tpls[name].push_back(templates[i]);
        }
        TplMap patch_tpls;
        templates = patch_typer->templates();
        for (unsigned i = 0; i < templates.size(); ++i) {
            std::string name = templates[i]->system()->residue(0).name;
            if (patch_tpls.find(name) == patch_tpls.end())
                patch_tpls[name] = std::vector<TemplatedSystemPtr>();
            patch_tpls[name].push_back(templates[i]);
        }
        /* If patch has template(s) of the same name as template(s) in src,
         * replace those in src with those in patch. Take union of remaining
         * templates. */
        for (TplMap::iterator p_iter = patch_tpls.begin();
                p_iter != patch_tpls.end(); ++p_iter) {
            TplMap::iterator s_iter = src_tpls.find(p_iter->first);
            if (s_iter != src_tpls.end()) {
                if (s_iter->second != p_iter->second) {
                    if (append_only) {
                        VIPARR_FAIL("Merge error: Cannot overwrite template " +
                                p_iter->first + " in append-only mode");
                    }
                    if (verbose) {
                        VIPARR_OUT << "...Overwriting template(s) "
                            << p_iter->first << std::endl;
                    }
                    for (unsigned i = 0; i < s_iter->second.size(); ++i)
                        src_typer->delTemplate(s_iter->second[i]);
                    for (unsigned i = 0; i < p_iter->second.size(); ++i)
                        src_typer->addTemplate(p_iter->second[i]);
                }
            } else {
                if (verbose) {
                    VIPARR_OUT << "...Adding template(s) " << p_iter->first
                        << std::endl;
                }
                for (unsigned i = 0; i < p_iter->second.size(); ++i)
                    src_typer->addTemplate(p_iter->second[i]);
            }
        }
    }

    std::list<msys::Id> MergeParams(const std::list<msys::Id>& src_rows,
            const std::list<msys::Id>& patch_rows, msys::ParamTablePtr table,
            bool append_only, bool verbose) {
        if (table->propIndex("type") == msys::BadId)
            VIPARR_FAIL("Merging params requires table to have a 'type'"
                   " column");
        typedef std::map<std::string, msys::IdList> TypeMap;
        TypeMap src_types;
        for (msys::Id row : src_rows) {
            std::string type = table->value(row, "type").asString();
            if (src_types.find(type) == src_types.end())
                src_types[type] = msys::IdList();
            src_types[type].push_back(row);
        }
        /* Add rows from patch with 'type' not found in src*/
        TypeMap patch_types;
        std::list<msys::Id> final_ids;
        for (msys::Id row : patch_rows) {
            std::string type = table->value(row, "type").asString();
            if (patch_types.find(type) == patch_types.end())
                patch_types[type] = msys::IdList();
            patch_types[type].push_back(row);
            if (src_types.find(type) == src_types.end()) {
                if (verbose)
                    VIPARR_OUT << "...Adding type '" << type << "'" << std::endl;
                final_ids.push_back(row);
            }
        }
        /* Add rows of src */
        std::set<std::string> done;
        for (msys::Id row : src_rows) {
            std::string type = table->value(row, "type").asString();
            TypeMap::iterator p_iter = patch_types.find(type);
            if (p_iter != patch_types.end()) {
                /* src 'type' found in patch -- add the version in patch */
                if (done.find(type) != done.end())
                    continue;
                /* check if params are different; if so, inform the user */
                bool all_same = true;
                TypeMap::iterator s_iter = src_types.find(type);
                if (s_iter->second.size() != p_iter->second.size())
                    all_same = false;
                else for (unsigned j = 0; j < p_iter->second.size(); ++j) {
                    bool j_matched = false;
                    for (unsigned k = 0; k < s_iter->second.size(); ++k) {
                        bool jk_match = true;
                        for (msys::Id prop = 0;
                                prop < table->propCount(); ++prop) {
                            if (table->propName(prop) == "memo") continue;
                            if (table->value(p_iter->second[j], prop)
                                    != table->value(s_iter->second[k], prop)) {
                                jk_match = false;
                                break;
                            }
                        }
                        if (jk_match) {
                            j_matched = true;
                            break;
                        }
                    }
                    if (!j_matched) {
                        all_same = false;
                        break;
                    }
                }
                if (!all_same) {
                    if (append_only)
                        VIPARR_FAIL("Merge error: Cannot overwrite type "
                                + type + " in append_only mode");
                    if (verbose)
                        VIPARR_OUT << "...Overwriting type '" << type << "'"
                            << std::endl;
                }
                final_ids.insert(final_ids.end(), p_iter->second.begin(),
                        p_iter->second.end());
                done.insert(type);
            } else {
                /* src 'type' not found in patch -- add from src */
                final_ids.push_back(row);
            }
        }
        return final_ids;
    }
}}
