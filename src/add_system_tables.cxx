#include "add_system_tables.hxx"
#include "append_params.hxx"
#include "ff.hxx"
#include <msys/override.hxx>

using namespace desres;
using namespace desres::viparr;

void desres::viparr::AddSystemTables(msys::SystemPtr sys,
        const std::map<std::string, std::list<msys::Id> >& share_params,
        const std::map<std::string, std::string>& name_conversions) {

    std::vector<std::string> tables = sys->tableNames();
    for (unsigned i = 0; i < tables.size(); ++i) {
        if (sys->table(tables[i])->params() == msys::ParamTablePtr())
            continue;
        msys::TermTablePtr term_table = sys->table(tables[i]);
        msys::IdList terms = term_table->terms();

        /* Convert term table name to global param table name */
        std::string name = tables[i];
        for (std::map<std::string, std::string>::const_iterator iter
                = name_conversions.begin(); iter != name_conversions.end();
                ++iter) {
            size_t loc = name.find(iter->first);
            if (loc != std::string::npos)
                name.replace(loc, iter->first.size(), iter->second);
        }

        msys::IdList rowIDs;
        msys::IdList new_param_ids;
        if (!Forcefield::HasParamTable(name) || term_table->params()
                != Forcefield::ParamTable(name)) {
            /* Append parameters of system to global table */
            std::map<std::string, std::list<msys::Id> >::const_iterator
                share_iter = share_params.find(name);
            if (share_iter == share_params.end())
                rowIDs = append_params::AppendParams<msys::IdList>(
                        term_table->params(), name);
            else
                rowIDs = append_params::AppendParams<msys::IdList>(
                        term_table->params(), name, share_iter->second);
            /* Get new term table param IDs */
            new_param_ids = msys::IdList(term_table->maxTermId(), msys::BadId);
            for (unsigned j = 0; j < terms.size(); ++j)
                if (term_table->param(terms[j]) != msys::BadId)
                    new_param_ids[terms[j]] = rowIDs.at(
                            term_table->param(terms[j]));
        }

        /* Convert name of override table */
        std::string overrides_name = tables[i];
        if (overrides_name == "nonbonded" || overrides_name == "vdw1")
            overrides_name = "vdw2";
        else
            overrides_name += "_overrides";

        std::vector<msys::IdPair> new_override_pairs;
        msys::IdList new_override_params;
        if (term_table->overrides()->count() > 0 &&
                (new_param_ids.size() > 0 ||
                 !Forcefield::HasParamTable(overrides_name) ||
                 term_table->overrides()->params()
                 != Forcefield::ParamTable(overrides_name))) {
            /* Append parameters of system's override table to global table */
            msys::IdList or_rowIDs;
            std::map<std::string, std::list<msys::Id> >::const_iterator
                share_iter = share_params.find(overrides_name);
            if (share_iter == share_params.end())
                or_rowIDs = append_params::AppendParams<msys::IdList>(
                        term_table->overrides()->params(), overrides_name);
            else
                or_rowIDs = append_params::AppendParams<msys::IdList>(
                        term_table->overrides()->params(), overrides_name,
                        share_iter->second);
            /* Get new override param and pair IDs */
            std::vector<msys::IdPair> pairs = term_table->overrides()->list();
            for (unsigned j = 0; j < pairs.size(); ++j) {
                if (rowIDs.size() > 0)
                    new_override_pairs.push_back(std::make_pair(rowIDs.at(
                                    pairs[j].first), rowIDs.at(pairs[j].second)));
                else
                    new_override_pairs.push_back(pairs[j]);
                new_override_params.push_back(or_rowIDs.at(
                            term_table->overrides()->get(pairs[j])));
            }
        }

        /* Reset term table params and overrides */
        if (new_param_ids.size() > 0) {
            term_table->resetParams(Forcefield::ParamTable(name));
            for (unsigned j = 0; j < terms.size(); ++j)
                term_table->setParam(terms[j], new_param_ids[terms[j]]);
        }
        if (new_override_pairs.size() > 0) {
            term_table->overrides()->resetParams(
                    Forcefield::ParamTable(overrides_name));
            for (unsigned j = 0; j < new_override_pairs.size(); ++j)
                term_table->overrides()->set(new_override_pairs[j],
                        new_override_params[j]);
        }
    }
}
