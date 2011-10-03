#ifndef desres_viparr_append_params_hxx
#define desres_viparr_append_params_hxx

#include "ff.hxx"
#include <msys/system.hxx>

namespace desres { namespace viparr { namespace append_params {

    typedef std::pair<bool, msys::Id> ParamToken;
    struct ParamComparator {
        msys::ParamTablePtr _dest;
        msys::ParamTablePtr _src;
        const msys::IdList& _d2s;
        ParamComparator(msys::ParamTablePtr dest, msys::ParamTablePtr src,
                const msys::IdList& dest_to_src);
        bool operator()(const ParamToken& pi, const ParamToken& pj) const;
    };

    /* Helper function to append parameters from a param table to a static
     * Forcefield table */
    template <class Container>
    Container AppendParams(msys::ParamTablePtr src_table,
            const std::string& dest_table_name,
            const std::list<msys::Id>& share_params=std::list<msys::Id>()) {
        /* If the param table is already a global param table, do
         * nothing */
        if (Forcefield::HasParamTable(dest_table_name) && src_table
                == Forcefield::ParamTable(dest_table_name)) {
            msys::IdList params = src_table->params();
            return Container(params.begin(), params.end());
        }

        /* Add new global param table, or check for consistency of columns */
        msys::ParamTablePtr ptable;
        msys::IdList prop_map;
        if (!Forcefield::HasParamTable(dest_table_name)) {
            ptable = msys::ParamTable::create();
            for (unsigned j = 0; j < src_table->propCount(); ++j) {
                ptable->addProp(src_table->propName(j), src_table->propType(j));
                prop_map.push_back(j);
            }
            Forcefield::AddParamTable(dest_table_name, ptable);
        } else {
            ptable = Forcefield::ParamTable(dest_table_name);
            for (unsigned j = 0; j < ptable->propCount(); ++j) {
                msys::Id index = src_table->propIndex(ptable->propName(j));
                if (index == msys::BadId) {
                    VIPARR_ERR << "WARNING: Property '"
                        << ptable->propName(j) << "' of existing table '"
                        << dest_table_name << "' not found in new table; setting "
                        << "to default values\n";
                    prop_map.push_back(msys::BadId);
                } else
                    prop_map.push_back(index);
            }
            for (unsigned j = 0; j < src_table->propCount(); ++j) {
                if (ptable->propIndex(src_table->propName(j)) == msys::BadId) {
                    VIPARR_ERR << "WARNING: Property '"
                        << src_table->propName(j) << "' of new table not found in "
                        << "existing table '" << dest_table_name
                        << "'; setting to default values\n";
                    ptable->addProp(src_table->propName(j), src_table->propType(j));
                    prop_map.push_back(j);
                }
            }
        }

        typedef std::map<ParamToken, msys::Id, ParamComparator> ParamMap;
        ParamComparator comp(ptable, src_table, prop_map);
        ParamMap cache(comp);
        /* Loop through dest params and cache unique rows */
        for (msys::Id param : share_params) {
            if (cache.find(ParamToken(true,param)) == cache.end())
                cache.insert(std::make_pair(ParamToken(true,param), param));
        }

        /* Merge param table rows */
        Container rows;
        for (unsigned j = 0; j < src_table->paramCount(); ++j) {
            if (share_params.size() > 0) {
                ParamMap::iterator iter = cache.find(ParamToken(false,j));
                if (iter != cache.end()) {
                    rows.push_back(iter->second);
                    continue;
                }
            }
            msys::Id paramid = ptable->addParam();
            for (unsigned k = 0; k < prop_map.size(); ++k) {
                if (prop_map[k] != msys::BadId) {
                    ptable->value(paramid, k) = src_table->value(j,
                            prop_map[k]);
                }
            }
            rows.push_back(paramid);
            if (share_params.size() > 0) {
                cache.insert(std::make_pair(ParamToken(true, paramid),
                            paramid));
            }
        }
        return rows;
    }

}}}

#endif
