#include "../base.hxx"
#include "../rules.hxx"
#include <msys/term_table.hxx>

using namespace desres;
using namespace desres::viparr;

namespace {

    inline int separation_lookup(msys::TermTablePtr table, msys::Id a0, msys::Id a1) {
        msys::IdList atoms(2);
        atoms[0] = a0;
        atoms[1] = a1;
        msys::IdList rows = table->findWithAll(atoms);
        if (rows.size() > 0)
            return table->propValue(rows[0], "separation").asInt();
        else
            /* Indicates separation > 4 */
            return 5;
    }

    inline msys::Id vdw_lookup(msys::TermTablePtr table, msys::Id atom,
            bool required) {
        msys::IdList rows = table->findExact(msys::IdList(1, atom));
        if (rows.size() > 0) return rows[0];
        if (required)
            VIPARR_ERR << "WARNING: Missing vdw term for atom " << atom
                << std::endl;
        return msys::BadId;
    }

    inline msys::TermTablePtr get_pairs_table(msys::SystemPtr sys) {
        std::map<std::string, Rules::VDWFunc>::const_iterator iter
            = Rules::VDWFuncRegistry().find(sys->nonbonded_info.vdw_funct);
        if (iter == Rules::VDWFuncRegistry().end()) {
            VIPARR_FAIL("Unsupported VDW function '" + sys->nonbonded_info.vdw_funct
                        + "'; make sure the VDW function was copied from the "
                        "forcefield to the system and is in the VDWFuncRegistry");
        }
        Rules::VDWFunc vdw_func = iter->second;
        std::string name = vdw_func.pair_table_name;

        if(sys->table(name) != msys::TermTablePtr())
            return sys->table(name);

        VIPARR_ERR << "Expected to find pairs table with name " << name <<
            ", but none such exists.\n";
        return msys::TermTablePtr();
                
}

    inline void update_pair_params(msys::TermTablePtr table, msys::SystemPtr sys,
            msys::Id ai, msys::Id aj,
            const std::map<std::string, double>& params) {
        msys::ParamTablePtr param_table = table->params();
        msys::IdList term(2);
        if (ai < aj) {
            term[0] = ai; term[1] = aj;
        } else {
            term[0] = aj; term[1] = ai;
        }
        msys::IdList rows = table->findWithAll(term);
        if (rows.size() == 0) {
            msys::Id param_row = param_table->addParam();
            rows.push_back(table->addTerm(term, param_row));
        }
        /* Update values */
        for (std::map<std::string, double>::const_iterator iter
                = params.begin(); iter != params.end(); ++iter) {
            msys::Id param_row = table->param(rows[0]);
            if (param_table->refcount(param_row) > 1) {
                param_row = param_table->duplicate(param_row);
                table->setParam(rows[0], param_row);
            }
            table->propValue(rows[0], iter->first) = iter->second;
        }
    }
}
