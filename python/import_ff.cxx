#include "wrap_obj.hxx"
#include "../src/importexport/import_ff.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {

    list import_templates(const std::string& path) {
        list L;
        std::vector<TemplatedSystemPtr> tpls = ImportTemplates(path);
        for (unsigned i = 0; i < tpls.size(); ++i)
            L.append(tpls[i]);
        return L;
    }

    list import_cmap(const std::string& path) {
        list L;
        std::vector<desres::msys::ParamTablePtr> cmaps = ImportCmap(path);
        for (unsigned i = 0; i < cmaps.size(); ++i)
            L.append(cmaps[i]);
        return L;
    }

    list import_params(const std::string& table_name, const std::string& path, const object& py_share_params, const std::string& nbfix_identifier) {
        list L(py_share_params);
        std::list<msys::Id> share_params;
        for (unsigned i = 0; i < len(L); ++i)
            share_params.push_back(extract<msys::Id>(L[i]));
        std::list<msys::Id> params = ImportParams(table_name, path, share_params, nbfix_identifier);
        L = list();
        for (msys::Id param : params)
            L.append(param);
        return L;
    }

    ForcefieldPtr import_forcefield(const std::string& dir, bool require_rules, const object& py_share_params) {
        dict D(py_share_params);
        list L(D.keys());
        std::map<std::string, std::list<msys::Id> > share_params;
        for (unsigned i = 0; i < len(L); ++i) {
            list py_vals(D[L[i]]);
            std::list<msys::Id> vals;
            for (unsigned j = 0; j < len(py_vals); ++j)
                vals.push_back(extract<msys::Id>(py_vals[j]));
            std::string key = extract<std::string>(L[i]);
            share_params.insert(std::make_pair(key, vals));
        }
        return ImportForcefield(dir, require_rules, share_params);
    }

}

namespace desres { namespace viparr {

    void export_imports() {
        def("ImportForcefield", import_forcefield);
        def("ImportRules", &ImportRules);
        def("ImportTemplates", import_templates);
        def("ImportCmap", import_cmap);
        def("ImportParams", import_params);
    }

}}

