#include "wrap_obj.hxx"
#include "../src/importexport/export_ff.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {

    void export_templates(const object& py_tpls, const std::string& path) {
        list L(py_tpls);
        std::vector<TemplatedSystemPtr> tpls(len(L));
        for (unsigned i = 0; i < tpls.size(); ++i)
            tpls[i] = extract<TemplatedSystemPtr>(L[i]);
        ExportTemplates(tpls, path);
    }

    void export_cmap(const object& py_cmaps, const std::string& path) {
        list L(py_cmaps);
        std::vector<msys::ParamTablePtr> cmaps(len(L));
        for (unsigned i = 0; i < cmaps.size(); ++i)
            cmaps[i] = extract<msys::ParamTablePtr>(L[i]);
        ExportCmap(cmaps, path);
    }

    void export_params(msys::ParamTablePtr table, const object& py_rows, const std::string& path) {
        list L(py_rows);
        std::list<msys::Id> rows;
        for (unsigned i = 0; i < len(L); ++i)
            rows.push_back(extract<msys::Id>(L[i]));
        ExportParams(table, rows, path);
    }

}

namespace desres { namespace viparr {

    void export_exports() {
        def("ExportForcefield", &ExportForcefield);
        def("ExportRules", &ExportRules);
        def("ExportTemplates", export_templates);
        def("ExportCmap", export_cmap);
        def("ExportParams", export_params);
    }

}}
