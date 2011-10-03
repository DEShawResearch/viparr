#include "wrap_obj.hxx"
#include "../src/merge_ff.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {

    list merge_params(const object& py_src_rows, const object& py_patch_rows, msys::ParamTablePtr table, bool append_only, bool verbose) {
        list L(py_src_rows);
        std::list<msys::Id> src_rows;
        for (unsigned i = 0; i < len(L); ++i)
            src_rows.push_back(extract<msys::Id>(py_src_rows[i]));
        L = list(py_patch_rows);
        std::list<msys::Id> patch_rows;
        for (unsigned i = 0; i < len(L); ++i)
            patch_rows.push_back(extract<msys::Id>(py_patch_rows[i]));
        std::list<msys::Id> merged = MergeParams(src_rows, patch_rows, table, append_only, verbose);
        L = list();
        for (msys::Id param : merged)
            L.append(param);
        return L;
    }
}

namespace desres { namespace viparr {

    void export_merges() {
        def("MergeForcefields", &MergeForcefields);
        def("MergeRules", &MergeRules);
        def("MergeTemplates", &MergeTemplates);
        def("MergeParams", merge_params);
    }

}}
