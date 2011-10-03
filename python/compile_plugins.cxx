#include "wrap_obj.hxx"
#include "../src/postprocess/compile_plugins.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {
    void compile_plugins(const msys::SystemPtr sys, const object& obj) {
        list L(obj);
        std::set<std::string> plugin_names;
        for(unsigned i = 0; i < len(L); ++i) {
            std::string name = extract<std::string>(L[i]);
            plugin_names.insert(name);
        }
        CompilePlugins(sys, plugin_names);
    }
}

namespace desres { namespace viparr {
    void export_compile_plugins() {
        def("CompilePlugins", compile_plugins);
        def("AddPairsTable", AddPairsTable);
        def("ApplyNBFix", ApplyNBFix);
        def("CleanupSystem", CleanupSystem);
    }
}}
