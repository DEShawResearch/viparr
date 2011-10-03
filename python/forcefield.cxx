#include "wrap_obj.hxx"
#include "../src/ff.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {

    list all_param_tables() {
        list L;
        std::vector<std::string> names = Forcefield::AllParamTables();
        for (unsigned i = 0; i < names.size(); ++i)
            L.append(names[i]);
        return L;
    }

    list filter_params_string(const std::string& table_name, const object& container, const std::string& key, const std::string& value) {
        list L(container);
        std::list<msys::Id> params;
        for (unsigned i = 0; i < len(L); ++i)
            params.push_back(extract<msys::Id>(L[i]));
        std::list<msys::Id> filtered = Forcefield::FilterParams<std::list<msys::Id> >(table_name, params, key, value);
        list L_out;
        for (msys::Id param : filtered)
            L_out.append(param);
        return L_out;
    }

    list filter_params_int(const std::string& table_name, const object& container, const std::string& key, int value) {
        list L(container);
        std::list<msys::Id> params;
        for (unsigned i = 0; i < len(L); ++i)
            params.push_back(extract<msys::Id>(L[i]));
        std::list<msys::Id> filtered = Forcefield::FilterParams<std::list<msys::Id> >(table_name, params, key, value);
        list L_out;
        for (msys::Id param : filtered)
            L_out.append(param);
        return L_out;
    }

    list filter_params_float(const std::string& table_name, const object& container, const std::string& key, double value) {
        list L(container);
        std::list<msys::Id> params;
        for (unsigned i = 0; i < len(L); ++i)
            params.push_back(extract<msys::Id>(L[i]));
        std::list<msys::Id> filtered = Forcefield::FilterParams<std::list<msys::Id> >(table_name, params, key, value);
        list L_out;
        for (msys::Id param : filtered)
            L_out.append(param);
        return L_out;
    }

    dict get_plugin_registry() {
        dict D;
        for (std::map<std::string, Forcefield::PluginPtr>::const_iterator iter = Forcefield::PluginRegistry().begin(); iter != Forcefield::PluginRegistry().end(); ++iter)
            D[iter->first] = iter->second;
        return D;
    }

    void set_plugin_registry(const object& obj) {
        dict D(obj);
        list L(D.keys());
        Forcefield::PluginRegistry().clear();
        for (unsigned i = 0; i < len(L); ++i) {
            std::string key = extract<std::string>(L[i]);
            Forcefield::PluginPtr plugin = extract<Forcefield::PluginPtr>(D[L[i]]);
            Forcefield::PluginRegistry().insert(std::make_pair(key, plugin));
        }
    }

    list row_ids(const Forcefield& ff, const std::string& name) {
        const std::list<msys::Id>& rowIDs = ff.rowIDs(name);
        list L;
        for (msys::Id row : rowIDs)
            L.append(row);
        return L;
    }

    void del_params(Forcefield& ff, const std::string& name, const object& obj) {
        list L(obj);
        std::list<msys::Id> params;
        for (unsigned i = 0; i < len(L); ++i)
            params.push_back(extract<msys::Id>(L[i]));
        ff.delParams(name, params);
    }

    void append_params(Forcefield& ff, const std::string& name, const object& obj) {
        list L(obj);
        std::list<msys::Id> params;
        for (unsigned i = 0; i < len(L); ++i)
            params.push_back(extract<msys::Id>(L[i]));
        ff.appendParams(name, params);
    }

    list param_tables(const Forcefield& ff) {
        list L;
        std::vector<std::string> tables = ff.paramTables();
        for (unsigned i = 0; i < tables.size(); ++i)
            L.append(tables[i]);
        return L;
    }

#if 0
    dict set_mode(Forcefield& ff, const std::string& mode) {
        std::map<std::string, msys::IdList> changed = ff.setMode(mode);
        dict D;
        for (std::map<std::string, msys::IdList>::iterator iter = changed.begin(); iter != changed.end(); ++iter) {
            list L;
            for (unsigned i = 0; i < iter->second.size(); ++i)
                L.append(iter->second[i]);
            D[iter->first] = L;
        }
        return D;
    }
#endif

    list cmap_tables(const Forcefield& ff) {
        list L;
        const std::vector<msys::ParamTablePtr>& cmaps = ff.cmapTables();
        for (unsigned i = 0; i < cmaps.size(); ++i)
            L.append(cmaps[i]);
        return L;
    }

    /* Python implementation of Plugin interface */
    class PluginPy : public Forcefield::Plugin {
        public:
            static Forcefield::PluginPtr create(const object& py_match,
                                                const object& py_compile) {
                return Forcefield::PluginPtr(new PluginPy(py_match, py_compile));
            }
            virtual void match(TemplatedSystemPtr tsys, ForcefieldPtr ff) const {
                try {
                    _py_match(tsys, ff);
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error calling Python plugin: "
                            + std::string(e.what()));
                }
            }

            virtual void compile(msys::SystemPtr sys) const {
                try {
                    _py_compile(sys);
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error in compile step of Python plugin: "
                                + std::string(e.what()));
                }
            }
        private:
            PluginPy(const object& py_match, const object& py_compile) : _py_match(py_match),
                                                                         _py_compile(py_compile)
                     { }
            object _py_match;
            object _py_compile;
    };
    list get_prerequisites(const Forcefield::Plugin& plugin) {
        list L;
        for (unsigned i = 0; i < plugin.prerequisites.size(); ++i)
            L.append(plugin.prerequisites[i]);
        return L;
    }
    void set_prerequisites(Forcefield::Plugin& plugin, const object& obj) {
        list L(obj);
        plugin.prerequisites.resize(len(L));
        for (unsigned i = 0; i < plugin.prerequisites.size(); ++i)
            plugin.prerequisites[i] = extract<std::string>(L[i]);
    }

    void clear_plugins() {
        Forcefield::PluginRegistry().clear();
    }
}

namespace desres { namespace viparr {

    void export_forcefield() {

        scope forcefield_class(class_<Forcefield, ForcefieldPtr>("Forcefield", no_init)
            .def("__eq__", _eq<ForcefieldPtr>)
            .def("__ne__", _ne<ForcefieldPtr>)
            .def("__hash__", _hash<ForcefieldPtr>)
            .def("HasParamTable", &Forcefield::HasParamTable)
            .staticmethod("HasParamTable")
            .def("ParamTable", &Forcefield::ParamTable)
            .staticmethod("ParamTable")
            .def("AddParamTable", &Forcefield::AddParamTable)
            .staticmethod("AddParamTable")
            .def("ClearParamTables", &Forcefield::ClearParamTables)
            .staticmethod("ClearParamTables")
            .def("AllParamTables", all_param_tables)
            .staticmethod("AllParamTables")
            .def("FilterParamsString", filter_params_string)
            .staticmethod("FilterParamsString")
            .def("FilterParamsInt", filter_params_int)
            .staticmethod("FilterParamsInt")
            .def("FilterParamsFloat", filter_params_float)
            .staticmethod("FilterParamsFloat")
            .add_static_property("PluginRegistry", get_plugin_registry, set_plugin_registry)
            .def("__init__", make_constructor(&Forcefield::create))
            .def("Copy", &Forcefield::copy)
            .staticmethod("Copy")
            .def_readwrite("name", &Forcefield::name)
            .def("rules", &Forcefield::rules)
            .def("resetRules", &Forcefield::resetRules)
            .def("typer", &Forcefield::typer)
            .def("resetTyper", &Forcefield::resetTyper)
            .def("rowIDs", row_ids)
            .def("delParam", &Forcefield::delParam)
            .def("delParams", del_params)
            .def("clearParams", &Forcefield::clearParams)
            .def("appendParam", &Forcefield::appendParam)
            .def("appendParams", append_params)
            .def("replaceParam", &Forcefield::replaceParam)
            .def("paramTables", param_tables)
#if 0
            .def("getModeParam", &Forcefield::getModeParam)
            .def("createModeParam", &Forcefield::createModeParam)
            .def("setMode", set_mode)
#endif
            .def("cmapTable", &Forcefield::cmapTable)
            .def("cmapTables", cmap_tables)
            .def("addCmapTable", &Forcefield::addCmapTable)
            .def("ClearPlugins", clear_plugins)
            .staticmethod("ClearPlugins")
        );

        class_<Forcefield::Plugin, Forcefield::PluginPtr, boost::noncopyable>("Plugin", no_init)
            .def("__eq__", _eq<Forcefield::PluginPtr>)
            .def("__ne__", _ne<Forcefield::PluginPtr>)
            .def("__hash__", _hash<Forcefield::PluginPtr>)
            .def("__init__", make_constructor(&PluginPy::create))
            .def("match", &Forcefield::Plugin::match)
            .def("compile", &Forcefield::Plugin::compile)
            .def("getPrerequisites", get_prerequisites)
            .def("setPrerequisites", set_prerequisites)
            ;
    }

}}

