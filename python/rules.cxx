#include "wrap_obj.hxx"
#include "../src/base.hxx"
#include "../src/rules.hxx"

using namespace desres::viparr;

namespace {

    dict get_func_registry() {
        dict D;
        for (std::map<std::string, Rules::VDWFunc>::const_iterator iter = Rules::VDWFuncRegistry().begin(); iter != Rules::VDWFuncRegistry().end(); ++iter)
            D[iter->first] = iter->second;
        return D;
    }
    void set_func_registry(const object& obj) {
        dict D(obj);
        list L(D.keys());
        Rules::VDWFuncRegistry().clear();
        for (unsigned i = 0; i < len(L); ++i) {
            std::string key = extract<std::string>(L[i]);
            Rules::VDWFunc func = extract<Rules::VDWFunc>(D[L[i]]);
            Rules::VDWFuncRegistry().insert(make_pair(key, func));
        }
    }

    dict get_rule_registry() {
        dict D;
        for (std::map<std::string, Rules::VDWCombRulePtr>::const_iterator iter = Rules::VDWCombRuleRegistry().begin(); iter != Rules::VDWCombRuleRegistry().end(); ++iter)
            D[iter->first] = iter->second;
        return D;
    }
    void set_rule_registry(const object& obj) {
        dict D(obj);
        list L(D.keys());
        Rules::VDWCombRuleRegistry().clear();
        for (unsigned i = 0; i < len(L); ++i) {
            std::string key = extract<std::string>(L[i]);
            Rules::VDWCombRulePtr rule = extract<Rules::VDWCombRulePtr>(D[L[i]]);
            Rules::VDWCombRuleRegistry().insert(make_pair(key, rule));
        }
    }

    list get_info(const Rules& rules) {
        list L;
        const std::vector<std::string>& info = rules.info;
        for (unsigned i = 0; i < info.size(); ++i)
            L.append(info[i]);
        return L;
    }
    void set_info(Rules& rules, const object& obj) {
        list L(obj);
        rules.info.resize(len(L));
        for (unsigned i = 0; i < rules.info.size(); ++i)
            rules.info[i] = extract<std::string>(L[i]);
    }

    list get_plugins(const Rules& rules) {
        list L;
        const std::vector<std::string>& plugins = rules.plugins;
        for (unsigned i = 0; i < plugins.size(); ++i)
            L.append(plugins[i]);
        return L;
    }
    void set_plugins(Rules& rules, const object& obj) {
        list L(obj);
        rules.plugins.resize(len(L));
        for (unsigned i = 0; i < rules.plugins.size(); ++i)
            rules.plugins[i] = extract<std::string>(L[i]);
    }

    void set_exclusions(Rules& rules, unsigned exclusions, const object& es, const object& lj) {
        list Les(es);
        std::vector<double> es_scale(len(Les));
        for (unsigned i = 0; i < es_scale.size(); ++i)
            es_scale[i] = extract<double>(Les[i]);
        list Llj(lj);
        std::vector<double> lj_scale(len(Llj));
        for (unsigned i = 0; i < lj_scale.size(); ++i)
            lj_scale[i] = extract<double>(Llj[i]);
        rules.setExclusions(exclusions, es_scale, lj_scale);
    }

    list get_param_names(const Rules::VDWFunc& func) {
        list L;
        const std::vector<std::string>& param_names = func.param_names;
        for (unsigned i = 0; i < param_names.size(); ++i)
            L.append(param_names[i]);
        return L;
    }
    void set_param_names(Rules::VDWFunc& func, const object& obj) {
        list L(obj);
        func.param_names.resize(len(L));
        for (unsigned i = 0; i < func.param_names.size(); ++i)
            func.param_names[i] = extract<std::string>(L[i]);
    }

    list get_pair_names(const Rules::VDWFunc& func) {
        list L;
        const std::vector<std::string>& pair_names = func.pair_param_names;
        for (unsigned i = 0; i < pair_names.size(); ++i)
            L.append(pair_names[i]);
        return L;
    }
    void set_pair_names(Rules::VDWFunc& func, const object& obj) {
        list L(obj);
        func.pair_param_names.resize(len(L));
        for (unsigned i = 0; i < func.pair_param_names.size(); ++i)
            func.pair_param_names[i] = extract<std::string>(L[i]);
    }

    list get_supported_rules(const Rules::VDWFunc& func) {
        list L;
        const std::vector<std::string>& rules = func.supported_rules;
        for (unsigned i = 0; i < rules.size(); ++i)
            L.append(rules[i]);
        return L;
    }
    void set_supported_rules(Rules::VDWFunc& func, const object& obj) {
        list L(obj);
        func.supported_rules.resize(len(L));
        for (unsigned i = 0; i < func.supported_rules.size(); ++i)
            func.supported_rules[i] = extract<std::string>(L[i]);
    }

    list call_rule(const Rules::VDWCombRule& rule, const object& obji, const object& objj, double lfact) {
        list Li(obji);
        std::vector<double> vi(len(Li));
        for (unsigned i = 0; i < vi.size(); ++i)
            vi[i] = extract<double>(Li[i]);
        list Lj(objj);
        std::vector<double> vj(len(Lj));
        for (unsigned i = 0; i < vj.size(); ++i)
            vj[i] = extract<double>(Lj[i]);
        std::vector<double> output = rule(vi, vj, lfact);
        list Loutput;
        for (unsigned i = 0; i < output.size(); ++i)
            Loutput.append(output[i]);
        return Loutput;
    }

    class VDWCombRulePy : public Rules::VDWCombRule {
        public:
            static Rules::VDWCombRulePtr create(const object& py_func) {
                return Rules::VDWCombRulePtr(new VDWCombRulePy(py_func));
            }
            virtual std::vector<double> operator()(const std::vector<double>& vi, const std::vector<double>& vj, double lfact) const {
                list Li;
                for (unsigned i = 0; i < vi.size(); ++i)
                    Li.append(vi[i]);
                list Lj;
                for (unsigned j = 0; j < vj.size(); ++j)
                    Lj.append(vj[j]);
                try {
                    list Loutput(_py_func(Li, Lj, lfact));
                    std::vector<double> output(len(Loutput));
                    for (unsigned i = 0; i < output.size(); ++i)
                        output[i] = extract<double>(Loutput[i]);
                    return output;
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error calling Python VDWCombRule: "
                            + std::string(e.what()));
                }
            }
            virtual bool operator==(const VDWCombRule& other) const {
                try {
                    const VDWCombRulePy& other_py
                        = dynamic_cast<const VDWCombRulePy&>(other);
                    return (_py_func == other_py._py_func);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
        private:
            VDWCombRulePy(const object& py_func) : _py_func(py_func) { }
            object _py_func;
    };

    void clear_vdw_func_registry() {
        Rules::VDWFuncRegistry().clear();
    }
    void clear_vdw_comb_rule_registry() {
        Rules::VDWCombRuleRegistry().clear();
    }
}

namespace desres { namespace viparr {

    void export_rules() {

        scope rules_class(class_<Rules, RulesPtr>("Rules", no_init)
            .def("__eq__", _eq<RulesPtr>)
            .def("__ne__", _ne<RulesPtr>)
            .def("__hash__", _hash<RulesPtr>)
            .def("__init__", make_constructor(&Rules::create))
            .add_static_property("VDWFuncRegistry", get_func_registry, set_func_registry)
            .add_static_property("VDWCombRuleRegistry", get_rule_registry, set_rule_registry)
            .def("MergeVDW", &Rules::MergeVDW)
            .staticmethod("MergeVDW")

            .def("ClearVdwFuncRegistry", clear_vdw_func_registry)
            .staticmethod("ClearVdwFuncRegistry")

            .def("ClearVdwCombRuleRegistry", clear_vdw_comb_rule_registry)
            .staticmethod("ClearVdwCombRuleRegistry")

            .add_property("info", get_info, set_info)
            .def_readwrite("vdw_func", &Rules::vdw_func)
            .def_readwrite("vdw_comb_rule", &Rules::vdw_comb_rule)
            .add_property("plugins", get_plugins, set_plugins)
            .def_readwrite("fatal", &Rules::fatal)
            .def_readwrite("nbfix_identifier", &Rules::nbfix_identifier)
            .def("setExclusions", set_exclusions)
            .def_readonly("exclusions", &Rules::exclusions)
            .def("es_scale", &Rules::es_scale)
            .def("lj_scale", &Rules::lj_scale)
        );

        class_<Rules::VDWFunc>("VDWFunc")
            .def_readwrite("vdw_table_name", &Rules::VDWFunc::vdw_table_name)
            .add_property("param_names", get_param_names, set_param_names)
            .def_readwrite("pair_table_name", &Rules::VDWFunc::pair_table_name)
            .add_property("pair_param_names", get_pair_names, set_pair_names)
            .add_property("supported_rules", get_supported_rules, set_supported_rules)
            .def(self == self)
            ;

        class_<Rules::VDWCombRule, Rules::VDWCombRulePtr, boost::noncopyable>("VDWCombRule", no_init)
            .def("__init__", make_constructor(&VDWCombRulePy::create))
            .def(self == self)
            .def("__call__", call_rule)
            ;
    }

}}
