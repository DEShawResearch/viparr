#ifndef desres_viparr_rules_hxx
#define desres_viparr_rules_hxx

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <boost/assign.hpp>

namespace desres { namespace viparr {

    /* Representation of rules forcefield file, containing the forcefield's
     * nonbonded info and list of plugins.
     * Also contains a static registry of supported VDW functions and VDW
     * combine rules. */
    class Rules {
            std::vector<double> _es_scale;
            std::vector<double> _lj_scale;

        public:

            /* Class to hold information about a VDW function */
            struct VDWFunc {
                std::string vdw_table_name; /* Name of nonbonded table */
                std::vector<std::string> param_names; /* Nonbonded params */
                std::string pair_table_name; /* Name of nonbonded pairs table */
                std::vector<std::string> pair_param_names; /* Pairs params */
                
                /* Combine rules supported by this vdw function */
                std::vector<std::string> supported_rules;

                bool operator==(const VDWFunc& other) const;
            };

            /* Static registry of all supported VDWFuncs */
            static std::map<std::string, VDWFunc>& VDWFuncRegistry();

            /* Wrapper class for a C++ function pointer or Python function
             * representing a VDW combine rule. The function must take two
             * param value lists of equal size and a scaling factor and return
             * the param list for the pair. */
            struct VDWCombRule {
                virtual std::vector<double> operator()(const
                        std::vector<double>& vi, const std::vector<double>& vj,
                        double lfact) const = 0;
                virtual bool operator==(const VDWCombRule& other) const = 0;

                virtual ~VDWCombRule() { }
            };
            typedef std::shared_ptr<VDWCombRule> VDWCombRulePtr;

            /* Static registry of all supported VDWCombRules */
            static std::map<std::string, VDWCombRulePtr>& VDWCombRuleRegistry();

            /* Call these constructors in a static initializer to register a vdw_func
             * or vdw_comb_rule on program initialization; registry is populated in
             * rules.cxx. */
            struct RegisterVDWFunc {
                explicit RegisterVDWFunc(const std::string& name,
                        const Rules::VDWFunc& func);
            };
            struct RegisterVDWCombRule {
                explicit RegisterVDWCombRule(const std::string& name,
                        std::vector<double> (*c_rule)(const std::vector<double>&,
                            const std::vector<double>&, double));
            };

            static std::shared_ptr<Rules> create() {
                return std::make_shared<Rules>();
            }

            /* Merge vdw_func or vdw_comb_rule names. If one of the values is
             * empty, the merged value is the nonempty value. Otherwise, if they
             * disagree, raises an error. */
            static std::string MergeVDW(const std::string& old_val,
                    const std::string& new_val);

            std::vector<std::string> info;
            std::string vdw_func;
            std::string vdw_comb_rule;
            std::vector<std::string> plugins;

            bool fatal = true;
            std::string nbfix_identifier;

            /* Set the exclusion rule and es and lj scales. es_scale[0]
             * is the scaling factor for 1-2 separation, es_scale[1] for 1-3,
             * etc., and likewise for lj_scale. exclusions must equal the
             * maximum separation for which a scaling factor is defined in both
             * scale. */
            void setExclusions(unsigned exclusions,
                    const std::vector<double>& es_scale,
                    const std::vector<double>& lj_scale);

            /* Return the exclusion rule */
            unsigned exclusions() const;

            /* Return the scale factor for 1-i separation */
            double es_scale(unsigned i) const;
            double lj_scale(unsigned i) const;
    };
    typedef std::shared_ptr<Rules> RulesPtr;

}}

#endif
