#include "base.hxx"
#include "rules.hxx"
#include <math.h>
#include <sstream>

using namespace desres::viparr;

namespace {

    class VDWCombRuleC : public Rules::VDWCombRule {
        public:
            VDWCombRuleC(std::vector<double> (*c_rule) (const
                        std::vector<double>&, const std::vector<double>&,
                        double)) : _c_rule(c_rule) {
                if (c_rule == NULL)
                    VIPARR_FAIL("Cannot create VDWCombRule with NULL pointer");
            }
            virtual bool operator==(const Rules::VDWCombRule& other) const {
                try {
                    const VDWCombRuleC& other_c
                        = dynamic_cast<const VDWCombRuleC&>(other);
                    return (_c_rule == other_c._c_rule);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
            virtual std::vector<double> operator()(const std::vector<double>&
                    vi, const std::vector<double>& vj, double lfact) const {
                return _c_rule(vi, vj, lfact);
            }
        private:
            std::vector<double> (*_c_rule) (const std::vector<double>&,
                    const std::vector<double>&, double);
    };

    std::vector<double> convert_sig_eps(double sij, double eij, double lfact)
    {
        std::vector<double> p(2);
        p[0] = lfact * pow( sij, 12) * eij * 4.0; /* aij */
        p[1] = lfact * pow( sij,  6) * eij * 4.0; /* bij */
        return p;
    }
    std::vector<double> geometric(const std::vector<double>& vi,
            const std::vector<double>& vj, double lfact) {
        double sij = sqrt( vi[0] * vj[0]);
        double eij = sqrt( vi[1] * vj[1]);
        return convert_sig_eps(sij, eij, lfact);
    }
    std::vector<double> arith_geom(const std::vector<double>& vi,
            const std::vector<double>& vj, double lfact) {
        double sij = 0.5 * (vi[0] + vj[0]);
        double eij =  sqrt( vi[1] * vj[1]);
        return convert_sig_eps(sij, eij, lfact);
    }
    std::vector<double> lb_geometric(const std::vector<double>& vi,
            const std::vector<double>& vj, double lfact) {
        double Ai, Bi, Ci, Aj, Bj, Cj;
        /* Convert type i, alpha==0 signifies that parameters are A,B */
        if (vi[0]==0) {
            Ai = vi[0];
            Bi = vi[1];
            Ci = 0;
        } else {
            Ai = (6.*vi[1]*exp(vi[0]))/(vi[0] - 6.);
            Bi = vi[2] / vi[0];
            Ci = (pow(vi[0],7)*vi[1])/(vi[0] - 6.);
        }
        /* Convert type j, alpha==0 signifies that parameters are A,B */
        if (vj[0]==0) {
            Aj = vj[0];
            Bj = vj[1];
            Cj = 0;
        } else {
            Aj = (6.*vj[1]*exp(vj[0]))/(vj[0] - 6.);
            Bj = vj[2] / vj[0];
            Cj = (pow(vj[0],7)*vj[1])/(vj[0] - 6.);
        }
        std::vector<double> p(3);
        p[0] = lfact * sqrt(Ai*Aj); /* aij */
        p[1] =         0.5*(Bi+Bj); /* bij */
        p[2] = lfact * sqrt(Ci*Cj); /* cij */
        return p;
    }
    std::vector<double> none(const std::vector<double>& vi,
            const std::vector<double>& vj, double lfact) {
        VIPARR_FAIL("Cannot call 'none' combine rule");
        return std::vector<double>();
    }
}

namespace desres { namespace viparr {

    bool Rules::VDWFunc::operator==(const Rules::VDWFunc& other) const {
        return (vdw_table_name == other.vdw_table_name
                && param_names == other.param_names
                && pair_table_name == other.pair_table_name
                && pair_param_names == other.pair_param_names
                && supported_rules == other.supported_rules);
    }

    void Rules::setExclusions(unsigned exclusions,
            const std::vector<double>& es_scale,
            const std::vector<double>& lj_scale) {
        if (es_scale.size() != lj_scale.size() ||
                (exclusions == 0 && es_scale.size() != 0) ||
                (exclusions != es_scale.size() + 1)) {
            std::stringstream msg;
            msg << "Incompatible inputs: exclusion rule = " << exclusions
                << ", es_scale size = " << es_scale.size()
                << ", lj_scale size = " << lj_scale.size();
            VIPARR_FAIL(msg.str());
        }
        _es_scale = es_scale;
        _lj_scale = lj_scale;
    }

    unsigned Rules::exclusions() const {
        return _es_scale.size() + 1;
    }

    double Rules::es_scale(unsigned i) const {
        if (i < 2 || i > exclusions()) {
            std::stringstream msg;
            msg << "es_scale not defined for i = " << i
                << ", exclusion rule = " << exclusions();
            VIPARR_FAIL(msg.str());
        }
        return _es_scale[i-2];
    }

    double Rules::lj_scale(unsigned i) const {
        if (i < 2 || i > exclusions()) {
            std::stringstream msg;
            msg << "lj_scale not defined for i = " << i
                << ", exclusion rule = " << exclusions();
            VIPARR_FAIL(msg.str());
        }
        return _lj_scale[i-2];
    }

    std::string Rules::MergeVDW(const std::string& old_val,
            const std::string& new_val) {
        std::string val = old_val;
        if (old_val.size() == 0)
            val = new_val;
        else if (new_val.size() == 0) { }
        else if (old_val != new_val)
            VIPARR_FAIL("Incompatible VDW functions or combine rules: "
                    + old_val + ", " + new_val);
        return val;
    }

    Rules::RegisterVDWFunc::RegisterVDWFunc(const std::string& name,
            const Rules::VDWFunc& func) {
        Rules::VDWFuncRegistry().insert(std::make_pair(name, func));
    }

    Rules::RegisterVDWCombRule::RegisterVDWCombRule(const std::string& name,
            std::vector<double> (*c_rule)(const std::vector<double>&,
                const std::vector<double>&, double)) {
        Rules::VDWCombRuleRegistry().insert(std::make_pair(name,
                    Rules::VDWCombRulePtr(new VDWCombRuleC(c_rule))));
    }

    std::map<std::string, Rules::VDWFunc>& Rules::VDWFuncRegistry() {
        static std::map<std::string, Rules::VDWFunc> registry;
        return registry;
    }

    std::map<std::string, Rules::VDWCombRulePtr>& Rules::VDWCombRuleRegistry() {
        static std::map<std::string, Rules::VDWCombRulePtr> registry;
        return registry;
    }
}}

static Rules::VDWFunc LJ12_6_SIG_EPSILON = {
    "vdw_12_6",
    {"sigma","epsilon"},
    "pair_12_6_es",
    {"aij","bij"},
    {"geometric","arithmetic/geometric"},
};
static Rules::VDWFunc EXP_6X = {
    "vdw_exp_6",
    {"alpha","epsilon","rmin"},
    "pair_exp_6_es",
    {"aij","bij","cij"},
    {"lb/geometric"},
};
static Rules::VDWFunc NONE = {
    "none",
    std::vector<std::string>(),
    "none",
    std::vector<std::string>(),
    {"none"},
};
static Rules::VDWFunc FROM_TABLES = {
    "from_tables",
    std::vector<std::string>(),
    "from_tables",
    std::vector<std::string>(),
    {"from_tables"},
};
static Rules::RegisterVDWFunc _1("lj12_6_sig_epsilon", LJ12_6_SIG_EPSILON);
static Rules::RegisterVDWFunc _2("exp_6x", EXP_6X);
static Rules::RegisterVDWFunc _4("none", NONE);
static Rules::RegisterVDWFunc _11("from_tables", FROM_TABLES);

static Rules::RegisterVDWCombRule _5("geometric", geometric);
static Rules::RegisterVDWCombRule _6("arithmetic/geometric", arith_geom);
static Rules::RegisterVDWCombRule _7("lb/geometric", lb_geometric);
static Rules::RegisterVDWCombRule _9("none", none);
