#include "import_ff.hxx"
#include <msys/fastjson/fastjson.hxx>
#include <boost/filesystem.hpp>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace dfj = desres::msys::fastjson;
namespace bfs = boost::filesystem;

using namespace desres;
using namespace desres::viparr;

namespace {

    std::vector<std::string> import_plugins(const dfj::Json& plugins) {
        std::vector<std::string> list;
        if (plugins.valid()) {
            if (plugins.kind() != dfj::Json::Array)
                VIPARR_FAIL("'plugins' field must be of array type");
            for (int i = 0; i < plugins.size(); ++i) {
                const dfj::Json& entry = plugins.elem(i);
                const char *pname = NULL;
                if (entry.kind() == dfj::Json::Array) {
                    if (entry.size() != 2
                            || entry.elem(0).kind() != dfj::Json::String
                            || entry.elem(1).kind() != dfj::Json::Int) {
                        std::stringstream ss;
                        ss << "Plugin " << i << " must be of string type";
                        VIPARR_FAIL(ss.str());
                    }
                    pname = entry.elem(0).as_string();
                    if (!strcmp(pname, "mass") && entry.elem(1).as_int()==1)
                        pname = "mass2";
                } else if (entry.kind() == dfj::Json::String)
                    pname = entry.as_string();
                else {
                    std::stringstream ss;
                    ss << "Plugin " << i << " must be of string type";
                    VIPARR_FAIL(ss.str());
                }
                list.push_back(pname);
            }
        }
        return list;
    }

    std::vector<double> import_scale(const dfj::Json& arr) {
        std::vector<double> v;
        if (arr.valid()) {
            try {
                for (int i = 0; i < arr.size(); i++)
                    v.push_back(arr.elem(i).as_float());
            } catch (std::exception& e) {
                VIPARR_FAIL("'es_scale' and 'lj_scale' values must be of "
                        "float type");
            }
        }
        return v;
    }
}

namespace desres { namespace viparr {

    RulesPtr ImportRules(const std::string& path) {

        if (!bfs::exists(path))
            VIPARR_FAIL("File not found");
        dfj::Json js;
        try {
            dfj::parse_json(path.c_str(), js);
        } catch (std::exception& e) {
            VIPARR_FAIL("Misformatted '" + path + "' file: " + e.what());
        }
        if (js.kind() != dfj::Json::Object)
            VIPARR_FAIL("'" + path + "' must be of object type");

        RulesPtr rules = Rules::create();

        const dfj::Json& info = js.get("info");
        if (info.valid())
            for (int i = 0; i < info.size(); ++i)
                rules->info.push_back(info.elem(i).as_string());

        const dfj::Json& exclusions = js.get("exclusions");
        const dfj::Json& es_scale = js.get("es_scale");
        const dfj::Json& lj_scale = js.get("lj_scale");
        if (exclusions.valid() && exclusions.kind() != dfj::Json::Int)
            VIPARR_FAIL("'exclusions' field must be of int type");
        if (es_scale.valid() && es_scale.kind() != dfj::Json::Array)
            VIPARR_FAIL("'es_scale' field must be of array type");
        if (lj_scale.valid() && lj_scale.kind() != dfj::Json::Array)
            VIPARR_FAIL("'lj_scale' field must be of array type");

        /* Determines exclusion rule based on the given value or the number of
         * provided es_scale/lj_scale values. Default set to 4. */
        unsigned excl_rule = exclusions.valid() ? exclusions.as_int() :
            es_scale.valid() ? es_scale.size()+1   :
            lj_scale.valid() ? lj_scale.size()+1   :
            4;

        /* If neither es nor lj scales are defined, set default values to 0.
         * Otherwise, check for compatibility between exclusions and scale
         * sizes. */
        std::vector<double> es = import_scale(es_scale);
        std::vector<double> lj = import_scale(lj_scale);
        if (lj.size() > 0 || es.size() > 0) {
            if (excl_rule != lj.size() + 1 || excl_rule != es.size() + 1)
                VIPARR_FAIL("Incorrect number of es or lj scaling factors for "
                        "exclusion rule");
        } else {
            es = std::vector<double>(excl_rule - 1, 0);
            lj = std::vector<double>(excl_rule - 1, 0);
        }
        rules->setExclusions(excl_rule, es, lj);

        /* Store vdw_func and vdw_comb_rule as lowercase */
        std::string vdw_func;
        std::string vdw_comb_rule;
        try {
            vdw_func = js.get("vdw_func").as_string("");
            vdw_comb_rule = js.get("vdw_comb_rule").as_string("");
        } catch (std::exception& e) {
            VIPARR_FAIL("'vdw_func' and 'vdw_comb_rule' fields must be of "
                    "string type");
        }
        std::transform(vdw_func.begin(), vdw_func.end(), vdw_func.begin(),
                tolower);
        std::transform(vdw_comb_rule.begin(), vdw_comb_rule.end(),
                vdw_comb_rule.begin(), tolower);
        rules->vdw_func = vdw_func;
        rules->vdw_comb_rule = vdw_comb_rule;

        try {
            rules->fatal = js.get("fatal").as_bool(1);
        } catch (std::exception& e) {
            VIPARR_FAIL("'fatal' field must be of boolean type");
        }
        try {
            rules->nbfix_identifier = js.get("nbfix_identifier").as_string("");
        } catch (std::exception& e) {
            VIPARR_FAIL("'NBFixIdentifier' field must be of string type");
        }
        rules->plugins = import_plugins(js.get("plugins"));
        return rules;
    }

}}
