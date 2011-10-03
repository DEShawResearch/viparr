#include "export_ff.hxx"
#include <msys/fastjson/fastjson.hxx>
#include <boost/filesystem.hpp>

namespace dfj = desres::msys::fastjson;
namespace bfs = boost::filesystem;

namespace desres { namespace viparr {

    void ExportRules(RulesPtr rules, const std::string& path) {
        if (bfs::exists(path))
            VIPARR_FAIL("File already exists; cannot overwrite");
        dfj::Json jinfo;
        jinfo.to_array();
        dfj::Json tmp;
        for (unsigned i = 0; i < rules->info.size(); ++i)
            jinfo.append(tmp.to_string(rules->info[i].c_str()));
        dfj::Json jes_scale;
        dfj::Json jlj_scale;
        jes_scale.to_array();
        jlj_scale.to_array();
        for (unsigned i = 2; i <= rules->exclusions(); ++i) {
            jes_scale.append(tmp.to_float(rules->es_scale(i)));
            jlj_scale.append(tmp.to_float(rules->lj_scale(i)));
        }
        dfj::Json jplugins;
        jplugins.to_array();
        for (unsigned i = 0; i < rules->plugins.size(); ++i)
            jplugins.append(tmp.to_string(rules->plugins[i].c_str()));
        dfj::Json js;
        js.to_object();
        js.append("info", jinfo);
        js.append("vdw_func", tmp.to_string(rules->vdw_func.c_str()));
        js.append("vdw_comb_rule", tmp.to_string(rules->vdw_comb_rule.c_str()));
        js.append("exclusions", tmp.to_int(rules->exclusions()));
        js.append("es_scale", jes_scale);
        js.append("lj_scale", jlj_scale);
        if (rules->nbfix_identifier != "") {
            js.append("nbfix_identifier",
                    tmp.to_string(rules->nbfix_identifier.c_str()));
        }
        if (!rules->fatal)
            js.append("fatal", tmp.to_bool(false));
        if (jplugins.size() > 0)
            js.append("plugins", jplugins);
        dfj::print_json(path.c_str(), js);
    }

}}
