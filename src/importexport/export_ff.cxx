#include "../base.hxx"
#include "export_ff.hxx"
#include <boost/filesystem.hpp>
#include <msys/fastjson/fastjson.hxx>

namespace dfj = desres::msys::fastjson;
namespace bfs = boost::filesystem;

namespace desres { namespace viparr {

    void ExportForcefield(ForcefieldPtr ff, const std::string& dir) {
        bfs::create_directories(dir);
        if (!bfs::is_empty(dir))
            VIPARR_FAIL("Cannot export forcefield to nonempty directory "
                    + dir);

        /* Export rules */
        bfs::path rules_path = bfs::path(dir) / "rules";
        try {
            ExportRules(ff->rules(), rules_path.string());
        } catch(std::exception& e) {
            VIPARR_FAIL("Error writing " + rules_path.string()
                    + ": " + e.what());
        }
        for (msys::Id param : ff->rowIDs("vdw1")) {
            if (Forcefield::ParamTable("vdw1")->value(param,
                        "nbfix_identifier").asString()
                    != ff->rules()->nbfix_identifier)
                VIPARR_FAIL("Cannot export forcefield: nbfix_identifier values"
                        " in vdw1 table do not match value in rules");
        }
        for (msys::Id param : ff->rowIDs("vdw2")) {
            if (Forcefield::ParamTable("vdw2")->value(param,
                        "nbfix_identifier").asString()
                    != ff->rules()->nbfix_identifier)
                VIPARR_FAIL("Cannot export forcefield: nbfix_identifier values"
                        " in vdw2 table do not match value in rules");
        }

            /* TemplateTyper -- export templates */
            bfs::path tpls_path = bfs::path(dir) / "templates";
            try {
                ExportTemplates(ff->typer()->templates(), tpls_path.string());
            } catch(std::exception& e) {
                VIPARR_FAIL("Error writing " + tpls_path.string()
                        + ": " + e.what());
            }

        /* Export cmaps */
        if (ff->cmapTables().size() > 0) {
            bfs::path cmaps_path = bfs::path(dir) / "cmap";
            try {
                ExportCmap(ff->cmapTables(), cmaps_path.string());
            } catch(std::exception& e) {
                VIPARR_FAIL("Error writing " + cmaps_path.string()
                        + ": " + e.what());
            }
        }

        /* Export param tables */
        std::vector<std::string> tables = ff->paramTables();
        for (unsigned i = 0; i < tables.size(); ++i) {
            const std::list<msys::Id>& rows = ff->rowIDs(tables[i]);
            std::string name = tables[i];
            bfs::path param_path = bfs::path(dir) / tables[i];
            try {
                ExportParams(Forcefield::ParamTable(name), rows,
                        param_path.string());
            } catch(std::exception& e) {
                VIPARR_FAIL("Error writing " + param_path.string()
                        + ": " + e.what());
            }
        }
    }

}}
