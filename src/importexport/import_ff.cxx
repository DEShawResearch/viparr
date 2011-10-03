#include "import_ff.hxx"
#include <msys/fastjson/parse.hxx>
#include <boost/filesystem.hpp>

namespace dfj = desres::msys::fastjson;
namespace bfs = boost::filesystem;

namespace desres { namespace viparr {

    ForcefieldPtr ImportForcefield(const std::string& dir, bool require_rules,
            const std::map<std::string, std::list<msys::Id> >& share_params) {

        if (!bfs::exists(dir))
            VIPARR_FAIL("Forcefield directory not found");
        if (!bfs::is_directory(dir))
            VIPARR_FAIL("Forcefield directory path is not a directory");

        /* Import rules */
        bfs::path rules_path = bfs::path(dir) / "rules";
        RulesPtr rules;
        if (bfs::exists(rules_path)) {
            try {
                rules = ImportRules(rules_path.string());
            } catch (std::exception& e) {
                VIPARR_FAIL("Error loading " + rules_path.string() + ": "
                        + e.what());
            }
        } else if (require_rules)
            VIPARR_FAIL("Forcefield directory is missing rules file");
        else 
            rules = Rules::create();

        /* Check for atomtypes.def to create appropriate typer */
        TemplateTyperPtr typer;
            /* Create empty TemplateTyper */
            typer = TemplateTyper::create();

        /* Create forcefield */
        ForcefieldPtr ff = Forcefield::create(rules, typer);

        /* Process template, param, and cmap files */
        bfs::directory_iterator iter(dir), end;
        for (; iter != end; ++iter) {
            std::string path = iter->path().string();
            std::string name = iter->path().filename().string();

            /* Ignore temp files, subdirectories, README, rules, .def files */
            if (name.substr(name.size()-1) == "~") continue;
            if (name.substr(0,1) == ".") continue;
            if (bfs::is_directory(path)) continue;
            if (name == "README") continue;
            if (name == "rules") continue;
            if (name.size() >= 4)
              if(name.substr(name.size()-4) == ".def")
                continue;

            /* Import templates */
            if (name.substr(0,9) == "templates") {
                std::vector<TemplatedSystemPtr> templates;
                try {
                    templates = ImportTemplates(path);
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error loading " + path + ": " + e.what());
                }
                for (unsigned i = 0; i < templates.size(); ++i)
                    typer->addTemplate(templates[i]);
                continue;
            }

            /* Import cmap tables */
            if (name == "cmap") {
                std::vector<msys::ParamTablePtr> cmaps;
                try {
                    cmaps = ImportCmap(path);
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error loading " + path + ": " + e.what());
                }
                for (unsigned i = 0; i < cmaps.size(); ++i)
                    ff->addCmapTable(cmaps[i]);
                continue;
            }

            /* Treat remaining files as param tables. */
            std::string nbfix_identifier;
            if (name == "vdw1")
                nbfix_identifier = rules->nbfix_identifier;
            else if (name == "vdw2") {
                if (rules->nbfix_identifier == "")
                    VIPARR_FAIL("Cannot have vdw2 table without "
                            "nbfix_identifier in rules file");
                nbfix_identifier = rules->nbfix_identifier;
            } else
                nbfix_identifier = "";

            std::list<msys::Id> rows;
            std::map<std::string, std::list<msys::Id> >::const_iterator
                share_iter = share_params.find(name);
            try {
                if (share_iter == share_params.end())
                    rows = ImportParams(name, path, std::list<msys::Id>(),
                            nbfix_identifier);
                else
                    rows = ImportParams(name, path, share_iter->second,
                            nbfix_identifier);
                ff->appendParams(name, rows);
            } catch (std::exception& e) {
              //VIPARR_FAIL("Error loading " + path + ": " + e.what());
              VIPARR_ERR << "\n\a!!!WARNING: failed to load " << path << ": " << e.what();
              VIPARR_ERR << "\n\aContinuing anyway, but BEWARE this missing component.\n\n";
            }
        }

        ff->name = dir;
        return ff;
    }
}}
