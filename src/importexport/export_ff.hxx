#ifndef desres_viparr_export_ff_hxx
#define desres_viparr_export_ff_hxx

#include "../ff.hxx"
#include <string>

namespace desres { namespace viparr {

    /* Export entire forcefield to an empty directory---Rules to "rules"
     * templates to "templates", SmartsTree to "atomtypes.def",
     * SmartsExclusions to "exclusions.def", SmartsImpropers to "impropers.def",
     * virtuals and drudes tables to "virtuals.def" and "drudes.def", cmap
     * tables to "cmap", and remaining param tables to name in the forcefield.
     * Templates are not exported with a SmartsTyper. */
    void ExportForcefield(ForcefieldPtr ff, const std::string& dir);

    /* Export individual forcefield files. "path" must specify a non-existent
     * file; exports will not overwrite existing files. Rules are exported
     * with default values. All templates are exported to a single file.
     * Atom type trees are exported with extended types. Exported param table
     * must have a "type" column, and all preceding columns are assumed to be 
     * param names. See import_ff.hxx for documentation of file formats. */
    void ExportRules(RulesPtr rules, const std::string& path);
    void ExportTemplates(const std::vector<TemplatedSystemPtr>& typer,
            const std::string& path);
    void ExportCmap(const std::vector<msys::ParamTablePtr>& cmap_tables,
            const std::string& path);
    void ExportParams(msys::ParamTablePtr table,
            const std::list<msys::Id>& rows, const std::string& path);

}}

#endif
