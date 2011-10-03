#ifndef desres_viparr_merge_forcefields_hxx
#define desres_viparr_merge_forcefields_hxx

#include "ff.hxx"
#include <string>

namespace desres { namespace viparr {

    /* Updates a source forcefield in place using a patch forcefield, merging
     * the rules, plugins, templates, param tables, scored SMARTS, SMARTS
     * exclusions, and SMARTS impropers. Templates in the patch forcefield 
     * overwrite those in the source forcefield of the same name, and param
     * table rows in the patch forcefield (including rows of param tables for 
     * drude and virtual definitions) overwrite those in the source forcefield 
     * of the same type string. Cmap tables and a SMARTS atomtype tree in the 
     * patch forcefield, if present, overwrite all cmap tables and the SMARTS 
     * atomtype tree in the source forcefield. */
    void MergeForcefields(ForcefieldPtr src, ForcefieldPtr patch,
            bool append_only=false, bool verbose=true);

    void MergeRules(RulesPtr src_rules, RulesPtr patch_rules,
            bool verbose=true);
    void MergeTemplates(TemplateTyperPtr src_typer, TemplateTyperPtr
            patch_typer, bool append_only=false, bool verbose=true);
    std::list<msys::Id> MergeParams(const std::list<msys::Id>& src_rows,
            const std::list<msys::Id>& patch_rows, msys::ParamTablePtr table,
            bool append_only=false, bool verbose=true);

}}

#endif
