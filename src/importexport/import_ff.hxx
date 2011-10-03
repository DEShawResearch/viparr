#ifndef desres_viparr_import_ff_hxx
#define desres_viparr_import_ff_hxx

#include "../ff.hxx"

namespace desres { namespace viparr {

    /* Import entire forcefield from forcefield directory. Rules is imported
     * from "rules" and templates from any files of name "templates...". If
     * "atomtypes.def" exists, then a SmartsTree is imported from
     * "atomtypes.def" and (optional) "extendtypes.def", SmartsExclusions are
     * imported from (optional) "exclusions.def", SmartsImpropers are imported
     * from (optional) "impropers.def", and "virtuals.def" and "drudes.def"
     * param files (optional) are imported, and the Forcefield is
     * constructed with a SmartsTyper. Otherwise, the Forcefield is constructed
     * with a TemplateTyper. "cmap" is treated as a cmap param file, and all
     * other files are treated as regular param files. Each param file is
     * imported into the static Forcefield::SharedParamTables map, with the list
     * of rows to the corresponding table stored with the Forcefield. If
     * share_params is true, a check is performed to see if an identical row
     * exists before importing, and if so, a new row is not created.
     * The rowID map key is always the file name. The SharedParamTables key is
     * the file name. */
    ForcefieldPtr ImportForcefield(const std::string& dir,
            bool require_rules=true,
            const std::map<std::string, std::list<msys::Id> >&
            share_params=std::map<std::string, std::list<msys::Id> >());

    /* Format:
     * { "info" : [string, ... , string]
     *   "vdw_func" : string
     *   "vdw_comb_rule" : string
     *   "exclusions" : int
     *   "es_scale" : [float, ... , float]
     *   "lj_scale" : [float, ... , float]
     *   "plugins" : [string, ... , string]
     *   "extendBondTypeEnvironment" : bool
     * }
     *
     * All fields are optional. exclusions n corresponds to exclusions up to
     * 1-n separation, default value set to 4. es_scale[0] corresponds to 1-2
     * scaling, es_scale[1] to 1-3 scaling, etc. Likewise for lj_scale. If
     * es_scale and lj_scale are both unspecified, default values for all scale
     * factors are set to 0. Otherwise, exclusions must match maximum separation
     * distance defined for both scales. */
    RulesPtr ImportRules(const std::string& path);

    /* Format:
     * { template_name: {
     *   "atoms" : [
     *     [atom_name, atomic_number, charge, [btype, nbtype]],
     *                        OR
     *     [atom_name, atomic_number, charge, [type]],
     *     ...
     *     ],
     *   "bonds" : [
     *     [atom_name, atom_name],
     *     ...
     *     ],
     *   "exclusions" : [
     *     [atom_name, atom_name],
     *     ...
     *     ],
     *   "impropers" : [
     *     [atom_name, atom_name, atom_name, atom_name],
     *     ...
     *     ],
     *   "cmap" : [
     *     [atom_name, ... , atom_name]
     *     ...
     *     ],
     *   "pseudos" : [
     *     [pseudo_name, charge, [btype, nbtype], field, site0_name,
     *       ... , siten_name, pset]
     *                        OR
     *     [pseudo_name, charge, [type], field, site0_name, ... , siten_name,
     *       pset],
     *     ...
     *     ]
     *   }
     *   ...
     * }
     * Bonds, impropers, etc. can reference external atom names not in "atoms",
     * but the reference names must be consistent. Atom types cannot contain
     * whitespace characters. Pseudo field "virtual_..." corresponds to param
     * table "virtuals_...".
     * */
    std::vector<TemplatedSystemPtr> ImportTemplates(const std::string& path);

    /* Format:
     * [ [ [phi, psi, energy],
     *     ...
     *     [phi, psi, energy] ],
     *   ...
     * ]
     * */
    std::vector<msys::ParamTablePtr> ImportCmap(const std::string& path);

    /* Import a parameter file. If table_name exists in Forcefield::ParamTables,
     * append new params to this table. Otherwise, create a new param table and
     * add it to Forcefield::ParamTables with key table_name. If share_params
     * is true, check if rows already exist in the shared tables, and if so,
     * do not create a new row. Adds columns named "type" and "memo". "type"
     * column is a ' '-concatenation of all of the "type" fields. If table is
     * "vdw1" or "vdw2", adds column named "nbfix_identifier" and sets value
     * for all imported rows to the given value. A list of added rows is
     * returned.
     *
     * Format:
     * [
     *   { "type": [string, ... , string],
     *     "params": {string: float, ... , string: float},
     *     "memo": string
     *   },
     *   ...
     * ]
     *
     * "type" fields cannot have whitespace characters. "params" keys cannot be
     * called "type", "memo", or "nbfix_identifier". A "params" key of "cmap"
     * is reserved for cmap ID; field of value 1 in "cmap" column is stored as a
     * string "cmap1", etc.
     * */
    std::list<msys::Id> ImportParams(const std::string& table_name,
            const std::string& path, const std::list<msys::Id>&
            share_params=std::list<msys::Id>(),
            const std::string& nbfix_identifier="");

}}

#endif
