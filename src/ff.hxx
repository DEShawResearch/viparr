#ifndef desres_viparr_ff_hxx
#define desres_viparr_ff_hxx

#include "base.hxx"
#include "rules.hxx"
#include "template_typer.hxx"
#include <msys/append.hxx>
#include <list>
#include <map>

namespace desres { namespace viparr {

    class Forcefield;
    typedef std::shared_ptr<Forcefield> ForcefieldPtr;

    /* Representation of a forcefield, consisting of
     * (1) a Rules object, containing the nonbonded info and list of plugins
     * (2) a TemplateTyper, SmartsTyper, or ScoredSmartsTyper object, containing
     *       information to assign atom types, extra exclusions, impropers,
     *       cmaps, and pseudo particles to the system
     * (3) a map from param table name to list of row IDs, indicating which
     *       static param tables and which rows in those param tables are from
     *       this forcefield
     * (4) an optional list of cmap param tables
     *
     * In addition, the Forcefield class has a static collection of all of the
     * param tables, with static accessors and modifiers, and a static
     * name-to-plugin registry of supported plugins. */
    class Forcefield {

        private:
            /* Static shared param tables */
            static std::map<std::string, msys::ParamTablePtr> SharedParamTables;
            static std::list<msys::Id> _empty_idlist;

            RulesPtr _rules;
            TemplateTyperPtr _typer;
            std::map<std::string, std::list<msys::Id> > _row_ids_map;
            std::vector<msys::ParamTablePtr> _cmaps;

            /* Private constructor; must construct using create() function */
            Forcefield(RulesPtr rules, TemplateTyperPtr typer) :
                _rules(rules), _typer(typer) { }

        public:
            /* Functions to access and modify static shared param tables */
            static bool HasParamTable(const std::string& name);
            static msys::ParamTablePtr ParamTable(const std::string& name);
            static void AddParamTable(const std::string& name,
                    msys::ParamTablePtr table);
            static std::vector<std::string> AllParamTables();
            /* Clears all shared tables; invalidates all existing Forcefield
             * objects */
            static void ClearParamTables() { SharedParamTables.clear(); }

            /* Functions to search within static shared param tables, by
             * taking a given list of params and returning only the params
             * matching the given key-value pair */
            template <class Container>
            static Container FilterParams(const std::string& table_name,
                    const Container& params, const std::string& key,
                    const std::string& value) {
                msys::ParamTablePtr table = ParamTable(table_name);
                if (table == msys::ParamTablePtr())
                    VIPARR_FAIL("Table '" + table_name + "' not found");
                msys::Id prop_index = table->propIndex(key);
                if (prop_index == msys::BadId)
                    VIPARR_FAIL("Parameter name '" + key + "' not found in "
                            + "table '" + table_name + "'");
                if (table->propType(prop_index) != msys::StringType)
                    VIPARR_FAIL("Column '" + key + "' of table '" + table_name
                            + "' does not have string type");
                Container filtered;
                for (msys::Id param : params) {
                    if (param >= table->paramCount()) {
                        std::stringstream msg;
                        msg << "Parameter " << param << " of table '"
                            << table_name << "' does not exist";
                        VIPARR_FAIL(msg.str());
                    }
                    if (table->value(param, prop_index).asString() == value)
                        filtered.push_back(param);
                }
                return filtered;
            }
            template <class Container>
            static Container FilterParams(const std::string& table_name,
                    const Container& params, const std::string& key,
                    int value) {
                msys::ParamTablePtr table = ParamTable(table_name);
                if (table == msys::ParamTablePtr())
                    VIPARR_FAIL("Table '" + table_name + "' not found");
                msys::Id prop_index = table->propIndex(key);
                if (prop_index == msys::BadId)
                    VIPARR_FAIL("Parameter name '" + key + "' not found in "
                            + "table '" + table_name + "'");
                if (table->propType(prop_index) != msys::IntType)
                    VIPARR_FAIL("Column '" + key + "' of table '" + table_name
                            + "' does not have int type");
                Container filtered;
                for (msys::Id param : params) {
                    if (param >= table->paramCount()) {
                        std::stringstream msg;
                        msg << "Parameter " << param << " of table '"
                            << table_name << "' does not exist";
                        VIPARR_FAIL(msg.str());
                    }
                    if (table->value(param, prop_index).asInt() == value)
                        filtered.push_back(param);
                }
                return filtered;
            }
            template <class Container>
            static Container FilterParams(const std::string& table_name,
                    const Container& params, const std::string& key,
                    double value) {
                msys::ParamTablePtr table = ParamTable(table_name);
                if (table == msys::ParamTablePtr())
                    VIPARR_FAIL("Table '" + table_name + "' not found");
                msys::Id prop_index = table->propIndex(key);
                if (prop_index == msys::BadId)
                    VIPARR_FAIL("Parameter name '" + key + "' not found in "
                            + "table '" + table_name + "'");
                if (table->propType(prop_index) != msys::FloatType)
                    VIPARR_FAIL("Column '" + key + "' of table '" + table_name
                            + "' does not have float type");
                Container filtered;
                for (msys::Id param : params) {
                    if (param >= table->paramCount()) {
                        std::stringstream msg;
                        msg << "Parameter " << param << " of table '"
                            << table_name << "' does not exist";
                        VIPARR_FAIL(msg.str());
                    }
                    if (table->value(param, prop_index).asFloat() == value)
                        filtered.push_back(param);
                }
                return filtered;
            }

            /* Plugin class is a wrapper for either C++ function pointers
             * or Python functions; the functions must be of the form
             *      void f(TemplatedSystemPtr sys, ForcefieldPtr ff)
             * General use-case is a plugin that matches the typed atoms in the
             * system to a particular param table */
            struct Plugin {
                /* Function to call when plugin is first applied */
                virtual void match(TemplatedSystemPtr sys,
                        ForcefieldPtr ff) const = 0;
                /* Function to call after match is called, subsequently
                   every time one of the tables in its dependencies is
                   modified. */
                virtual void compile(msys::SystemPtr sys) const = 0;

                /* The list of names of ParamTables that should trigger
                   recompilation upon modification.
                 */
                std::vector<std::string> dependencies;
                /* A list of Plugin names which, if also applied in this 
                 * forcefield, must be applied before this one */
                std::vector<std::string> prerequisites;

                virtual ~Plugin() { }
            };
            typedef std::shared_ptr<Plugin> PluginPtr;

            /* Static registry of all supported plugins */
            static std::map<std::string, PluginPtr>& PluginRegistry();

            /* Call this constructor in a static initializer to register a
             * plugin on program initialization; registry is populated in
             * plugins/....cxx */
            struct RegisterPlugin {
                explicit RegisterPlugin(const std::string& name,
                        void (*c_match)(TemplatedSystemPtr, ForcefieldPtr),
                        const std::vector<std::string>& prerequisites
                                        =std::vector<std::string>(),
                        void (*c_compile)(msys::SystemPtr)=NULL,
                        const std::vector<std::string>& dependencies=std::vector<std::string>());
            };
            /* Call this constructor in a static initializer to register an
             * additional Plugin-Plugin prerequisite on program
             * initialization */
            struct RegisterPluginPrerequisite {
                explicit RegisterPluginPrerequisite(const std::string&
                        first_plugin, const std::string& second_plugin);
            };

            /* Create forcefield with given rules and typer, with no pattern
             * tables and no cmap tables */
            static ForcefieldPtr create(RulesPtr rules,
                    TemplateTyperPtr typer) {
                return ForcefieldPtr(new Forcefield(rules, typer));
            }

            /* Create a shallow copy of a forcefield sharing the same rules and
             * typer objects */
            static ForcefieldPtr copy(ForcefieldPtr ff) {
                return ForcefieldPtr(new Forcefield(*ff.get()));
            }

            /* Can store the file-path from which the forcefield was imported */
            std::string name;

            RulesPtr rules() { return _rules; }
            void resetRules(RulesPtr rules) { _rules = rules; }
            TemplateTyperPtr typer() { return _typer; }
            void resetTyper(TemplateTyperPtr typer) { _typer = typer; }

            /* Access which rows of a param table are from this forcefield.
             * Returns an empty list if this forcefield does not
             * have the given param table */
            const std::list<msys::Id>& rowIDs(const std::string& name) const;

            /* Delete a parameter, list of parameters, or entire parameter
             * table from this forcefield. */
            void delParam(const std::string& name, msys::Id param);
            void delParams(const std::string& name,
                    const std::list<msys::Id>& params);
            void clearParams(const std::string& name);

            /* Add parameters to this forcefield, or replace an existing
             * parameter */
            void appendParam(const std::string& name, msys::Id param);
            void appendParams(const std::string& name,
                    const std::list<msys::Id>& params);
            void replaceParam(const std::string& name, msys::Id old_param,
                    msys::Id new_param);

#if 0
            /* Mode manipulation, for forcefields with modes */
            msys::Id getModeParam(const std::string& name, msys::Id rowID,
                    const std::string& mode);
            msys::Id createModeParam(const std::string& name, msys::Id rowID,
                    const std::string& mode);
            bool setMode(const std::string& name, msys::Id rowID,
                    const std::string& mode);
            std::map<std::string, msys::IdList> setMode(const
                    std::string& mode);
#endif

            std::vector<std::string> paramTables() const;

            /* Functions to access and modify cmap tables; the cmapTable
             * accessor uses a 1-based index */
            msys::ParamTablePtr cmapTable(unsigned cmap) const;
            const std::vector<msys::ParamTablePtr>& cmapTables() const {
                return _cmaps; }
            void addCmapTable(msys::ParamTablePtr cmap_table) { 
                _cmaps.push_back(cmap_table); }
            void delCmapTables() { _cmaps.clear(); }
    };

}}

#endif
