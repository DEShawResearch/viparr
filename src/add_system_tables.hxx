#ifndef desres_viparr_add_system_tables_hxx
#define desres_viparr_add_system_tables_hxx

#include <msys/system.hxx>
#include <list>

namespace desres { namespace viparr {

    const std::map<std::string, std::string> DefaultTableNameConversions = {
        {"nonbonded", "vdw1"},
        {"virtual_", "virtuals_"},
    };

    /* Merge parameters in a chemical system into the static Forcefield tables,
     * and point term tables in the system to the static Forcefield tables.
     * Parameter tables in the system are merged into the Forcefield tables of
     * the same name, with the exceptions of the replacements specified by
     * name_conversions. */
    void AddSystemTables(msys::SystemPtr sys,
            const std::map<std::string, std::list<msys::Id> >&
            share_params=std::map<std::string, std::list<msys::Id> >(),
            const std::map<std::string, std::string>& name_conversions
            =DefaultTableNameConversions);

}}

#endif
