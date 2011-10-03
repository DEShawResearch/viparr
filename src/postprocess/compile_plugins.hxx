#ifndef desres_viparr_compile_plugins
#define desres_viparr_compile_plugins

#include <msys/system.hxx>
#include <msys/term_table.hxx>

namespace desres { namespace viparr {

    void CompilePlugins(msys::SystemPtr sys, std::set<std::string> plugins);

    msys::TermTablePtr AddPairsTable(msys::SystemPtr sys);
    void ApplyNBFix(msys::SystemPtr sys);
    /* This removes all temporary tables added by viparr (marked by
     * msys::NO_CATEGORY) as well as tables with no terms */
    void CleanupSystem(msys::SystemPtr sys);

}}

#endif
