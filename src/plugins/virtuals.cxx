#include "../ff.hxx"
#include "../templated_system.hxx"

using namespace desres;
using namespace desres::viparr;

/* This plugin is an amalgamation of the virtuals_regular and virtuals_shift,
 * plugins. Its purpose is to support current in-use forcefields. */
static void match_virtuals(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    Forcefield::PluginRegistry()["virtuals_regular"]->match(sys, ff);
}

static void compile_virtuals(msys::SystemPtr sys) {
}
    

static Forcefield::RegisterPlugin _("virtuals", match_virtuals,
                                    std::vector<std::string>(),
                                    compile_virtuals);
