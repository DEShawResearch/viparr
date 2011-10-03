#include "../ff.hxx"

using namespace desres;
using namespace desres::viparr;

/* This plugin is an amalgamation of the charges_formal, charges_bci,
 * exclusions_and_scaled_pairs, pairs_es_balanced, and
 * pairs_lj_scaled_14 plugins. Its purpose is to support current in-use
 * forcefields. */
static void match_exclusions(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    Forcefield::PluginRegistry()["exclusions_and_scaled_pairs"]->match(sys, ff);
    if (ff->rowIDs("vdw1_14").size() > 0)
        Forcefield::PluginRegistry()["pairs_lj_scaled_14"]->match(sys, ff);
}

static void compile_exclusions(msys::SystemPtr sys) {
    
    Forcefield::PluginRegistry()["charges_formal"]->compile(sys);

    /* exclusions_and_scaled_pairs MUST be called after
     * pairs_es_balanced to match viparr3's (presumably correct)
     * output */
    Forcefield::PluginRegistry()["exclusions_and_scaled_pairs"]->compile(sys);

    Forcefield::PluginRegistry()["pairs_lj_scaled_14"]->compile(sys);
}

static Forcefield::RegisterPlugin _("exclusions", match_exclusions,
                                    std::vector<std::string>(),
                                    compile_exclusions);
