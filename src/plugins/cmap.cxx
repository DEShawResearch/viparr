#include "add_nbody_table.hxx"

using namespace desres;
using namespace desres::viparr;

/*
  DESRESCode#3431 Note that while cmaps live in a std::vector and
  are therefore 0-indexed, there are referred to in the
  torsion_torsion_cmap table 1-indexed. The conversion happens in
  Forcefield::cmapTable.
*/
static void apply_cmap(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    /* Match forward and reverse permutations with no bonds */
    std::vector<PermutationPtr> perms;
    perms.push_back(Permutation::Identity);
    perms.push_back(Permutation::Reverse);
    AddNbodyTable(sys, ff, "torsiontorsion_cmap", "cmap", 8, sys->cmaps(),
            SystemToPattern::BType, TypeToPattern::Default, perms, true,
            msys::BOND);

    msys::TermTablePtr table = sys->system()->table("torsiontorsion_cmap");
    /* Add cmap param tables as auxiliary tables to system */
    msys::IdList terms = table->terms();
    std::map<std::string, msys::Id> row_map; // old cmapid to new param row
    for (unsigned i = 0, n = terms.size(); i < n; ++i) {
        std::string name = table->propValue(terms[i], "cmapid").asString();
        if (row_map.find(name) == row_map.end()) {
            /* This is the first time we encounter this cmapid for this ff */
            int cmap = atoi(name.substr(4).c_str());
            if (cmap > (int)(ff->cmapTables().size()))
                VIPARR_FAIL("Missing cmap table " + name.substr(4));
            std::string new_name = name;
            if (sys->system()->auxTable(name) != msys::ParamTablePtr()) {
                /* This cmapid clashes with something already in the system;
                 * rename it */
                int new_cmap = 0;
                do {
                    ++new_cmap;
                    std::stringstream ss;
                    ss << "cmap" << new_cmap;
                    new_name = ss.str();
                } while (sys->system()->auxTable(new_name)
                        != msys::ParamTablePtr());
            }
            /* Create a new row in the torsiontorsion table for this system,
             * with possibly a different cmapid from the forcefield */
            msys::Id new_p = Forcefield::ParamTable(
                    "torsiontorsion_cmap")->duplicate(table->param(terms[i]));
            Forcefield::ParamTable("torsiontorsion_cmap")->value(new_p,
                    "cmapid") = new_name;
            row_map[name] = new_p;
            sys->system()->addAuxTable(new_name, ff->cmapTable(cmap));
        }
        table->setParam(terms[i], row_map[name]);
    }
}

static Forcefield::RegisterPlugin _("cmap", apply_cmap);
