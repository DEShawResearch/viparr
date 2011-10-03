#include "add_nbody_table.hxx"

void desres::viparr::AddNbodyTable(TemplatedSystemPtr sys, ForcefieldPtr ff,
        const std::string& table_name, const std::string& plugin_name,
        int natoms, const std::vector<msys::IdList>& nbodies,
        SystemToPatternPtr sys_to_pattern, TypeToPatternPtr type_to_pattern,
        const std::vector<PermutationPtr>& perms, bool required,
        msys::Category category) {
    if (ff->rowIDs(table_name).size() == 0)
        VIPARR_FAIL("Must have '" + table_name + "' table for '" +
                plugin_name + "' plugin");
    msys::TermTablePtr table = sys->system()->addTable(table_name, natoms,
            Forcefield::ParamTable(table_name));
    table->category = category;
    ParameterMatcherPtr matcher = ParameterMatcher::create(ff, table_name,
            sys_to_pattern, type_to_pattern, perms);
    unsigned old_size = table->termCount();
    for (unsigned i = 0, n = nbodies.size(); i < n; ++i) {
        const msys::IdList& term = nbodies[i];
        msys::Id row = matcher->match(sys, term);
        if (row == msys::BadId) {
            if (!required) continue;
            Pattern patt = (*sys_to_pattern)(sys, term);
            if (!matcher->match_bonds)
                patt.bonds.clear();
            std::stringstream msg;
            msg << "No match found for table '" << table_name << "', pattern "
                << patt.print() << ", atoms (" << term[0];
            for (int j = 1; j < natoms; ++j)
                msg << "," << term[j];
            msg << ")";
            if (!ff->rules()->fatal) {
                VIPARR_ERR << "WARNING: " << msg.str() << std::endl;
                continue;
            }
            else
                VIPARR_FAIL(msg.str());
        }
        table->addTerm(term, row);
    }
    /* Double-check that if all nbodies are required to be matched, we have
     * indeed added a match for each nbody */
    if (table->termCount() - old_size != nbodies.size()
            && required && ff->rules()->fatal)
        VIPARR_FAIL("VIPARR BUG--incorrect number of terms added");
}
