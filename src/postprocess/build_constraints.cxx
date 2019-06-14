#include "build_constraints.hxx"
#include "optimize_vsitedefs.hxx"
#include "../add_system_tables.hxx"
#include "../base.hxx"
#include "../ff.hxx"
#include <sstream>
#include <vector>
#include <map>

using namespace desres;
using namespace desres::viparr;
using desres::msys::Id;
using desres::msys::IdList;

namespace {

    /* Fast IdSet with known max ID */
    struct IdSet {
        std::vector<bool> contains;
        unsigned size;

        IdSet(unsigned n) {
            contains = std::vector<bool>(n, false);
            size = 0;
        }
        bool insert(unsigned i) {
            if (contains[i]) return false;
            contains[i] = true;
            ++size;
            return true;
        }
        bool insert(msys::IdList::const_iterator start,
                msys::IdList::const_iterator end) {
            bool success = true;
            for (msys::IdList::const_iterator iter = start; iter != end; ++iter)
                success &= insert(*iter);
            return success;
        }
    };

    /* alist of the form (O, H1, H2) with H1's ID < H2's ID */
    void build_hoh_constraint(msys::SystemPtr sys, const IdList& alist,
            IdSet& constrained_bonds, IdSet& constrained_angles) {

        /* Add tables if they do not exist */
        if (!Forcefield::HasParamTable("constraint_hoh")) {
            msys::ParamTablePtr params = msys::ParamTable::create();
            params->addProp("theta", msys::FloatType);
            params->addProp("r1", msys::FloatType);
            params->addProp("r2", msys::FloatType);
            Forcefield::AddParamTable("constraint_hoh", params);
        }
        msys::ParamTablePtr params = Forcefield::ParamTable("constraint_hoh");
        msys::TermTablePtr constrained = sys->table("constraint_hoh");
        if (constrained == msys::TermTablePtr()) {
            constrained = sys->addTable("constraint_hoh", 3, params);
            constrained->category = msys::CONSTRAINT;
        }
        if (constrained->params() != params)
            VIPARR_FAIL("VIPARR bug: constraint term table does not point to "
                    "global param table");

        /* Create tuple of angle harm and stretch harm params */
        msys::TermTablePtr angle_harm = sys->table("angle_harm");
        msys::TermTablePtr stretch_harm = sys->table("stretch_harm");
        if (angle_harm == msys::TermTablePtr()
                || stretch_harm == msys::TermTablePtr())
            VIPARR_FAIL("Must have stretch harm and angle harm tables before "
                    "adding hoh constraints");
        std::vector<double> param_values;
        IdList angles = angle_harm->findWithAll(alist);
        if (angles.size() == 0) {
            std::stringstream msg;
            msg << "Cannot build constraint: No angle_harm term found for "
                << "water (" << alist[1] << ", " << alist[0] << ", " << alist[2]
                << ")";
            VIPARR_FAIL(msg.str());
        }
        if (angles.size() > 1) {
            std::stringstream msg;
            msg << "Cannot build constraint: Multiple angle_harm terms found "
                << "for water (" << alist[1] << ", " << alist[0] << ", "
                << alist[2] << ")";
            VIPARR_FAIL(msg.str());
        }
        constrained_angles.insert(angles[0]);
        Id param_id = angle_harm->param(angles[0]);
        param_values.push_back(angle_harm->params()->value(param_id,
                    "theta0").asFloat());
        IdList bond(2);
        bond[0] = alist[0];
        for (int i = 1; i <= 2; ++i) {
            bond[1] = alist[i];
            IdList bonds = stretch_harm->findWithAll(bond);
            if (bonds.size() == 0) {
                std::stringstream msg;
                msg << "Cannot build constraint: No stretch_harm term found "
                    << "for water bond (" << bond[0] << ", " << bond[1] << ")";
                VIPARR_FAIL(msg.str());
            }
            if (bonds.size() > 1) {
                std::stringstream msg;
                msg << "Cannot build constraint: Multiple stretch_harm terms "
                    << "found for water bond (" << bond[0] << ", " << bond[1]
                    << ")";
                VIPARR_FAIL(msg.str());
            }
            constrained_bonds.insert(bonds[0]);
            Id param_id = stretch_harm->param(bonds[0]);
            param_values.push_back(stretch_harm->params()->value(
                        param_id, "r0").asFloat());
        }

        /* If param value exists, use existing param. Otherwise, add new
         * param */
        Id param = msys::BadId;
        IdList matches = params->findFloat(params->propIndex("theta"),
                param_values[0]);
        for (unsigned i = 0; i < matches.size(); ++i) {
            if (params->value(matches[i], "r1") == param_values[1]
                    && params->value(matches[i], "r2") == param_values[2]) {
                param = matches[i];
                break;
            }
        }
        if (param == msys::BadId) {
            param = params->addParam();
            params->value(param, "theta") = param_values[0];
            params->value(param, "r1") = param_values[1];
            params->value(param, "r2") = param_values[2];
        }
        constrained->addTerm(alist, param);
    }

    /* alist of the form (A, H1, H2, ..., Hn) */
    void build_ahn_constraint(msys::SystemPtr sys, const IdList& alist,
            IdSet& constrained_bonds) {

        unsigned n = alist.size() - 1;
        std::stringstream name;
        name << "constraint_ah" << n;

        /* Add tables if they do not exist */
        if (!Forcefield::HasParamTable(name.str())) {
            msys::ParamTablePtr params = msys::ParamTable::create();
            for (unsigned i = 1; i <= n; ++i) {
                std::stringstream prop_name;
                prop_name << "r" << i;
                params->addProp(prop_name.str(), msys::FloatType);
            }
            Forcefield::AddParamTable(name.str(), params);
        }
        msys::ParamTablePtr params = Forcefield::ParamTable(name.str());
        msys::TermTablePtr constrained = sys->table(name.str());
        if (constrained == msys::TermTablePtr()) {
            constrained = sys->addTable(name.str(), n+1, params);
            constrained->category = msys::CONSTRAINT;
        }
        if (constrained->params() != params)
            VIPARR_FAIL("VIPARR bug: constraint term table does not point to "
                    "global param table");

        msys::TermTablePtr stretch_harm = sys->table("stretch_harm");
        if (stretch_harm == msys::TermTablePtr())
            VIPARR_FAIL("Must have stretch harm table before adding ahn "
                    "constraints");
        std::vector<double> param_values;
        IdList bond(2);
        bond[0] = alist[0];
        for (unsigned i = 1; i <= n; ++i) {
            bond[1] = alist[i];
            IdList bonds = stretch_harm->findWithAll(bond);
            if (bonds.size() == 0) {
                std::stringstream msg;
                msg << "Cannot build constraint: No stretch_harm term found "
                    << "for bond (" << bond[0] << ", " << bond[1] << ")";
                VIPARR_FAIL(msg.str());
            }
            if (bonds.size() > 1) {
                std::stringstream msg;
                msg << "Cannot build constraint: Multiple stretch_harm terms "
                    << "found for bond (" << bond[0] << ", " << bond[1] << ")";
                VIPARR_FAIL(msg.str());
            }
            constrained_bonds.insert(bonds[0]);
            Id param_id = stretch_harm->param(bonds[0]);
            param_values.push_back(stretch_harm->params()->value(
                        param_id, "r0").asFloat());
        }

        /* If param value exists, use existing param. Otherwise, add new
         * param */
        Id param = msys::BadId;
        IdList matches = params->findFloat(params->propIndex("r1"),
                param_values[0]);
        for (unsigned i = 0; i < matches.size(); ++i) {
            bool match = true;
            for (unsigned j = 1; j <= n; ++j) {
                std::stringstream prop_name;
                prop_name << "r" << j;
                if (params->value(matches[i], prop_name.str())
                        != param_values[j-1]) {
                    match = false;
                    break;
                }
            }
            if (match) {
                param = matches[i];
                break;
            }
        }
        if (param == msys::BadId) {
            param = params->addParam();
            for (unsigned j = 1; j <= n; ++j) {
                std::stringstream prop_name;
                prop_name << "r" << j;
                params->value(param, prop_name.str()) = param_values[j-1];
            }
        }
        constrained->addTerm(alist, param);
    }
}
 
namespace desres { namespace viparr {

    void BuildConstraints(msys::SystemPtr sys, const msys::IdList& atoms,
            bool keep, const std::set<std::string>& exclude, 
            bool optimize_vsite_defs, bool verbose) {

        AddSystemTables(sys);

        /* Remove terms with selected atoms from all constraint tables */
        std::vector<std::string> tables = sys->tableNames();
        for (unsigned i = 0; i < tables.size(); ++i) {
            msys::TermTablePtr table = sys->table(tables[i]);
            if (table->category == msys::CONSTRAINT) {
                msys::IdList terms = table->findWithAny(atoms);
                for (unsigned j = 0; j < terms.size(); ++j)
                    table->delTerm(terms[j]);
            }
        }

        /* Create set of existing constrained atoms, to ensure constraints do
         * not overlap */
        IdSet constrained_atoms(sys->maxAtomId());
        tables = sys->tableNames();
        for (unsigned i = 0; i < tables.size(); ++i) {
            msys::TermTablePtr table = sys->table(tables[i]);
            if (table->category != msys::CONSTRAINT)
                continue;
            IdList terms = table->terms();
            for (unsigned j = 0; j < terms.size(); ++j) {
                IdList term_atoms = table->atoms(terms[j]);
                if (!constrained_atoms.insert(term_atoms.begin(),
                            term_atoms.end()))
                    VIPARR_FAIL("Existing constraints in system overlap");
            }
        }
        
        /* We store constrained bonds and angles and update their
         * "constrained" values later */
        IdSet constrained_bonds(sys->table("stretch_harm")
                == msys::TermTablePtr() ? 0
                : sys->table("stretch_harm")->maxTermId());
        IdSet constrained_angles(sys->table("angle_harm")
                == msys::TermTablePtr() ? 0
                : sys->table("angle_harm")->maxTermId());

        /* Store boolean array of input atom IDs, to reset constrained values
         * of stretch and angle terms for these atoms later. This seems to be
         * much faster than using msys::TermTable::findWithAny. */
        IdSet atomsel(sys->maxAtomId());

        IdSet processed(sys->maxAtomId());
        for (unsigned i = 0; i < atoms.size(); ++i) {
            Id a = atoms[i];
            atomsel.insert(a);
            if (sys->atom(a).atomic_number == 0) continue;
            if (sys->atom(a).atomic_number == 1) {
                /* Switch a with bonded heavy atom */
                Id b = msys::BadId;
                IdList bonded = sys->bondedAtoms(a);
                for (unsigned j = 0; j < bonded.size(); ++j) {
                    if (sys->atom(bonded[j]).atomic_number > 1) {
                        b = bonded[j];
                        break;
                    }
                }
                if (b == msys::BadId) continue;
                a = b;
            }
            if (processed.contains[a]) continue;
            processed.insert(a);
            IdList alist(1, a);
            IdList bonded = sys->bondedAtoms(a);
            unsigned n_real_bonded=0;
            for (unsigned j = 0; j < bonded.size(); ++j){
                if (sys->atom(bonded[j]).atomic_number == 1)
                    alist.push_back(bonded[j]);
                if (sys->atom(bonded[j]).atomic_number > 0)
                    n_real_bonded++;
            }
            unsigned n = alist.size() - 1;
            if (n == 0) // No bonded H
                continue;

            if (sys->atom(a).atomic_number == 8 && n == 2 && n_real_bonded == 2 ) {
                /* HOH constraint */
                if (exclude.find("hoh") == exclude.end()) {
                    /* Ensure constraints do not overlap */
                    if (!constrained_atoms.insert(alist.begin(),
                                alist.end())) {
                        std::stringstream msg;
                        msg << "Constraint with heavy atom " << a
                            << " would overlap other constraints";
                        VIPARR_FAIL(msg.str());
                    }
                    /* Add constraint */
                    if (alist[1] > alist[2]) {
                        Id tmp = alist[1]; alist[1] = alist[2]; alist[2] = tmp;
                    }
                    build_hoh_constraint(sys, alist, constrained_bonds,
                            constrained_angles);
                }
            } else {
                /* AHn constraint */
                std::stringstream name;
                name << "ah" << n;
                if (exclude.find(name.str()) == exclude.end()) {
                    /* Ensure constraints do not overlap */
                    if (!constrained_atoms.insert(alist.begin(),
                                alist.end())) {
                        std::stringstream msg;
                        msg << "Constraint with heavy atom " << a
                            << " would overlap other constraints";
                        VIPARR_FAIL(msg.str());
                    }
                    /* Add constraint */
                    build_ahn_constraint(sys, alist, constrained_bonds);
                }
            }
        }

        /* Update "constrained" values for stretch_harm and angle_harm 
         * tables */
        if (constrained_bonds.size > 0) {
            msys::TermTablePtr stretch_harm = sys->table("stretch_harm");
            IdList terms = stretch_harm->terms();
            for (unsigned i = 0; i < terms.size(); ++i) {
                if (atomsel.contains[stretch_harm->atoms(terms[i])[0]]
                        || atomsel.contains[stretch_harm->atoms(terms[i])[1]]) {
                    stretch_harm->termPropValue(terms[i],"constrained") = 0;
                    if (!keep && constrained_bonds.contains[terms[i]])
                        stretch_harm->termPropValue(terms[i],"constrained") = 1;
                }
            }
        }
        if (constrained_angles.size > 0) {
            msys::TermTablePtr angle_harm = sys->table("angle_harm");
            IdList terms = angle_harm->terms();
            for (unsigned i = 0; i < terms.size(); ++i) {
                if (atomsel.contains[angle_harm->atoms(terms[i])[0]]
                        || atomsel.contains[angle_harm->atoms(terms[i])[1]]
                        || atomsel.contains[angle_harm->atoms(terms[i])[2]]) {
                    angle_harm->termPropValue(terms[i], "constrained") = 0;
                    if (!keep && constrained_angles.contains[terms[i]])
                        angle_harm->termPropValue(terms[i], "constrained") = 1;
                }
            }
        }
       
        /* We have added constraints.
         * See if we can use more efficient virtual site routines 
         */
        if (optimize_vsite_defs)
            OptimizeVsiteDefs(sys, verbose);
    }
}}
