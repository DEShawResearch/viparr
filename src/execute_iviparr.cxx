#include "execute_iviparr.hxx"
#include <viparr/version.hxx> /* Auto-generated in SConscript */
#include "base.hxx"
#include "ff.hxx"
#include "parameter_matcher.hxx"
#include "pattern.hxx"
#include "plugins/pairs_helper.hxx"
#include "importexport/import_ff.hxx"
#include "importexport/export_ff.hxx"
#include <msys/atomsel.hxx>
#include <msys/analyze.hxx>
#include <msys/clone.hxx>
#include <msys/dms.hxx>
#include <msys/io.hxx>
#include <msys/elements.hxx>
#include <msys/override.hxx>

using namespace desres;
using desres::msys::Id;
using desres::msys::IdList;
using namespace desres::viparr;

namespace {

// FIXME - unused, should it be?
//Id get_center_atom(msys::SystemPtr sys, const IdList& improper) {
//    for (unsigned i = 0; i < improper.size(); ++i) {
//        bool center = true;
//        for (unsigned j = 0; j < improper.size(); ++j) {
//            if (i == j) continue;
//            if (sys->findBond(improper[i], improper[j]) == msys::BadId) {
//                center = false;
//                break;
//            }
//        }
//        if (center) return improper[i];
//    }
//    VIPARR_FAIL("Improper term does not have correct bond geometry");
//}

bool is_dihedral(msys::SystemPtr sys, const IdList& dihedral) {
    if (dihedral.size() != 4)
        VIPARR_FAIL("Dihedral term does not have 4 atoms");
    if (sys->findBond(dihedral[0], dihedral[1]) != msys::BadId
            && sys->findBond(dihedral[1], dihedral[2]) != msys::BadId
            && sys->findBond(dihedral[2], dihedral[3]) != msys::BadId)
        return true;
    return false;
}

unsigned get_separation(msys::SystemPtr sys, const IdList& atoms) {
    /* If atom is pseudo, replace with bonded non-pseudo atom */
    IdList base_atoms = atoms;
    for (unsigned idx = 0; idx < 2; ++idx) {
        if (sys->atom(base_atoms[idx]).atomic_number == 0) {
            IdList bonded = sys->bondedAtoms(base_atoms[idx]);
            for (unsigned i = 0; i < bonded.size(); ++i) {
                if (sys->atom(bonded[i]).atomic_number > 0) {
                    base_atoms[idx] = bonded[i];
                    break;
                }
            }
        }
    }
    if (base_atoms[0] == base_atoms[1])
        return 1;
    if (sys->findBond(base_atoms[0], base_atoms[1]) != msys::BadId)
        return 2;
    IdList bonded = sys->bondedAtoms(base_atoms[0]);
    for (unsigned i = 0; i < bonded.size(); ++i)
        if (sys->findBond(bonded[i], base_atoms[1]) != msys::BadId)
            return 3;
    IdList second_bonded;
    for (unsigned i = 0; i < bonded.size(); ++i) {
        IdList tmp = sys->bondedAtoms(bonded[i]);
        second_bonded.insert(second_bonded.end(), tmp.begin(), tmp.end());
    }
    for (unsigned i = 0; i < second_bonded.size(); ++i)
        if (sys->findBond(second_bonded[i], base_atoms[1]) != msys::BadId)
            return 4;
    /* Return 5 if >4 */
    return 5;
}

TemplatedSystemPtr create_template(TemplatedSystemPtr tsys, const IdList& atoms,
        RulesPtr rules) {
    msys::SystemPtr sys = tsys->system();
    /* Construct new template */
    msys::SystemPtr tpl_sys = msys::Clone(sys, atoms);
    std::stringstream resname;
    resname << sys->residue(sys->atom(atoms[0]).residue).name
        << sys->residue(sys->atom(atoms[0]).residue).resid;
    tpl_sys->residue(0).name = resname.str();
    TemplatedSystemPtr tpl = TemplatedSystem::create(tpl_sys);
    IdList sys_to_tpl(sys->maxAtomId(), msys::BadId);
    /* Set atom names and types */
    for (unsigned j = 0; j < atoms.size(); ++j) {
        sys_to_tpl[atoms[j]] = j;
        std::stringstream type;
        type << msys::AbbreviationForElement(
                sys->atom(atoms[j]).atomic_number) << atoms[j];
        tpl_sys->atom(j).name = type.str();
        if (sys->atom(atoms[j]).atomic_number != 0) {
            tpl->setTypes(j, type.str(), type.str());
            tsys->setTypes(atoms[j], type.str(), type.str());
        }
        else {
            tpl->setTypes(j, type.str(), type.str(), type.str());
            tsys->setTypes(atoms[j], type.str(), type.str(), type.str());
        }
    }
    /* Add external bonded atoms and bonds */
    int extern_id = 1;
    for (unsigned j = 0; j < atoms.size(); ++j) {
        IdList bonded = sys->bondedAtoms(atoms[j]);
        for (unsigned k = 0; k < bonded.size(); ++k) {
            if (sys_to_tpl[bonded[k]] == msys::BadId) {
                Id external = tpl_sys->addAtom(tpl_sys->atom(0).residue);
                sys_to_tpl[bonded[k]] = external;
                std::stringstream type;
                type << "$" << extern_id;
                ++extern_id;
                tpl_sys->atom(external).name = type.str();
                tpl_sys->atom(external).atomic_number = -1;
                tpl_sys->addBond(sys_to_tpl[atoms[j]], external);
            }
        }
    }
    /* Add impropers; put center atom last */
    if (sys->table("dihedral_trig") != msys::TermTablePtr()) {
        IdList term_list = sys->table("dihedral_trig")->findWithAny(atoms);
        for (unsigned i = 0; i < term_list.size(); ++i) {
            IdList term = sys->table("dihedral_trig")->atoms(term_list[i]);
            if (is_dihedral(sys, term)) continue;
            bool in_template = true;
            for (unsigned j = 0; j < 4; ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                tpl->addImproper(term);
        }
    }
    if (sys->table("improper_harm") != msys::TermTablePtr()) {
        IdList term_list = sys->table("improper_harm")->findWithAny(atoms);
        for (unsigned i = 0; i < term_list.size(); ++i) {
            IdList term = sys->table("improper_harm")->atoms(term_list[i]);
            bool in_template = true;
            for (unsigned j = 0; j < 4; ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                tpl->addImproper(term);
        }
    }
    if (sys->table("improper_anharm") != msys::TermTablePtr()) {
        IdList term_list = sys->table("improper_anharm")->findWithAny(atoms);
        for (unsigned i = 0; i < term_list.size(); ++i) {
            IdList term = sys->table("improper_anharm")->atoms(term_list[i]);
            bool in_template = true;
            for (unsigned j = 0; j < 4; ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                tpl->addImproper(term);
        }
    }
    /* Add exclusions */
    if (sys->table("exclusion") != msys::TermTablePtr()) {
        IdList terms = sys->table("exclusion")->findWithAny(atoms);
        for (unsigned i = 0; i < terms.size(); ++i) {
            IdList term = sys->table("exclusion")->atoms(terms[i]);
            /* Pseudo exclusions are automatically generated by plug-in */
            if (sys->atom(term[0]).atomic_number == 0
                    || sys->atom(term[1]).atomic_number == 0) continue;
            if (get_separation(sys, term) <= rules->exclusions()) continue;
            term[0] = sys_to_tpl[term[0]];
            term[1] = sys_to_tpl[term[1]];
            if (term[0] != msys::BadId && term[1] != msys::BadId)
                tpl->addExclusion(term);
        }
    }
    /* Add cmap */
    if (sys->table("torsiontorsion_cmap") != msys::TermTablePtr()) {
        IdList terms = sys->table("torsiontorsion_cmap")->findWithAny(atoms);
        for (unsigned i = 0; i < terms.size(); ++i) {
            IdList term = sys->table("torsiontorsion_cmap")->atoms(terms[i]);
            bool in_template = true;
            for (unsigned j = 0; j < 8; ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                tpl->addCmap(term);
        }
    }
    /* Add pseudo sites */
    std::vector<std::string> tables = sys->tableNames();
    for (unsigned i = 0; i < tables.size(); ++i) {
        if (tables[i].compare(0, 8, "virtual_") != 0
                && tables[i].compare(0, 6, "drude_") != 0)
            continue;
        IdList terms = sys->table(tables[i])->findWithAny(atoms);
        for (unsigned j = 0; j < terms.size(); ++j) {
            IdList term = sys->table(tables[i])->atoms(terms[j]);
            if (sys_to_tpl[term[0]] == msys::BadId)
                continue;
            for (unsigned k = 0; k < term.size(); ++k)
                term[k] = sys_to_tpl[term[k]];
            tpl->addPseudoSites(tables[i], term);
        }
    }
    return tpl;
}

void check_template(TemplatedSystemPtr tsys, const IdList& atoms,
        RulesPtr rules, TemplatedSystemPtr tpl, IdList& tpl_to_sys) {
    msys::SystemPtr sys = tsys->system();
    msys::SystemPtr tpl_sys = tpl->system();
    IdList sys_to_tpl(sys->maxAtomId(), msys::BadId);
    /* Map atom names and types for non-pseudos */
    for (unsigned i = 0; i < tpl_to_sys.size(); ++i) {
        if (tpl_to_sys[i] != msys::BadId) {
            if (tpl_sys->atom(i).charge != sys->atom(tpl_to_sys[i]).charge) {
                std::stringstream msg;
                msg << "Atom " << tpl_to_sys[i] << " has charge different "
                    << "from previous atom of the same type '"
                    << tpl->btype(i) << "'";
                VIPARR_FAIL(msg.str());
            }
            sys_to_tpl[tpl_to_sys[i]] = i;
            if (tpl_sys->atom(i).atomic_number > 0) {
                sys->atom(tpl_to_sys[i]).name = tpl_sys->atom(i).name;
                tsys->setTypes(tpl_to_sys[i], tpl->btype(i), tpl->nbtype(i));
            }
        }
    }
    /* Check and map external and pseudo atoms */
    bool new_match = true;
    /* This outer loop allows us to handle pseudos attached to other pseudos */
    while (new_match) {
        new_match = false;
        for (unsigned i = 0; i < tpl_to_sys.size(); ++i) {
            if (tpl_to_sys[i] == msys::BadId)
                continue;
            /* The following assumes that host atoms of pseudos cannot be
             * externally bonded atoms */
            if (tpl_sys->atom(i).atomic_number == -1)
                continue;
            IdList bonded = tpl_sys->bondedAtoms(i);
            IdList sys_bonded = sys->bondedAtoms(tpl_to_sys[i]);
            for (unsigned j = 0; j < bonded.size(); ++j) {
                if (tpl_to_sys[bonded[j]] != msys::BadId)
                    continue;
                if (tpl_sys->atom(bonded[j]).atomic_number > 0)
                    VIPARR_FAIL("VIPARR bug--non-pseudo atom not mapped by "
                            "template match");
                unsigned k = 0;
                for ( ; k < sys_bonded.size(); ++k) {
                    if (sys_to_tpl[sys_bonded[k]] != msys::BadId)
                        continue;
                    if (sys->atom(sys_bonded[k]).atomic_number == 0
                            && tpl_sys->atom(bonded[j]).atomic_number == -1)
                        continue;
                    if (sys->atom(sys_bonded[k]).atomic_number > 0
                            && tpl_sys->atom(bonded[j]).atomic_number == 0)
                        continue;
                    /* This assumes that if corresponding atoms in identical
                     * residues have multiple attached pseudos, the pseudos are
                     * specified in the same atom ID order. Otherwise a separate
                     * graph match must be implemented here to match the
                     * attached pseudo graphs */
                    sys_to_tpl[sys_bonded[k]] = bonded[j];
                    tpl_to_sys[bonded[j]] = sys_bonded[k];
                    if (sys->atom(sys_bonded[k]).atomic_number == 0) {
                        if (tpl_sys->atom(bonded[j]).charge
                                != sys->atom(sys_bonded[k]).charge) {
                            std::stringstream msg;
                            msg << "Pseudo " << sys_bonded[k] << " has charge different "
                                << "from previous atom of the same type '"
                                << tpl->btype(bonded[j]) << "'";
                            VIPARR_FAIL(msg.str());
                        }
                        sys->atom(sys_bonded[k]).name
                            = tpl_sys->atom(bonded[j]).name;
                        tsys->setTypes(sys_bonded[k], tpl->btype(bonded[j]),
                                tpl->nbtype(bonded[j]), tpl->pset(bonded[j]));
                    }
                    new_match = true;
                    break;
                }
                if (k == sys_bonded.size()) {
                    std::stringstream msg;
                    msg << "Atom " << tpl_to_sys[i] << " has fewer attached "
                        << "pseudos or bonds to a different residue than "
                        << "previous atom of same type '"
                        << tpl->btype(i) << "'";
                    VIPARR_FAIL(msg.str());
                }
            }
            for (unsigned j = 0; j < sys_bonded.size(); ++j) {
                if (sys_to_tpl[sys_bonded[j]] == msys::BadId) {
                    std::stringstream msg;
                    msg << "Atom " << tpl_to_sys[i] << " has more attached "
                        << "pseudos or bonds to a different residue than "
                        << "previous atom of same type '"
                        << tpl->btype(i) << "'";
                    VIPARR_FAIL(msg.str());
                }
            }
        }
    }

    /* Shared error message */
    std::stringstream msg;
    msg << "Residues " << sys->residue(sys->atom(atoms[0]).residue).name
        << sys->residue(sys->atom(atoms[0]).residue).resid << " and "
        << tpl_sys->residue(0).name
        << " have matching topologies but different ";

    /* Check impropers */
    std::set<IdList> impropers;
    if (sys->table("dihedral_trig") != msys::TermTablePtr()) {
        IdList term_list = sys->table("dihedral_trig")->findWithAny(atoms);
        for (unsigned i = 0; i < term_list.size(); ++i) {
            IdList term = sys->table("dihedral_trig")->atoms(term_list[i]);
            if (is_dihedral(sys, term)) continue;
            bool in_template = true;
            for (unsigned j = 0; j < 4; ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                impropers.insert(term);
        }
    }
    if (sys->table("improper_harm") != msys::TermTablePtr()) {
        IdList term_list = sys->table("improper_harm")->findWithAny(atoms);
        for (unsigned i = 0; i < term_list.size(); ++i) {
            IdList term = sys->table("improper_harm")->atoms(term_list[i]);
            bool in_template = true;
            for (unsigned j = 0; j < 4; ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                impropers.insert(term);
        }
    }
    if (sys->table("improper_anharm") != msys::TermTablePtr()) {
        IdList term_list = sys->table("improper_anharm")->findWithAny(atoms);
        for (unsigned i = 0; i < term_list.size(); ++i) {
            IdList term = sys->table("improper_anharm")->atoms(term_list[i]);
            bool in_template = true;
            for (unsigned j = 0; j < 4; ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                impropers.insert(term);
        }
    }
    if (impropers.size() != tpl->impropers().size()) {
        msg << "impropers";
        VIPARR_FAIL(msg.str());
    }
    for (unsigned i = 0; i < tpl->impropers().size(); ++i)
        if (impropers.find(tpl->impropers()[i]) == impropers.end()) {
            msg << "impropers";
            VIPARR_FAIL(msg.str());
        }
    /* Check exclusions */
    std::set<IdList> exclusions;
    if (sys->table("exclusion") != msys::TermTablePtr()) {
        IdList terms = sys->table("exclusion")->findWithAny(atoms);
        for (unsigned i = 0; i < terms.size(); ++i) {
            IdList term = sys->table("exclusion")->atoms(terms[i]);
            if (get_separation(sys, term) <= rules->exclusions()) continue;
            term[0] = sys_to_tpl[term[0]];
            term[1] = sys_to_tpl[term[1]];
            if (term[0] != msys::BadId && term[1] != msys::BadId
                    && tpl_sys->atom(term[0]).atomic_number != 0
                    && tpl_sys->atom(term[1]).atomic_number != 0)
                exclusions.insert(term);
        }
    }
    if (exclusions.size() != tpl->exclusions().size()) {
        msg << "exclusions";
        VIPARR_FAIL(msg.str());
    }
    for (unsigned i = 0; i < tpl->exclusions().size(); ++i) {
        IdList excl = tpl->exclusions()[i];
        if (exclusions.find(excl) == exclusions.end()) {
            std::reverse(excl.begin(), excl.end());
            if (exclusions.find(excl) == exclusions.end()) {
                msg << "exclusions";
                VIPARR_FAIL(msg.str());
            }
        }
    }
    /* Check cmap */
    std::set<IdList> cmaps;
    if (sys->table("torsiontorsion_cmap") != msys::TermTablePtr()) {
        IdList terms = sys->table("torsiontorsion_cmap")->findWithAny(atoms);
        for (unsigned i = 0; i < terms.size(); ++i) {
            IdList term = sys->table("torsiontorsion_cmap")->atoms(terms[i]);
            bool in_template = true;
            for (unsigned j = 0; j < term.size(); ++j) {
                term[j] = sys_to_tpl[term[j]];
                if (term[j] == msys::BadId)
                    in_template = false;
            }
            if (in_template)
                cmaps.insert(term);
        }
    }
    if (cmaps.size() != tpl->cmaps().size()) {
        msg << "cmaps";
        VIPARR_FAIL(msg.str());
    }
    for (unsigned i = 0; i < tpl->cmaps().size(); ++i) {
        IdList cmap = tpl->cmaps()[i];
        if (cmaps.find(cmap) == cmaps.end()) {
            std::reverse(cmap.begin(), cmap.end());
            if (cmaps.find(cmap) == cmaps.end()) {
                msg << "cmaps";
                VIPARR_FAIL(msg.str());
            }
        }
    }
    /* Check pseudo sites */
    std::vector<std::string> tables = sys->tableNames();
    for (unsigned i = 0; i < tables.size(); ++i) {
        if (tables[i].compare(0, 8, "virtual_") != 0
                )
            continue;
        std::set<IdList> pseudos;
        IdList terms = sys->table(tables[i])->findWithAny(atoms);
        for (unsigned j = 0; j < terms.size(); ++j) {
            IdList term = sys->table(tables[i])->atoms(terms[j]);
            if (sys_to_tpl[term[0]] == msys::BadId)
                continue;
            for (unsigned k = 0; k < term.size(); ++k)
                term[k] = sys_to_tpl[term[k]];
            pseudos.insert(term);
        }
        TemplatedSystem::PseudoType pseudo_type;
        unsigned j = 0;
        for ( ; j < tpl->pseudoTypes().size(); ++j) {
            if (tpl->pseudoTypes()[j].name == tables[i]) {
                pseudo_type = tpl->pseudoTypes()[j];
                break;
            }
        }
        if (j == tpl->pseudoTypes().size()) {
            if (pseudos.size() == 0)
                continue;
            else {
                msg << "pseudo types";
                VIPARR_FAIL(msg.str());
            }
        }
        if (pseudos.size() != pseudo_type.sites_list.size()) {
            msg << "pseudo-site lists for pseudo type " << tables[i];
            VIPARR_FAIL(msg.str());
        }
        for (unsigned j = 0; j < pseudo_type.sites_list.size(); ++j)
            if (pseudos.find(pseudo_type.sites_list[j]) == pseudos.end()) {
                msg << "pseudo-site lists for pseudo type " << tables[i];
                VIPARR_FAIL(msg.str());
            }
    }
}

struct PluginInfo {
    std::string ff_file;
    SystemToPatternPtr stp;
    std::vector<PermutationPtr> perms;
    PluginInfo() { }
    PluginInfo(const std::string& _ff_file, SystemToPatternPtr _stp,
            const std::vector<PermutationPtr>& _perms)
        : ff_file(_ff_file), stp(_stp), perms(_perms) { }
};

void add_parameters(TemplatedSystemPtr tsys, const IdList& atoms,
        const PluginInfo& info, ForcefieldPtr ff) {
    msys::SystemPtr sys = tsys->system();
    msys::TermTablePtr term_table;
    IdList terms;
    msys::ParamTablePtr param_table;
    if (info.ff_file == "mass") {
        terms = atoms;
        if (Forcefield::HasParamTable("mass"))
            param_table = Forcefield::ParamTable("mass");
        else {
            param_table = msys::ParamTable::create();
            Forcefield::AddParamTable("mass", param_table);
        }
        param_table->addProp("type", msys::StringType);
        param_table->addProp("amu", msys::FloatType);
    } else {
        if (info.ff_file.compare(0, 9, "virtuals_") == 0)
            term_table = sys->table("virtual_"
                    + info.ff_file.substr(9));
        else if (info.ff_file == "improper_trig")
            term_table = sys->table("dihedral_trig");
        else if (info.ff_file == "ureybradley_harm")
            term_table = sys->table("stretch_harm");
        else if (info.ff_file == "vdw1")
            term_table = sys->table("nonbonded");
        else
            term_table = sys->table(info.ff_file);
        if (term_table == msys::TermTablePtr())
            return;
        if (term_table->params() == msys::ParamTablePtr())
            VIPARR_FAIL("Term table '" + term_table->name() + "' has no "
                    "param table");
        terms = term_table->findWithAny(atoms);
        if (Forcefield::HasParamTable(info.ff_file))
            param_table = Forcefield::ParamTable(info.ff_file);
        else {
            param_table = msys::ParamTable::create();
            Forcefield::AddParamTable(info.ff_file, param_table);
        }
        param_table->addProp("type", msys::StringType);
        for (unsigned i = 0; i < term_table->params()->propCount(); ++i)
            param_table->addProp(term_table->params()->propName(i),
                    term_table->params()->propType(i));
        if (info.ff_file == "vdw1")
            param_table->addProp("nbfix_identifier", msys::StringType);
    }
    if (info.ff_file == "pseudopol_fermi") {
        if (terms.size() > 0) {
            Id param = param_table->addParam();
            /* 'type' is not actually used for pseudopol_fermi */
            for (unsigned k = 0; k < term_table->params()->propCount(); ++k)
                param_table->value(param, term_table->params()->propName(k))
                    = term_table->propValue(terms[0], k);
            param_table->value(param, "type") = "DUMMY";
            VIPARR_OUT << "Creating forcefield table '" << info.ff_file << "'"
                << std::endl;
            ff->appendParam("pseudopol_fermi", param);
            return;
        }
    }
    ff->clearParams(info.ff_file);
    /* Store which params are in this ff table */
    std::vector<bool> in_ff;
    if (info.ff_file == "mass")
        in_ff.resize(sys->maxAtomId(), false);
    else
        in_ff.resize(term_table->params()->paramCount(), false);
    for (unsigned i = 0; i < terms.size(); ++i) {
        IdList term_atoms;
        if (info.ff_file == "mass")
            term_atoms = IdList(1, terms[i]);
        else
            term_atoms = term_table->atoms(terms[i]);
        /* The following disambiguation of dihedrals/impropers and
         * stretch/ureybradley may not work if the molecule has cycles of
         * length <= 4 atoms */
        if (info.ff_file == "improper_trig" && is_dihedral(sys, term_atoms))
            continue;
        if (info.ff_file == "dihedral_trig" && !is_dihedral(sys, term_atoms))
            continue;
        if (info.ff_file == "ureybradley_harm"
                && sys->findBond(term_atoms[0], term_atoms[1]) != msys::BadId)
            continue;
        if (info.ff_file == "stretch_harm" &&
                sys->findBond(term_atoms[0], term_atoms[1]) == msys::BadId)
            continue;
        if (info.ff_file == "ureybradley_harm") {
            IdList bonded = sys->bondedAtoms(term_atoms[0]);
            for (unsigned j = 0; j < bonded.size(); ++j)
                if (sys->findBond(bonded[j], term_atoms[1]) != msys::BadId) {
                    term_atoms.insert(term_atoms.begin() + 1, bonded[j]);
                    break;
                }
            if (term_atoms.size() == 2) {
                std::stringstream ss;
                ss << "Middle atom for ureybradley term " << term_atoms[0]
                    << " " << term_atoms[1] << " not found";
                VIPARR_FAIL(ss.str());
            }
        }
        Pattern pattern = (*info.stp)(tsys, term_atoms);
        bool found = false;
        std::stringstream type;
        for (unsigned j = 0; j < info.perms.size(); ++j) {
            /* Convert permuted pattern into type string */
            Pattern dup = (*info.perms[j])(pattern);
            type.str("");
            type << dup.atoms[0];
            for (unsigned atom_ind = 1; atom_ind < dup.atoms.size(); ++atom_ind)
                type << " " << dup.atoms[atom_ind];
            for (unsigned k = 0; k < dup.flags.size(); ++k)
                type << " " << dup.flags[k];
            /* Check if type string already exists */
            IdList param_matches_vec = param_table->findString(
                    param_table->propIndex("type"), type.str());
            std::list<Id> param_matches(param_matches_vec.begin(),
                    param_matches_vec.end());
            for (std::list<Id>::iterator iter = param_matches.begin();
                    iter != param_matches.end(); ) {
                if (*iter >= in_ff.size() || !in_ff[*iter])
                    iter = param_matches.erase(iter);
                else
                    ++iter;
            }
            if (param_matches.size() > 0) {
                /* Check that existing parameters agree with new ones */
                found = true;
                if (info.ff_file == "mass") {
                    if (sys->atom(terms[i]).mass
                            != param_table->value(param_matches.front(),
                                "amu").asFloat()) {
                        std::stringstream msg;
                        msg << "Atom " << i << " has mass different from "
                            << "previous atom of the same type '"
                            << type.str() << "'";
                        VIPARR_FAIL(msg.str());
                    }
                } else {
                    for (unsigned k = 0; k < term_table->params()->propCount();
                            ++k) {
                        if (term_table->params()->propName(k) == "type"
                                || term_table->params()->propName(k) == "memo")
                            continue;
                        if (param_table->value(param_matches.front(),
                                    term_table->params()->propName(k))
                                != term_table->propValue(terms[i], k)) {
                            std::stringstream msg;
                            msg << "Term " << i << " in '" << info.ff_file
                                << "' table has params different from "
                                << "previous term of the same type '"
                                << type.str() << "'";
                            VIPARR_FAIL(msg.str());
                        }
                    }
                    if (info.ff_file != "dihedral_trig"
                            && param_matches.size() != 1)
                        VIPARR_FAIL("Atom type " + type.str()
                                + " cannot match multiple "
                                + info.ff_file + " params");
                    else if (info.ff_file == "dihedral_trig") {
                        /* Assume that multiple dihedral trig params must match
                         * previous group of dihedral trig params in the same
                         * order */
                        IdList next_term_atoms;
                        if (i+1 < terms.size())
                            next_term_atoms = term_table->atoms(terms[i+1]);
                        std::list<Id>::iterator iter = param_matches.begin();
                        while (term_atoms.size() == next_term_atoms.size()
                                && std::equal(term_atoms.begin(),
                                    term_atoms.end(),
                                    next_term_atoms.begin())) {
                            ++i;
                            ++iter;
                            if (iter == param_matches.end()) {
                                std::stringstream msg;
                                msg << "Term " << i << " in 'dihedral_trig' "
                                    << "table has more duplicate params than "
                                    << "previous term of the same type '"
                                    << type.str() << "'";
                                VIPARR_FAIL(msg.str());
                            }
                            for (unsigned k = 0;
                                    k < term_table->params()->propCount();
                                    ++k) {
                                if (term_table->params()->propName(k) == "type"
                                        || term_table->params()->propName(k) ==
                                        "memo")
                                    continue;
                                if (param_table->value(*iter,
                                            term_table->params()->propName(k))
                                        != term_table->propValue(terms[i], k)) {
                                    std::stringstream msg;
                                    msg << "Term " << i << " in 'dihedral_trig'"
                                        << " table has params different from "
                                        << "previous term of the same type '"
                                        << type.str() << "'";
                                    VIPARR_FAIL(msg.str());
                                }
                            }
                            if (i+1 < terms.size())
                                next_term_atoms = term_table->atoms(terms[i+1]);
                            else
                                next_term_atoms.clear();
                        }
                        ++iter;
                        if (iter != param_matches.end()) {
                            std::stringstream msg;
                            msg << "Term " << i << " in 'dihedral_trig' "
                                << "table has fewer duplicate params than "
                                << "previous term of the same type '"
                                << type.str() << "'";
                            VIPARR_FAIL(msg.str());
                        }
                    }
                }
            }
        }
        if (!found) {
            if (ff->rowIDs(info.ff_file).size() == 0)
                VIPARR_OUT << "Creating forcefield table '" << info.ff_file
                    << "'" << std::endl;
            /* Add new param to ff */
            Id param = param_table->addParam();
            if (info.ff_file == "mass")
                param_table->value(param, "amu")
                    = sys->atom(terms[i]).mass;
            else {
                for (unsigned j = 0; j < term_table->params()->propCount(); ++j)
                    param_table->value(param, term_table->params()->propName(j))
                        = term_table->propValue(terms[i], j);
            }
            /* Use type string of last checked permutation */
            param_table->value(param, "type") = type.str();
            if (info.ff_file == "vdw1")
                param_table->value(param, "nbfix_identifier")
                    = ff->rules()->nbfix_identifier;
            ff->appendParam(info.ff_file, param);
            if (param >= in_ff.size())
                in_ff.resize(param + 100, false);
            in_ff[param] = true;
            if (info.ff_file == "dihedral_trig") {
                /* Add possible multiple rows for dihedral trig terms */
                IdList next_term_atoms;
                if (i+1 < terms.size())
                    next_term_atoms = term_table->atoms(terms[i+1]);
                while (term_atoms.size() == next_term_atoms.size()
                        && std::equal(term_atoms.begin(), term_atoms.end(),
                        next_term_atoms.begin())) {
                    ++i;
                    Id param = param_table->addParam();
                    for (unsigned j = 0; j < term_table->params()->propCount();
                            ++j)
                        param_table->value(param,
                                term_table->params()->propName(j))
                            = term_table->propValue(terms[i], j);
                    param_table->value(param, "type") = type.str();
                    ff->appendParam(info.ff_file, param);
                    if (param >= in_ff.size())
                        in_ff.resize(param + 100, false);
                    in_ff[param] = true;
                    if (i+1 < terms.size())
                        next_term_atoms = term_table->atoms(terms[i+1]);
                    else
                        next_term_atoms.clear();
                }
            }
        }
    }
}

void add_vdw2_parameters(TemplatedSystemPtr tsys, const IdList& atoms,
        ForcefieldPtr ff) {
    msys::SystemPtr sys = tsys->system();
    msys::TermTablePtr vdw1_term_table = sys->table("nonbonded");
    if (vdw1_term_table == msys::TermTablePtr())
        vdw1_term_table = sys->table("vdw1");
    if (vdw1_term_table == msys::TermTablePtr()) {
        VIPARR_ERR << "WARNING: nonbonded table not found" << std::endl;
        return;
    }
    msys::OverrideTablePtr overrides = vdw1_term_table->overrides();
    if (overrides == msys::OverrideTablePtr()) return;

    VIPARR_OUT << "Creating forcefield table 'vdw2'" << std::endl;
    std::set<Id> atom_set(atoms.begin(), atoms.end());
    std::vector<std::set<std::string> > vdw1_param_to_types(
            vdw1_term_table->params()->paramCount());
    IdList vdw1_terms = vdw1_term_table->terms();
    for (unsigned i = 0; i < vdw1_terms.size(); ++i) {
        Id atom = vdw1_term_table->atoms(vdw1_terms[i])[0];
        if (atom_set.find(atom) != atom_set.end()) {
            vdw1_param_to_types[vdw1_term_table->param(vdw1_terms[i])].insert(
                    tsys->nbtype(atom));
        }
    }
    msys::ParamTablePtr vdw2_param_table;
    if (Forcefield::HasParamTable("vdw2"))
        vdw2_param_table = Forcefield::ParamTable("vdw2");
    else {
        vdw2_param_table = msys::ParamTable::create();
        Forcefield::AddParamTable("vdw2", vdw2_param_table);
    }
    vdw2_param_table->addProp("type", msys::StringType);
    for (unsigned i = 0; i < overrides->params()->propCount(); ++i)
        vdw2_param_table->addProp(overrides->params()->propName(i),
                overrides->params()->propType(i));
    vdw2_param_table->addProp("nbfix_identifier", msys::StringType);
    std::vector<msys::IdPair> override_list = overrides->list();
    for (unsigned i = 0; i < override_list.size(); ++i) {
        std::set<std::string> t1 = vdw1_param_to_types[override_list[i].first];
        std::set<std::string> t2 = vdw1_param_to_types[override_list[i].second];
        std::set<std::pair<std::string, std::string> > included;
        for (std::set<std::string>::iterator iter1 = t1.begin();
                iter1 != t1.end(); ++iter1) {
            for (std::set<std::string>::iterator iter2 = t2.begin();
                    iter2 != t2.end(); ++iter2) {
                if (included.find(std::make_pair(*iter1, *iter2))
                        != included.end() ||
                        included.find(std::make_pair(*iter1, *iter2))
                        != included.end())
                    continue;
                included.insert(std::make_pair(*iter1, *iter2));
                std::string type = (*iter1) + " " + (*iter2);
                Id param = vdw2_param_table->addParam();
                for (unsigned j = 0; j < overrides->params()->propCount(); ++j) {
                    vdw2_param_table->value(param,
                            overrides->params()->propName(j))
                        = overrides->params()->value(
                                overrides->get(override_list[i]), j);
                }
                vdw2_param_table->value(param, "type") = type;
                vdw2_param_table->value(param, "nbfix_identifier")
                    = ff->rules()->nbfix_identifier;
                ff->appendParam("vdw2", param);
            }
        }
    }
}


void generate_scaled_pairs_if_necessary(TemplatedSystemPtr tsys, const IdList& atoms,
        ForcefieldPtr ff) {
    const Rules::VDWFunc& func = Rules::VDWFuncRegistry().find(
            ff->rules()->vdw_func)->second;
    const Rules::VDWCombRulePtr& vdw_rule = Rules::VDWCombRuleRegistry().find(
            ff->rules()->vdw_comb_rule)->second;

    std::vector<std::string> tables = tsys->system()->tableNames();
    for (unsigned i = 0; i < tables.size(); ++i) {
        if (tables[i].substr(0,4) == "pair"
                && tables[i] != func.pair_table_name)
            VIPARR_FAIL("System has pairs table '" + tables[i]
                    + "' which does not agree with VDW functional form '"
                    + ff->rules()->vdw_func + "'");
    }

    msys::TermTablePtr pair_table = tsys->system()->table(func.pair_table_name);
    if (pair_table == msys::TermTablePtr()) return;
    IdList terms = pair_table->findWithAny(atoms);
    if (terms.size() == 0) return;

    msys::SystemPtr sys = tsys->system();
    msys::TermTablePtr vdw_term_table = sys->table("nonbonded");
    if (vdw_term_table == msys::TermTablePtr())
        vdw_term_table = sys->table("vdw1");
    if (vdw_term_table == msys::TermTablePtr()) {
        VIPARR_FAIL("nonbonded table not found");
    }
    const std::vector<std::string>& vdw_props = func.param_names;

    std::vector<std::string> pair_params = func.pair_param_names;
    pair_params.push_back("qij");
    VIPARR_OUT << "Creating forcefield table 'scaled_pair_overrides'"
        << std::endl;
    ff->clearParams("scaled_pair_overrides");
    ff->rules()->plugins.push_back("scaled_pair_overrides");
    msys::ParamTablePtr scaled_pair_overrides = msys::ParamTable::create();
    scaled_pair_overrides->addProp("type", msys::StringType);
    scaled_pair_overrides->addProp("separation", msys::IntType);
    for (unsigned p = 0; p < pair_params.size(); ++p) {
        if (pair_table->params()->propIndex(pair_params[p]) == msys::BadId)
            VIPARR_FAIL("Pairs table is missing property '" + pair_params[p]
                    + "'");
        scaled_pair_overrides->addProp(pair_params[p], msys::FloatType);
    }
    Forcefield::AddParamTable("scaled_pair_overrides", scaled_pair_overrides);
    bool added_scaled_pair_terms = false;
    for (unsigned i = 0; i < terms.size(); ++i) {
        IdList term_atoms = pair_table->atoms(terms[i]);
        unsigned sep = get_separation(tsys->system(), term_atoms);
        if (sep < 2 || sep > ff->rules()->exclusions()) {
            std::stringstream msg;
            msg << "Atom pair (" << term_atoms[0] << "," << term_atoms[1]
                << ") in pairs table has bond separation " << sep
                << "; this must be between 2 and the exclusion rule "
                << "of the input forcefield";
            VIPARR_FAIL(msg.str());
        }
        /* See if a scaled pair parameter for this type pair already exists */
        Id existing = msys::BadId;
        std::string type = tsys->nbtype(term_atoms[0]) + " "
            + tsys->nbtype(term_atoms[1]);
        IdList param_matches_vec = scaled_pair_overrides->findString(
                scaled_pair_overrides->propIndex("type"), type);
        std::list<Id> param_matches(param_matches_vec.begin(),
                param_matches_vec.end());
        if (tsys->nbtype(term_atoms[0]) != tsys->nbtype(term_atoms[1])) {
            type = tsys->nbtype(term_atoms[1]) + " "
                + tsys->nbtype(term_atoms[0]);
            param_matches_vec = scaled_pair_overrides->findString(
                    scaled_pair_overrides->propIndex("type"), type);
            param_matches.insert(param_matches.end(),
                    param_matches_vec.begin(), param_matches_vec.end());
        }
        for (std::list<Id>::iterator iter = param_matches.begin();
                iter != param_matches.end(); ) {
            if (scaled_pair_overrides->value(*iter, "separation") != int(sep))
                iter = param_matches.erase(iter);
            else
                ++iter;
        }
        if (param_matches.size() == 1)
            existing = param_matches.front();
        else if (param_matches.size() > 1)
            VIPARR_FAIL("VIPARR BUG--multiple scaled pairs for same types");
        if (existing != msys::BadId) {
            /* Check that parameters are same as previous scaled pair */
            for (unsigned p = 0; p < pair_params.size(); ++p) {
                if (fabs(scaled_pair_overrides->value(existing,
                            pair_params[p]).asFloat()
                        - pair_table->propValue(terms[i],
                            pair_params[p]).asFloat()) > 1e-8) {
                    std::stringstream msg;
                    msg << "Atoms " << term_atoms[0] << " and " << term_atoms[1]
                        << " at separation " << sep << " have param '"
                        << pair_params[p] << "' value "
                        << pair_table->propValue(terms[i],
                                pair_params[p]).asFloat()
                        << "; previous pair of same type '"
                        << scaled_pair_overrides->value(existing,
                                "type").asString() << "' has value "
                        << scaled_pair_overrides->value(existing,
                                pair_params[p]).asFloat();
                    VIPARR_FAIL(msg.str());
                }
            }
        } else {
            msys::Id row0 = vdw_lookup(vdw_term_table, term_atoms[0], true);
            msys::Id row1 = vdw_lookup(vdw_term_table, term_atoms[1], true);
            if (row0 == msys::BadId || row1 == msys::BadId)
                VIPARR_FAIL("couldnt find row in nonbonded table atom(s)");
            std::vector<double> vi(vdw_props.size());
            std::vector<double> vj(vdw_props.size());
            for (unsigned p = 0; p < vdw_props.size(); ++p){
                vi[p] = vdw_term_table->propValue(row0, vdw_props[p]).asFloat();
                vj[p] = vdw_term_table->propValue(row1, vdw_props[p]).asFloat();
            }
            std::vector<double> pcomb = (*vdw_rule)(vi, vj, ff->rules()->lj_scale(sep));
            if (pcomb.size() != func.pair_param_names.size())
                VIPARR_FAIL("not enough combined parameters returned");
            double qij = ff->rules()->es_scale(sep) * tsys->system()->atom(term_atoms[0]).charge
            * tsys->system()->atom(term_atoms[1]).charge;
            pcomb.push_back(qij);
            bool equivalent = true;
            for (unsigned p = 0; p < pair_params.size(); ++p){
                double pval = pair_table->propValue(terms[i], pair_params[p]).asFloat();
                double magdelta = fabs( pcomb[p] - pval);
                double magref = (1e-8 + 1e-05 * fabs(pval));
                bool isclose = magdelta <= magref;
                if( not isclose) {
                    VIPARR_OUT << "Atoms " << term_atoms[0] << " and " << term_atoms[1]
                        << " at separation " << sep << " have param '"
                        << pair_params[p] << "' value "
                        << pval << " and simple value "
                        << pcomb[p] << " magdelta is " << magdelta
                        << " tolcheck is " << isclose << std::endl;
                    equivalent = false;
                    break;
                }
            }
            if(not equivalent){
                Id param = scaled_pair_overrides->addParam();
                std::string type = tsys->nbtype(term_atoms[0]) + " "
                            + tsys->nbtype(term_atoms[1]);
                scaled_pair_overrides->value(param, "type") = type;
                scaled_pair_overrides->value(param, "separation") = sep;
                for (unsigned p = 0; p < pair_params.size(); ++p)
                    scaled_pair_overrides->value(param, pair_params[p])
                           = pair_table->propValue(terms[i], pair_params[p]);
                ff->appendParam("scaled_pair_overrides", param);
                added_scaled_pair_terms = true;
            }
        }
    }
    if ( not added_scaled_pair_terms){
        VIPARR_OUT << "  Forcefield table 'scaled_pair_overrides' was unnecessary... Removing"
        << std::endl;

        ff->rules()->plugins.pop_back();
    }
}

} // namespace

ForcefieldPtr desres::viparr::ExecuteIviparr(msys::SystemPtr sys, const msys::IdList& atoms,
        RulesPtr rules, bool templateonly) {

    VIPARR_OUT << "Applying iViparr to " << atoms.size() << " atoms"
        << std::endl;

    /* Check compatibility of vdw_func and vdw_comb_rule */
    std::transform(sys->nonbonded_info.vdw_funct.begin(),
            sys->nonbonded_info.vdw_funct.end(),
            sys->nonbonded_info.vdw_funct.begin(),
            tolower);
    std::transform(sys->nonbonded_info.vdw_rule.begin(),
            sys->nonbonded_info.vdw_rule.end(),
            sys->nonbonded_info.vdw_rule.begin(),
            tolower);
    if (rules->vdw_func != "" && Rules::VDWFuncRegistry().find(rules->vdw_func)
            == Rules::VDWFuncRegistry().end())
        VIPARR_FAIL("Unrecognized vdw_func in input forcefield");
    std::string sys_vdw_func;
    for (std::map<std::string, Rules::VDWFunc>::const_iterator iter
            = Rules::VDWFuncRegistry().begin(); iter
            != Rules::VDWFuncRegistry().end(); ++iter)
        if (iter->second.vdw_table_name == sys->nonbonded_info.vdw_funct)
            sys_vdw_func = iter->first;
    if (sys_vdw_func == "")
        VIPARR_FAIL("Unrecognized vdw_funct in input system");
    if (rules->vdw_func == "")
        rules->vdw_func = sys_vdw_func;
    else if (rules->vdw_func != sys_vdw_func)
        VIPARR_FAIL("VDW functional form in input system and forcefield "
                "do not agree");
    if (rules->vdw_comb_rule != ""
            && Rules::VDWCombRuleRegistry().find(rules->vdw_comb_rule)
            == Rules::VDWCombRuleRegistry().end())
        VIPARR_FAIL("Unrecognized vdw_comb_rule in input forcefield");
    if (Rules::VDWCombRuleRegistry().find(sys->nonbonded_info.vdw_rule)
            == Rules::VDWCombRuleRegistry().end())
        VIPARR_FAIL("Unrecognized vdw_rule in input system");
    if (rules->vdw_comb_rule == "")
        rules->vdw_comb_rule = sys->nonbonded_info.vdw_rule;
    else if (rules->vdw_comb_rule != sys->nonbonded_info.vdw_rule)
        VIPARR_FAIL("VDW combine rules in input system and forcefield "
                "do not agree");

    std::vector<bool> in_selection(sys->maxAtomId(), false);
    for (unsigned i = 0; i < atoms.size(); ++i)
        in_selection[atoms[i]] = true;
    TemplatedSystemPtr tsys = TemplatedSystem::create(sys);
    std::vector<IdList> fragments;
    sys->updateFragids(&fragments);

    /* Create templates */
    VIPARR_OUT << "Creating forcefield templates" << std::endl;
    TemplateTyperPtr typer = TemplateTyper::create();
    for (unsigned i = 0; i < fragments.size(); ++i) {
        /* Split into residues */
        std::map<Id, IdList> residues;
        bool process = in_selection[fragments[i][0]];
        for (unsigned j = 0; j < fragments[i].size(); ++j) {
            if (process ^ in_selection[fragments[i][j]])
                VIPARR_FAIL("Atom selection contains an incomplete fragment");
            Id resid = sys->atom(fragments[i][j]).residue;
            std::map<Id, IdList>::iterator iter
                = residues.insert(std::make_pair(resid, IdList())).first;
            iter->second.push_back(fragments[i][j]);
        }
        if (!process) continue;
        for (std::map<Id, IdList>::iterator iter = residues.begin();
                iter != residues.end(); ++iter) {
            IdList tpl_to_sys;
            std::stringstream why_not;
            TemplatedSystemPtr tpl = typer->findMatch(tsys, iter->second,
                    "", 0, tpl_to_sys, why_not);
            if (tpl == TemplatedSystemPtr())
                /* Residue does not match any existing templates; create new
                 * template */
                typer->addTemplate(create_template(tsys, iter->second, rules));
            else
                /* Residue matches existing template in structure; check that
                 * all other template properties match */
                check_template(tsys, iter->second, rules, tpl, tpl_to_sys);
        }
    }

    /* If pseudopol_fermi table is present, need to rename certain atom types
     * O, H, C, and NH1 */
    msys::TermTablePtr pseudopol_fermi = sys->table("pseudopol_fermi");
    if (pseudopol_fermi != msys::TermTablePtr()
            && pseudopol_fermi->termCount() > 0) {
        IdList terms = pseudopol_fermi->findWithAny(atoms);
        std::set<std::string> O_types;
        std::set<std::string> H_types;
        std::set<std::string> C_types;
        std::set<std::string> NH1_types;
        for (unsigned i = 0; i < terms.size(); ++i) {
            IdList term_atoms = pseudopol_fermi->atoms(terms[i]);
            O_types.insert(tsys->btype(term_atoms[0]));
            O_types.insert(tsys->btype(term_atoms[2]));
            H_types.insert(tsys->btype(term_atoms[1]));
            H_types.insert(tsys->btype(term_atoms[3]));
            for (unsigned idx = 0; idx < 4; ++idx) {
                IdList bonded = sys->bondedAtoms(term_atoms[idx]);
                for (unsigned j = 0; j < bonded.size(); ++j) {
                    if (sys->atom(bonded[j]).atomic_number == 0)
                        continue;
                    if (idx % 2 == 0)
                        C_types.insert(tsys->btype(bonded[j]));
                    else
                        NH1_types.insert(tsys->btype(bonded[j]));
                    break;
                }
            }
        }
        for (unsigned i = 0; i < atoms.size(); ++i) {
            if (O_types.find(tsys->btype(atoms[i])) != O_types.end())
                tsys->setTypes(atoms[i], "O", tsys->nbtype(atoms[i]));
            else if (H_types.find(tsys->btype(atoms[i])) != H_types.end())
                tsys->setTypes(atoms[i], "H", tsys->nbtype(atoms[i]));
            else if (C_types.find(tsys->btype(atoms[i])) != C_types.end())
                tsys->setTypes(atoms[i], "C", tsys->nbtype(atoms[i]));
            else if (NH1_types.find(tsys->btype(atoms[i])) != NH1_types.end())
                tsys->setTypes(atoms[i], "NH1", tsys->nbtype(atoms[i]));
        }
        std::vector<TemplatedSystemPtr> templates = typer->templates();
        for (unsigned j = 0; j < templates.size(); ++j) {
            TemplatedSystemPtr tpl = templates[j];
            IdList tpl_atoms = tpl->system()->atoms();
            for (unsigned i = 0; i < tpl_atoms.size(); ++i) {
                if (O_types.find(tpl->btype(tpl_atoms[i])) != O_types.end())
                    tpl->setTypes(tpl_atoms[i], "O", tpl->nbtype(tpl_atoms[i]));
                else if (H_types.find(tpl->btype(atoms[i])) != H_types.end())
                    tpl->setTypes(tpl_atoms[i], "H", tpl->nbtype(tpl_atoms[i]));
                else if (C_types.find(tpl->btype(atoms[i])) != C_types.end())
                    tpl->setTypes(tpl_atoms[i], "C", tpl->nbtype(tpl_atoms[i]));
                else if (NH1_types.find(tpl->btype(atoms[i]))
                        != NH1_types.end())
                    tpl->setTypes(tpl_atoms[i], "NH1",
                            tpl->nbtype(tpl_atoms[i]));
            }
        }
    }

    ForcefieldPtr ff = Forcefield::create(rules, typer);
    if (templateonly)
        return ff;

    /* Add parameters based on plugins in rules.
     * NOTE: Inverse plugin behavior is hard-coded here and in the
     * add_parameters function; this may need to be changed if the behavior of a
     * plugin is changed. */
    for (std::vector<std::string>::iterator iter = rules->plugins.begin();
           iter != rules->plugins.end(); ) {
        /* Ignore these plugins */
        if (*iter == "exclusions"
                || *iter == "exclusions_without_pairs"
                || *iter == "pairs_es_scaled"
                || *iter == "pairs_lj_scaled") {
            ++iter;
            continue;
        }
        /* The functionality of these plugins will be handled by hard-coded
         * template charges and the scaled_pair_overrides plugin */
        if (*iter == "charges_bci"
                || *iter == "charges_formal"
                || *iter == "pairs_es_balanced"
                || *iter == "pairs_lj_scaled_14"
                || *iter == "scaled_pair_overrides") {
            iter = rules->plugins.erase(iter);
            continue;
        }
        /* Handle these plugins */
        std::vector<PluginInfo> infos;
        std::vector<PermutationPtr> fwbw_perm;
        fwbw_perm.push_back(Permutation::Identity);
        fwbw_perm.push_back(Permutation::Reverse);
        std::vector<PermutationPtr> id_perm(1, Permutation::Identity);
        if (*iter == "angles")
            infos.push_back(PluginInfo("angle_harm", SystemToPattern::Bonded,
                        fwbw_perm));
        else if (*iter == "bonds")
            infos.push_back(PluginInfo("stretch_harm", SystemToPattern::Bonded,
                        fwbw_perm));
        else if (*iter == "cmap")
            infos.push_back(PluginInfo("torsiontorsion_cmap",
                    SystemToPattern::BType, fwbw_perm));
        else if (*iter == "impropers") {
            infos.push_back(PluginInfo("improper_harm", SystemToPattern::BType,
                        fwbw_perm));
            infos.push_back(PluginInfo("improper_anharm",
                        SystemToPattern::BType, id_perm));
            infos.push_back(PluginInfo("improper_trig", SystemToPattern::BType,
                        fwbw_perm));
        }
        else if (*iter == "mass")
            infos.push_back(PluginInfo("mass", SystemToPattern::BType,
                        id_perm));
        else if (*iter == "mass2")
            infos.push_back(PluginInfo("mass", SystemToPattern::NBType,
                        id_perm));
        else if (*iter == "propers" || *iter == "propers_allowmissing")
            infos.push_back(PluginInfo("dihedral_trig",
                        SystemToPattern::Bonded, fwbw_perm));
        else if (*iter == "pseudopol_fermi")
            infos.push_back(PluginInfo("pseudopol_fermi", SystemToPatternPtr(),
                        std::vector<PermutationPtr>()));
        else if (*iter == "ureybradley")
            infos.push_back(PluginInfo("ureybradley_harm",
                        SystemToPattern::Bonded, fwbw_perm));
        else if (*iter == "virtuals" || *iter == "virtuals_regular"
                ) {
            std::vector<std::string> tables = sys->tableNames();
            for (unsigned j = 0; j < tables.size(); ++j) {
                if (tables[j].compare(0, 8, "virtual_") == 0) {
                    std::string vname = "virtuals_";
                    vname += tables[j].substr(8);
                    if (tables[j] == "virtual_midpoint")
                        VIPARR_FAIL("iviparr currently does not support "
                                "virtual_midpoint pseudos");
                    else if (tables[j] == "virtual_fdat3")
                        infos.push_back(PluginInfo(vname,
                                    SystemToPattern::PseudoBondToSecond,
                                    id_perm));
                    else
                        infos.push_back(PluginInfo(vname,
                                    SystemToPattern::PseudoBType, id_perm));
                }
            }
        } else if (*iter == "vdw1"
                )
            infos.push_back(PluginInfo("vdw1", SystemToPattern::NBType,
                        fwbw_perm));
        else if (*iter == "vdw2") {
            add_vdw2_parameters(tsys, atoms, ff);
            if (ff->rowIDs("vdw2").size() == 0) {
                VIPARR_ERR << "WARNING: Override table for plugin vdw2"
                    << " not found; removing plugin from output forcefield"
                    << std::endl;
                iter = rules->plugins.erase(iter);
            } else
                ++iter;
            continue;
        } else
            VIPARR_FAIL("iviparr does not recognize plugin: " + *iter);
        unsigned total_params = 0;
        for (unsigned j = 0; j < infos.size(); ++j) {
            add_parameters(tsys, atoms, infos[j], ff);
            total_params += ff->rowIDs(infos[j].ff_file).size();
        }
        if (total_params == 0) {
            VIPARR_ERR << "WARNING: Table for plugin " << *iter
                << " not found; removing plugin from output forcefield"
                << std::endl;
            iter = rules->plugins.erase(iter);
        } else
            ++iter;
    }

    /* Add cmap tables */
    std::vector<std::string> aux_tables = sys->auxTableNames();
    std::vector<msys::ParamTablePtr> cmap_tables;
    for (unsigned i = 0; i < aux_tables.size(); ++i) {
        if (aux_tables[i].compare(0, 4, "cmap") != 0)
            continue;
        if (cmap_tables.size() == 0)
            VIPARR_OUT << "Creating cmap tables" << std::endl;
        unsigned cmap_id = atoi(aux_tables[i].substr(4).c_str());
        if (cmap_id > cmap_tables.size())
            cmap_tables.resize(cmap_id);
        cmap_tables[cmap_id - 1] = sys->auxTable(aux_tables[i]);
    }
    for (unsigned i = 0; i < cmap_tables.size(); ++i)
        ff->addCmapTable(cmap_tables[i]);

    /* Generate scaled_pair_overrides table from pairs table if they
     * dont conform to the forcefields ES and LJ scale factors. This is to
     * preserve a user's hand-modifications of pair terms in their input DMS. */
    generate_scaled_pairs_if_necessary(tsys, atoms, ff);
    return ff;
}
