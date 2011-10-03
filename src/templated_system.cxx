#include "base.hxx"
#include "templated_system.hxx"
#include <sstream>
#include <msys/clone.hxx>
#include <msys/param_table.hxx>

using namespace desres::viparr;

TemplatedSystem::TemplatedSystem() :
    _sys(msys::System::create()), _atom_count(0), _max_atom_id(0),
    _bond_count(0), _max_bond_id(0), _type_table(msys::ParamTable::create()),
    _arom_table(msys::ParamTable::create()) {

    _type_table->addProp("btype", msys::StringType);
    _type_table->addProp("nbtype", msys::StringType);
    _type_table->addProp("pset", msys::StringType);
    _arom_table->addProp("aromatic", msys::IntType);
}

TemplatedSystem::TemplatedSystem(msys::SystemPtr sys) :
    _sys(sys), _atom_count(sys->atomCount()),
    _max_atom_id(sys->maxAtomId()), _bond_count(sys->bondCount()),
    _max_bond_id(sys->maxBondId()), _type_table(msys::ParamTable::create()),
    _arom_table(msys::ParamTable::create()) {

    _type_table->addProp("btype", msys::StringType);
    _type_table->addProp("nbtype", msys::StringType);
    _type_table->addProp("pset", msys::StringType);
    _arom_table->addProp("aromatic", msys::IntType);
    for (unsigned i = 0; i < _max_atom_id; ++i) {
        _type_table->addParam();
    }
    for (unsigned i = 0; i < _max_bond_id; ++i)
        _arom_table->addParam();
}

std::shared_ptr<TemplatedSystem> TemplatedSystem::create() {
    return std::make_shared<TemplatedSystem>();
}

std::shared_ptr<TemplatedSystem> TemplatedSystem::create(
        msys::SystemPtr sys) {
    return std::make_shared<TemplatedSystem>(sys);
}

TemplatedSystemPtr TemplatedSystem::clone(const IdList& atoms) {
    /* Clone contained system */
    msys::SystemPtr clone = msys::Clone(_sys, atoms);
    TemplatedSystemPtr tclone = create(clone);

    /* Copy atom types, create atom index map */
    IdList atoms_map(_sys->maxAtomId(), msys::BadId);
    for (unsigned i = 0; i < atoms.size(); ++i) {
        tclone->setTypes(i, this->btype(atoms[i]), this->nbtype(atoms[i]),
                this->pset(atoms[i]));
        atoms_map[atoms[i]] = i;
    }
    
    /* Copy bond aromaticity */
    IdList bonds = _sys->bonds();
    for (unsigned i = 0; i < bonds.size(); ++i) {
        Id ai = _sys->bond(bonds[i]).i;
        Id aj = _sys->bond(bonds[i]).j;
        if (atoms_map[ai] != msys::BadId && atoms_map[aj] != msys::BadId)
            tclone->setAromatic(clone->findBond(atoms_map[ai], atoms_map[aj]),
                    aromatic(_sys->findBond(ai, aj)));
    }

    /* Copy relevant items in lists of typed atoms and tuples */
    for (unsigned i = 0; i < _typed_atoms.size(); ++i) {
        if (atoms_map[_typed_atoms[i][0]] != msys::BadId)
            tclone->addTypedAtom(atoms_map[_typed_atoms[i][0]]);
    }
    const TupleList* tuple_lists[7] = {&_non_pseudo_bonds, &_pseudo_bonds,
        &_angles, &_dihedrals, &_exclusions, &_impropers, &_cmaps};
    void (TemplatedSystem::*adders[7])(const IdList&) = {
        &TemplatedSystem::addNonPseudoBond,
        &TemplatedSystem::addPseudoBond,
        &TemplatedSystem::addAngle,
        &TemplatedSystem::addDihedral,
        &TemplatedSystem::addExclusion,
        &TemplatedSystem::addImproper,
        &TemplatedSystem::addCmap
    };
    for (unsigned i = 0; i < 7; ++i) {
        const TupleList* tuple_list = tuple_lists[i];
        for (unsigned j = 0; j < tuple_list->size(); ++j) {
            IdList new_tuple((*tuple_list)[j].size(), msys::BadId);
            bool copy = true;
            for (unsigned c = 0; c < (*tuple_list)[j].size(); ++c) {
                if (atoms_map[(*tuple_list)[j][c]] != msys::BadId)
                    new_tuple[c] = atoms_map[(*tuple_list)[j][c]];
                else {
                    copy = false;
                    break;
                }
            }
            if (copy)
                ((*tclone.get()).*(adders[i]))(new_tuple);
        }
    }

    /* Copy pseudos */
    for (unsigned i = 0; i < _pseudo_types.size(); ++i) {
        const PseudoType& ptype = _pseudo_types[i];
        for (unsigned j = 0; j < ptype.sites_list.size(); ++j) {
            IdList new_tuple(ptype.nsites, msys::BadId);
            bool copy = true;
            for (unsigned c = 0; c < ptype.nsites; ++c) {
                if (atoms_map[ptype.sites_list[j][c]] != msys::BadId)
                    new_tuple[c] = atoms_map[ptype.sites_list[j][c]];
                else {
                    copy = false;
                    break;
                }
            }
            if (copy)
                tclone->addPseudoSites(ptype.name, new_tuple);
        }
    }
    return tclone;
}

const std::string& TemplatedSystem::hash() {
    updateSystem();
    if (_hash == "")
        _hash = msys::Graph::hash(_sys, _sys->atoms());
    return _hash;
}

desres::msys::GraphPtr TemplatedSystem::graph() {
    updateSystem();
    if (_graph == msys::GraphPtr())
        _graph = msys::Graph::create(_sys, _sys->atoms());
    return _graph;
}

std::string TemplatedSystem::btype(Id atom) const {
    if (atom >= _type_table->paramCount())
        return "";
    else
        return _type_table->value(atom, "btype").asString();
}

std::string TemplatedSystem::nbtype(Id atom) const {
    if (atom >= _type_table->paramCount())
        return "";
    else
        return _type_table->value(atom, "nbtype").asString();
}

std::string TemplatedSystem::pset(Id atom) const {
    if (atom >= _type_table->paramCount())
        return "";
    else
        return _type_table->value(atom, "pset").asString();
}

bool TemplatedSystem::aromatic(Id bond) const {
    if (bond >= _arom_table->paramCount())
        return false;
    else
        return _arom_table->value(bond, "aromatic").asInt();
}

void TemplatedSystem::setTypes(Id atom, const std::string& btype,
        const std::string& nbtype, const std::string& pset) {
    updateSystem();
    _type_table->value(atom, "btype") = btype;
    _type_table->value(atom, "nbtype") = nbtype;
    _type_table->value(atom, "pset") = pset;
}

void TemplatedSystem::setAromatic(Id bond, bool arom) {
    updateSystem();
    _arom_table->value(bond, "aromatic") = arom;
}

void TemplatedSystem::updateSystem() {
    if (_max_atom_id != _type_table->paramCount()) {
        std::stringstream msg;
        msg << "Inconsistent system: max atom ID = " << _max_atom_id
            << ", type table nrows = " << _type_table->paramCount();
        VIPARR_FAIL(msg.str());
    }
    if (_max_bond_id != _arom_table->paramCount()) {
        std::stringstream msg;
        msg << "Inconsistent system: max bond ID = " << _max_bond_id
            << ", aromaticity table nrows = " << _arom_table->paramCount();
        VIPARR_FAIL(msg.str());
    }
    if (_atom_count != _sys->atomCount()
            || _max_atom_id != _sys->maxAtomId()
            || _bond_count != _sys->bondCount()
            || _max_bond_id != _sys->maxBondId()) {
        /* System topology has changed: reset _graph and _hash, and add rows
         * to _type_table and _arom_table if necessary */
        _hash = "";
        _graph.reset();
        for (unsigned i = _max_atom_id; i < _sys->maxAtomId(); ++i) {
            _type_table->addParam();
        }
        for (unsigned i = _max_bond_id; i < _sys->maxBondId(); ++i)
            _arom_table->addParam();
        _atom_count = _sys->atomCount();
        _max_atom_id = _sys->maxAtomId();
        _bond_count = _sys->bondCount();
        _max_bond_id = _sys->maxBondId();
    }
}

void TemplatedSystem::addTypedAtom(Id atom) {
    _typed_atoms.push_back(msys::IdList(1, atom));
}

void TemplatedSystem::addNonPseudoBond(const IdList& atoms) {
    assert(atoms.size() == 2);
    _non_pseudo_bonds.push_back(atoms);
}

void TemplatedSystem::addPseudoBond(const IdList& atoms) {
    assert(atoms.size() == 2);
    _pseudo_bonds.push_back(atoms);
}

void TemplatedSystem::addAngle(const IdList& atoms) {
    assert(atoms.size() == 3);
    _angles.push_back(atoms);
}

void TemplatedSystem::addDihedral(const IdList& atoms) {
    assert(atoms.size() == 4);
    _dihedrals.push_back(atoms);
}

void TemplatedSystem::addExclusion(const IdList& atoms) {
    assert(atoms.size() == 2);
    _exclusions.push_back(atoms);
}

void TemplatedSystem::addImproper(const IdList& atoms) {
    assert(atoms.size() == 4);
    _impropers.push_back(atoms);
}

void TemplatedSystem::addCmap(const IdList& atoms) {
    assert(atoms.size() == 8);
    _cmaps.push_back(atoms);
}

unsigned TemplatedSystem::addPseudoType(const std::string& type,
        unsigned nsites) {

    for (unsigned i = 0; i < _pseudo_types.size(); ++i) {
        if (_pseudo_types[i].name == type) { // Type already exists
            if (_pseudo_types[i].nsites != nsites) {
                std::stringstream msg;
                msg << "Pseudo type " << type << " already exists with "
                    << _pseudo_types[i].nsites << " sites";
                VIPARR_FAIL(msg.str());
            } else
                return i;
        }
    }
    /* Add new type */
    PseudoType p = {type, nsites, TupleList()};
    _pseudo_types.push_back(p);
    return _pseudo_types.size() - 1;
}

void TemplatedSystem::addPseudoSites(const std::string& type,
        const IdList& atoms) {

    unsigned i = addPseudoType(type, atoms.size());
    _pseudo_types[i].sites_list.push_back(atoms);
}

void TemplatedSystem::removeTypedAtom(Id atom) {
    auto it = std::find(_typed_atoms.begin(),
			_typed_atoms.end(), msys::IdList(1, atom));
    if (it == _typed_atoms.end()) {
        std::string s = "{" + std::to_string(atom) + "}";
	VIPARR_FAIL("TemplatedSystem::removeTypedAtom: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_typed_atoms.erase(it);
    }
}

void TemplatedSystem::removeNonPseudoBond(const IdList& atoms) {
    assert(atoms.size() == 2);
    auto it = std::find(_non_pseudo_bonds.begin(), 
			_non_pseudo_bonds.end(), atoms);
    if (it == _non_pseudo_bonds.end()) {
        std::string s = "{";
	s += std::to_string(atoms[0]) + ", ";
	s += std::to_string(atoms[1]) + "}";
	VIPARR_FAIL("TemplatedSystem::removeNonPseudoBond: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_non_pseudo_bonds.erase(it);
    }
}

void TemplatedSystem::removePseudoBond(const IdList& atoms) {
    assert(atoms.size() == 2);
    auto it = std::find(_pseudo_bonds.begin(), _pseudo_bonds.end(), atoms);
    if (it == _pseudo_bonds.end()) {
        std::string s = "{";
	s += std::to_string(atoms[0]) + ", ";
	s += std::to_string(atoms[1]) + "} ";
	VIPARR_FAIL("TemplatedSystem::removePseudoBonds: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_pseudo_bonds.erase(it);
    }
}

void TemplatedSystem::removeAngle(const IdList& atoms) {
    assert(atoms.size() == 3);
    auto it = std::find(_angles.begin(), _angles.end(), atoms);
    if (it == _angles.end()) {
        std::string s = "{";
	for (unsigned int i = 0; i < 2; i++)
	    s += std::to_string(atoms[i]) + ", ";
	s += std::to_string(atoms[2]) + "}";
	VIPARR_FAIL("TemplatedSystem::removeAngle: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_angles.erase(it);
    }
}

void TemplatedSystem::removeDihedral(const IdList& atoms) {
    assert(atoms.size() == 4);
    auto it = std::find(_dihedrals.begin(), _dihedrals.end(), atoms);
    if (it == _dihedrals.end()) {
        std::string s = "{";
	for (unsigned int i = 0; i < 3; i++)
	    s += std::to_string(atoms[i]) + ", ";
	s += std::to_string(atoms[3]) + "}";
	VIPARR_FAIL("TemplatedSystem::removeDihedral: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_dihedrals.erase(it);
    }
}

void TemplatedSystem::removeExclusion(const IdList& atoms) {
    assert(atoms.size() == 2);
    auto it = std::find(_exclusions.begin(), _exclusions.end(), atoms);
    if (it == _exclusions.end()) {
 	IdList reversed;
        reversed.push_back(atoms.back());
	reversed.push_back(atoms.front());
        it = std::find(_exclusions.begin(), _exclusions.end(), reversed);
    }
    if (it == _exclusions.end()) {
        std::string s = "{";
	s += std::to_string(atoms[0]) + ", ";
	s += std::to_string(atoms[1]) + "}";
	VIPARR_FAIL("TemplatedSystem::removeExclusion: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_exclusions.erase(it);
    }
}

void TemplatedSystem::removeImproper(const IdList& atoms) {
    assert(atoms.size() == 4);
    auto it = std::find(_impropers.begin(), _impropers.end(), atoms);
    if (it == _impropers.end()) {
        std::string s = "{";
	for (unsigned int i = 0; i < 3; i++)
	    s += std::to_string(atoms[i]) + ", ";
	s += std::to_string(atoms[3]) + "}";
	VIPARR_FAIL("TemplatedSystem::removeImproper: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_impropers.erase(it);
    }
}

void TemplatedSystem::removeCmap(const IdList& atoms) {
    assert(atoms.size() == 8);
    auto it = std::find(_cmaps.begin(), _cmaps.end(), atoms);
    if (it == _cmaps.end()) {
        std::string s = "{";
	for (unsigned int i = 0; i < 7; i++)
	    s += std::to_string(atoms[i]) + ", ";
	s += std::to_string(atoms[7]) + "}";
	VIPARR_FAIL("TemplatedSystem::removeCmap: " + s + 
		    " was not found in the TemplatedSystem.");
    }
    else {
	_cmaps.erase(it);
    }
}

