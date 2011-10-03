#include "wrap_obj.hxx"
#include "../src/templated_system.hxx"
#include "../src/base.hxx"

using namespace desres::viparr;
using desres::msys::Id;
using desres::msys::IdList;

namespace {

    TemplatedSystemPtr create1() {
        return TemplatedSystem::create();
    }

    TemplatedSystemPtr create2(desres::msys::SystemPtr sys) {
        return TemplatedSystem::create(sys);
    }

    TemplatedSystemPtr tsys_clone(TemplatedSystem& sys, const object& atoms) {
        list L(atoms);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        return sys.clone(ids);
    }

    void add_non_pseudo_bond(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addNonPseudoBond(ids);
    }

    void add_pseudo_bond(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addPseudoBond(ids);
    }

    void add_angle(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addAngle(ids);
    }

    void add_dihedral(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addDihedral(ids);
    }

    void add_exclusion(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addExclusion(ids);
    }

    void add_improper(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addImproper(ids);
    }

    void add_cmap(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addCmap(ids);
    }

    void add_pseudo_sites(TemplatedSystem& sys, const std::string& type, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.addPseudoSites(type, ids);
    }

    void remove_non_pseudo_bond(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.removeNonPseudoBond(ids);
    }

    void remove_pseudo_bond(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.removePseudoBond(ids);
    }

    void remove_angle(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.removeAngle(ids);
    }

    void remove_dihedral(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.removeDihedral(ids);
    }

    void remove_exclusion(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.removeExclusion(ids);
    }

    void remove_improper(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.removeImproper(ids);
    }

    void remove_cmap(TemplatedSystem& sys, const object& elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i = 0; i < ids.size(); ++i)
            ids[i] = extract<Id>(L[i]);
        sys.removeCmap(ids);
    }

    list typed_atoms(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& atoms = sys.typedAtoms();
        for (unsigned i = 0; i < atoms.size(); ++i) {
            list elem;
            elem.append(atoms[i][0]);
            L.append(elem);
        }
        return L;
    }

    list non_pseudo_bonds(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& bonds = sys.nonPseudoBonds();
        for (unsigned i = 0; i < bonds.size(); ++i) {
            list elem;
            elem.append(bonds[i][0]);
            elem.append(bonds[i][1]);
            L.append(elem);
        }
        return L;
    }

    list pseudo_bonds(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& bonds = sys.pseudoBonds();
        for (unsigned i = 0; i < bonds.size(); ++i) {
            list elem;
            elem.append(bonds[i][0]);
            elem.append(bonds[i][1]);
            L.append(elem);
        }
        return L;
    }

    list angles(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& angles = sys.angles();
        for (unsigned i = 0; i < angles.size(); ++i) {
            list elem;
            elem.append(angles[i][0]);
            elem.append(angles[i][1]);
            elem.append(angles[i][2]);
            L.append(elem);
        }
        return L;
    }

    list dihedrals(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& dihedrals = sys.dihedrals();
        for (unsigned i = 0; i < dihedrals.size(); ++i) {
            list elem;
            elem.append(dihedrals[i][0]);
            elem.append(dihedrals[i][1]);
            elem.append(dihedrals[i][2]);
            elem.append(dihedrals[i][3]);
            L.append(elem);
        }
        return L;
    }

    list impropers(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& impropers = sys.impropers();
        for (unsigned i = 0; i < impropers.size(); ++i) {
            list elem;
            elem.append(impropers[i][0]);
            elem.append(impropers[i][1]);
            elem.append(impropers[i][2]);
            elem.append(impropers[i][3]);
            L.append(elem);
        }
        return L;
    }

    list exclusions(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& exclusions = sys.exclusions();
        for (unsigned i = 0; i < exclusions.size(); ++i) {
            list elem;
            elem.append(exclusions[i][0]);
            elem.append(exclusions[i][1]);
            L.append(elem);
        }
        return L;
    }

    list cmaps(TemplatedSystem& sys) {
        list L;
        const std::vector<IdList>& cmaps = sys.cmaps();
        for (unsigned i = 0; i < cmaps.size(); ++i) {
            list elem;
            for (unsigned j = 0; j < 8; ++j)
                elem.append(cmaps[i][j]);
            L.append(elem);
        }
        return L;
    }

    list pseudo_types(TemplatedSystem& sys) {
        list L;
        const std::vector<TemplatedSystem::PseudoType>& types = sys.pseudoTypes();
        for (unsigned i = 0; i < types.size(); ++i)
            L.append(types[i]);
        return L;
    }

    list get_sites_list(TemplatedSystem::PseudoType& ptype) {
        list L;
        for (unsigned i = 0; i < ptype.sites_list.size(); ++i) {
            list elem;
            for (unsigned j = 0; j < ptype.sites_list[i].size(); ++j)
                elem.append(ptype.sites_list[i][j]);
            L.append(elem);
        }
        return L;
    }

    void set_sites_list(TemplatedSystem::PseudoType& ptype, const object& elems) {
        ptype.sites_list.clear();
        list L(elems);
        for (unsigned i = 0; i < len(L); ++i) {
            ptype.sites_list.push_back(IdList());
            list elem(L[i]);
            for (unsigned j = 0; j < len(elem); ++j)
                ptype.sites_list[i].push_back(extract<Id>(elem[j]));
        }
    }

}

namespace desres { namespace viparr {

    void export_templated_system() {

        scope templated_system_class(class_<TemplatedSystem, TemplatedSystemPtr>("TemplatedSystem", no_init)
            .def("__eq__", _eq<TemplatedSystemPtr>)
            .def("__ne__", _ne<TemplatedSystemPtr>)
            .def("__hash__", _hash<TemplatedSystemPtr>)
            .def("__init__", make_constructor(create1))
            .def("__init__", make_constructor(create2))
            .def("system", &TemplatedSystem::system)
            .def("clone", tsys_clone)
            .def("btype", &TemplatedSystem::btype)
            .def("nbtype", &TemplatedSystem::nbtype)
            .def("pset", &TemplatedSystem::pset)
            .def("setTypes", &TemplatedSystem::setTypes)
            .def("hash", &TemplatedSystem::hash, return_const())
            .def("graph", &TemplatedSystem::graph)
            .def("addTypedAtom", &TemplatedSystem::addTypedAtom)
            .def("addNonPseudoBond", add_non_pseudo_bond)
            .def("addPseudoBond", add_pseudo_bond)
            .def("addAngle", add_angle)
            .def("addDihedral", add_dihedral)
            .def("addExclusion", add_exclusion)
            .def("addImproper", add_improper)
            .def("addCmap", add_cmap)
            .def("addPseudoType", &TemplatedSystem::addPseudoType)
            .def("addPseudoSites", add_pseudo_sites)
            .def("removeTypedAtom", &TemplatedSystem::removeTypedAtom)
            .def("removeNonPseudoBond", remove_non_pseudo_bond)
            .def("removePseudoBond", remove_pseudo_bond)
            .def("removeAngle", remove_angle)
            .def("removeDihedral", remove_dihedral)
            .def("removeExclusion", remove_exclusion)
            .def("removeImproper", remove_improper)
            .def("removeCmap", remove_cmap)
            .def("typedAtoms", typed_atoms)
            .def("nonPseudoBonds", non_pseudo_bonds)
            .def("pseudoBonds", pseudo_bonds)
            .def("angles", angles)
            .def("dihedrals", dihedrals)
            .def("exclusions", exclusions)
            .def("impropers", impropers)
            .def("cmaps", cmaps)
            .def("pseudoTypes", pseudo_types)
        );

        class_<TemplatedSystem::PseudoType>("PseudoType")
            .def_readonly("name", &TemplatedSystem::PseudoType::name)
            .def_readonly("nsites", &TemplatedSystem::PseudoType::nsites)
            .add_property("sites_list", get_sites_list, set_sites_list)
            ;
    }

}}
