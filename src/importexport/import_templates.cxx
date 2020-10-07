#include "import_ff.hxx"
#include <msys/fastjson/parse.hxx>
#include <sstream>

namespace dfj = desres::msys::fastjson;

using namespace desres;
using namespace desres::viparr;

namespace {

    typedef std::map<std::string, msys::Id> AtomMap;
    
    void extract_types(const dfj::Json& types, TemplatedSystemPtr sys,
            msys::Id atm, const std::string& pset="") {

        std::string btype;
        std::string nbtype;

        if (types.kind() != dfj::Json::Array
                || types.size() < 1
                || types.elem(0).kind() != dfj::Json::String
                || (types.size() == 2
                    && types.elem(1).kind() != dfj::Json::String)
                || types.size() > 2)
            VIPARR_FAIL("Type must be of the form [type] or [btype, nbtype]");

        /* Extract types */
        btype = types.elem(0).as_string();
        nbtype = (types.size() == 1 ? btype : types.elem(1).as_string());

        /* Check types */
        if (btype.size() == 0 || nbtype.size() == 0)
            VIPARR_FAIL("Type has 0 length");
        if (btype.find_first_of(" \t\r\n\v\f") != std::string::npos) {
            VIPARR_FAIL("Type '" + btype
                    + "' contains embedded white space");
        }
        if (nbtype.find_first_of(" \t\r\n\v\f") != std::string::npos) {
            VIPARR_FAIL("Type '" + nbtype
                    + "' contains embedded white space");
        }

        sys->setTypes(atm, btype, nbtype, pset);
    }

    void build_atoms(const dfj::Json& atoms, TemplatedSystemPtr sys,
            AtomMap& map) {
        if (!atoms)
            VIPARR_FAIL("No atoms in template");
        if (atoms.kind() != dfj::Json::Array)
            VIPARR_FAIL("atoms list must be of array type");
        sys->system()->addAtomProp("memo", msys::StringType);
        for (int i=0; i<atoms.size(); i++) {
            const dfj::Json& arr = atoms.elem(i);
            if ((arr.size() != 4 && arr.size() != 5)
                    || arr.elem(0).kind() != dfj::Json::String
                    || arr.elem(1).kind() != dfj::Json::Int
                    || arr.elem(2).kind() != dfj::Json::Float
                    || arr.elem(3).kind() != dfj::Json::Array
                    || (arr.size() == 5
                        && arr.elem(4).kind() != dfj::Json::String)) {
                std::stringstream msg;
                msg << "Atom " << i << " must be of the form"
                    << "[name, atomic_number, charge, [type(s)] (, memo)]";
                VIPARR_FAIL(msg.str());
            }
            AtomMap::const_iterator iter
                = map.find(arr.elem(0).as_string());
            if (iter != map.end()) {
                VIPARR_FAIL("Duplicate atom name '"
                        + std::string(arr.elem(0).as_string()) + "'");
            }
            msys::Id id = sys->system()->addAtom(0);
            sys->system()->atom(id).name = arr.elem(0).as_string();
            sys->system()->atom(id).atomic_number = arr.elem(1).as_int();
            sys->system()->atom(id).charge = arr.elem(2).as_float();
            try {
                extract_types(arr.elem(3), sys, id);
            } catch (std::exception& e) {
                std::stringstream msg;
                msg << "Atom " << i << " has misformatted type: "
                    << e.what();
                VIPARR_FAIL(msg.str());
            }
            if (arr.size() == 5)
                sys->system()->atomPropValue(id, "memo")
                    = arr.elem(4).as_string();
            map[sys->system()->atom(id).name] = id;
        }
    }

    /* The viparr files use ext atom names ($1, $2, etc) to indicate
     * that a template has a bond to an atom outside the template. When
     * such names are encountered, we create a ext atom in the template,
     * which in particular will have atomic number -1 but still be bonded
     * to a real template atom. We are strict about the interpretation of
     * the ext atom name, e.g., $1 in the bonds part of the template means
     * the same atom as $1 in cmap declarations. */
    void build_bonds(const dfj::Json& bonds, TemplatedSystemPtr sys,
            AtomMap& map) {
        if (!bonds) return;
        if (bonds.kind() != dfj::Json::Array)
            VIPARR_FAIL("bonds list must be of array type");
        for (int bi = 0; bi < bonds.size(); bi++) {
            const dfj::Json& arr = bonds.elem(bi);
            if (arr.size() != 2
                    || arr.elem(0).kind() != dfj::Json::String
                    || arr.elem(1).kind() != dfj::Json::String) {
                std::stringstream msg;
                msg << "Bond " << bi << " must be of the form "
                    "[atom_name_1, atom_name_2]";
                VIPARR_FAIL(msg.str());
            }
            const char * name1 = arr.elem(0).as_string();
            const char * name2 = arr.elem(1).as_string();
            AtomMap::iterator i = map.find(name1), j = map.find(name2),
                e = map.end();
            if (i != e && j != e) { // Internal bond
                sys->system()->addBond(i->second,j->second);
            } else if (i != e) { // j is external
                msys::Id ext = sys->system()->addAtom(0);
                sys->system()->atom(ext).name = name2;
                sys->system()->atom(ext).atomic_number = -1;
                sys->system()->addBond(i->second,ext);
                map[name2] = ext;
            } else if (j != e) { // i is external
                msys::Id ext = sys->system()->addAtom(0);
                sys->system()->atom(ext).name = name1;
                sys->system()->atom(ext).atomic_number = -1;
                sys->system()->addBond(j->second,ext);
                map[name1] = ext;
            }
        }
    }

    void build_tuples(const dfj::Json& js, const std::string& type,
            TemplatedSystemPtr sys, AtomMap& map ) {
        const dfj::Json& defs = js.get(type.c_str());
        if (!defs) return;
        if (defs.kind() != dfj::Json::Array)
            VIPARR_FAIL(type + " list must be of array type");
        for (int i = 0; i < defs.size(); i++) {
            const dfj::Json& arr = defs.elem(i);
            if (arr.kind() != dfj::Json::Array) {
                std::stringstream msg;
                msg << "Item " << i << " of " << type << " list must be of the "
                    "form [atom_name_1,...,atom_name_n]";
                VIPARR_FAIL(msg.str());
            }
            int N = arr.size();
            msys::IdList atoms(N);
            for (int j=0; j<N; j++) {
                if (arr.elem(j).kind() != dfj::Json::String) {
                    std::stringstream msg;
                    msg << "Item " << i << " of " << type << " list must be of "
                        "the form [atom_name_1,...,atom_name_n]";
                    VIPARR_FAIL(msg.str());
                }
                AtomMap::iterator iter=map.find(arr.elem(j).as_string());
                if (iter==map.end()) {
                    std::stringstream msg;
                    msg << "Item " << i << " of " << type
                        << " list includes atom '"
                        << arr.elem(j).as_string() << "' not in template";
                    VIPARR_FAIL(msg.str());
                }
                else
                    atoms[j] = iter->second;
            }
            if (type == "exclusions") {
                if (N != 2) {
                    std::stringstream msg;
                    msg << "Item " << i << " of " << type
                        << " list must have exactly 2 atoms";
                    VIPARR_FAIL(msg.str());
                }
                sys->addExclusion(atoms);
            }
            else if (type == "impropers") {
                if (N != 4) {
                    std::stringstream msg;
                    msg << "Item " << i << " of " << type
                        << " list must have exactly 4 atoms";
                    VIPARR_FAIL(msg.str());
                }
                sys->addImproper(atoms);
            }
            else if (type == "cmap") {
                if (N != 8) {
                    std::stringstream msg;
                    msg << "Item " << i << " of " << type
                        << " list must have exactly 8 atoms";
                    VIPARR_FAIL(msg.str());
                }
                sys->addCmap(atoms);
            }
        }
    }

    /* Format is
     * [ name, charge, [types], field, site1, site2, ..., siteN, pset ]
     * Pseudo name, charge, type, and pset info are stored as atom
     * properties. Site info is stored in the pseudo_sites table of
     * the template. Pseudos are bonded to the first site atom, and, for
     * virtual_midpoint, also to the second site atom. */
    void build_pseudos(const dfj::Json& defs, TemplatedSystemPtr sys,
            AtomMap& map) {
        if (!defs) return;
        if (defs.kind() != dfj::Json::Array)
            VIPARR_FAIL("pseudos list must be of array type");
        for (int i = 0; i < defs.size(); i++) {
            const dfj::Json& arr = defs.elem(i);
            if (arr.kind() != dfj::Json::Array
                    || arr.size() < 6
                    || arr.elem(0).kind() != dfj::Json::String
                    || arr.elem(1).kind() != dfj::Json::Float
                    || arr.elem(2).kind() != dfj::Json::Array
                    || arr.elem(3).kind() != dfj::Json::String
                    || arr.elem(arr.size() - 1).kind() != dfj::Json::String) {
                std::stringstream msg;
                msg << "Item " << i << " of pseudos list must be of the form "
                    << "[name,charge,[types],field,site_1,...,site_n,pset]";
                VIPARR_FAIL(msg.str());
            }
            msys::Id id = sys->system()->addAtom(0);
            sys->system()->atom(id).name = arr.elem(0).as_string();
            sys->system()->atom(id).charge = arr.elem(1).as_float();
            sys->system()->atom(id).atomic_number = 0;
            std::string pset = arr.elem(arr.size()-1).as_string();
            try {
                extract_types(arr.elem(2), sys, id, pset);
            } catch (std::exception& e) {
                std::stringstream msg;
                msg << "Pseudo " << i << " has misformatted type: "
                    << e.what();
                VIPARR_FAIL(msg.str());
            }
            int N = arr.size()-5;
            std::string pseudo_type = arr.elem(3).as_string();
            sys->addPseudoType(pseudo_type, N+1);
            msys::IdList atoms(N+1);
            atoms[0] = id;
            for (int j=0; j<N; j++) {
                if (arr.elem(j+4).kind() != dfj::Json::String) { 
                    std::stringstream msg;
                    msg << "Item " <<i<< " of pseudos list must be of the form "
                        << "[name,charge,[types],field,site_1,...,site_n,pset]";
                    VIPARR_FAIL(msg.str());
                }
                AtomMap::iterator iter=map.find(arr.elem(j+4).as_string());
                if (iter==map.end()) {
                    std::stringstream msg;
                    msg << "Virtual site " << arr.elem(j+4).as_string()
                        << " of pseudo " << i << " is not in template";
                    VIPARR_FAIL(msg.str());
                }
                else
                    atoms[j+1] = iter->second;
            }
            sys->addPseudoSites(pseudo_type, atoms);
            sys->system()->addBond(atoms[0], atoms[1]);
            if (pseudo_type == "virtual_midpoint")
                sys->system()->addBond(atoms[0], atoms[2]);
        }
    }

    TemplatedSystemPtr import_template(const dfj::Json& js) {
        TemplatedSystemPtr sys = TemplatedSystem::create();
        sys->system()->addResidue(sys->system()->addChain());

        /* Create a lookup table for atoms by name */
        AtomMap atommap;
        build_atoms(js.get("atoms"), sys, atommap);
        build_bonds(js.get("bonds"), sys, atommap);
        msys::MultiIdList fragments;
        unsigned num_frags = sys->system()->updateFragids(&fragments);
        if (num_frags > 1) {
            std::stringstream msg;
            msg << "Atoms are disconnected" << std::endl;
            for (unsigned i = 0; i < num_frags; ++i) {
                msg << "Fragment " << i << ": ";
                for (unsigned j = 0; j < fragments[i].size(); ++j)
                    msg << sys->system()->atom(fragments[i][j]).name << " ";
                msg << std::endl;
            }
            VIPARR_FAIL(msg.str());
        }
        build_tuples(js, "exclusions", sys, atommap);
        build_tuples(js, "impropers", sys, atommap);
        build_tuples(js, "cmap", sys, atommap);
        build_pseudos(js.get("pseudos"), sys, atommap);
        return sys;
    }
}

namespace desres { namespace viparr {

    std::vector<TemplatedSystemPtr> ImportTemplates(const std::string& path) {

        if (!fs::exists(path))
            VIPARR_FAIL("File not found");
        dfj::Json js;
        try {
            dfj::parse_json(path.c_str(), js);
        } catch (std::exception& e) {
            VIPARR_FAIL("Misformatted '" + path + "' file: " + e.what());
        }
        if (js.kind() != dfj::Json::Object)
            VIPARR_FAIL("Template file must be of object type");

        std::vector<TemplatedSystemPtr> templates;
        for (int i = 0; i < js.size(); i++) {
            try {
                TemplatedSystemPtr sys = import_template(js.elem(i));
                sys->system()->residue(0).name = js.key(i);
                templates.push_back(sys);
            } catch(std::runtime_error& e) {
                VIPARR_FAIL("Template " + std::string(js.key(i)) + ": "
                        + e.what());
            }
        }
        return templates;
    }

}}
