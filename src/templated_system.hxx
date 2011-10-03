#ifndef desres_viparr_templated_system_hxx
#define desres_viparr_templated_system_hxx

#include <msys/graph.hxx>
#include <msys/system.hxx>

namespace desres { namespace viparr {

    /* A wrapper class for msys::System to hold template information, as a
     * vehicle for communication between the TemplateTyper and parameter
     * matcher plug-ins. Can hold atom btype, nbtype, and pset information,
     * bond aromaticity, and lists of atoms, bonds, angles, dihedrals,
     * exclusions, impropers, cmaps, and pseudo-sites needed for
     * parameter-matching. TemplatedSystems are also used to store templates;
     * it can hold the msys::Graph and msys::Graph::hash information for the
     * template to cache the msys::Graph computation. */
    class TemplatedSystem {

        typedef msys::Id Id;
        typedef msys::IdList IdList;
        typedef std::vector<IdList> TupleList;

        public:
            TemplatedSystem();
            TemplatedSystem(msys::SystemPtr sys);

            struct PseudoType {
                std::string name;
                unsigned nsites;
                TupleList sites_list;
            };

            /* Create empty system */
            static std::shared_ptr<TemplatedSystem> create();

            /* Wrap a given msys::System */
            static std::shared_ptr<TemplatedSystem> create(
                    msys::SystemPtr sys);

            /* Access the contained msys::System */
            msys::SystemPtr system() { return _sys; }

            /* Clone the contained msys::System and copy relevant items */
            std::shared_ptr<TemplatedSystem> clone(const IdList& atoms);

            /* On the first call, graph() and hash() compute the msys::Graph and
             * msys::Graph::hash for all atoms in the system. Values are stored
             * in _graph and _hash and are returned without computation
             * for subsequent calls. _graph and _hash are reset when atoms or
             * bonds are added to or deleted from the system. */
            msys::GraphPtr graph();
            const std::string& hash();

            /* Return or set atom type properties */
            std::string btype(Id atom) const;
            std::string nbtype(Id atom) const;
            std::string pset(Id atom) const;
            void setTypes(Id atom, const std::string& btype,
                    const std::string& nbtype, const std::string& pset="");

            /* Return or set bond aromaticity */
            bool aromatic(Id bond) const;
            void setAromatic(Id bond, bool arom);

            /* Add element to the corresponding list */
            void addTypedAtom(Id atom);
            void addNonPseudoBond(const IdList& atoms);
            void addPseudoBond(const IdList& atoms);
            void addAngle(const IdList& atoms);
            void addDihedral(const IdList& atoms);
            void addExclusion(const IdList& atoms);
            void addImproper(const IdList& atoms);
            void addCmap(const IdList& atoms);

            /* Add a new pseudo type with empty site-tuple list if type does
             * not already exist; returns the index in _pseudo_types of
             * the given type */
            unsigned addPseudoType(const std::string& type,
                    unsigned nsites);

            /* Add a site-tuple list to the given pseudo type */
            void addPseudoSites(const std::string& type,
                    const IdList& atoms);

            /* Remove element from the corresponding list */
	    void removeTypedAtom(Id atom);
            void removeNonPseudoBond(const IdList& atoms);
            void removePseudoBond(const IdList& atoms);
            void removeAngle(const IdList& atoms);
            void removeDihedral(const IdList& atoms);
            void removeExclusion(const IdList& atoms);
            void removeImproper(const IdList& atoms);
            void removeCmap(const IdList& atoms);

            /* Return corresponding list */
            const TupleList& typedAtoms() const { return _typed_atoms; }
            const TupleList& nonPseudoBonds() const {
                return _non_pseudo_bonds; }
            const TupleList& pseudoBonds() const { return _pseudo_bonds; }
            const TupleList& angles() const { return _angles; }
            const TupleList& dihedrals() const { return _dihedrals; }
            const TupleList& exclusions() const { return _exclusions; }
            const TupleList& impropers() const { return _impropers; }
            const TupleList& cmaps() const { return _cmaps; }
            const std::vector<PseudoType>& pseudoTypes() const {
                return _pseudo_types; }


        private:
            /* Internal system, available for public access through
             * system() method */
            msys::SystemPtr _sys;

            /* Keep track of atoms and bonds in the system, to add rows to
             * _type_table and update _graph and _hash when the system topology
             * has changed */
            Id _atom_count;
            Id _max_atom_id;
            Id _bond_count;
            Id _max_bond_id;
            void updateSystem();

            /* msys::Graph and msys::Graph::hash for all atoms in the system */
            msys::GraphPtr _graph;
            std::string _hash;

            /* Tables containing atom btypes, nbtypes, and psets and bond
             * aromaticity */
            msys::ParamTablePtr _type_table;
            msys::ParamTablePtr _arom_table;

            /* Lists of all typed atoms, non-pseudo bonds, pseudo bonds, angles,
             * dihedrals, exclusion pairs, impropers, and cmaps */
            TupleList _typed_atoms;
            TupleList _non_pseudo_bonds;
            TupleList _pseudo_bonds;
            TupleList _angles;
            TupleList _dihedrals;
            TupleList _exclusions;
            TupleList _impropers;
            TupleList _cmaps;

            /* List of pseudo types, the number of sites for each type
             * (including the pseudo particle itself), and a list of
             * site-tuples for each type, each site-tuple consisting of the
             * pseudo Id followed by the Ids of the other site atoms */
            std::vector<PseudoType> _pseudo_types;

            /* Mapping from atom IDs to row in nonbonded pairs param table,
             * for easy update of pairs params */
            std::map<std::pair<msys::Id, msys::Id>, msys::Id> _pair_to_row;
    };

    typedef std::shared_ptr<TemplatedSystem> TemplatedSystemPtr;
}}

#endif
