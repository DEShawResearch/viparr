#ifndef desres_viparr_template_typer_hxx
#define desres_viparr_template_typer_hxx

#include "templated_system.hxx"
#include <vector>
#include <map>
#include <string>

namespace desres { namespace viparr {

    typedef std::vector<TemplatedSystemPtr> TemplateList;
    typedef std::map<std::string, TemplateList> TemplateMap;

    /* Helper class for Forcefield. Performs atom typing using templates
     * and applies template information to system. Templates are stored
     * as a hash map, with the hash function provided by msys::Graph. */
    class TemplateTyper {

            typedef std::pair<msys::SystemPtr, uint64_t> FormulaKey;
            typedef std::map<FormulaKey, std::string> FormulaMap;
            mutable FormulaMap _formula_cache;

            std::string const& get_formula(msys::SystemPtr sys,
                                           msys::IdList const& atoms) const;

        public:

            /* Create a new typer with no templates */
            static std::shared_ptr<TemplateTyper> create() {
                return std::make_shared<TemplateTyper>();
            }

            /* Compute the graph hash for a template and add to the hash
             * map */
            void addTemplate(TemplatedSystemPtr tpl);

            /* Remove a template from the hash map */
            void delTemplate(TemplatedSystemPtr tpl);

            TemplateList findTemplateByName(const std::string& name) const;

            /* Return list of templates for a given formula hash */
            TemplateList findTemplateByHash(const std::string& hash) const;

            /* Match templates to a given fragment of a given system. Templates
             * are matched separately to the separate residues of sys->system();
             * matches are based on isomorphism of the bond-graph. If all
             * residues are successfully matched, returns true and stores
             * matching with template-to-system Id map in matches. Otherwise,
             * returns false and prints the error message to why_not. */
            virtual bool matchFragment(TemplatedSystemPtr sys,
                    const msys::IdList& fragment, std::vector<std::pair<
                    TemplatedSystemPtr, msys::IdList> >& matches,
                    std::ostream& why_not);

            /* Assign atom type and other template information from matched 
             * templates to the given system. Permanent changes to sys->system()
             * are the mapping of atom charges and bond aromaticity (if
             * available) from the templates and the addition of pseudo
             * particles and bonds. Other changes to sys not reflected in
             * sys->system() are the addition of btype, nbtype, and pset,
             * and hierarchy information for atoms in the fragment, and the
             * addition of typed atoms, non-pseudo bonds, angles, dihedrals,
             * impropers, exclusions, cmaps, and pseudo-sites for the fragment.
             * */
            void assignMatch(TemplatedSystemPtr sys, const std::vector<
                    std::pair<TemplatedSystemPtr, msys::IdList> >& matches,
                    bool rename_atoms=false, bool rename_residues=false) const;

            /* Helper function used within TemplateTyper::matchFragment and
             * SmartsTyper::matchFragment to match a residue to a single
             * template */
            TemplatedSystemPtr findMatch(TemplatedSystemPtr sys,
                    const msys::IdList& atoms, const std::string& resname,
                    msys::Id resid, msys::IdList& tmap,
                    std::ostream& why_not) const;

            /* Helper function used within SmartsTyper::matchFragment to
             * construct a template for a fragment using SMARTS definitions */
            virtual TemplatedSystemPtr buildTemplate(TemplatedSystemPtr sys,
                   const msys::IdList& fragment, std::ostream& why_not) {
                   return TemplatedSystemPtr();
            };

            std::vector<TemplatedSystemPtr> templates() const;

            virtual ~TemplateTyper() { }

        protected:
            /* Hash map of templates */
            TemplateMap _templates;
    };
    typedef std::shared_ptr<TemplateTyper> TemplateTyperPtr;

}}

#endif
