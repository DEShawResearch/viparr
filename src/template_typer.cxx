#include "base.hxx"
#include "template_typer.hxx"
#include "parameter_matcher.hxx"
#include "util/get_bonds_angles_dihedrals.hxx"
#include <msys/elements.hxx>
#include <stdio.h>
#include <sstream>
#include <msys/MsysThreeRoe.hpp>

using namespace desres::msys;

typedef std::vector<IdList> TupleList;

namespace desres { namespace viparr {

    void TemplateTyper::addTemplate(TemplatedSystemPtr tpl) {
      _templates[tpl->hash()].push_back(tpl);
    }

    void TemplateTyper::delTemplate(TemplatedSystemPtr tpl) {
      TemplateMap::iterator iter = _templates.find(tpl->hash());
      if (iter == _templates.end())
        VIPARR_FAIL("Cannot delete template: template not found");
      for (std::vector<TemplatedSystemPtr>::iterator viter
             = iter->second.begin(); viter != iter->second.end(); ++viter) {
        if (*viter == tpl) {
          iter->second.erase(viter);
          return;
        }
      }
      VIPARR_FAIL("Cannot delete template: template not found");
    }

    std::vector<TemplatedSystemPtr> TemplateTyper::templates() const {
      std::vector<TemplatedSystemPtr> vec;
      for (TemplateMap::const_iterator iter = _templates.begin();
           iter != _templates.end(); ++iter) {
        vec.insert(vec.end(), iter->second.begin(), iter->second.end());
      }
      return vec;
    }

    TemplateList TemplateTyper::findTemplateByHash(const std::string&
                                                   hash) const {
      TemplateMap::const_iterator iter = _templates.find(hash);
      if (iter == _templates.end()) return TemplateList();
      return iter->second;
    }

    TemplateList TemplateTyper::findTemplateByName(const std::string&
                                                   name) const {
      TemplateList tpls;
      for (TemplateMap::const_iterator iter = _templates.begin(); iter !=
             _templates.end(); ++iter) {
        for (unsigned i = 0; i < iter->second.size(); ++i) {
          if (iter->second[i]->system()->residue(0).name == name)
            tpls.push_back(iter->second[i]);
        }
      }
      return tpls;
    }

    bool TemplateTyper::matchFragment(TemplatedSystemPtr sys,
                                      const IdList& fragment, std::vector<std::pair<TemplatedSystemPtr,
                                      msys::IdList> >& matches, std::ostream& why_not) {

      /* Partition the system into residues and split fragment based
       * on residue */
      std::map<Id, IdList> residues;
      for (unsigned i=0; i<fragment.size(); i++) {
        Id resid = sys->system()->atom(fragment[i]).residue;
        std::map<Id, IdList>::iterator iter
          = residues.insert(std::make_pair(resid, IdList())).first;
        iter->second.push_back(fragment[i]);
      }

      /* Save template matches */
      matches.clear();
      matches.reserve(residues.size());

      /* Find a template for each residue. Note that since bonds don't
       * straddle forcefields, if any template matches, then the entire
       * fragment must be processed by this forcefield. */
      bool matched = true;
      for (std::map<Id, IdList>::const_iterator r=residues.begin();
           r!=residues.end(); ++r) {
        const IdList& atoms = r->second;
        std::string resname = sys->system()->residue(r->first).name;
        Id resid = r->first;
        //if (atoms.size() < sys->system()->atomCountForResidue(resid)) {
        //    VIPARR_ERR << "WARNING: Residue " << resid << ": " << resname
        //        << " straddles fragments; will match as separate residues"
        //        << std::endl;
        //}

        IdList tmap;
        TemplatedSystemPtr tpl = findMatch(sys, atoms, resname,
                                           resid, tmap, why_not);
        if (tpl == TemplatedSystemPtr()) {
          matched = false;
          break;
        } else
          matches.push_back(std::make_pair(tpl, tmap));
      }
      return matched;
    }

    std::string const& TemplateTyper::get_formula(
      msys::SystemPtr sys,
      const IdList& atoms) const {

      uint64_t h = ThreeRoe::Hash64(&atoms[0], atoms.size()*sizeof(atoms[0]));
      std::pair<FormulaMap::iterator,bool> r = _formula_cache.insert(
        FormulaMap::value_type(FormulaKey(sys, h),""));
      if (r.second) {
        int counts[128];
        memset(counts,0,sizeof(counts));
        for (Id i : atoms) ++counts[sys->atom(i).atomic_number];
        std::stringstream formula;
        for (unsigned j=1; j<128; ++j) {
          if (!counts[j]) continue;
          formula << msys::AbbreviationForElement(j);
          if (counts[j]>1)
            formula << counts[j];
        }
        r.first->second = formula.str();
      }
      return r.first->second;
    }

    /* Find template match for atoms in sys. Returns the matched template
     * and stores the matching in tmap, where atom i in tpl corresponds to 
     * atom tmap[i] in sys. If no match is found, outputs the error message
     * to why_not. */
    TemplatedSystemPtr TemplateTyper::findMatch(TemplatedSystemPtr sys,
                                                const IdList& atoms, const std::string& resname, Id resid,
                                                IdList& tmap, std::ostream& why_not) const {

      TemplateList candidates = findTemplateByHash(msys::Graph::hash(
                                                     sys->system(), atoms));
      if (!candidates.size()) { // No hash match
        std::string formula = get_formula(sys->system(), atoms);
          
        why_not << "has no template with matching formula for residue " 
                << resid << " (" << resname << ", " << formula << ")";
        for (TemplateMap::const_iterator iter = _templates.begin();
             iter != _templates.end(); ++iter)
          for (unsigned i = 0; i < iter->second.size(); ++i) {
            if (iter->second[i]->system()->residue(0).name == resname)
              why_not << "\n\ta template with name " << resname 
                      << " was found but has different "
                      << "chemical formula and/or terminal locations";

            if(get_formula(iter->second[i]->system(), iter->second[i]->system()->atoms()) == formula) {
              why_not << "\n\ta template (" << iter->second[i]->system()->residue(0).name;
              why_not << ") with same chemical formula but different terminal extensions was found. ";
              why_not << "\n\tDid you remember to include or exclude connections to external residues?";
            }
          }
        why_not << "." << std::endl;
        return TemplatedSystemPtr();
      }
      msys::GraphPtr target = msys::Graph::create(sys->system(), atoms);
      TemplatedSystemPtr tpl;
      std::vector<std::pair<Id, Id> > perm;
      for (TemplatedSystemPtr t : candidates) {
        std::vector<std::pair<Id, Id> > p;
        if (t->graph()->match(target, p)) {
          if (tpl == TemplatedSystemPtr()) {
            tpl = t;
            perm = p;
          } else {
            std::stringstream msg;
            msg << "Multiple templates " << t->system()->residue(0).name
                << " and " << tpl->system()->residue(0).name
                << " from a single forcefield match residue "
                << resid << " (" << resname << ")"  << std::endl;
            VIPARR_FAIL(msg.str());
          }
        }
      }
      if (tpl == TemplatedSystemPtr()) { // No graph match
        why_not << "has no template with matching topology for "
                << "residue " << resid << " (" << resname
                << "), but templates found with matching formula and "
                << "different bond topology:";
        for (unsigned i = 0; i < candidates.size(); ++i)
          why_not << " " << candidates[i]->system()->residue(0).name;
        why_not << "." << std::endl;
        return tpl;
      }
      /* Found match */
      tmap.clear();
      tmap.resize(tpl->system()->atomCount(), BadId);
      for (unsigned j=0; j<perm.size(); j++)
        tmap[perm[j].first] = perm[j].second;
      return tpl;
    }

    void TemplateTyper::assignMatch(TemplatedSystemPtr sys,
                                    const std::vector<std::pair<TemplatedSystemPtr, IdList> >& matches,
                                    bool rename_atoms, bool rename_residues) const {

      IdList assigned_atoms;
      for (unsigned match_i = 0, n = matches.size(); match_i < n; ++match_i) {
        TemplatedSystemPtr tpl = matches[match_i].first;
        IdList tmap = matches[match_i].second;

        /* Map atom types and charges, and add typed atoms */
        for (unsigned i=0; i<tmap.size(); i++) {
          Id j=tmap.at(i);
          if (bad(j)) continue; // Pseudo and external atoms not mapped
          sys->setTypes(j, tpl->btype(i), tpl->nbtype(i));
          sys->system()->atom(j).charge = tpl->system()->atom(i).charge;
          if (rename_atoms) {
            sys->system()->atom(j).name = tpl->system()->atom(i).name;
          }
          if (rename_residues) {
            Id res = sys->system()->atom(j).residue;
            sys->system()->residue(res).name
              = tpl->system()->residue(0).name;
          }
          sys->addTypedAtom(j);
          assigned_atoms.push_back(j);
        }

        /* Validate bonds, add aromaticity for internal bonds, add external
         * bonded atoms to tmap */
        std::set<Id> ambiguous_externals;
        for (unsigned i=0; i<tpl->system()->bondCount(); i++) {
          bond_t const& tbond = tpl->system()->bond(i);
          if (tpl->system()->atom(tbond.i).atomic_number == 0 || 
              tpl->system()->atom(tbond.j).atomic_number == 0) {
            /* Pseudo bond, ignore for now */
          }
          else if (tpl->system()->atom(tbond.i).atomic_number>0 &&
                   tpl->system()->atom(tbond.j).atomic_number>0) {
            /* Case of internal bond */
            Id ai = tmap[tbond.i];
            Id aj = tmap[tbond.j];
            Id bond = sys->system()->findBond(ai,aj);
            if (bad(bond))
              VIPARR_FAIL("Incorrect match; system is missing bond");
            /* Copy aromaticity */
            sys->setAromatic(bond, tpl->aromatic(i));
          } else {
            /* Case of external bond */
            Id t_in = tbond.i;
            Id t_ex = tbond.j;
            if (tpl->system()->atom(t_in).atomic_number == -1)
              std::swap(t_in,t_ex);

            Id s_in = tmap[t_in];
            Id s_ex = BadId;
            IdList const& s_bonds = sys->system()->bondsForAtom(s_in);
            for (unsigned j=0; j<s_bonds.size(); j++) {
              Id tmp = sys->system()->bond(s_bonds[j]).other(s_in);
              if (std::find(tmap.begin(), tmap.end(), tmp)
                  != tmap.end()) {
                continue;
              }
              if (s_ex != BadId)
                ambiguous_externals.insert(t_ex);
              s_ex = tmp;
            }
            if (bad(sys->system()->findBond(s_in, s_ex)))
              VIPARR_FAIL("Incorrect match; system is missing bond");
            tmap.at(t_ex) = s_ex;
          }
        }

        /* Add impropers, exclusions, cmaps */
        const TupleList& exclusions = tpl->exclusions();
        for (TupleList::const_iterator iter = exclusions.begin();
             iter != exclusions.end(); ++iter) {
          IdList tuple = *iter;
          for (unsigned j = 0; j < 2; ++j) {
            if (ambiguous_externals.find(tuple[j])
                != ambiguous_externals.end())
              VIPARR_FAIL("Exclusion in template "
                          << tpl->system()->residue(0).name
                          << " references an ambiguous externally "
                          "bonded atom");
            tuple[j] = tmap[tuple[j]];
          }
          sys->addExclusion(tuple);
        }
        const TupleList& impropers = tpl->impropers();
        for (TupleList::const_iterator iter = impropers.begin();
             iter != impropers.end(); ++iter) {
          IdList tuple = *iter;
          for (unsigned j = 0; j < 4; ++j) {
            if (ambiguous_externals.find(tuple[j])
                != ambiguous_externals.end())
              VIPARR_FAIL("Improper in template "
                          << tpl->system()->residue(0).name
                          << " references an ambiguous externally "
                          "bonded atom");
            tuple[j] = tmap[tuple[j]];
          }
          sys->addImproper(tuple);
        }
        const TupleList& cmaps = tpl->cmaps();
        for (TupleList::const_iterator iter = cmaps.begin();
             iter != cmaps.end(); ++iter) {
          IdList tuple = *iter;
          for (unsigned j = 0; j < 8; ++j) {
            if (ambiguous_externals.find(tuple[j])
                != ambiguous_externals.end())
              VIPARR_FAIL("Cmap in template "
                          << tpl->system()->residue(0).name
                          << " references an ambiguous externally "
                          "bonded atom");
            tuple[j] = tmap[tuple[j]];
          }
          sys->addCmap(tuple);
        }

        /* Reorder pseudo types with drudes at the end, to support drudes
         * attached to virtuals */
        std::vector<TemplatedSystem::PseudoType> pseudo_types;
        for (unsigned i = 0; i < tpl->pseudoTypes().size(); ++i)
          if (tpl->pseudoTypes()[i].name.substr(0,5) != "drude")
            pseudo_types.push_back(tpl->pseudoTypes()[i]);
        for (unsigned i = 0; i < tpl->pseudoTypes().size(); ++i)
          if (tpl->pseudoTypes()[i].name.substr(0,5) == "drude")
            pseudo_types.push_back(tpl->pseudoTypes()[i]);

        /* Add pseudo atoms, pseudo-sites, and pseudo-bonds */
        for (unsigned ind = 0; ind < pseudo_types.size(); ++ind) {
          const TupleList& sites_list = pseudo_types[ind].sites_list;
          for (unsigned i = 0; i < sites_list.size(); ++i) {
            IdList tuple = sites_list[i];
            for (unsigned j = 1; j < tuple.size(); ++j) {
              if (ambiguous_externals.find(tuple[j])
                  != ambiguous_externals.end())
                VIPARR_FAIL("Virtual site definition in template "
                            << tpl->system()->residue(0).name
                            << " references an ambiguous externally "
                            "bonded atom");
              if (tmap[tuple[j]] == msys::BadId)
                VIPARR_FAIL("VIPARR BUG: pseudo has a site atom "
                            "not yet processed");
              tuple[j] = tmap[tuple[j]];
            }
            Id tid = tuple[0];
            atom_t const& tpseudo = tpl->system()->atom(tid);
            Id parent = tuple[1];
            Id id = tuple[0] = sys->system()->addAtom(
              sys->system()->atom(parent).residue);
            tmap[tid] = id;
            sys->addPseudoSites(pseudo_types[ind].name, tuple);

            atom_t& pseudo = sys->system()->atom(id);
            pseudo.name=tpseudo.name;
            pseudo.x = sys->system()->atom(parent).x;
            pseudo.y = sys->system()->atom(parent).y;
            pseudo.z = sys->system()->atom(parent).z;
            pseudo.charge = tpseudo.charge;
            pseudo.formal_charge = 0;
            //pseudo.resonant_charge  = 0;
            sys->setTypes(id, tpl->btype(tid), tpl->nbtype(tid), 
                          tpl->pset(tid));
            sys->addTypedAtom(id);
            assigned_atoms.push_back(id);
            /* Add pseudo bonds */
            const IdList& bonded = tpl->system()->bondedAtoms(tid);
            for (unsigned j = 0; j < bonded.size(); ++j) {
              if (tmap[bonded[j]] == msys::BadId) {
                if (tpl->system()->atom(bonded[j]).atomic_number
                    != 0)
                  VIPARR_FAIL("VIPARR BUG: pseudo bonded to an "
                              "atom not yet processed");
                /* Bonded to another pseudo not yet added; add this
                 * bond later */
                continue;
              }
              Id pbond = sys->system()->addBond(id, tmap[bonded[j]]);
              sys->system()->bond(pbond).order = 1;
              //sys->system()->bond(pbond).resonant_order = 1;
              sys->setAromatic(pbond, false);
            }
          }
        }
      }

      /* Add bonds, angles, and dihedrals lists */
      std::vector<IdList> non_pseudo_bonds;
      std::vector<IdList> pseudo_bonds;
      std::vector<IdList> angles;
      std::vector<IdList> dihedrals;
      GetBondsAnglesDihedrals(sys->system(), assigned_atoms,
                              non_pseudo_bonds, pseudo_bonds, angles, dihedrals);
      for (unsigned i = 0; i < non_pseudo_bonds.size(); ++i)
        sys->addNonPseudoBond(non_pseudo_bonds[i]);
      for (unsigned i = 0; i < pseudo_bonds.size(); ++i)
        sys->addPseudoBond(pseudo_bonds[i]);
      for (unsigned i = 0; i < angles.size(); ++i)
        sys->addAngle(angles[i]);
      for (unsigned i = 0; i < dihedrals.size(); ++i)
        sys->addDihedral(dihedrals[i]);
    }
  }
}
