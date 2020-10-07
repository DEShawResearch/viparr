#include "export_ff.hxx"
#include <msys/fastjson/fastjson.hxx>
#include <iostream>
#include <fstream>
#include <msys/elements.hxx>

namespace dfj = desres::msys::fastjson;

using namespace desres;
using namespace desres::viparr;

namespace {

    void export_template(TemplatedSystemPtr tpl, std::ofstream& out,
            std::string end) {
        dfj::Json jatoms;
        jatoms.to_array();
        dfj::Json tmp;
        msys::IdList atoms = tpl->system()->atoms();
        std::map<std::string, int> name_counts;
        std::vector<std::string> names(tpl->system()->maxAtomId());
        for (unsigned i = 0; i < atoms.size(); ++i) {
            const msys::atom_t& atom = tpl->system()->atom(atoms[i]);
            std::stringstream name;
            if(atom.name == ""){
                name << msys::AbbreviationForElement(atom.atomic_number);
            }else{
                name << atom.name;
            }
            if (name_counts.find(name.str()) == name_counts.end())
                name_counts[name.str()] = 0;
            else {
                ++name_counts[name.str()];
                name << "_" << name_counts[name.str()];
            }
            names[atoms[i]] = name.str();
            if (atom.atomic_number < 1)
                continue;
            dfj::Json jatom;
            jatom.to_array();
            jatom.append(tmp.to_string(name.str().c_str()));
            jatom.append(tmp.to_int(atom.atomic_number));
            jatom.append(tmp.to_float(atom.charge));
            dfj::Json jtypes;
            jtypes.to_array();
            jtypes.append(tmp.to_string(tpl->btype(atoms[i]).c_str()));
            jtypes.append(tmp.to_string(tpl->nbtype(atoms[i]).c_str()));
            jatom.append(jtypes);
            if (tpl->system()->atomPropIndex("memo") != msys::BadId) {
                std::string memo = tpl->system()->atomPropValue(
                        atoms[i], "memo").asString();
                if (memo != "")
                    jatom.append(tmp.to_string(memo.c_str()));
            }
            jatoms.append(jatom);
        }

        dfj::Json jbonds;
        jbonds.to_array();
        msys::IdList bonds = tpl->system()->bonds();
        for (unsigned i = 0; i < bonds.size(); ++i) {
            msys::Id ai = tpl->system()->bond(bonds[i]).i;
            msys::Id aj = tpl->system()->bond(bonds[i]).j;
            if (tpl->system()->atom(ai).atomic_number == 0
                    || tpl->system()->atom(aj).atomic_number == 0)
                continue;
            dfj::Json jbond;
            jbond.to_array();
            jbond.append(tmp.to_string(names[ai].c_str()));
            jbond.append(tmp.to_string(names[aj].c_str()));
            jbonds.append(jbond);
        }

        dfj::Json jexcls;
        jexcls.to_array();
        const std::vector<msys::IdList>& excls = tpl->exclusions();
        for (unsigned i = 0; i < excls.size(); ++i) {
            dfj::Json jexcl;
            jexcl.to_array();
            for (unsigned j = 0; j < excls[i].size(); ++j)
                jexcl.append(tmp.to_string(names[excls[i][j]].c_str()));
            jexcls.append(jexcl);
        }

        dfj::Json jimprs;
        jimprs.to_array();
        const std::vector<msys::IdList>& imprs = tpl->impropers();
        for (unsigned i = 0; i < imprs.size(); ++i) {
            dfj::Json jimpr;
            jimpr.to_array();
            for (unsigned j = 0; j < imprs[i].size(); ++j)
                jimpr.append(tmp.to_string(names[imprs[i][j]].c_str()));
            jimprs.append(jimpr);
        }

        dfj::Json jcmaps;
        jcmaps.to_array();
        const std::vector<msys::IdList>& cmaps = tpl->cmaps();
        for (unsigned i = 0; i < cmaps.size(); ++i) {
            dfj::Json jcmap;
            jcmap.to_array();
            for (unsigned j = 0; j < cmaps[i].size(); ++j)
                jcmap.append(tmp.to_string(names[cmaps[i][j]].c_str()));
            jcmaps.append(jcmap);
        }

        dfj::Json jpseudos;
        jpseudos.to_array();
        const std::vector<TemplatedSystem::PseudoType>&
            pseudoTypes = tpl->pseudoTypes();
        for (unsigned i = 0; i < pseudoTypes.size(); ++i) {
            for (unsigned j = 0; j < pseudoTypes[i].sites_list.size(); ++j) {
                const msys::IdList& sites = pseudoTypes[i].sites_list[j];
                const msys::atom_t& pseudo = tpl->system()->atom(sites[0]);
                dfj::Json jpseudo;
                jpseudo.to_array();
                jpseudo.append(tmp.to_string(names[sites[0]].c_str()));
                jpseudo.append(tmp.to_float(pseudo.charge));
                dfj::Json jtypes;
                jtypes.to_array();
                jtypes.append(tmp.to_string(tpl->btype(sites[0]).c_str()));
                jtypes.append(tmp.to_string(tpl->nbtype(sites[0]).c_str()));
                jpseudo.append(jtypes);
                jpseudo.append(tmp.to_string(pseudoTypes[i].name.c_str()));
                for (unsigned k = 1; k < sites.size(); ++k)
                    jpseudo.append(tmp.to_string(names[sites[k]].c_str()));
                jpseudo.append(tmp.to_string(tpl->pset(sites[0]).c_str()));
                jpseudos.append(jpseudo);
            }
        }

        /* Pretty printing */
        std::string ind = "   ";
        std::string name = tpl->system()->residue(0).name;
        out << ind << "\""  << name << "\": {\n";
        std::vector<dfj::Json*> jsons;
        std::vector<std::string> labels;
        std::vector<std::string> endls;
        if (jatoms.size() > 0) {
            jsons.push_back(&jatoms);
            labels.push_back("atoms");
            endls.push_back("],\n");
        }
        if (jbonds.size() > 0) {
            jsons.push_back(&jbonds);
            labels.push_back("bonds");
            endls.push_back("],\n");
        }
        if (jexcls.size() > 0) {
            jsons.push_back(&jexcls);
            labels.push_back("exclusions");
            endls.push_back("],\n");
        }
        if (jimprs.size() > 0) {
            jsons.push_back(&jimprs);
            labels.push_back("impropers");
            endls.push_back("],\n");
        }
        if (jcmaps.size() > 0) {
            jsons.push_back(&jcmaps);
            labels.push_back("cmap");
            endls.push_back("],\n");
        }
        if (jpseudos.size() > 0) {
            jsons.push_back(&jpseudos);
            labels.push_back("pseudos");
            endls.push_back("],\n");
        }
        if (endls.size() > 0)
            endls[endls.size() - 1] = "]\n";
        for (unsigned i = 0; i < jsons.size(); ++i) {
            out << ind << ind << "\"" << labels[i] << "\": [\n";
            for (int j = 0; j < jsons[i]->size(); ++j) {
                out << ind << ind << ind;
                dfj::print_json(out, jsons[i]->elem(j), "", " ");
                if (j == jsons[i]->size() - 1)
                    out << "\n";
                else
                    out << ",\n";
            }
            out << ind << ind << endls[i];
        }
        out << ind << end;
    }

}

namespace desres { namespace viparr {

    void ExportTemplates(const std::vector<TemplatedSystemPtr>& tpls,
            const std::string& path) {
        if (fs::exists(path))
            VIPARR_FAIL("File already exists; cannot overwrite");
        std::ofstream out(path.c_str());
        out << "{\n";
        for (unsigned i = 0; i < tpls.size(); ++i) {
            if (i == tpls.size() - 1)
                export_template(tpls[i], out, "}\n");
            else
                export_template(tpls[i], out, "},\n");
        }
        out << "}";
    }
}}
