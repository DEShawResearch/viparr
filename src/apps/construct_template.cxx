#include "../base.hxx"
#include "../importexport/export_ff.hxx"
#include <msys/clone.hxx>
#include <msys/dms.hxx>
#include <msys/io.hxx>
#include <msys/atomsel.hxx>
#include <msys/elements.hxx>
#include <boost/program_options.hpp>

using namespace desres;
using namespace desres::viparr;

int main(int argc, char *argv[]) {

    /* Parse command line */
    namespace po = boost::program_options;
    po::options_description available_opts("Available options");
    available_opts.add_options()
        ("by-residue", "create a separate template for each residue")
        ;
    po::options_description hidden_opts("Hidden options");
    hidden_opts.add_options()
        ("help", "display help message")
        ("structure-file", po::value<std::string>(), "structure file")
        ("selection", po::value<std::string>(), "atomsel to templatize")
        ("out-template", po::value<std::string>(), "output template file")
        ;
    po::options_description all_opts("All options");
    all_opts.add(available_opts).add(hidden_opts);
    po::positional_options_description pos_opts;
    pos_opts.add("structure-file", 1).add("selection", 1).add("out-template", 1);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc,
                    argv).options(all_opts).positional(pos_opts).run(), vm);
        po::notify(vm);
    } catch (po::error) {
        std::cerr << "Error: Unrecognized arguments\n";
        std::cerr << "Usage:\n";
        std::cerr << "  " << argv[0]
            << " structure-file selection out-template [ options ]\n";
        std::cerr << available_opts << std::endl;
        return 1;
    }
    if (vm.count("help")) {
        std::cout << "Usage:\n";
        std::cout << "  " << argv[0]
            << " structure-file selection out-template [ options ]\n";
        std::cout << available_opts << std::endl;
        return 0;
    }
    std::string structure_file;
    std::string selection;
    std::string out_template;
    if (vm.count("structure-file") && vm.count("selection")
            && vm.count("out-template")) {
        structure_file = vm["structure-file"].as<std::string>();
        selection = vm["selection"].as<std::string>();
        out_template = vm["out-template"].as<std::string>();
    } else {
        std::cerr << "Error: Missing arguments\n";
        std::cerr << "Usage:\n";
        std::cerr << "  " << argv[0]
            << " structure-file selection out-template [ options ]\n";
        std::cerr << available_opts;
        return 1;
    }
    bool by_residue = vm.count("by-residue");

    /* Load system and get selection */
    msys::SystemPtr sys = msys::Load(structure_file, false);
    if (sys->atomCount() == 0)
        VIPARR_FAIL("Imported system is empty; make sure atomic numbers are "
                "present");

    // Cleanup if necessary
    for (msys::Id aid : sys->atoms()){
       msys::atom_t & atom=sys->atom(aid);
       if(atom.atomic_number==0 && atom.mass != 0.0){
          atom.atomic_number=msys::GuessAtomicNumber(atom.mass);
       }
    }

    msys::IdList atomsel = msys::Atomselect(sys, selection);
    std::vector<bool> in_selection(sys->maxAtomId(), false);
    for (unsigned i = 0; i < atomsel.size(); ++i)
        in_selection[atomsel[i]] = true;

    /* Split selected atoms into fragments and residues, get counts of
     * external bonds */
    std::vector<msys::IdList> fragments;
    unsigned nfrags = sys->updateFragids(&fragments);
    std::vector<msys::IdList> template_sels;
    std::vector<int> external_count(sys->maxAtomId(), 0);
    for (unsigned i = 0; i < nfrags; ++i) {
        if (by_residue) {
            std::map<msys::Id, msys::IdList> resmap;
            for (unsigned j = 0; j < fragments[i].size(); ++j) {
                if (in_selection[fragments[i][j]]) {
                    resmap[sys->atom(fragments[i][j]).residue].push_back(
                            fragments[i][j]);
                    msys::IdList bonded = sys->bondedAtoms(fragments[i][j]);
                    for (unsigned k = 0; k < bonded.size(); ++k) {
                        if (sys->atom(fragments[i][j]).residue
                                != sys->atom(bonded[k]).residue ||
                                !in_selection[bonded[k]])
                            external_count[fragments[i][j]] += 1;
                    }
                }
            }
            for (std::map<msys::Id, msys::IdList>::iterator iter
                    = resmap.begin(); iter != resmap.end(); ++iter)
                template_sels.push_back(iter->second);
        } else {
            msys::IdList tmp_sel;
            for (unsigned j = 0; j < fragments[i].size(); ++j) {
                if (in_selection[fragments[i][j]]) {
                    tmp_sel.push_back(fragments[i][j]);
                    msys::IdList bonded = sys->bondedAtoms(fragments[i][j]);
                    for (unsigned k = 0; k < bonded.size(); ++k)
                        if (!in_selection[bonded[k]])
                            external_count[fragments[i][j]] += 1;
                }
            }
            if (tmp_sel.size() > 0)
                template_sels.push_back(tmp_sel);
        }
    }

    /* Create templates */
    std::map<std::string, int> resname_count;
    TemplateTyperPtr typer = TemplateTyper::create();
    for (unsigned i = 0; i < template_sels.size(); ++i) {
        msys::SystemPtr clone = msys::Clone(sys, template_sels[i]);
        int external_ind = 1;
        /* Add externally bonded atoms */
        for (unsigned j = 0; j < template_sels[i].size(); ++j) {
            while (external_count[template_sels[i][j]] > 0) {
                msys::Id ext_atom = clone->addAtom(clone->atom(0).residue);
                clone->atom(ext_atom).atomic_number = -1;
                std::stringstream name;
                name << "$" << external_ind;
                clone->atom(ext_atom).name = name.str();
                clone->addBond(ext_atom, j);
                ++external_ind;
                external_count[template_sels[i][j]] -= 1;
            }
        }
        /* Rename duplicate atom and residue names */
        std::map<std::string, int> atomname_count;
        for(msys::Id j : clone->atoms()) {
            if(clone->atom(j).name == ""){
                clone->atom(j).name = msys::AbbreviationForElement(clone->atom(j).atomic_number);
            }
            if (atomname_count.find(clone->atom(j).name)
                    != atomname_count.end()) {
                atomname_count[clone->atom(j).name] += 1;
                std::stringstream name;
                name << clone->atom(j).name << "_"
                    << atomname_count[clone->atom(j).name];
                clone->atom(j).name = name.str();
            } else
                atomname_count[clone->atom(j).name] = 0;
        }
        if (resname_count.find(clone->residue(0).name) != resname_count.end()) {
            resname_count[clone->residue(0).name] += 1;
            std::stringstream name;
            name << clone->residue(0).name << "_"
                << resname_count[clone->residue(0).name];
            clone->residue(0).name = name.str();
        } else
            resname_count[clone->residue(0).name] = 0;

        TemplatedSystemPtr tpl = TemplatedSystem::create(clone);

        /* use real vdw types if present */
        msys::TermTablePtr vdw=clone->table("nonbonded");
        if (vdw == msys::TermTablePtr() || vdw->params()->propIndex("type") == msys::BadId){
            for(msys::Id j: clone->atoms()){
                tpl->setTypes(j, clone->atom(j).name, clone->atom(j).name, "");
            }
        }else{
            msys::Id typeIdx=vdw->params()->propIndex("type");
            for(msys::Id j: clone->atoms()){
                std::string vdwType=clone->atom(j).name;
                msys::IdList rows=vdw->findExact({j});
                if(rows.size()==0){
                    if(clone->atom(j).atomic_number == -1){
                        continue;
                    }
                    printf("Couldnt find a vdw row for atom %u in template %u of %lu. Original atom id= %u\n",j,i,template_sels.size(),template_sels[i][j]);
                    printf("Name=%s AtomicNumber= %d\n",vdwType.c_str(), clone->atom(j).atomic_number );
                    assert(false);
                }else if(rows.size()>1){
                    printf(" Found duplicate vdw rows for atom %u in template %u of %lu. Original atom id= %u. RowIds=",j,i,template_sels.size(),template_sels[i][j]);
                    for(msys::Id id: rows){
                        printf(" %u",id);
                    }
                    printf("\n");
                }
                vdwType=vdw->propValue(rows[0],typeIdx).asString();
                tpl->setTypes(j, vdwType, vdwType, "");
            }
        }

        /* add impropers if present */
        /* FIXME!! WILL NOT HANDLE IMPROPERS THAT STRADDLE RESIDUES IF --by-residue IS USED!!! */
        msys::TermTablePtr impropers=clone->table("improper_harm");
        if (impropers != msys::TermTablePtr()){
            for( msys::Id tid : impropers->terms()){
                tpl->addImproper(impropers->atoms(tid));
            }
        }

        msys::IdList tpl_to_sys;
        std::stringstream why_not;
        TemplatedSystemPtr foundTemplate = typer->findMatch(tpl, clone->atoms(), "", 0, tpl_to_sys, why_not);
        if (foundTemplate == TemplatedSystemPtr()){
            /* Residue does not match any existing templates; create new template */
            typer->addTemplate(tpl);
        }else{
            /* Residue matches existing template in structure; check that
             * all other template properties match */
            printf("FoundDuplicates: %s %s\n",tpl->system()->residue(0).name.c_str(),foundTemplate->system()->residue(0).name.c_str());
            //check_template(tsys, iter->second, rules, tpl, tpl_to_sys); 
        }
    }


    /* Export templates */
    ExportTemplates(typer->templates(), out_template);
}
