#include <viparr/version.hxx> /* Auto-generated in SConscript */
#include "../importexport/import_ff.hxx"
#include "stdlib.h"
#include <cstdio>
#include <boost/program_options.hpp>
#include <string>
#include <set>
#include <msys/dms.hxx>
#include <msys/clone.hxx>

using namespace desres;
using namespace desres::viparr;
namespace po = boost::program_options;

static const char* helpmsg = 
"\n"
"Constructs DMS system files from templates. This will construct one DMS file\n"
"for each specified template and export it to `tplname.dms`. The DMS file will\n"
"not contain forcefield information.\n"
"\n";

static void show_help(const char* prog, po::options_description const& opts) {
    std::cerr << "Usage:\n";
    std::cerr << "  " << prog
              << " template_file [ template_name ... ]\n";
    std::cerr << helpmsg;
    std::cerr << opts << std::endl;
}

int main(int argc, char *argv[]) {

    /* Parse command line */
    po::options_description hidden_opts("Hidden options");
    hidden_opts.add_options()
        ("help", "display help message")
        ("tplfile", po::value<std::string>(), "template file")
        ("tplname", po::value<std::vector<std::string> >(), "names of templates to convert")
        ;
    po::options_description all_opts("All options");
    all_opts.add(hidden_opts);
    po::positional_options_description pos_opts;
    pos_opts.add("tplfile", 1).add("tplname", -1);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc,
                                          argv).options(all_opts).positional(pos_opts).run(), vm);
        po::notify(vm);
    } catch (po::error) {
        std::cerr << "Error: Unrecognized arguments\n";
        show_help(argv[0], all_opts);
        return 1;
    }
    if (vm.count("help")) {
        show_help(argv[0], all_opts);
        return 0;
    }
    std::string tplfile;
    std::set<std::string> tplname;
    if(vm.count("tplfile")) {
        tplfile = vm["tplfile"].as<std::string>();
        std::vector<std::string> tmp;
        if(vm.count("tplname")) tmp=vm["tplname"].as<std::vector<std::string> >();
        std::copy( tmp.begin(), tmp.end(), std::inserter( tplname, tplname.end() ) );
    } else {
        std::cerr << "Error: Missing arguments\n";
        show_help(argv[0], all_opts);
        return 1;
    }
    
    msys::Provenance prov = msys::Provenance::fromArgs(argc, argv);
    prov.version = "viparr/";
    prov.version += VIPARR_VERSION;

    std::vector<TemplatedSystemPtr> templates;
    templates = ImportTemplates(tplfile);

   /* Create Fake vdw Param Table/matcher based on nbtypes */
    std::map<std::string, msys::Id> vdw_idx;
    msys::ParamTablePtr param_table = msys::ParamTable::create();
    param_table->addProp("type", msys::StringType);

    /* Import templates and export dms files */
    std::set<std::string> exported;
    for (unsigned i = 0; i < templates.size(); ++i) {
        TemplatedSystemPtr tpl=templates[i];
        msys::SystemPtr sys = tpl->system();
        std::string name=sys->residue(0).name;
        if (tplname.size()==0 || tplname.count(name)){

            /* Fill in some information that is otherwise missing on direct export */
            msys::TermTablePtr table = 
                tpl->system()->addTable("nonbonded", 1, param_table);
            table->category = msys::NONBONDED;
            bool skipExport=false;
            for (msys::Id j=0; j<sys->atomCount(); j++) {
                std::string aname=sys->atom(j).name;
                std::string type=tpl->nbtype(j);
                if(aname.size()==0 || aname[0]=='$' || type.size()==0 ){
                    printf("Will not export %s due to funny atom name or empty atom type\n",name.c_str());
                    skipExport=true;
                    break;
                }
                std::map<std::string, msys::Id>::iterator iter=vdw_idx.lower_bound(type);
                if(iter ==vdw_idx.end() || vdw_idx.key_comp()(type,iter->first)){
                    msys::Id idx=vdw_idx.size();
                    iter=vdw_idx.insert(iter, std::make_pair(type,idx) );
                    msys::Id paramid = param_table->addParam();
                    assert(idx==paramid);
                    param_table->value(paramid,"type")=type;
                }
                msys::Id row=iter->second;
                msys::IdList term(1, j);
                table->addTerm(term, row);

                /* FIXME: Need accesses to msys's mass table */
                sys->atom(j).mass=0.0;
            }     
            if(skipExport) continue;
            std::string ofile=name+".dms";
            printf("Exporting template %s\n", name.c_str());
            sys = msys::Clone(sys, sys->atoms()); /* Remove unused parameters */
            msys::ExportDMS(sys, ofile, prov);
            exported.insert(name);
        }           
    }
    printf("Exported %lu templates to dms files\n",exported.size());
    
    return 0;
}
