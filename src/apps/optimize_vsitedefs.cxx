#include "../base.hxx"
#include "../postprocess/optimize_vsitedefs.hxx"
#include <msys/dms.hxx>
#include <boost/program_options.hpp>
#include <set>
#include <string>
#include <vector>

using namespace desres;
using namespace desres::viparr;

int main(int argc, char *argv[]) {

    /* Parse command line */
    namespace po = boost::program_options;
    po::options_description available_opts("Available options");
    available_opts.add_options();
    po::options_description hidden_opts("Hidden options");
    hidden_opts.add_options()
        ("help", "display help message")
        ("input-dms", po::value<std::string>(), "dms input")
        ;
    po::options_description all_opts("All options");
    all_opts.add(available_opts).add(hidden_opts);
    po::positional_options_description pos_opts;
    pos_opts.add("input-dms", 1);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc,
                    argv).options(all_opts).positional(pos_opts).run(), vm);
        po::notify(vm);
    } catch (po::error) {
        VIPARR_ERR << "Error: Unrecognized arguments\n";
        VIPARR_ERR << "Usage:\n";
        VIPARR_ERR << "  " << argv[0] << " structure.dms [ options ]\n";
        VIPARR_ERR << available_opts << std::endl;
        return 1;
    }
    if (vm.count("help")) {
        std::cout << "Usage:\n";
        std::cout << "  " << argv[0] << " structure.dms [ options ]\n";
        std::cout << available_opts << std::endl;
        return 0;
    }
    std::string inputdms;
    if (vm.count("input-dms")) {
        inputdms = vm["input-dms"].as<std::string>();
    } else {
        VIPARR_ERR << "Error: Missing DMS structure file name\n";
        VIPARR_ERR << "Usage:\n";
        VIPARR_ERR << "  " << argv[0] << " structure.dms [ options ]\n";
        VIPARR_ERR << available_opts;
        return 1;
    }

    /* optimize vsite defs */
    msys::SystemPtr sys = msys::ImportDMS(inputdms, false);
    VIPARR_OUT << "Loaded " << sys->chainCount() << " chains, "
        << sys->residueCount() << " residues, "
        << sys->atomCount() << " atoms" << std::endl;
    VIPARR_OUT << "Optimizing Virtual Site Definitions" << std::endl;
    OptimizeVsiteDefs(sys);
    msys::ExportDMS(sys, inputdms, msys::Provenance::fromArgs(argc,argv));
} 
