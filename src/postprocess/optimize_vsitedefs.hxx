#ifndef desres_viparr_optimize_vsitedef_hxx
#define desres_viparr_optimize_vsitedef_hxx

#include <msys/system.hxx>

namespace desres { namespace viparr {

    /* Trys to turn lcXn virtual site definintions into lcX by
     * looking for definitons that completly depend on constrained bonds
     */
    void OptimizeVsiteDefs(msys::SystemPtr sys, bool verbose=true);

}}

#endif
