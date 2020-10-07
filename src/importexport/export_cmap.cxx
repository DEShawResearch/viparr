#include "export_ff.hxx"
#include <msys/fastjson/fastjson.hxx>

namespace dfj = desres::msys::fastjson;

namespace desres { namespace viparr {

    /*
      DESRESCode#3431 Note that while cmaps live in a std::vector and
      are therefore 0-indexed, there are referred to in the
      torsion_torsion_cmap table 1-indexed. The conversion happens in
      Forcefield::cmapTable.
     */
    void ExportCmap(const std::vector<msys::ParamTablePtr>& cmap_tables,
            const std::string& path) {
        if (fs::exists(path))
            VIPARR_FAIL("File already exists; cannot overwrite");
        dfj::Json jcmaps;
        jcmaps.to_array();
        for (unsigned i = 0; i < cmap_tables.size(); ++i) {
            dfj::Json jcmap;
            if (cmap_tables[i] == msys::ParamTablePtr()) {
                jcmap.to_null();
                jcmaps.append(jcmap);
                continue;
            }
            jcmap.to_array();
            unsigned nrows = cmap_tables[i]->paramCount();
            for (unsigned j = 0; j < nrows; ++j) {
                dfj::Json jrow;
                jrow.to_array();
                dfj::Json tmp;
                jrow.append(tmp.to_float(cmap_tables[i]->value(j,0).asFloat()));
                jrow.append(tmp.to_float(cmap_tables[i]->value(j,1).asFloat()));
                jrow.append(tmp.to_float(cmap_tables[i]->value(j,2).asFloat()));
                jcmap.append(jrow);
            }
            jcmaps.append(jcmap);
        }
        dfj::print_json(path.c_str(), jcmaps);
    }

}}


