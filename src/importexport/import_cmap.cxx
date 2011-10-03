#include "import_ff.hxx"
#include <msys/fastjson/parse.hxx>
#include <boost/filesystem.hpp>
#include <sstream>

namespace dfj = desres::msys::fastjson;
namespace bfs = boost::filesystem;

namespace desres { namespace viparr {
    /*
      DESRESCode#3431 Note that while cmaps live in a std::vector and
      are therefore 0-indexed, there are referred to in the
      torsion_torsion_cmap table 1-indexed. The conversion happens in
      Forcefield::cmapTable.
     */
    std::vector<msys::ParamTablePtr> ImportCmap(const std::string& path) {

        if (!bfs::exists(path))
            VIPARR_FAIL("File not found");
        dfj::Json js;
        try {
            dfj::parse_json(path.c_str(), js);
        } catch (std::exception& e) {
            VIPARR_FAIL("Misformatted '" + path + "' file: " + e.what());
        }
        if (js.kind() != dfj::Json::Array)
            VIPARR_FAIL("'" + path + "' must be of array type");

        std::vector<msys::ParamTablePtr> cmaps;
        for (int i = 0; i < js.size(); ++i) {
            const dfj::Json& cmap_js = js.elem(i);
            if (cmap_js.kind() == dfj::Json::Null) {
                cmaps.push_back(msys::ParamTablePtr());
                continue;
            }
            msys::ParamTablePtr cmap_table = msys::ParamTable::create();
            cmap_table->addProp("phi", msys::FloatType);
            cmap_table->addProp("psi", msys::FloatType);
            cmap_table->addProp("energy", msys::FloatType);
            if (cmap_js.kind() != dfj::Json::Array) {
                std::stringstream msg;
                msg << "cmap table " << i << " is not of array type";
                VIPARR_FAIL(msg.str());
            }
            for (int j = 0; j < cmap_js.size(); ++j) {
                const dfj::Json& elem = cmap_js.elem(j);
                if (elem.kind() != dfj::Json::Array || elem.size() != 3) {
                    std::stringstream msg;
                    msg << "Row " << j << " of cmap table " << i
                        << " must be a [phi, psi, energy] triple";
                    VIPARR_FAIL(msg.str());
                }
                msys::Id row = cmap_table->addParam();
                for (unsigned col = 0; col < 3; ++col)
                    cmap_table->value(row, col)
                        = elem.elem(col).as_float();
            }
            cmaps.push_back(cmap_table);
        }
        return cmaps;
    }

}}

