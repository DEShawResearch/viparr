#include "import_ff.hxx"
#include "../append_params.hxx"
#include <msys/fastjson/parse.hxx>
#include <boost/tokenizer.hpp>
#include <sstream>

namespace dfj = desres::msys::fastjson;
using namespace desres;

std::list<msys::Id> desres::viparr::ImportParams(const std::string& name,
        const std::string& path, const std::list<msys::Id>& share_params,
        const std::string& nbfix_identifier) {

    if (!fs::exists(path))
        VIPARR_FAIL("File not found");
    dfj::Json js;
    try {
        dfj::parse_json(path.c_str(), js);
    } catch (std::exception& e) {
        VIPARR_FAIL("Misformatted '" + path + "' file: " + e.what());
    }
    if (js.kind() != dfj::Json::Array)
        VIPARR_FAIL("Param file must be of array type");
    if (js.size() == 0)
        return std::list<msys::Id>();

    const dfj::Json& first = js.elem(0);
    if (first.kind() != dfj::Json::Object)
        VIPARR_FAIL("Row 0 must be an object with fields 'type' and 'params'");
    const dfj::Json& type = first.get("type");
    const dfj::Json& params = first.get("params");
    if (!type.valid() || !params.valid())
        VIPARR_FAIL("Row 0 must be an object with fields 'type' and "
                "'params'");
    if (type.kind() != dfj::Json::Array
            && type.kind() != dfj::Json::String)
        VIPARR_FAIL("Field 'type' of row 0 must be of array type");
    if (params.kind() != dfj::Json::Object)
        VIPARR_FAIL("Field 'params' of row 0 must be of object type");

#if 0
    /* Get set of all modes, for Alex-style forcefields */
    std::set<std::string> modes;
    for (int i = 0; i < js.size(); ++i) {
        const dfj::Json& js_row = js.elem(i);
        const dfj::Json& js_mode = js_row.get("mode");
        if (js_mode.valid() && std::string(js_mode.as_string()) != "") {
            if (js_mode.kind() != dfj::Json::String) {
                std::stringstream msg;
                msg << "Field 'mode' of row "
                    << i << " must be of string type";
                VIPARR_FAIL(msg.str());
            }
            std::string mode = js_mode.as_string();
            if (mode.find_first_of(" \t\r\n\v\f") != std::string::npos) {
                std::stringstream msg;
                msg << "'mode' value in row " << i
                    << " cannot contain whitespaces";
                VIPARR_FAIL(msg.str());
            }
            modes.insert(mode);
        }
    }
#endif

    unsigned npatterns = 0;
    if (type.kind() == dfj::Json::Array)
        npatterns = type.size();
    else {
        boost::char_separator<char> sep(" ");
        std::string type_str(type.as_string());
        boost::tokenizer<boost::char_separator<char> >
            boost_tokens(type_str, sep);
        std::vector<std::string> tokens(boost_tokens.begin(),
                boost_tokens.end());
        npatterns = tokens.size();
    }
    if (name == "vdw1" && npatterns != 1)
        VIPARR_FAIL("vdw1 table must have exactly one type");
    if (name == "vdw2" && npatterns != 2)
        VIPARR_FAIL("vdw2 table must have exactly two types");
    unsigned nparams = params.size();
    if (name == "vdw2" && nbfix_identifier == "")
        VIPARR_FAIL("Cannot import vdw2 table with empty nbfix_identifier");

    /* Create new param table */
    msys::ParamTablePtr param_table = msys::ParamTable::create();
    param_table->addProp("type", msys::StringType);
    for (unsigned i = 0; i < nparams; ++i) {
        std::string key = params.key(i);
        /* Cannot have a parameter called "type" or "memo" */
        if (key == "type" || key == "memo" || key == "nbfix_identifier"
                || key.substr(0,5) == "mode_") {
            VIPARR_FAIL("Parameter name '" + key + "' is reserved, "
                    " please rename");
        }
        msys::ValueType type;
        switch (params.elem(i).kind()) {
            case dfj::Json::Float:
            case dfj::Json::Int:
                type = msys::FloatType; break;
            case dfj::Json::String: type = msys::StringType; break;
            default: VIPARR_FAIL("Parameter " + key
                             + " must be of float, int, or string type");
        }
        /* Store cmap values as "cmap1", "cmap2", etc. */
        if (key == "cmapid")
            type = msys::StringType;
        param_table->addProp(key, type);
    }
    /* If vdw1 or vdw2, add "nbfix_identifier" column */
    if (name == "vdw1" || name == "vdw2")
        param_table->addProp("nbfix_identifier", msys::StringType);

#if 0
    /* If there are modes, add "mode_..." columns */
    for (std::set<std::string>::iterator iter = modes.begin();
            iter != modes.end(); ++iter)
        param_table->addProp("mode_" + (*iter), msys::IntType);
#endif
    param_table->addProp("memo", msys::StringType);

    /* Add type patterns and params */
    for (unsigned i = 0, n = js.size(); i < n; ++i) {
        const dfj::Json& row = js.elem(i);
        if (row.kind() != dfj::Json::Object) {
            std::stringstream msg;
            msg << "Row " << i << " must be an object with fields 'type' "
                "and 'params'";
            VIPARR_FAIL(msg.str());
        }
        const dfj::Json& type = row.get("type");
        const dfj::Json& params = row.get("params");
        if (!type.valid() || !params.valid()) {
            std::stringstream msg;
            msg << "Row " << i << " must be an object with fields 'type' "
                "and 'params'";
            VIPARR_FAIL(msg.str());
        }
        /* Create concatenated type string */
        std::string type_string;
        const dfj::Json& mode = row.get("mode");
        if (mode.valid() && std::string(mode.as_string()) != "")
            type_string = "__mode__" + std::string(mode.as_string()) + " ";
        if (type.kind() == dfj::Json::Array) {
            std::vector<std::string> type_vec(npatterns);
            try {
                for (unsigned j = 0; j < npatterns; ++j)
                    type_vec[j] = type.elem(j).as_string();
            } catch (std::exception& e) {
                std::stringstream msg;
                msg << "'type' field of row " << i
                    << " must be an array of " << npatterns << " strings";
                VIPARR_FAIL(msg.str());
            }
            type_string += type_vec[0];
            if (type_vec[0].find_first_of(" \t\r\n\v\f") != std::string::npos) {
                std::stringstream msg;
                msg << "'type' strings in row " << i
                    << " cannot contain whitespaces";
                VIPARR_FAIL(msg.str());
            }
            for (unsigned j = 1; j < npatterns; ++j) {
                if (type_vec[j].find_first_of(" \t\r\n\v\f")
                        != std::string::npos) {
                    std::stringstream msg;
                    msg << "'type' strings in row " << i
                        << " cannot contain whitespaces";
                    VIPARR_FAIL(msg.str());
                }
                type_string += " " + type_vec[j];
            }
        } else {
            boost::char_separator<char> sep(" ");
            std::string type_str(type.as_string());
            boost::tokenizer<boost::char_separator<char> >
                boost_tokens(type_str, sep);
            std::vector<std::string> tokens(boost_tokens.begin(),
                    boost_tokens.end());
            if (tokens.size() != npatterns) {
                std::stringstream msg;
                msg << "'type' field of row " << i
                    << " must have " << npatterns << " patterns";
                VIPARR_FAIL(msg.str());
            }
            type_string += type_str;
        }
        /* Add rows to param table */
        msys::Id paramid = param_table->addParam();
        param_table->value(paramid, "type") = type_string;
        if (name == "vdw1" || name == "vdw2")
            param_table->value(paramid, "nbfix_identifier")
                = nbfix_identifier;
        try {
            const dfj::Json& memo = row.get("memo");
            param_table->value(paramid, "memo") = memo.as_string("");
        } catch (std::exception& e) {
            param_table->value(paramid, "memo") = "";
        }
        for (unsigned j = 0; j < nparams; ++j) {
            std::string key = params.key(j);
            msys::Id ind = param_table->propIndex(key);
            if (ind == msys::BadId) {
                std::stringstream msg;
                msg << "Unrecognized param '" << key << "' in row " << i;
                VIPARR_FAIL(msg.str());
            }
            msys::ValueRef ref = param_table->value(paramid, ind);
            if (param_table->propName(ind) == "cmapid") {
                std::stringstream cmap;
                cmap << "cmap" << params.elem(j).as_int();
                ref = cmap.str();
            } else {
                switch (param_table->propType(ind)) {
                    case msys::FloatType:
                        ref = params.elem(j).as_float(); break;
                    case msys::IntType:
                        ref = params.elem(j).as_int(); break;
                    case msys::StringType:
                        ref = params.elem(j).as_string(); break;
                    default: VIPARR_FAIL("Parameter " + key +
                                     " must be float, int, or string type");
                } 
            }
        }
    }

    /* Append new param table to existing table, if present */
    std::list<msys::Id> rows = append_params::AppendParams<
        std::list<msys::Id> >(param_table, name, share_params);

#if 0
    /* For each added type with modes, create a new row in the param table with
     * pointers to the mode rows */
    msys::ParamTablePtr ptable = Forcefield::ParamTable(name);
    std::map<std::string, msys::Id> type_map;
    std::list<msys::Id> new_rows;
    for (unsigned i = 0; i < rows.size(); ++i) {
        std::string type = ptable->value(rows[i], "type").asString();
        if (type.substr(0,8) != "__mode__") {
            new_rows.push_back(rows[i]);
            continue;
        }
        std::string raw_type = type.substr(type.find_first_of(" ") + 1);
        std::map<std::string, msys::Id>::iterator iter
            = type_map.find(raw_type);
        if (iter == type_map.end()) {
            msys::Id param = ptable->addParam();
            /* Set starting values of new row to those of the first encountered
             * mode */
            for (unsigned j = 0; j < ptable->propCount(); ++j) {
                if (ptable->propName(j) == "type")
                    ptable->value(param, j) = raw_type;
                else if (ptable->propName(j).substr(0,5) != "mode_")
                    ptable->value(param, j) = ptable->value(rows[i], j);
            }
            /* Add new row to returned rows */
            new_rows.push_back(param);
            iter = type_map.insert(std::make_pair(raw_type, param)).first;
        }
        std::string mode_name = type.substr(8, type.find_first_of(" ")-8);
        /* Set the 'mode_...' column to 1 + rowID of the mode row; we add 1
         * because if new rows are added to the param table without modes,
         * their 'mode_...' values are set by default to 0 */
        ptable->value(iter->second, "mode_" + mode_name) = rows[i] + 1;
    }
    return new_rows;
#endif
    return rows;
}
