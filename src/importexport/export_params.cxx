#include "../base.hxx"
#include "export_ff.hxx"
#include <msys/fastjson/fastjson.hxx>
#include <iostream>
#include <fstream>
#include <sstream>

namespace dfj = desres::msys::fastjson;

namespace desres { namespace viparr {

    void ExportParams(msys::ParamTablePtr table, const std::list<msys::Id>&
            rows, const std::string& path) {

        if (fs::exists(path))
            VIPARR_FAIL("File already exists; cannot overwrite");
        msys::Id typeCol=table->propIndex("type");
        if (msys::bad(typeCol))
            VIPARR_FAIL("Param table is missing 'type' column, "
                    "cannot be exported");

        /* Replace rows pointing to mode rows with the mode rows */
        msys::IdList mode_cols;
        for (unsigned i = 0; i < table->propCount(); ++i)
            if (table->propName(i).substr(0,5) == "mode_")
                mode_cols.push_back(i);
        msys::IdList new_rows;
        for (msys::Id row : rows) {
            bool has_modes = false;
            for (unsigned j = 0; j < mode_cols.size(); ++j) {
                msys::Id mode_idx = table->value(row, j).asInt();
                if (mode_idx > 0) {
                    /* Mode index is 1 + rowID of mode row */
                    if (table->value(mode_idx - 1, typeCol).asString(
                                ).substr(0,8) != "__mode__") {
                        std::stringstream msg;
                        msg << "Row " << mode_idx << " of param table"
                            << " should be of type '__mode__...'";
                        VIPARR_FAIL(msg.str());
                    }
                    new_rows.push_back(mode_idx - 1);
                    has_modes = true;
                }
            }
            if (!has_modes)
                new_rows.push_back(row);
        }

        std::vector<std::string> endl(new_rows.size(),",\n");
        if (new_rows.size() > 0)
            endl[new_rows.size()-1]="\n";

        std::ofstream out(path.c_str());
        out << "[\n";
        dfj::Json tmp;
        for (unsigned i = 0; i < new_rows.size(); ++i) {
            msys::Id row = new_rows[i];
            dfj::Json jrow;
            jrow.to_object();
            dfj::Json jmode; /* Default is invalid */
            std::string type = table->value(row,typeCol).asString();
            if (type.substr(0,8) == "__mode__") {
                /* If first token is '__mode__...', save the mode name */
                std::size_t space = type.find(" ");
                jmode.to_string(type.substr(8, space).c_str());
                type = type.substr(space + 1);
            }
            dfj::Json jtype;
            jtype.to_string(type.c_str());
            jrow.append("type", jtype);
            /* all columns that are not "type", "memo", or "nbfix_identifier"
             * are treated as params */
            dfj::Json jparams;
            jparams.to_object();
            for (msys::Id col = 0; col < table->propCount(); ++col) {
                std::string key = table->propName(col);
                if (key == "type" || key == "memo"
                        || key == "nbfix_identifier")
                    continue;
                msys::ValueRef val = table->value(row, col);
                switch (table->propType(col)) {
                    case msys::IntType:
                        jparams.append(key.c_str(),tmp.to_int(val.asInt()));
                        break;
                    case msys::FloatType:
                        jparams.append(key.c_str(),tmp.to_float(val.asFloat()));
                        break;
                    case msys::StringType:
                        if (key != "cmapid")
                            jparams.append(key.c_str(),
                                    tmp.to_string(val.asString().c_str()));
                        else {
                            int id = atoi(val.asString().c_str() + 4);
                            jparams.append(key.c_str(), tmp.to_int(id));
                        }
                        break;
                    default: VIPARR_FAIL("Bad parameter value type");
                } 
            }
            if (jparams.size() > 0)
                jrow.append("params", jparams);
            if (jmode.valid())
                jrow.append("mode", jmode);
            if (table->propIndex("memo") != msys::BadId) {
                dfj::Json jmemo;
                jmemo.to_string(table->value(row, "memo").asString().c_str());
                jrow.append("memo", jmemo);
            }
            out<<"   ";
            dfj::print_json(out, jrow, "", " ");
            out<<endl[i];
        }
        out << "]\n";
        out.close();
    }
}}
