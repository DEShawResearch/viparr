#include "compile_plugins.hxx"
#include "../base.hxx"
#include "../ff.hxx"

#include <algorithm>
#include <string>

desres::msys::TermTablePtr desres::viparr::AddPairsTable(msys::SystemPtr sys) {
  auto iter = Rules::VDWFuncRegistry().find(sys->nonbonded_info.vdw_funct);
  if (iter == Rules::VDWFuncRegistry().end()) {
    std::string msg = "Unsupported VDW function '" + sys->nonbonded_info.vdw_funct
      + "'; make sure the VDW function was copied from the "
      "forcefield to the system and is in the VDWFuncRegistry. Available functions:";
    for(auto s : Rules::VDWFuncRegistry())
      msg += " " + s.first;
    VIPARR_FAIL(msg);
  }
  Rules::VDWFunc vdw_func = iter->second;
  std::string name = vdw_func.pair_table_name;

  /* DESRESCode#1810 */
  auto lowercase = [](std::string s) {
    std::string answer;
    answer.reserve(s.length());
    std::transform(s.begin(), s.end(), answer.begin(), ::tolower);
    return answer;
  };
  assert(lowercase(name) != "none");

  msys::ParamTablePtr pairs_params;
  if (!Forcefield::HasParamTable(name)) {
    pairs_params = msys::ParamTable::create();
    Forcefield::AddParamTable(name, pairs_params);
    for (unsigned i = 0; i < vdw_func.pair_param_names.size(); ++i)
      pairs_params->addProp(vdw_func.pair_param_names[i],
                            msys::FloatType);
    pairs_params->addProp("qij", msys::FloatType);
  } else
    pairs_params = Forcefield::ParamTable(name);
  pairs_params->addProp("type", msys::StringType);
  pairs_params->addProp("memo", msys::StringType);
  msys::TermTablePtr pairs = sys->addTable(name, 2, pairs_params);
  pairs->category = msys::BOND;
  return pairs;
}

void desres::viparr::CompilePlugins(msys::SystemPtr sys,
                                    std::set<std::string> plugins) {
  ApplyNBFix(sys);
  AddPairsTable(sys);

  for (std::string plugin_name : plugins) {
    std::map<std::string, Forcefield::PluginPtr>::const_iterator
      iter = Forcefield::PluginRegistry().find(plugin_name);
    if(iter == Forcefield::PluginRegistry().end())
      VIPARR_FAIL("Asked to compile unknown plugin " + plugin_name);
    iter->second->compile(sys);
  }
}

void desres::viparr::CleanupSystem(msys::SystemPtr sys) {

  /* Change VDW functional form name */
  if(Rules::VDWFuncRegistry().find(sys->nonbonded_info.vdw_funct) !=
     Rules::VDWFuncRegistry().end()) {
  
    sys->nonbonded_info.vdw_funct = Rules::VDWFuncRegistry().find(
      sys->nonbonded_info.vdw_funct)->second.vdw_table_name;
  }
  
    if (sys->nonbonded_info.vdw_funct == "none"
        || sys->nonbonded_info.vdw_rule == "none")
    VIPARR_FAIL("Cannot prepare system with 'none' vdw_funct or "
                "vdw_rule for export");

  /* Remove NO_CATEGORY tables and tables with no terms */
  std::vector<std::string> tables = sys->tableNames();
  for (unsigned i = 0; i < tables.size(); ++i) {
    if (sys->table(tables[i])->category == msys::NO_CATEGORY
        || sys->table(tables[i])->termCount() == 0)
      sys->delTable(tables[i]);
  }
}
