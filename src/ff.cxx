#include "ff.hxx"
#include "base.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {
  /* Implementation of Plugin interface in C++ */
  class PluginC : public Forcefield::Plugin {
  public:
    PluginC(void (*c_match)(TemplatedSystemPtr, ForcefieldPtr),
            const std::vector<std::string>& prerequisites,
            void (*c_compile)(msys::SystemPtr),
            const std::vector<std::string>& dependencies)
      : _c_match(c_match), _c_compile(c_compile) {
      /*
      if (c_match == NULL)
        VIPARR_FAIL("Plugin match function cannot be NULL");
      */
      this->prerequisites = prerequisites;
      this->dependencies = dependencies;
    }
    virtual void match(TemplatedSystemPtr sys, ForcefieldPtr ff) const {
      if(_c_match == NULL) return;
      _c_match(sys, ff);
    }
    virtual void compile(msys::SystemPtr sys) const {
      if(_c_compile == NULL) return;
      _c_compile(sys);
    }
  private:
    void (*_c_match)(TemplatedSystemPtr, ForcefieldPtr);
    void (*_c_compile)(msys::SystemPtr);
  };
}

namespace desres { namespace viparr {

    typedef std::map<std::string, msys::ParamTablePtr> ParamTableMap;

    ParamTableMap Forcefield::SharedParamTables = ParamTableMap();
    std::list<msys::Id> Forcefield::_empty_idlist = std::list<msys::Id>();

    bool Forcefield::HasParamTable(const std::string& name) {
      ParamTableMap::const_iterator iter = SharedParamTables.find(name);
      return (iter != SharedParamTables.end());
    }

    msys::ParamTablePtr Forcefield::ParamTable(const std::string& name) {
      ParamTableMap::const_iterator iter = SharedParamTables.find(name);
      if (iter == SharedParamTables.end())
        VIPARR_FAIL("Param table " + name + " not found");
      return iter->second;
    }

    void Forcefield::AddParamTable(const std::string& name,
                                   msys::ParamTablePtr table) {
      ParamTableMap::const_iterator iter = SharedParamTables.find(name);
      if (iter != SharedParamTables.end())
        VIPARR_FAIL("Cannot add table: Param table " + name
                    + " already exists");
      SharedParamTables.insert(std::make_pair(name, table));
    }

    std::vector<std::string> Forcefield::AllParamTables() {
      std::vector<std::string> names;
      for (ParamTableMap::const_iterator iter = SharedParamTables.begin();
           iter != SharedParamTables.end(); ++iter)
        names.push_back(iter->first);
      return names;
    }

    const std::list<msys::Id>& Forcefield::rowIDs(const std::string&
                                                  name) const {
      std::map<std::string, std::list<msys::Id> >::const_iterator iter
        = _row_ids_map.find(name);
      return (iter == _row_ids_map.end() ? _empty_idlist : iter->second);
    }

    void Forcefield::delParam(const std::string& name, msys::Id param) {
      delParams(name, std::list<msys::Id>(1, param));
    }

    void Forcefield::delParams(const std::string& name,
                               const std::list<msys::Id>& params) {
      std::set<msys::Id> param_set(params.begin(), params.end());
      std::map<std::string, std::list<msys::Id> >::iterator iter
        = _row_ids_map.find(name);
      if (iter == _row_ids_map.end())
        return;
      std::list<msys::Id>& all_params = iter->second;
      for (std::list<msys::Id>::iterator list_iter = all_params.begin();
           list_iter != all_params.end(); ) {
        if (param_set.find(*list_iter) != param_set.end())
          list_iter = all_params.erase(list_iter);
        else
          ++list_iter;
      }
    }

    void Forcefield::clearParams(const std::string& name) {
      std::map<std::string, std::list<msys::Id> >::iterator iter
        = _row_ids_map.find(name);
      if (iter == _row_ids_map.end())
        return;
      iter->second.clear();
    }

    void Forcefield::appendParam(const std::string& name, msys::Id param) {
      appendParams(name, std::list<msys::Id>(1, param));
    }

    void Forcefield::appendParams(const std::string& name,
                                  const std::list<msys::Id>& params) {
      if (params.size() == 0)
        return;
      std::map<std::string, msys::ParamTablePtr>::iterator static_map_iter
        = SharedParamTables.find(name);
      if (static_map_iter == SharedParamTables.end())
        VIPARR_FAIL("Shared param table '" + name + "' does not exist");
      std::map<std::string, std::list<msys::Id> >::iterator row_ids_iter
        = _row_ids_map.find(name);
      if (row_ids_iter == _row_ids_map.end())
        row_ids_iter = _row_ids_map.insert(std::make_pair(name,
                                                          std::list<msys::Id>())).first;
      for (msys::Id param : params) {
        if (param >= static_map_iter->second->paramCount()) {
          std::stringstream msg;
          msg << "Parameter " << param << " of shared param table '"
              << name << "' does not exist";
          VIPARR_FAIL(msg.str());
        }
        row_ids_iter->second.push_back(param);
      }
    }

    void Forcefield::replaceParam(const std::string& name, msys::Id old_param,
                                  msys::Id new_param) {
      ParamTableMap::iterator static_iter = SharedParamTables.find(name);
      if (static_iter == SharedParamTables.end()
          || new_param >= static_iter->second->paramCount()) {
        std::stringstream msg;
        msg << "Parameter " << new_param << " of shared param table '"
            << name << "' does not exist";
        VIPARR_FAIL(msg.str());
      }
      std::map<std::string, std::list<msys::Id> >::iterator map_iter
        = _row_ids_map.find(name);
      if (map_iter == _row_ids_map.end())
        return;
      for (std::list<msys::Id>::iterator iter = map_iter->second.begin();
           iter != map_iter->second.end(); ++iter) {
        if (*iter != old_param)
          continue;
        *iter = new_param;
        return;
      }
    }

    std::vector<std::string> Forcefield::paramTables() const {
      std::vector<std::string> tables;
      for (std::map<std::string, std::list<msys::Id> >::const_iterator iter
             = _row_ids_map.begin(); iter != _row_ids_map.end(); ++iter) {
        if (iter->second.size() != 0)
          tables.push_back(iter->first);
      }
      return tables;
    }

    msys::ParamTablePtr Forcefield::cmapTable(unsigned cmap) const {
      if (cmap > _cmaps.size() || cmap <= 0)
        VIPARR_FAIL("Invalid cmap id");
      return _cmaps[cmap-1];
    }

    std::map<std::string, Forcefield::PluginPtr>& Forcefield::PluginRegistry() {
      static std::map<std::string, Forcefield::PluginPtr> registry;
      return registry;
    }

    Forcefield::RegisterPlugin::RegisterPlugin(const std::string& name,
                                               void (*c_match)(TemplatedSystemPtr, ForcefieldPtr),
                                               const std::vector<std::string>& prerequisites,
                                               void (*c_compile)(msys::SystemPtr),
                                               const std::vector<std::string>& dependencies) {
      Forcefield::PluginRegistry().insert(std::make_pair(name,
                                                         Forcefield::PluginPtr(new PluginC(c_match,
                                                                                           prerequisites,
                                                                                           c_compile, dependencies))));
    }

    Forcefield::RegisterPluginPrerequisite::RegisterPluginPrerequisite(const std::string& first_plugin,
                                                                       const std::string& second_plugin) {
      std::map<std::string, Forcefield::PluginPtr>::iterator iter
        = Forcefield::PluginRegistry().find(second_plugin);
      if (iter == Forcefield::PluginRegistry().end()) {
        VIPARR_FAIL("Must register plugin '" + second_plugin
                    + "' before registering its prerequisites");
      }
      iter->second->prerequisites.push_back(first_plugin);
    }

  }}
