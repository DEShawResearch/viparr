#include "optimize_vsitedefs.hxx"
#include "../base.hxx"

#include <cstdio>
#include <sstream>
#include <vector>
#include <map>

#include <boost/assign/list_inserter.hpp>

#include <msys/term_table.hxx>

using namespace desres;
using desres::msys::Id;
using desres::msys::IdList;

namespace {

    /* Cache of param table values, to avoid addition of duplicate rows */
    typedef std::map<std::pair<std::string, std::vector<double> >, Id> ParamCache;

    /* True if we converted the term */
    bool convert_vdef_term(msys::SystemPtr sys, std::string const& newName,
                           msys::TermTablePtr vtable, Id vterm,
                           std::vector<msys::TermTablePtr> &ctables,
                           ParamCache &cache){
        
        /* sites includes the Id of the virtual site */
        IdList sites=vtable->atoms(vterm);
 
#if 0       
        printf("Optimizing %s: Nparams=%u Term=%u Ids=",
               vtable->name().c_str(),
               vtable->termPropCount(), vterm);
        for (Id idx=0; idx<sites.size();++idx){
            printf(" %u",sites[idx]);
        }
        printf("\n");
#endif

        /* See if we can convert this lcXn term into its non-normalized version.
         * Requirements: 
         * 1) # of atoms in constraint is >= # of atoms in vdef
         * 2) Central (eg first) lcXn site is central (eg first) atom in the constraint def
         * 3) All non-central sites in vdef are found in the constraint
         */
        for (msys::TermTablePtr ctable : ctables){
            /* Condition #1 */
            if (ctable->atomCount() < sites.size()-1) continue;
            
            //printf("Contraint table %s passes check 1\n",ctable->name().c_str());

            for(const auto cterm : ctable->findWithAny(sites)) {
                auto atoms_size = ctable->atomCount();
                
                /* condition #2 */
                if(sites[1] != ctable->atom(cterm, 0)) continue;
                //printf("  Contraint term %u passes check 2\n",cterm);

                IdList idxMap;
                for (Id idx=2; idx<sites.size();++idx){
                    for(Id jdx=1;jdx<atoms_size;++jdx){
                        if(sites[idx] !=  ctable->atom(cterm, jdx))continue;
                        idxMap.push_back(jdx);
                        
                    }
                }
                /* condition #3 */
                if(idxMap.size()<sites.size()-2) continue;

                //printf("    Contraint term %u passes check 3\n",cterm);

                msys::TermTablePtr lcxTable = sys->table(newName);
                msys::ParamTablePtr oldParams=vtable->params();
                if (lcxTable == msys::TermTablePtr()) {
                    //printf("Creating new table: %s\n   propnames=",newName.c_str());
                    msys::ParamTablePtr params = msys::ParamTable::create();
                    /* Clone Schema... Should be made a low-level function with msys */
                    for(Id prop=0;prop<oldParams->propCount();++prop){
                        //printf(" %s",oldParams->propName(prop).c_str());
                        params->addProp(oldParams->propName(prop), oldParams->propType(prop));
                    }
                    //printf("\n");
                    lcxTable = sys->addTable(newName, sites.size(), params);
                    lcxTable->category = msys::VIRTUAL;
                }

                /* Rewrite vterm using constraint information */
                std::vector<double> newParams;
                for (Id idx=0; idx<idxMap.size();++idx){
                    /* get vsite param */
                    std::stringstream prop_name;
                    prop_name << "c" << idx+1;
                    double oldP=vtable->propValue(vterm, prop_name.str()).asFloat();
                    
                    /* get constraint param */
                    prop_name.str("");
                    prop_name << "r" << idxMap[idx];
                    double dist=ctable->propValue(cterm, prop_name.str()).asFloat();

                    newParams.push_back(oldP/dist);
                }

                /* If param value exists, use existing param. Otherwise, add new
                 * param */
                //printf("Cache Size: %lu\n",cache.size());
                ParamCache::iterator cache_iter = cache.find(ParamCache::key_type(newName, newParams));
                Id param;
                if (cache_iter != cache.end())
                    param = cache_iter->second;
                else {
                    param = lcxTable->params()->addParam();
                    for (Id idx = 0; idx <newParams.size(); ++idx){
                        //printf("Adding new parameter: %u -> %f\n",idx,newParams[idx]);
                        lcxTable->params()->value(param, idx) = newParams[idx];
                    }
                    cache.insert(ParamCache::value_type(ParamCache::key_type(newName,newParams),param));
                }
                lcxTable->addTerm(sites, param);
                        
                return true;
                
            }
            
        }
        /* couldnt do anything */
        return false;
    }
}


namespace desres { namespace viparr {
    
    void OptimizeVsiteDefs(msys::SystemPtr sys, bool verbose){
        
        std::map<std::string,std::string> toOptimize;
        boost::assign::insert(toOptimize)
            ( "virtual_lc2n", "virtual_lc2" ) 
            ( "virtual_lc3n", "virtual_lc3" ) 
            ( "virtual_lc4n", "virtual_lc4" );
        
        std::set<std::string> supportedCons;
        supportedCons.insert("constraint_hoh");
        std::stringstream ss;
        for(unsigned i=1;i<9;++i){
            ss.str("");
            ss<<"constraint_ah"<<i;
            supportedCons.insert(ss.str());
        }

        std::vector<msys::TermTablePtr> vtables;
        std::vector<msys::TermTablePtr> ctables;
        bool hasVsites=false;
        /* Found out if there are any constraints */
        std::vector<std::string> names = sys->tableNames();
        for (unsigned i = 0; i < names.size(); ++i) {
            msys::TermTablePtr table = sys->table(names[i]);
            if (table->category == msys::CONSTRAINT){
                if(supportedCons.count(names[i])){
                    ctables.push_back(table);
                }else{
                    VIPARR_ERR
                        << "WARNING: Skipping unsupported constraint type "
                        << names[i] << " in OptimizeVsiteDefs" << std::endl;
                }
            }else if(table->category == msys::VIRTUAL){
                hasVsites=true;
                if (toOptimize.count(names[i]))
                    vtables.push_back(table);
            }
        }
        
        if(vtables.size()==0){
            if(hasVsites && verbose)
                VIPARR_OUT << "No virtual site definitions need to be optimized"
                    << std::endl;
            return;
        }else if( ctables.size()==0 ){
            if (verbose)
                VIPARR_OUT << "No virtual site definitions can be optimized "
                    << "(no constraint terms found)" << std::endl;
            return;
        }
        
        ParamCache cache;
        for ( msys::TermTablePtr vtable : vtables){
            Id total=0;
            Id converted=0;
            std::string newType=toOptimize.find(vtable->name())->second;
            for (Id vterm : vtable->terms()){
                ++total;
                if(convert_vdef_term(sys, newType, vtable, vterm, ctables,
                            cache)){
                    vtable->delTerm(vterm);
                    ++converted;
                }
            }
            if (verbose)
                VIPARR_OUT << "Converted " << converted << " of " << total
                    << " " << vtable->name() << " terms to " << newType
                    << std::endl;
            if (vtable->termCount()==0){
                sys->delTable(vtable->name());
            }
        }
        
    }
    
}}
