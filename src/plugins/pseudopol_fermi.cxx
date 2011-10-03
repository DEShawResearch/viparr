#include "../parameter_matcher.hxx"
#include <algorithm>

using namespace desres;
using namespace desres::viparr;
using msys::Id;
using msys::IdList;

static const int natoms = 4;
static const std::string tablename = "pseudopol_fermi";

/* Carbonyl oxygens are those that have btype O, and are bonded to a C 
 * with btype C, which has a total of three bonds. */
static bool is_carbonyl(TemplatedSystemPtr sys, Id id) {
    msys::SystemPtr mol = sys->system();

    if (mol->atom(id).atomic_number!=8 ||
        mol->bondCountForAtom(id)!=1 ||
        sys->btype(id)!="O") return false;

    Id c = mol->bondedAtoms(id).at(0);

    if (mol->atom(c).atomic_number!=6 ||
        mol->bondCountForAtom(c)!=3 ||
        sys->btype(c)!="C") return false;

    return true;
}

/* Amide hydrogens are those with type H, bonded to an N with type NH1. */
static bool is_amide(TemplatedSystemPtr sys, Id id) {
    msys::SystemPtr mol = sys->system();

    if (mol->atom(id).atomic_number!=1 ||
        mol->bondCountForAtom(id)!=1 ||
        sys->btype(id)!="H") return false;

    Id n = mol->bondedAtoms(id).at(0);

    if (mol->atom(n).atomic_number!=7 ||
        mol->bondCountForAtom(n)!=3 ||
        sys->btype(n)!="NH1") return false;

    return true;
}

namespace {
    struct Entry {
        msys::Id O;
        msys::Id H;

        Entry() : O(msys::BadId), H(msys::BadId) {}
        Entry(msys::Id o, msys::Id h) : O(o), H(h) {}
    };
}

static void pseudopol_fermi( TemplatedSystemPtr sys, ForcefieldPtr ff) {

    if (ff->rowIDs(tablename).size() == 0)
        VIPARR_FAIL("Must have '" + tablename + "' table for '" + tablename +
                "' plugin");
    msys::ParamTablePtr params = ff->ParamTable(tablename);

    /* term table */
    msys::SystemPtr mol = sys->system();
    msys::TermTablePtr table = mol->addTable(tablename, natoms, params);
    table->category = msys::BOND;

    /* make sure there is precisely one parameter entry */
    if (ff->rowIDs(tablename).size() != 1)
        VIPARR_FAIL("Expected exactly one entry pseudopol_fermi");

    /* mapping from residue to carbonyl, amide atoms */
    typedef std::map<Id, Entry> ResMap;
    ResMap resmap;

    /* Scan for carbonyls and amides */
    static const Entry defval;
    for (IdList atm : sys->typedAtoms()) {
        Id res = mol->atom(atm[0]).residue;
        if (mol->atomCountForResidue(res)<6) continue;
        /* create an entry for this residue if we don't have one already */
        ResMap::iterator p = resmap.insert(std::make_pair(res,defval)).first;
        if (is_carbonyl(sys,atm[0]))   p->second.O = atm[0];
        else if (is_amide(sys,atm[0])) p->second.H = atm[0];
    }

    /* locate O in residue n and n+3; H in residue n+4 and n+7.  */
    IdList ids(4);
    for (ResMap::const_iterator r0=resmap.begin(); r0!=resmap.end(); ++r0) {
        ResMap::const_iterator r3=resmap.find(r0->first+3);
        ResMap::const_iterator r4=resmap.find(r0->first+4);
        ResMap::const_iterator r7=resmap.find(r0->first+7);
        if (r3==resmap.end() ||
            r4==resmap.end() ||
            r7==resmap.end()) continue;

        ids[0] = r0->second.O;
        ids[1] = r4->second.H;
        ids[2] = r3->second.O;
        ids[3] = r7->second.H;
        if (std::find(ids.begin(), ids.end(), (Id)msys::BadId)
                != ids.end()) continue;

        /* make sure atoms are in the same fragment */
        const Id f0 = mol->atom(ids[0]).fragid;
        const Id f1 = mol->atom(ids[1]).fragid;
        const Id f2 = mol->atom(ids[2]).fragid;
        const Id f3 = mol->atom(ids[3]).fragid;
        if (f0!=f1 || f0!=f2 || f0!=f3) continue;

        /* we have a match - generate a term */
        table->addTerm(ids, ff->rowIDs(tablename).front());
    }
}

static Forcefield::RegisterPlugin _(tablename, pseudopol_fermi);
