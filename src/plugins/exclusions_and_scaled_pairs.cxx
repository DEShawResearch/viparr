#include "../parameter_matcher.hxx"
#include "pairs_helper.hxx"

using namespace desres;
using namespace desres::viparr;
using desres::msys::Id;
using desres::msys::IdList;

static void add_exclusions(msys::SystemPtr sys, Id atm1, Id atm2, Id param) {
    /* atm1 and its pseudos */
    IdList atoms1;
    atoms1.push_back(atm1);
    IdList bonded = sys->bondedAtoms(atm1);
    for (unsigned i = 0; i < bonded.size(); ++i)
        if (sys->atom(bonded[i]).atomic_number == 0)
            atoms1.push_back(bonded[i]);
    /* atm2 and its pseudos */
    IdList atoms2;
    atoms2.push_back(atm2);
    bonded = sys->bondedAtoms(atm2);
    for (unsigned i = 0; i < bonded.size(); ++i)
        if (sys->atom(bonded[i]).atomic_number == 0)
            atoms2.push_back(bonded[i]);
    /* insert exclusions */
    msys::TermTablePtr table = sys->table("exclusion");
    for (unsigned i = 0; i < atoms1.size(); ++i) {
        for (unsigned j = 0; j < atoms2.size(); ++j) {
            if (atoms1[i] == atoms2[j]) continue;
            IdList pair(2);
            if (atoms1[i] < atoms2[j]) {
                pair[0] = atoms1[i];
                pair[1] = atoms2[j];
            } else {
                pair[0] = atoms2[j];
                pair[1] = atoms1[i];
            }
            if (table->findWithAll(pair).size() == 0)
                table->addTerm(pair, param);
        }
    }
}

static void apply_exclusions(TemplatedSystemPtr sys, ForcefieldPtr ff) {
    msys::TermTablePtr table = sys->system()->table("exclusion");
    if (table == msys::TermTablePtr()) {
        table = sys->system()->addTable("exclusion", 2,
                msys::ParamTable::create());
        table->category = msys::EXCLUSION;
    }
    if (table->params() == msys::ParamTablePtr())
        table->resetParams(msys::ParamTable::create());
    table->params()->addProp("separation", msys::IntType);
    table->params()->addProp("es_scale", msys::FloatType);
    table->params()->addProp("lj_scale", msys::FloatType);

    /* Extra exclusions provided by the templates or SMARTS are added first;
     * no scaled pairs or BCI-balanced pairs will be generated for these
     * exclusions */
    const std::vector<IdList>& exclusions = sys->exclusions();
    for (unsigned i = 0, n = exclusions.size(); i < n; ++i) {
        Id p = table->params()->addParam();
        table->params()->value(p, "separation") = -1;
        table->params()->value(p, "es_scale") = 0;
        table->params()->value(p, "lj_scale") = 0;
        Id atm1 = exclusions[i][0];
        Id atm2 = exclusions[i][1];
        add_exclusions(sys->system(), atm1, atm2, p);
    }

    /* Generate all the tuple exclusions, up to the exclusion rule; these may
     * be added to the pairs table as scaled or BCI-balanced pairs */
    std::vector<IdList> alltuples[4];
    for (IdList atom : sys->typedAtoms())
        if (sys->system()->atom(atom[0]).atomic_number > 0)
            alltuples[0].push_back(atom);
    alltuples[1] = sys->nonPseudoBonds();
    alltuples[2] = sys->angles();
    alltuples[3] = sys->dihedrals();
    int rule = ff->rules()->exclusions();
    for (int i = 1; i <= rule; ++i) {
        Id p = table->params()->addParam();
        table->params()->value(p, "separation") = i;
        table->params()->value(p, "es_scale") = (i == 1
                ? 0 : ff->rules()->es_scale(i));
        table->params()->value(p, "lj_scale") = (i == 1
                ? 0 : ff->rules()->lj_scale(i));
        const std::vector<IdList>& tuples = alltuples[i-1];
        for (unsigned j = 0, n = tuples.size(); j < n; ++j) {
            Id atm1 = tuples[j][0];
            Id atm2 = tuples[j][i-1];
            add_exclusions(sys->system(), atm1, atm2, p);
        }
    }
}

static void compile_pairs_es_scaled(msys::SystemPtr sys) {
    msys::TermTablePtr pairs = get_pairs_table(sys);
    if(pairs == msys::TermTablePtr())
        VIPARR_FAIL("Unable to find pairs table.\n");
    
    msys::TermTablePtr exclusions = sys->table("exclusion");
    if (exclusions == msys::TermTablePtr()) {
        VIPARR_ERR << "WARNING: Missing exclusion table" << std::endl;
        return;
    }
    if (exclusions->params() == msys::ParamTablePtr()
            || exclusions->params()->propIndex("es_scale") == msys::BadId)
        return;
    msys::IdList terms = exclusions->terms();
    for (msys::Id term : terms) {
        if (exclusions->param(term) == msys::BadId) continue;
        double es_scale = exclusions->propValue(term, "es_scale").asFloat();
        if (es_scale == 0) continue;
        msys::IdList atoms = exclusions->atoms(term);
        /* Compute qij and add/update pair */
        double qij = es_scale * sys->atom(atoms[0]).charge
            * sys->atom(atoms[1]).charge;
        std::map<std::string, double> params;
        params.insert(std::make_pair("qij", qij));
        update_pair_params(pairs, sys, atoms[0], atoms[1], params);
    }
}

static Forcefield::RegisterPlugin _("exclusions_and_scaled_pairs",
                                    apply_exclusions,
                                    std::vector<std::string>(),
                                    compile_pairs_es_scaled);
