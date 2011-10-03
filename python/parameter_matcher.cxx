#include "wrap_obj.hxx"
#include "../src/parameter_matcher.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {

    ParameterMatcherPtr construct1(ForcefieldPtr ff, const std::string& table, SystemToPatternPtr sys_to_pattern,
            TypeToPatternPtr type_to_pattern, const object& obj) {
        list L(obj);
        std::vector<PermutationPtr> perms;
        for (unsigned i = 0; i < len(L); ++i)
            perms.push_back(extract<PermutationPtr>(L[i]));
        return ParameterMatcher::create(ff, table, sys_to_pattern, type_to_pattern, perms);
    }

    ParameterMatcherPtr construct2(msys::ParamTablePtr param_table, const object& py_rows, SystemToPatternPtr sys_to_pattern,
            TypeToPatternPtr type_to_pattern, const object& py_perms) {
        list L_rows(py_rows);
        list L_perms(py_perms);
        std::list<msys::Id> rows;
        for (unsigned i = 0; i < len(L_rows); ++i)
            rows.push_back(extract<msys::Id>(L_rows[i]));
        std::vector<PermutationPtr> perms(len(L_perms));
        for (unsigned i = 0; i < perms.size(); ++i)
            perms[i] = extract<PermutationPtr>(L_perms[i]);
        return ParameterMatcher::create(param_table, rows, sys_to_pattern, type_to_pattern, perms);
    }

    tuple match(ParameterMatcher& matcher, TemplatedSystemPtr sys, const object& obj, bool allow_repeat) {
        list L(obj);
        msys::IdList term(len(L));
        for (unsigned i = 0; i < term.size(); ++i)
            term[i] = extract<msys::Id>(L[i]);
        PermutationPtr perm;
        msys::Id row_id = matcher.match(sys, term, &perm, allow_repeat);
        return boost::python::make_tuple(row_id, perm);
    }

    void write_multiple(const ParameterMatcher& matcher, msys::Id row_id, const object& obj, msys::TermTablePtr table) {
        list L(obj);
        msys::IdList term(len(L));
        for (unsigned i = 0; i < term.size(); ++i)
            term[i] = extract<msys::Id>(L[i]);
        matcher.writeMultiple(row_id, term, table);
    }

    list row_ids(const ParameterMatcher& matcher) {
        list L;
        const std::list<msys::Id>& rowIDs = matcher.rowIDs();
        for (msys::Id row : rowIDs)
            L.append(row);
        return L;
    }
}

namespace desres { namespace viparr {

    void export_parameter_matcher() {

        class_<ParameterMatcher, ParameterMatcherPtr>("ParameterMatcher", no_init)
            .def("__eq__", _eq<ParameterMatcherPtr>)
            .def("__ne__", _ne<ParameterMatcherPtr>)
            .def("__hash__", _hash<ParameterMatcherPtr>)
            .def("__init__", make_constructor(construct1))
            .def("__init__", make_constructor(construct2))
            .def("match", match)
            .def("writeMultiple", write_multiple)
            .def("paramTable", &ParameterMatcher::paramTable)
            .def("rowIDs", row_ids)
            .def("sysToPattern", &ParameterMatcher::sysToPattern)
            .def("typeToPattern", &ParameterMatcher::typeToPattern)
            ;
    }
}}
