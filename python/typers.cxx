#include "wrap_obj.hxx"
#include "../src/template_typer.hxx"
#include "../src/base.hxx"

using namespace desres;
using namespace desres::viparr;
using desres::msys::Id;
using desres::msys::IdList;

namespace {

    list find_template(TemplateTyper& typer, const std::string& name) {
        list L;
        std::vector<TemplatedSystemPtr> tpls = typer.findTemplateByName(name);
        for (unsigned i = 0; i < tpls.size(); ++i)
            L.append(tpls[i]);
        return L;
    }

    tuple match_fragment(TemplateTyper& typer, TemplatedSystemPtr sys, const object& py_atoms) {
        list L(py_atoms);
        IdList atoms(len(L));
        for (unsigned i = 0; i < atoms.size(); ++i)
            atoms[i] = extract<Id>(L[i]);
        std::vector<std::pair<TemplatedSystemPtr, IdList> > matches;
        std::stringstream why_not;
        bool matched = typer.matchFragment(sys, atoms, matches, why_not);
        L = list();
        if (matched) {
            for (unsigned i = 0; i < matches.size(); ++i) {
                list IDs;
                for (unsigned j = 0; j < matches[i].second.size(); ++j)
                    IDs.append(matches[i].second[j]);
                L.append(boost::python::make_tuple(matches[i].first, IDs));
            }
            return boost::python::make_tuple(L, "");
        }
        return boost::python::make_tuple(L, why_not.str());
    }

    void assign_match(const TemplateTyper& typer, TemplatedSystemPtr sys, const object& py_matches, bool rename_atoms, bool rename_residues) {
        list L(py_matches);
        std::vector<std::pair<TemplatedSystemPtr, IdList> > matches(len(L));
        for (unsigned i = 0; i < matches.size(); ++i) {
            if (len(L[i]) != 2)
                VIPARR_FAIL("Error assigning match: Incorrect argument type for matches");
            TemplatedSystemPtr tpl = extract<TemplatedSystemPtr>(L[i][0]);
            list py_IDs(L[i][1]);
            IdList IDs(len(py_IDs));
            for (unsigned j = 0; j < IDs.size(); ++j)
                IDs[j] = extract<Id>(py_IDs[j]);
            matches[i] = std::make_pair(tpl, IDs);
        }
        typer.assignMatch(sys, matches, rename_atoms, rename_residues);
    }

    list templates(const TemplateTyper& typer) {
        list L;
        std::vector<TemplatedSystemPtr> templates = typer.templates();
        for (unsigned i = 0; i < templates.size(); ++i)
            L.append(templates[i]);
        return L;
    }

}

namespace desres { namespace viparr {

    void export_typers() {

        class_<TemplateTyper, TemplateTyperPtr>("TemplateTyper", no_init)
            .def("__eq__", _eq<TemplateTyperPtr>)
            .def("__ne__", _ne<TemplateTyperPtr>)
            .def("__hash__", _hash<TemplateTyperPtr>)
            .def("__init__", make_constructor(&TemplateTyper::create))
            .def("addTemplate", &TemplateTyper::addTemplate)
            .def("delTemplate", &TemplateTyper::delTemplate)
            .def("findTemplate", find_template)
            .def("matchFragment", match_fragment)
            .def("assignMatch", assign_match)
            .def("templates", templates)
            ;
    }

}}
