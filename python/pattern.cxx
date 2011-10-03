#include "wrap_obj.hxx"
#include "../src/base.hxx"
#include "../src/pattern.hxx"

using namespace desres;
using namespace desres::viparr;

namespace {

    list get_bond_strings() {
        list L;
        for (std::set<std::string>::const_iterator iter = TypeToPattern::BondStrings.begin(); iter != TypeToPattern::BondStrings.end(); ++iter)
            L.append(*iter);
        return L;
    }

    void set_bond_strings(const object& obj) {
        list L(obj);
        TypeToPattern::BondStrings.clear();
        for (unsigned i = 0; i < len(L); ++i)
            TypeToPattern::BondStrings.insert(extract<std::string>(L[i]));
    }

    list get_atoms(Pattern& p) {
        list L;
        for (unsigned i = 0; i < p.atoms.size(); ++i)
            L.append(p.atoms[i]);
        return L;
    }

    void set_atoms(Pattern& p, const object& obj) {
        list L(obj);
        p.atoms.resize(len(L));
        for (unsigned i = 0; i < p.atoms.size(); ++i)
            p.atoms[i] = extract<std::string>(L[i]);
    }

    list get_bonds(Pattern& p) {
        list L;
        for (unsigned i = 0; i < p.bonds.size(); ++i)
            L.append(p.bonds[i]);
        return L;
    }

    void set_bonds(Pattern& p, const object& obj) {
        list L(obj);
        p.bonds.resize(len(L));
        for (unsigned i = 0; i < p.bonds.size(); ++i)
            p.bonds[i] = extract<std::string>(L[i]);
    }

    list get_flags(Pattern& p) {
        list L;
        for (unsigned i = 0; i < p.flags.size(); ++i)
            L.append(p.flags[i]);
        return L;
    }

    void set_flags(Pattern& p, const object& obj) {
        list L(obj);
        p.flags.resize(len(L));
        for (unsigned i = 0; i < p.flags.size(); ++i)
            p.flags[i] = extract<std::string>(L[i]);
    }

    Pattern sys_to_pattern_apply(const SystemToPattern& sys_to_pattern, TemplatedSystemPtr sys, const object& obj) {
        list L(obj);
        msys::IdList atoms(len(L));
        for (unsigned i = 0; i < atoms.size(); ++i)
            atoms[i] = extract<msys::Id>(L[i]);
        return sys_to_pattern(sys, atoms);
    }

    list perm_impr_list() {
        list L;
        L.append(Permutation::Identity);
        L.append(Permutation::Improper[1]);
        L.append(Permutation::Improper[2]);
        L.append(Permutation::Improper[3]);
        L.append(Permutation::Improper[4]);
        L.append(Permutation::Improper[5]);
        return L;
    }

    class SystemToPatternPy : public SystemToPattern {
        public:
            static SystemToPatternPtr create(const object& py_sys_to_pattern) {
                return SystemToPatternPtr(new SystemToPatternPy(py_sys_to_pattern));
            }
            virtual Pattern operator()(TemplatedSystemPtr sys, const msys::IdList& atoms) const {
                list L;
                for (unsigned i = 0; i < atoms.size(); ++i)
                    L.append(atoms[i]);
                try {
                    return extract<Pattern>(_py_sys_to_pattern(sys, L));
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error calling Python SystemToPattern: "
                            + std::string(e.what()));
                }
            }
            virtual bool operator==(const SystemToPattern& other) const {
                try {
                    const SystemToPatternPy& other_py
                        = dynamic_cast<const SystemToPatternPy&>(other);
                    return (_py_sys_to_pattern == other_py._py_sys_to_pattern);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
        private:
            SystemToPatternPy(const object& py_sys_to_pattern) : _py_sys_to_pattern(py_sys_to_pattern) { }
            object _py_sys_to_pattern;
    };

    class TypeToPatternPy : public TypeToPattern {
        public:
            static TypeToPatternPtr create(const object& py_type_to_pattern) {
                return TypeToPatternPtr(new TypeToPatternPy(py_type_to_pattern));
            }
            virtual Pattern operator()(const std::string& type) const {
                try {
                    return extract<Pattern>(_py_type_to_pattern(type));
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error calling Python TypeToPattern: "
                            + std::string(e.what()));
                }
            }
            virtual bool operator==(const TypeToPattern& other) const {
                try {
                    const TypeToPatternPy& other_py
                        = dynamic_cast<const TypeToPatternPy&>(other);
                    return (_py_type_to_pattern == other_py._py_type_to_pattern);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
        private:
            TypeToPatternPy(const object& py_type_to_pattern) : _py_type_to_pattern(py_type_to_pattern) { }
            object _py_type_to_pattern;
    };

    class PermutationPy : public Permutation {
        public:
            static PermutationPtr create(const object& py_perm) {
                return PermutationPtr(new PermutationPy(py_perm));
            }
            static PermutationPtr create_empty() {
                return PermutationPtr(new PermutationPy());
            }
            virtual Pattern operator()(const Pattern& pattern) const {
                try {
                    return extract<Pattern>(_py_perm(pattern));
                } catch (std::exception& e) {
                    VIPARR_FAIL("Error calling Python Permutation: "
                            + std::string(e.what()));
                }
            }
            bool operator==(const Permutation& other) const {
                try {
                    const PermutationPy& other_py
                        = dynamic_cast<const PermutationPy&>(other);
                    return (_py_perm == other_py._py_perm);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
        private:
            PermutationPy(const object& py_perm=object()) : _py_perm(py_perm) { }
            object _py_perm;
    };
}

namespace desres { namespace viparr {

    void export_pattern() {
        class_<Pattern>("Pattern")
            .add_property("atoms", get_atoms, set_atoms)
            .add_property("bonds", get_bonds, set_bonds)
            .add_property("flags", get_flags, set_flags)
            .def(self < self)
            .def(self == self)
            .def(self != self)
            .def("__str__", &Pattern::print)
            ;
    }

    void export_system_to_pattern() {
        scope system_to_pattern_class(class_<SystemToPattern, SystemToPatternPtr, boost::noncopyable>("SystemToPattern", no_init)
            .def("__init__", make_constructor(&SystemToPatternPy::create))
            .def("__call__", sys_to_pattern_apply)
            .def(self == self)
            .def(self != self)
            );
        system_to_pattern_class.attr("NBType") = SystemToPattern::NBType;
        system_to_pattern_class.attr("BType") = SystemToPattern::BType;
        system_to_pattern_class.attr("Bonded") = SystemToPattern::Bonded;
        system_to_pattern_class.attr("BondToFirst") = SystemToPattern::BondToFirst;
        system_to_pattern_class.attr("PseudoBType") = SystemToPattern::PseudoBType;
        system_to_pattern_class.attr("PseudoBondToFirst") = SystemToPattern::PseudoBondToFirst;
        system_to_pattern_class.attr("PseudoBondToSecond") = SystemToPattern::PseudoBondToSecond;
    }

    void export_type_to_pattern() {
        scope type_to_pattern_class(class_<TypeToPattern, TypeToPatternPtr, boost::noncopyable>("TypeToPattern", no_init)
            .def("__init__", make_constructor(&TypeToPatternPy::create))
            .def("__call__", &TypeToPattern::operator())
            .def(self == self)
            .def(self != self)
            .add_static_property("BondStrings", get_bond_strings, set_bond_strings)
            );
        type_to_pattern_class.attr("Default") = TypeToPattern::Default;
        type_to_pattern_class.attr("Pseudo") = TypeToPattern::Pseudo;
    }

    void export_permutation() {
        scope permutation_class(class_<Permutation, PermutationPtr, boost::noncopyable>("Permutation", no_init)
            .def("__init__", make_constructor(&PermutationPy::create))
            .def("__init__", make_constructor(&PermutationPy::create_empty))
            .def("__call__", &Permutation::operator())
            .def(self == self)
            .def(self != self)
            );
        permutation_class.attr("Identity") = Permutation::Identity;
        permutation_class.attr("Reverse") = Permutation::Reverse;
        permutation_class.attr("Improper") = perm_impr_list();
        permutation_class.attr("Null") = Permutation::Null;
    }

}}
