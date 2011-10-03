#include "base.hxx"
#include "pattern.hxx"
#include <boost/assign/list_of.hpp>
#include <boost/tokenizer.hpp>
#include <sstream>

using namespace desres;
using namespace desres::viparr;
namespace {

    Pattern sp_nbtype(TemplatedSystemPtr sys, const msys::IdList& atoms) {
        Pattern p;
        p.atoms.resize(atoms.size());
        for (unsigned i = 0; i < atoms.size(); ++i)
            p.atoms[i] = sys->nbtype(atoms[i]);
        return p;
    }

    Pattern sp_btype(TemplatedSystemPtr sys, const msys::IdList& atoms) {
        Pattern p;
        p.atoms.resize(atoms.size());
        for (unsigned i = 0; i < atoms.size(); ++i)
            p.atoms[i] = sys->btype(atoms[i]);
        return p;
    }

    Pattern sp_bonded(TemplatedSystemPtr sys, const msys::IdList& atoms) {
        Pattern p;
        p.atoms.resize(atoms.size());
        p.bonds.resize(atoms.size() - 1);
        p.atoms[0] = sys->btype(atoms[0]);
        for (unsigned i = 1; i < atoms.size(); ++i) {
            p.atoms[i] = sys->btype(atoms[i]);
            std::string bond = Pattern::GetBondString(sys, atoms[i-1],
                    atoms[i]);
            p.bonds[i-1] = bond;
        }
        return p;
    }

    Pattern sp_bond_to_first(TemplatedSystemPtr sys,
            const msys::IdList& atoms) {
        Pattern p;
        p.atoms.resize(atoms.size());
        p.bonds.resize(atoms.size() - 1);
        p.atoms[0] = sys->btype(atoms[0]);
        for (unsigned i = 1; i < atoms.size(); ++i) {
            p.atoms[i] = sys->btype(atoms[i]);
            std::string bond = Pattern::GetBondString(sys, atoms[0], atoms[i]);
            p.bonds[i-1] = bond;
        }
        return p;
    }

    Pattern sp_pseudo_btype(TemplatedSystemPtr sys,
            const msys::IdList& atoms) {
        Pattern p;
        p.atoms.resize(atoms.size() - 1);
        p.flags.resize(1, sys->pset(atoms[0]));
        p.atoms[0] = sys->btype(atoms[1]);
        for (unsigned i = 2; i < atoms.size(); ++i)
            p.atoms[i-1] = sys->btype(atoms[i]);
        return p;
    }

    Pattern sp_pseudo_bond_to_first(TemplatedSystemPtr sys,
            const msys::IdList& atoms) {
        Pattern p;
        p.atoms.resize(atoms.size() - 1);
        p.bonds.resize(atoms.size() - 2);
        p.flags.resize(1, sys->pset(atoms[0]));
        p.atoms[0] = sys->btype(atoms[1]);
        for (unsigned i = 2; i < atoms.size(); ++i) {
            p.atoms[i-1] = sys->btype(atoms[i]);
            std::string bond = Pattern::GetBondString(sys, atoms[1], atoms[i]);
            p.bonds[i-2] = bond;
        }
        return p;
    }

    Pattern sp_pseudo_bond_to_second(TemplatedSystemPtr sys,
            const msys::IdList& atoms) {
        Pattern p;
        p.atoms.resize(atoms.size() - 1);
        p.bonds.resize(atoms.size() - 2);
        p.flags.resize(1, sys->pset(atoms[0]));
        p.atoms[0] = sys->btype(atoms[1]);
        p.atoms[1] = sys->btype(atoms[2]);
        p.bonds[0] = Pattern::GetBondString(sys, atoms[1], atoms[2]);
        for (unsigned i = 3; i < atoms.size(); ++i) {
            p.atoms[i-1] = sys->btype(atoms[i]);
            std::string bond = Pattern::GetBondString(sys, atoms[2], atoms[i]);
            p.bonds[i-2] = bond;
        }
        return p;
    }

    void tokenize(const std::string& type, std::vector<std::string>& tokens) {
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char> > boost_tokens(type, sep);
        tokens = std::vector<std::string>(boost_tokens.begin(),
                boost_tokens.end());
    }

    Pattern tp_default(const std::string& type) {
        std::vector<std::string> tokens;
        tokenize(type, tokens);
        Pattern p;
        for (unsigned i = 0; i < tokens.size(); ++i) {
            if (TypeToPattern::BondStrings.find(tokens[i])
                    == TypeToPattern::BondStrings.end())
                p.atoms.push_back(tokens[i]);
            else
                p.bonds.push_back(tokens[i]);
        }
        return p;
    }

    Pattern tp_pseudo(const std::string& type) {
        std::vector<std::string> tokens;
        tokenize(type, tokens);
        Pattern p;
        for (unsigned i = 0; i < tokens.size() - 1; ++i) {
            if (TypeToPattern::BondStrings.find(tokens[i])
                    == TypeToPattern::BondStrings.end())
                p.atoms.push_back(tokens[i]);
            else
                p.bonds.push_back(tokens[i]);
        }
        p.flags.push_back(tokens[tokens.size() - 1]);
        return p;
    }

    Pattern perm_identity(const Pattern& pattern) { return pattern; }

    Pattern perm_reverse(const Pattern& pattern) {
        Pattern p = pattern;
        std::reverse(p.atoms.begin(), p.atoms.end());
        std::reverse(p.bonds.begin(), p.bonds.end());
        return p;
    }

    /* Atoms 0132, bonds 021 */
    Pattern perm_improper1(const Pattern& pattern) {
        if (pattern.atoms.size() != 4 || pattern.bonds.size() != 3)
            VIPARR_FAIL("Improper Permutation must be used on patterns with "
                    "4 atoms and 3 bonds");
        Pattern p = pattern;
        std::swap(p.atoms[2], p.atoms[3]);
        std::swap(p.bonds[1], p.bonds[2]);
        return p;
    }

    /* Atoms 0213, bonds 102 */
    Pattern perm_improper2(const Pattern& pattern) {
        if (pattern.atoms.size() != 4 || pattern.bonds.size() != 3)
            VIPARR_FAIL("Improper Permutation must be used on patterns with "
                    "4 atoms and 3 bonds");
        Pattern p = pattern;
        std::swap(p.atoms[1], p.atoms[2]);
        std::swap(p.bonds[0], p.bonds[1]);
        return p;
    }

    /* Atoms 0231, bonds 120 */
    Pattern perm_improper3(const Pattern& pattern) {
        if (pattern.atoms.size() != 4 || pattern.bonds.size() != 3)
            VIPARR_FAIL("Improper Permutation must be used on patterns with "
                    "4 atoms and 3 bonds");
        Pattern p = pattern;
        std::swap(p.atoms[1], p.atoms[3]);
        std::swap(p.bonds[0], p.bonds[2]);
        std::swap(p.atoms[1], p.atoms[2]);
        std::swap(p.bonds[0], p.bonds[1]);
        return p;
    }

    /* Atoms 0312, bonds 201 */
    Pattern perm_improper4(const Pattern& pattern) {
        if (pattern.atoms.size() != 4 || pattern.bonds.size() != 3)
            VIPARR_FAIL("Improper Permutation must be used on patterns with "
                    "4 atoms and 3 bonds");
        Pattern p = pattern;
        std::swap(p.atoms[1], p.atoms[2]);
        std::swap(p.bonds[0], p.bonds[1]);
        std::swap(p.atoms[1], p.atoms[3]);
        std::swap(p.bonds[0], p.bonds[2]);
        return p;
    }

    /* Atoms 0321, bonds 210 */
    Pattern perm_improper5(const Pattern& pattern) {
        if (pattern.atoms.size() != 4 || pattern.bonds.size() != 3)
            VIPARR_FAIL("Improper Permutation must be used on patterns with "
                    "4 atoms and 3 bonds");
        Pattern p = pattern;
        std::swap(p.atoms[1], p.atoms[3]);
        std::swap(p.bonds[0], p.bonds[2]);
        return p;
    }

    class SystemToPatternC : public SystemToPattern {
        public:
            SystemToPatternC(Pattern (*c_sys_to_pattern)(TemplatedSystemPtr,
                        const msys::IdList&))
                : _c_sys_to_pattern(c_sys_to_pattern) {
                    if (c_sys_to_pattern == NULL)
                        VIPARR_FAIL("SystemToPattern function cannot be NULL");
            }
            virtual Pattern operator()(TemplatedSystemPtr sys,
                    const msys::IdList& atoms) const {
                return _c_sys_to_pattern(sys, atoms);
            }
            virtual bool operator==(const SystemToPattern& other) const {
                try {
                    const SystemToPatternC& other_c
                        = dynamic_cast<const SystemToPatternC&>(other);
                    return (_c_sys_to_pattern == other_c._c_sys_to_pattern);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
        private:
            Pattern (*_c_sys_to_pattern)(TemplatedSystemPtr,
                    const msys::IdList&);
    };

    class TypeToPatternC : public TypeToPattern {
        public:
            TypeToPatternC(Pattern (*c_type_to_pattern)(const std::string&))
                : _c_type_to_pattern(c_type_to_pattern) {
                    if (c_type_to_pattern == NULL)
                        VIPARR_FAIL("TypeToPattern function cannot be NULL");
            }
            virtual Pattern operator()(const std::string& type) const {
                return _c_type_to_pattern(type);
            }
            virtual bool operator==(const TypeToPattern& other) const {
                try {
                    const TypeToPatternC& other_c
                        = dynamic_cast<const TypeToPatternC&>(other);
                    return (_c_type_to_pattern == other_c._c_type_to_pattern);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
        private:
            Pattern (*_c_type_to_pattern)(const std::string&);
    };

    class PermutationC : public Permutation {
        public:
            PermutationC(Pattern (*c_perm)(const Pattern&)=NULL)
                : _c_perm(c_perm) { }
            virtual Pattern operator()(const Pattern& pattern) const {
                if (_c_perm == NULL)
                    VIPARR_FAIL("Permutation function is NULL");
                return _c_perm(pattern);
            }
            bool operator==(const Permutation& other) const {
                try {
                    const PermutationC& other_c
                        = dynamic_cast<const PermutationC&>(other);
                    return (_c_perm == other_c._c_perm);
                } catch (std::bad_cast& error) {
                    return false;
                }
            }
        private:
            Pattern (*_c_perm)(const Pattern&);
    };
}

namespace desres { namespace viparr {

    std::string Pattern::GetBondString(TemplatedSystemPtr sys, msys::Id ai,
            msys::Id aj) {
        msys::Id bond = sys->system()->findBond(ai, aj);
        if (bond == msys::BadId) {
            std::stringstream msg;
            msg << "GetBondString error: bond between atom "
                << ai << " and atom " << aj << " does not exist";
            VIPARR_FAIL(msg.str());
        }
        if (sys->aromatic(bond))
            return ":";
        int order = sys->system()->bond(bond).order;
        switch (order) {
            case 1: return "-";
            case 2: return "=";
            case 3: return "#";
            default: std::stringstream msg;
                     msg << "GetBondString error: Cannot handle bond order "
                         << order << " between atoms " << ai << " and " << aj;
                     VIPARR_FAIL(msg.str());
        }
        return "";
    }

    bool Pattern::operator<(const Pattern& rhs) const {
        return atoms < rhs.atoms
            || (atoms == rhs.atoms && (bonds < rhs.bonds
                        || (bonds == rhs.bonds && flags < rhs.flags)));
    }

    bool Pattern::operator==(const Pattern& rhs) const {
        return atoms == rhs.atoms && bonds == rhs.bonds && flags == rhs.flags;
    }

    bool Pattern::operator!=(const Pattern& rhs) const {
        return atoms != rhs.atoms || bonds != rhs.bonds || flags != rhs.flags;
    }

    std::string Pattern::print() const {
        std::stringstream out;
        out << "(" << atoms[0];
        for (unsigned i = 1; i < atoms.size(); ++i)
            out << ", " << atoms[i];
        for (unsigned i = 0; i < bonds.size(); ++i)
            out << ", " << bonds[i];
        for (unsigned i = 0; i < flags.size(); ++i)
            out << ", " << flags[i];
        out << ")";
        return out.str();
    }

    std::ostream& operator<<(std::ostream& stream, const Pattern& pattern) {
        return (stream << pattern.print());
    }

    SystemToPatternPtr SystemToPattern::NBType(new SystemToPatternC(sp_nbtype));
    SystemToPatternPtr SystemToPattern::BType(new SystemToPatternC(sp_btype));
    SystemToPatternPtr SystemToPattern::Bonded(new SystemToPatternC(sp_bonded));
    SystemToPatternPtr
        SystemToPattern::BondToFirst(new SystemToPatternC(sp_bond_to_first));
    SystemToPatternPtr
        SystemToPattern::PseudoBType(new SystemToPatternC(sp_pseudo_btype));
    SystemToPatternPtr SystemToPattern::PseudoBondToFirst(
            new SystemToPatternC(sp_pseudo_bond_to_first));
    SystemToPatternPtr SystemToPattern::PseudoBondToSecond(
            new SystemToPatternC(sp_pseudo_bond_to_second));

    std::set<std::string> TypeToPattern::BondStrings
        = boost::assign::list_of("-")("=")("#")(":")("~");

    TypeToPatternPtr TypeToPattern::Default(new TypeToPatternC(tp_default));
    TypeToPatternPtr TypeToPattern::Pseudo(new TypeToPatternC(tp_pseudo));

    PermutationPtr Permutation::Identity(new PermutationC(perm_identity));
    PermutationPtr Permutation::Reverse(new PermutationC(perm_reverse));
    PermutationPtr Permutation::Null(new PermutationC());
    PermutationPtr Permutation::Improper[6] = {
        PermutationPtr(new PermutationC(perm_identity)),
        PermutationPtr(new PermutationC(perm_improper1)),
        PermutationPtr(new PermutationC(perm_improper2)),
        PermutationPtr(new PermutationC(perm_improper3)),
        PermutationPtr(new PermutationC(perm_improper4)),
        PermutationPtr(new PermutationC(perm_improper5))
    };
}}
