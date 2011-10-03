#ifndef desres_viparr_pattern_hxx
#define desres_viparr_pattern_hxx

#include "templated_system.hxx"
#include <string>
#include <vector>
#include <set>
#include <iostream>

namespace desres { namespace viparr {

    /* Structure that holds atom types, bond types, and additional string flags
     * to match. (Flags currently only used for pseudos, to match "pset".)
     * Bond type conventions:
     *   "-" : order-1 bond
     *   "=" : order-2 bond
     *   "#" : order-3 bond
     *   ":" : aromatic bond
     * */
    struct Pattern {
        static std::string GetBondString(TemplatedSystemPtr sys, msys::Id ai,
                msys::Id aj);

        std::vector<std::string> atoms;
        std::vector<std::string> bonds;
        std::vector<std::string> flags;

        bool operator<(const Pattern &rhs) const;
        bool operator==(const Pattern& rhs) const;
        bool operator!=(const Pattern& rhs) const;

        std::string print() const;
    };
    std::ostream& operator<<(std::ostream& stream, const Pattern& pattern);

    /* Wrapper class for a C++ function pointer or Python function that takes
     * a system and a list of atoms and constructs a Pattern object from them.
     * Contains pre-defined instances of SystemToPattern as static members;
     * these use the above bond type conventions. */
    struct SystemToPattern {
        virtual Pattern operator()(TemplatedSystemPtr sys,
                const msys::IdList& atoms) const = 0;
        virtual bool operator==(const SystemToPattern& other) const = 0;
        bool operator!=(const SystemToPattern& other) const {
            return (!operator==(other));
        }
        virtual ~SystemToPattern() { }

        /* atoms: nbtype(a_0), ... , nbtype(a_n)
         * bonds: None
         * flags: None */
        static std::shared_ptr<SystemToPattern> NBType;

        /* atoms: btype(a_0), ... , btype(a_n)
         * bonds: None
         * flags: None */
        static std::shared_ptr<SystemToPattern> BType;

        /* atoms: btype(a_0), ... , btype(a_n)
         * bonds: bond(a_0,a_1), bond(a_1,a_2), ... , bond(a_{n-1},a_n)
         * flags: None */
        static std::shared_ptr<SystemToPattern> Bonded;

        /* atoms: btype(a_0), ... , btype(a_n)
         * bonds: bond(a_0,a_1), bond(a_0,a_2), ... , bond(a_0,a_n)
         * flags: None */
        static std::shared_ptr<SystemToPattern> BondToFirst;

        /* atoms: btype(a_1), ... , btype(a_n)
         * bonds: None
         * flags: pset(a_0) */
        static std::shared_ptr<SystemToPattern> PseudoBType;

        /* atoms: btype(a_1), ... , btype(a_n)
         * bonds: bond(a_1,a_2), bond(a_1,a_3), ... , bond(a_1,a_n)
         * flags: pset(a_0) */
        static std::shared_ptr<SystemToPattern> PseudoBondToFirst;

        /* atoms: btype(a_1), ... , btype(a_n)
         * bonds: bond(a_1,a_2), bond(a_3,a_2), ... , bond(a_n,a_2)
         * flags: pset(a_0) */
        static std::shared_ptr<SystemToPattern> PseudoBondToSecond;
    };
    typedef std::shared_ptr<SystemToPattern> SystemToPatternPtr;

    /* Wrapper class for a C++ function pointer or Python function that takes
     * a "type" string from a pattern table and constructs a Pattern object it.
     * Contains pre-defined instances of TypeToPattern as static members. */
    struct TypeToPattern {
        virtual Pattern operator()(const std::string&) const = 0;
        virtual bool operator==(const TypeToPattern& other) const = 0;
        bool operator!=(const TypeToPattern& other) const {
            return (!operator==(other));
        }
        virtual ~TypeToPattern() { }

        /* Strings in this set are treated by Default and Pseudo as bond
         * types */
        static std::set<std::string> BondStrings;

        /* Splits string by ' ', adds all tokens in BondStrings to bonds,
         * adds all other tokens to atoms */
        static std::shared_ptr<TypeToPattern> Default;

        /* Splits string by ' ', adds all tokens in BondStrings to bonds,
         * adds last token to flags, adds all other tokens to atoms */
        static std::shared_ptr<TypeToPattern> Pseudo;
    };
    typedef std::shared_ptr<TypeToPattern> TypeToPatternPtr;

    /* Wrapper class for a C++ function pointer or Python function that takes
     * a Pattern and returns a copy with permuted attributes. Contains
     * pre-defined instances of Permutation as static members. */
    struct Permutation {
        virtual Pattern operator()(const Pattern& pattern) const = 0;
        virtual bool operator==(const Permutation& other) const = 0;
        bool operator!=(const Permutation& other) const {
            return (!operator==(other));
        }
        virtual ~Permutation() { }

        /* Returns the same Pattern back */
        static std::shared_ptr<Permutation> Identity;

        /* Reverses the order of atoms and the order of bonds, keeps flags
         * fixed */
        static std::shared_ptr<Permutation> Reverse;

        /* An array of 6 permutations that handle patterns with exactly 4
         * atoms and 3 bonds by assuming a bond-to-first structure, fixing
         * the first (center) atom, and going through all six permutations
         * of the remaining three atoms and corresponding bonds */
        static std::shared_ptr<Permutation> Improper[6];

        /* Invalid Permutation, returned by ParameterMatcher::match if no
         * match is found */
        static std::shared_ptr<Permutation> Null;
    };
    typedef std::shared_ptr<Permutation> PermutationPtr;

}}

#endif
