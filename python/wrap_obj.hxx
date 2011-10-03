#include <boost/python.hpp>

using namespace boost::python;

namespace {

    typedef with_custodian_and_ward_postcall<0,1> return_obj;
    typedef return_internal_reference<> return_ptr;
    typedef return_value_policy<copy_const_reference> return_const;

    template <class Obj>
        bool _eq(const Obj& self, const Obj& other) { return self==other; }

    template <class Obj>
        bool _ne(const Obj& self, const Obj& other) { return self!=other; }

    template <class Obj>
        unsigned long _hash(const Obj& obj) { return reinterpret_cast<unsigned long>(obj.get()); }

}
