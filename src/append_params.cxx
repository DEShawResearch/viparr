#include "append_params.hxx"

using namespace desres;

/* Compare two parameters, possibly from different tables dest and src,
 * lexicographically in the order of properties in dest. dest must contain
 * the properties of src as a subset; any additional properties are assumed
 * to be the default empty value in src. */
desres::viparr::append_params::ParamComparator::ParamComparator(
        msys::ParamTablePtr dest, msys::ParamTablePtr src, const msys::IdList&
        dest_to_src) : _dest(dest), _src(src), _d2s(dest_to_src)
{}

bool desres::viparr::append_params::ParamComparator::operator() (
        const ParamToken& pi, const ParamToken& pj) const {
    msys::Value _empty_val; _empty_val.s=nullptr;
    if (pi.first && pj.first)
        /* Both in dest */
        return _dest->compare(pi.second, pj.second) < 0;
    else if (!pi.first && !pj.first) {
        /* Both in src */
        for (unsigned i = 0; i < _d2s.size(); ++i) {
            if (_d2s[i] == msys::BadId)
                continue;
            int c = _src->value(pi.second, _d2s[i]).compare(
                    _src->value(pj.second, _d2s[i]));
            if (c < 0) return true; // pi < pj
            if (c > 0) return false; // pi > pj
        }
        return false; // pi == pj
    } else if (pi.first && !pj.first) {
        /* pi in dest, pj in src */
        for (unsigned i = 0; i < _d2s.size(); ++i) {
            int c;
            if (_d2s[i] != msys::BadId)
                c = _dest->value(pi.second, i).compare(
                        _src->value(pj.second, _d2s[i]));
            else
                c = _dest->value(pi.second, i).compare(
                        msys::ValueRef(_dest->propType(i),
                            _empty_val));
            if (c < 0) return true; // pi < pj
            if (c > 0) return false; // pi > pj
        }
        return false; // pi == pj
    } else {
        /* pi in src, pj in dest */
        for (unsigned i = 0; i < _d2s.size(); ++i) {
            int c;
            if (_d2s[i] != msys::BadId)
                c = _src->value(pi.second, _d2s[i]).compare(
                        _dest->value(pj.second, i));
            else
                c = msys::ValueRef(_dest->propType(i),
                        _empty_val).compare(_dest->value(pj.second,
                                i));
            if (c < 0) return true; // pi < pj
            if (c > 0) return false; // pi > pj
        }
        return false; // pi == pj
    }
}
