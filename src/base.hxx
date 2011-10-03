#ifndef desres_viparr_base_hxx
#define desres_viparr_base_hxx

#include <iostream>
#include <stdexcept>

#ifdef _MSC_VER
#define VIPARR_LOC __FILE__ << ":" << __LINE__ << "\n" << __FUNCSIG__
#else
#define VIPARR_LOC __FILE__ << ":" << __LINE__ << "\n" << __PRETTY_FUNCTION__
#endif

#define VIPARR_OUT std::cout
#define VIPARR_ERR std::cerr
#define VIPARR_FAIL(args) do { \
    std::stringstream _viparr_tmp_ss; \
    _viparr_tmp_ss << args << "\nlocation: " << VIPARR_LOC; \
    throw std::runtime_error(_viparr_tmp_ss.str()); \
} while (0)

#endif
