#include "fix_masses.hxx"
#include "../base.hxx"
#include <algorithm>

void desres::viparr::FixMasses(msys::SystemPtr sys,
        const msys::IdList& atoms, bool verbose) {

    std::map<int, std::map<double, msys::IdList> > atoms_by_anum_mass;
    std::map<int, std::vector<double> > masses_by_anum;
    for (unsigned i = 0; i < atoms.size(); ++i) {
        int anum = sys->atom(atoms[i]).atomic_number;
        double mass = sys->atom(atoms[i]).mass;
        atoms_by_anum_mass[anum][mass].push_back(atoms[i]);
        masses_by_anum[anum].push_back(mass);
    }
    for (std::map<int, std::map<double, msys::IdList> >::iterator iter
            = atoms_by_anum_mass.begin(); iter != atoms_by_anum_mass.end();
            ++iter) {
        if (iter->first == 0) continue;
        if (iter->second.size() == 1) continue;
        if (verbose)
            VIPARR_OUT << "Found " << iter->second.size()
                << " distinct masses for atomic number " << iter->first
                << std::endl;
        std::vector<double>& masses = masses_by_anum[iter->first];
        std::nth_element(masses.begin(), masses.begin() + masses.size() / 2,
                masses.end());
        for (std::map<double, msys::IdList>::iterator m_iter
                = iter->second.begin(); m_iter != iter->second.end();
                ++m_iter) {
            if (verbose)
                VIPARR_OUT << "   " << m_iter->first << " ("
                    << m_iter->second.size() << ")" << std::endl;
            for (unsigned i = 0; i < m_iter->second.size(); ++i)
                sys->atom(m_iter->second[i]).mass = masses[masses.size() / 2];
        }
        if (verbose)
            VIPARR_OUT << "   Replaced with median value "
                << masses[masses.size() / 2] << std::endl;
    }
}
