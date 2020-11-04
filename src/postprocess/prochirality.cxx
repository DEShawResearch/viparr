#include "prochirality.hxx"
#include <algorithm>
#include <array>
#include <msys/atomsel.hxx>
#include <msys/system.hxx>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

namespace {
using namespace desres;
enum Chirality
{
    Chirality_R,
    Chirality_S,
};

struct vec3f_t {
  float x, y, z;

  vec3f_t() = default;
  vec3f_t(float x_, float y_, float z_) : x(x_), y(y_), z(z_) { }
};

float dot(vec3f_t a, vec3f_t b) {
  return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

vec3f_t cross(vec3f_t a, vec3f_t b) {
  return vec3f_t(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  );
}

vec3f_t operator-(vec3f_t a, vec3f_t b) {
  return vec3f_t(
    a.x - b.x,
    a.y - b.y,
    a.z - b.z
  );
}

using fliplist_t = std::vector<std::pair<std::string, std::string>>;
using priorities_t = std::array<std::string, 4>;
const std::vector<std::tuple<std::string, priorities_t, fliplist_t>> IUPAC_DATA = {
    { "resname VAL and name CB", {"HB", "CG2", "CG1", "CA"}, { { "CG1", "CG2" }, { "HG11", "HG21" }, { "HG12", "HG22" }, { "HG13", "HG23" } } },
    { "resname TYR and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname TRP and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname SER and name CB", {"HB2", "HB3", "CA", "OG"}, { { "HB2", "HB3" } } },
    { "resname PHE and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname MET and name CB", {"HB2", "HB3", "CA", "CG"}, { { "HB2", "HB3" } } },
    { "resname MET and name CG", {"HG2", "HG3", "CB", "SD"}, { { "HG2", "HG3" } } },
    { "resname LYS and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname LYS and name CG", {"HG3", "HG2", "CD", "CB"}, { { "HG2", "HG3" } } },
    { "resname LYS and name CD", {"HD2", "HD3", "CG", "CE"}, { { "HD2", "HD3" } } },
    { "resname LYS and name CE", {"HE2", "HE3", "CD", "NZ"}, { { "HE2", "HE3" } } },
    { "resname LEU and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname LEU and name CG", {"HG", "CD2", "CD1", "CB"}, { { "CD1", "CD2" }, { "HD11", "HD21" }, { "HD12", "HD22" }, { "HD13", "HD23" } } },
    { "resname ILE and name CG1", {"HG13", "HG12", "CD1", "CB"}, { { "HG12", "HG13" } } },
    { "resname HIS and name CB", {"HB2", "HB3", "CA", "CG"}, { { "HB2", "HB3" } } },
    { "resname GLU and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname GLU and name CG", {"HG2", "HG3", "CB", "CD"}, { { "HG2", "HG3" } } },
    { "resname GLN and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname GLN and name CG", {"HG2", "HG3", "CB", "CD"}, { { "HG2", "HG3" } } },
    { "resname CYS and name CB", {"HB2", "HB3", "CA", "SG"}, { { "HB2", "HB3" } } },
    { "resname ASP and name CB", {"HB2", "HB3", "CA", "CG"}, { { "HB2", "HB3" } } },
    { "resname ASN and name CB", {"HB2", "HB3", "CA", "CG"}, { { "HB2", "HB3" } } },
    { "resname ARG and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname ARG and name CG", {"HG2", "HG3", "CB", "CD"}, { { "HG2", "HG3" } } },
    { "resname ARG and name CD", {"HD2", "HD3", "CG", "NE"}, { { "HD2", "HD3" } } },
    { "resname GLY and name CA", {"HA3", "HA2", "C", "N"}, { { "HA2", "HA3" } } },
    { "resname PRO and name CB", {"HB3", "HB2", "CG", "CA"}, { { "HB2", "HB3" } } },
    { "resname PRO and name CG", {"HG2", "HG3", "CB", "CD"}, { { "HG2", "HG3" } } },
    { "resname PRO and name CD", {"HD2", "HD3", "CG", "N"}, { { "HD2", "HD3" } } }
};


template<typename T>
int32_t
inversionCount(const T& array)
{
    int inv_count = 0;
    for (size_t i = 0; i < array.size() - 1; i++)
        for (size_t j = i + 1; j < array.size(); j++)
            if (array[i] > array[j])
                inv_count++;

    return inv_count;
}

Chirality
classifyTetrahedralCenter(msys::SystemPtr s, msys::Id centerAtom, const priorities_t& priorityNames)
{
    std::array<vec3f_t, 4> r;
    std::array<int32_t, 4> priorities;
    msys::IdList bonds = s->bondsForAtom(centerAtom);

    size_t bondIterationK = 0;
    auto neighborK = 0;
    for (bondIterationK = 0; bondIterationK < bonds.size(); bondIterationK++) {
        auto bond = s->bondFAST(bonds[bondIterationK]);
        auto neighborAtomId = bond.i != centerAtom ? bond.i : bond.j;
        auto n = s->atomFAST(neighborAtomId);
        if (n.atomic_number > 0) {
            if (neighborK >= 4) {
                throw std::runtime_error("More than a 4-coordinate center");
            }
            auto offset = std::find(priorityNames.begin(), priorityNames.end(), std::string(n.name));
            if (offset == priorityNames.end()) {
                throw std::runtime_error("Didn't find expected atoms around " + std::to_string(centerAtom));
            }
            priorities[neighborK] = std::distance(priorityNames.begin(), offset);
            r[neighborK] = vec3f_t(n.x, n.y, n.z);
            neighborK++;
        }
    }

    double volume = dot(r[1] - r[0], cross(r[2] - r[0], r[3] - r[0]));
    auto order_parity = inversionCount(priorities) % 2 == 0 ? -1 : 1;
    auto coord_parity = volume < 0 ? 1 : -1;
    auto parity = order_parity * coord_parity;
    return parity == 1 ? Chirality_R : Chirality_S;
}
}


namespace desres { namespace viparr {

msys::IdList
FixProchiralProteinAtomNames(msys::SystemPtr s, bool checkOnly)
{
    auto performSwap = [&](msys::Id centerAtom, const fliplist_t& fliplist) {
        auto residueId = s->atomFAST(centerAtom).residue;

        std::vector<std::pair<msys::Id, msys::Id>> flipIds(fliplist.size(), { msys::BadId, msys::BadId });

        for (auto atomId : s->atomsForResidue(residueId)) {
            auto const& name = s->atomFAST(atomId).name;
            for (size_t i = 0; i < fliplist.size(); i++) {
                if (name == fliplist[i].first) {
                    if (flipIds[i].first != msys::BadId) {
                        throw std::runtime_error("Residue contains two atoms with same name");
                    }
                    flipIds[i].first = atomId;
                }
                if (name == fliplist[i].second) {
                    if (flipIds[i].second != msys::BadId) {
                        throw std::runtime_error("Residue contains two atoms with same name");
                    }
                    flipIds[i].second = atomId;
                }
            }
        }

        for (auto id1_id2_pair : flipIds) {
            auto id1 = id1_id2_pair.first;
            auto id2 = id1_id2_pair.second;
            if (id1 == msys::BadId || id2 == msys::BadId) {
                throw std::runtime_error("Could not find atoms to flip residueId={}" + std::to_string(residueId));
            }
            auto& a1 = s->atomFAST(id1);
            auto& a2 = s->atomFAST(id2);
            auto const name1 = a1.name;
            auto const name2 = a2.name;
            a1.name = name2;
            a2.name = name1;
        }
    };

    msys::IdList incorrectCenters;
    for (const auto& item : IUPAC_DATA) {
        const auto centerSelection = std::get<0>(item);
        const auto priorities = std::get<1>(item);
        const auto flipList = std::get<2>(item);

        auto ids = msys::Atomselect(s, centerSelection);
        for (auto id : ids) {
            if (classifyTetrahedralCenter(s, id, priorities) != Chirality_R) {
                incorrectCenters.push_back(id);
                if (!checkOnly) {
                    performSwap(id, flipList);
                }
            }
        }
    }
    return incorrectCenters;
}

}}