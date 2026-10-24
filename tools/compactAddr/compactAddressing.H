/*---------------------------------------------------------------------------*\
    compactAddressing.H: compact form of an ldu level's addressing, built once
    and cached per lduAddressing (shared by tools/compactAddr and
    tools/fusedVectorSolve).

      off[face] = upper - lower  (uint16, faces in owner order)
      cnt[cell] = faces owned by the cell  (uint8)

    compactOf() returns nullptr when a level does not fit (an offset above
    65535, more than 255 owned faces, or faces not in owner order), so callers
    keep their int32 path. STEP_REUSE_OFF containing "compact" disables it.
\*---------------------------------------------------------------------------*/

#ifndef compactAddressing_H
#define compactAddressing_H

#include "lduAddressing.H"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <vector>

namespace Foam
{
namespace compactAddr
{

inline bool off()
{
    static const bool o = []
    {
        const char* e = std::getenv("STEP_REUSE_OFF");
        return e && std::strstr(e, "compact");
    }();
    return o;
}

struct Compact
{
    long nCells = -1;
    long nFaces = -1;
    const label* uBegin = nullptr;
    long check = 0;
    bool ok = false;
    std::vector<std::uint16_t> off;   // upper - lower per face
    std::vector<std::uint8_t> cnt;    // faces per owner cell
};

inline long sample(const labelUList& u)
{
    long s = 0;
    const long n = u.size();
    for (long i = 0; i < n && i < 64; ++i) s += u[i];
    for (long i = (n > 64 ? n - 64 : 0); i < n; ++i) s += 3*long(u[i]);
    return s;
}

// Compact addressing of an ldu level, or nullptr when it does not fit
inline const Compact* compactOf(const lduAddressing& a)
{
    static std::unordered_map<const lduAddressing*, Compact> cache;

    const labelUList& l = a.lowerAddr();
    const labelUList& u = a.upperAddr();
    const long nCells = a.size();
    const long nFaces = l.size();

    Compact& c = cache[&a];
    if
    (
        c.nCells == nCells && c.nFaces == nFaces
     && c.uBegin == u.cdata() && c.check == sample(u)
    )
    {
        return c.ok ? &c : nullptr;
    }

    c.nCells = nCells;
    c.nFaces = nFaces;
    c.uBegin = u.cdata();
    c.check = sample(u);
    c.ok = false;
    c.off.assign(nFaces, 0);
    c.cnt.assign(nCells, 0);

    const labelUList& os = a.ownerStartAddr();
    for (long cell = 0; cell < nCells; ++cell)
    {
        const long n = os[cell + 1] - os[cell];
        if (n < 0 || n > 255)
        {
            return nullptr;
        }
        c.cnt[cell] = std::uint8_t(n);
        for (long f = os[cell]; f < os[cell + 1]; ++f)
        {
            const long d = long(u[f]) - long(l[f]);
            if (l[f] != cell || d <= 0 || d > 65535)
            {
                return nullptr;
            }
            c.off[f] = std::uint16_t(d);
        }
    }
    if (os[nCells] != nFaces)
    {
        return nullptr;
    }
    c.ok = true;
    return &c;
}

} // End namespace compactAddr
} // End namespace Foam

#endif

// ************************************************************************* //
