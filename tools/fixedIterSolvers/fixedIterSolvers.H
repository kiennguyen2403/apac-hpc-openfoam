/*---------------------------------------------------------------------------*\
    fixedIterSolvers: shared helpers for the LD_PRELOAD replacements in this
    library (smoothSolverFixedPath.C, solutionCache.C, pcgFixedPath.C).

    STEP_REUSE_OFF  comma list of features to switch off (same variable as
                    libstepReuse / libcompactAddr / libfusedVectorSolve):
                      fixed   smoothSolver fixed-iteration path
                      gcache  grad(p), grad(k), grad(omega) cache
                      pcg     PCG fixed-count path (pcgFixedPath.C; not built: no gain)
                    The stock code then runs, through the next definition of
                    the replaced symbol.
    STEP_REUSE_REPORT=1  one line per path on stderr (master) saying which
                    path ran
\*---------------------------------------------------------------------------*/

#ifndef fixedIterSolvers_H
#define fixedIterSolvers_H

#include <cstdlib>
#include <cstring>
#include <string>

namespace Foam
{
namespace fixedIterSolvers
{

inline bool off(const char* key)
{
    const char* e = std::getenv("STEP_REUSE_OFF");
    return e && std::strstr(e, key);
}

//- The next definition (normally the stock one) of a replaced symbol;
//  FatalError when it cannot be found
void* next(const char* mangledName);

//- True when STEP_REUSE_REPORT is set and this is the master rank
bool reporting();

//- Print "fixedIterSolvers: <msg>" once per message when reporting()
void report(const std::string& msg);

} // End namespace fixedIterSolvers
} // End namespace Foam

#endif
