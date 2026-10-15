// True when solverPerformance::checkConvergence can never return true for these
// controls (absolute tolerance <= 0 and relative tolerance not active) and it
// would print nothing. A fixed-iteration solve then needs no residual norm until
// the last iteration.
#ifndef fixedIterControl_H
#define fixedIterControl_H

#include "lduMatrix.H"

namespace Foam
{
inline bool convergenceImpossible
(
    const scalar tolerance,
    const scalar relTol,
    const int log
)
{
    return
        tolerance <= 0
     && relTol <= solverPerformance::small_
     && log < 2
     && lduMatrix::debug < 2
     && solverPerformance::debug < 2;
}
}

#endif
