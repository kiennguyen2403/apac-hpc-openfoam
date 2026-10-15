/*---------------------------------------------------------------------------*\
    smoothSolverFixed: smoothSolver that skips the per-sweep residual (a full
    matrix pass plus a global reduction) when convergence cannot be reached
    (tolerance 0, relTol 0), computing it once after the last sweep. Same
    sweeps, same psi and same reported residuals as smoothSolver; otherwise
    falls back to smoothSolver::solve.
\*---------------------------------------------------------------------------*/

#ifndef smoothSolverFixed_H
#define smoothSolverFixed_H

#include "smoothSolver.H"

namespace Foam
{

class smoothSolverFixed
:
    public smoothSolver
{
public:

    TypeName("smoothSolverFixed");

    smoothSolverFixed
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    virtual ~smoothSolverFixed() = default;

    virtual solverPerformance solve
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt=0
    ) const;
};

}

#endif
