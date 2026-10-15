/*---------------------------------------------------------------------------*\
    PCGFixed: copy of PCG (v2512) that skips the per-iteration residual norm
    (one global reduction) when convergence cannot be reached (tolerance 0,
    relTol 0), computing it once after the loop. Same iterates and same
    reported residuals as PCG.
\*---------------------------------------------------------------------------*/

#ifndef PCGFixed_H
#define PCGFixed_H

#include "lduMatrix.H"

namespace Foam
{

class PCGFixed
:
    public lduMatrix::solver
{
    mutable autoPtr<lduMatrix::preconditioner> preconPtr_;

    PCGFixed(const PCGFixed&) = delete;
    void operator=(const PCGFixed&) = delete;

public:

    TypeName("PCGFixed");

    PCGFixed
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    virtual ~PCGFixed() = default;

    virtual solverPerformance scalarSolve
    (
        solveScalarField& psi,
        const solveScalarField& source,
        const direction cmpt=0
    ) const;

    virtual solverPerformance solve
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt=0
    ) const;
};

}

#endif
