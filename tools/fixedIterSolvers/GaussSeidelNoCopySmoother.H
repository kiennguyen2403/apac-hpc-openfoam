/*---------------------------------------------------------------------------*\
    GaussSeidelNoCopySmoother ("GaussSeidelNoCopy"): GaussSeidel smoother
    without the full copy bPrime = source at the start of every sweep.

    In the stock sweep, bPrime[u] receives source[u], then coupled-interface
    terms (boundary cells only), then the lower-triangle terms of earlier
    cells, and is read when cell u is reached. Here only coupled boundary
    cells are copied up front; any other cell is initialised on its first
    lower-triangle update (bPrime[u] = source[u] - l*psi, i.e. the same
    operation as copying first and subtracting), and cells with no such
    update read source directly. Same operations in the same order, so
    results are bit-identical to GaussSeidel; one full-field read and write
    per sweep is saved.
\*---------------------------------------------------------------------------*/

#ifndef GaussSeidelNoCopySmoother_H
#define GaussSeidelNoCopySmoother_H

#include "lduMatrix.H"

namespace Foam
{

class GaussSeidelNoCopySmoother
:
    public lduMatrix::smoother
{
    // Per face: this face makes the first update of bPrime[upperAddr[face]]
    // (and that cell is not a coupled boundary cell)
    boolList firstWrite_;

    // Per cell: bPrime is never written for this cell, read source instead
    boolList readSource_;

    // Cells referenced by coupled interfaces (copied before the interface update)
    labelList boundaryCells_;

    void sweep
    (
        solveScalarField& psi,
        const solveScalarField& source,
        const direction cmpt,
        const label nSweeps
    ) const;

public:

    TypeName("GaussSeidelNoCopy");

    GaussSeidelNoCopySmoother
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    virtual void smooth
    (
        solveScalarField& psi,
        const scalarField& source,
        const direction cmpt,
        const label nSweeps
    ) const;

    virtual void scalarSmooth
    (
        solveScalarField& psi,
        const solveScalarField& source,
        const direction cmpt,
        const label nSweeps
    ) const;
};

}

#endif
