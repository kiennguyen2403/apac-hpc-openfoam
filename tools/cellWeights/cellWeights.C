/*---------------------------------------------------------------------------*\
    cellWeights: write a per-cell weight field for decomposePar (weightField),
    so that scotch/KaHIP balance estimated work instead of cell count.

    -mode faces:    w = 4*(internal faces of the cell) + 7*(wall faces of the cell)
                    (fitted per-rank cost: ~7.5 us per internal face, shared by
                    two cells, and ~6.5 us extra per wall face; job 181351948)
    -mode measured: w = faces-weight * f_r, where r is the cell's rank in an
                    existing decomposition (-dict) and f_r = (measured cost of r
                    / faces-weight of r), normalised by the median over ranks
                    (-costs file: "rank cost" per line; missing ranks = median)

    Weights are written as integers that OpenFOAM's conversion keeps exactly:
    cell 0 gets weight 1, so the minimum is 1 and int(w/min) = w; fractional
    weights (measured mode) are rounded with a golden-ratio dither, so sums
    over many cells are unbiased.

    Usage: cellWeights -mode faces|measured [-dict decomposeParDict -costs file]
                       [-name cellWeight]
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "decompositionMethod.H"
#include "IFstream.H"
#include "SortableList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("mode", "word", "faces | measured");
    argList::addOption("dict", "file", "decomposeParDict of the measured decomposition");
    argList::addOption("costs", "file", "per-rank measured cost: 'rank cost' lines");
    argList::addOption("name", "word", "field name (default cellWeight)");

    #include "setRootCase.H"
    #include "createTime.H"
    fvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    const word mode(args.get<word>("mode"));
    const word name(args.getOrDefault<word>("name", "cellWeight"));

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const label nCells = mesh.nCells();

    // faces weight (integer)
    scalarField fw(nCells, Zero);
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        fw[own[facei]] += 4;
        fw[nei[facei]] += 4;
    }
    label nWall = 0;
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (isA<wallPolyPatch>(pp))
        {
            forAll(pp, i) { fw[own[pp.start() + i]] += 7; ++nWall; }
        }
    }

    scalarField w(fw);

    if (mode == "measured")
    {
        IFstream dis(args.get<fileName>("dict"));
        dictionary dict(dis);
        autoPtr<decompositionMethod> method = decompositionMethod::New(dict);
        const labelList cellToProc(method->decompose(mesh, mesh.cellCentres()));
        const label n = method->nDomains();

        scalarField cost(n, -1);
        {
            IFstream cis(args.get<fileName>("costs"));
            std::string line;
            while (cis.good() && cis.getLine(line))
            {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream iss(line);
                label r; double c;
                if (iss >> r >> c && r >= 0 && r < n) cost[r] = c;
            }
        }
        scalarField fwRank(n, Zero);
        forAll(cellToProc, celli) fwRank[cellToProc[celli]] += fw[celli];

        // cost per unit faces-weight, normalised by the median over measured ranks
        DynamicList<scalar> ratios;
        forAll(cost, r) if (cost[r] > 0) ratios.push_back(cost[r]/fwRank[r]);
        SortableList<scalar> sorted(ratios);
        const scalar medRatio = sorted[sorted.size()/2];
        scalarField f(n, 1);
        label nMissing = 0;
        forAll(cost, r)
        {
            if (cost[r] > 0) f[r] = (cost[r]/fwRank[r])/medRatio;
            else ++nMissing;
        }
        Info<< "measured: ranks " << n << " without cost (median used) " << nMissing
            << ", rank factor min " << min(f) << " max " << max(f) << endl;

        const scalar phi = 0.6180339887498949;
        forAll(w, celli)
        {
            const scalar x = fw[celli]*f[cellToProc[celli]];
            const scalar u = (celli*phi) - std::floor(celli*phi);
            w[celli] = std::floor(x + u);
        }
    }
    else if (mode != "faces")
    {
        FatalErrorInFunction << "unknown -mode " << mode << exit(FatalError);
    }

    w[0] = 1;   // minimum weight 1: OpenFOAM's int(w/min) then keeps every weight exactly
    forAll(w, celli) if (w[celli] < 1) w[celli] = 1;

    Info<< "weights: min " << min(w) << " max " << max(w)
        << " mean " << average(w) << " sum " << sum(w)
        << " (wall faces " << nWall << ")" << endl;

    volScalarField weights
    (
        IOobject(name, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, 1),
        fvPatchFieldBase::calculatedType()
    );
    weights.primitiveFieldRef() = w;
    weights.writeObject(IOstreamOption(IOstreamOption::BINARY), true);
    Info<< "Wrote " << weights.objectPath() << nl << "End" << endl;
    return 0;
}
