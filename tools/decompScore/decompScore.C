/*---------------------------------------------------------------------------*\
    decompScore: read the mesh once, then for every candidate decomposeParDict
    run the decomposition method and write per-rank statistics (no processor
    directories are written). Used to screen decompositions with a per-rank
    cost model before any solver run.

    Each candidate runs in its own fork()ed child process. The child starts
    from the parent's state after the mesh load (mesh and cell centres shared
    copy-on-write) with the decomposition libraries untouched, as in a fresh
    decomposePar. Running candidates one after another in a single process
    does not reproduce decomposePar for scotch: state left by one scotch run
    changes the next result.

    Usage: decompScore -candidates <file: "name dictPath [A B C]" per line>
                       -outDir <dir> [-ranksPerNode 104] [-nParallel 1]
    Output: <outDir>/<name>.ranks with columns
            rank cells procFaces bndFaces neighbours offNodeFaces
            internalFaces wallFaces wallCells

    With "A B C" the candidate is decomposed with cell weights
        w = A*(internal faces of the cell) + B*(wall faces of the cell) + C,
    rounded with a golden-ratio dither, cell 0 set to 1 and every weight >= 1:
    the integer weights cellWeights -mode linear writes for decomposePar's
    weightField (identical formula), so the result is the decomposition
    decomposePar would make with that field.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "decompositionMethod.H"
#include "IFstream.H"
#include "OFstream.H"
#include "clockTime.H"
#include "wallPolyPatch.H"

#include <sys/wait.h>
#include <unistd.h>
#include <map>
#include <cmath>
#include <sstream>

using namespace Foam;

// Integer cell weights, as cellWeights -mode linear
static scalarField linearCellWeights
(
    const polyMesh& mesh,
    const scalar A,
    const scalar B,
    const scalar C
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const label nCells = mesh.nCells();

    scalarField nInt(nCells, Zero);
    scalarField nWall(nCells, Zero);
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        nInt[own[facei]] += 1;
        nInt[nei[facei]] += 1;
    }
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (isA<wallPolyPatch>(pp))
        {
            forAll(pp, i) { nWall[own[pp.start() + i]] += 1; }
        }
    }

    const scalar phi = 0.6180339887498949;
    scalarField w(nCells);
    forAll(w, celli)
    {
        const scalar x = A*nInt[celli] + B*nWall[celli] + C;
        const scalar u = (celli*phi) - std::floor(celli*phi);
        w[celli] = std::floor(x + u);
        if (w[celli] < 1) w[celli] = 1;
    }
    w[0] = 1;
    return w;
}

// Decompose with one candidate dictionary and write its per-rank table
static void scoreCandidate
(
    const polyMesh& mesh,
    const pointField& cc,
    const word& name,
    const fileName& dictPath,
    const fileName& outDir,
    const label ranksPerNode,
    const scalarField* weightCoeffs
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const label nInternal = mesh.nInternalFaces();
    const label nFaces = mesh.nFaces();

    clockTime timer;
    IFstream dis(dictPath);
    dictionary dict(dis);
    autoPtr<decompositionMethod> method = decompositionMethod::New(dict);
    labelList cellToProc;
    if (weightCoeffs)
    {
        const scalarField w
        (
            linearCellWeights
            (
                mesh, (*weightCoeffs)[0], (*weightCoeffs)[1], (*weightCoeffs)[2]
            )
        );
        Info<< "CANDIDATE_WEIGHTS " << name << " min " << min(w) << " max " << max(w)
            << " mean " << average(w) << " sum " << sum(w) << endl;
        cellToProc = method->decompose(mesh, cc, w);
    }
    else
    {
        cellToProc = method->decompose(mesh, cc);
    }
    const label n = method->nDomains();
    const scalar tDecomp = timer.timeIncrement();

    labelList nCells(n, Zero), nProc(n, Zero), nBnd(n, Zero), nOff(n, Zero);
    labelList nInt(n, Zero), nWall(n, Zero), nWallCells(n, Zero);
    labelList shared(n*n, Zero);
    for (const label proci : cellToProc) ++nCells[proci];
    for (label facei = 0; facei < nInternal; ++facei)
    {
        const label a = cellToProc[own[facei]];
        const label b = cellToProc[nei[facei]];
        if (a != b)
        {
            ++nProc[a]; ++nProc[b];
            ++shared[a*n + b]; ++shared[b*n + a];
            if (a/ranksPerNode != b/ranksPerNode) { ++nOff[a]; ++nOff[b]; }
        }
        else
        {
            ++nInt[a];
        }
    }
    for (label facei = nInternal; facei < nFaces; ++facei)
    {
        ++nBnd[cellToProc[own[facei]]];
    }
    boolList isWallCell(mesh.nCells(), false);
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (isA<wallPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                const label celli = own[pp.start() + i];
                ++nWall[cellToProc[celli]];
                isWallCell[celli] = true;
            }
        }
    }
    forAll(isWallCell, celli)
    {
        if (isWallCell[celli]) ++nWallCells[cellToProc[celli]];
    }

    label maxCells = 0, maxProc = 0, maxNbr = 0, totProc = 0;
    {
        OFstream os(outDir/name + ".ranks");
        os  << "# candidate " << name << " dict " << dictPath;
        if (weightCoeffs)
        {
            os  << " weights A " << (*weightCoeffs)[0] << " B " << (*weightCoeffs)[1]
                << " C " << (*weightCoeffs)[2];
        }
        os  << nl
            << "# rank cells procFaces bndFaces neighbours offNodeFaces"
               " internalFaces wallFaces wallCells" << nl;
        for (label proci = 0; proci < n; ++proci)
        {
            label nNbr = 0;
            for (label j = 0; j < n; ++j) if (shared[proci*n + j]) ++nNbr;
            os  << proci << ' ' << nCells[proci] << ' ' << nProc[proci] << ' '
                << nBnd[proci] << ' ' << nNbr << ' ' << nOff[proci] << ' '
                << nInt[proci] << ' ' << nWall[proci] << ' ' << nWallCells[proci] << nl;
            maxCells = max(maxCells, nCells[proci]);
            maxProc = max(maxProc, nProc[proci]);
            maxNbr = max(maxNbr, nNbr);
            totProc += nProc[proci];
        }
    }
    Info<< "CANDIDATE " << name << " nDomains " << n
        << " maxCells " << maxCells << " maxProcFaces " << maxProc
        << " maxNeighbours " << maxNbr << " totalProcFaces " << totProc/2
        << " decompose_s " << tDecomp << endl;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("candidates", "file", "List of 'name dictPath [A B C]' lines");
    argList::addOption("outDir", "dir", "Directory for <name>.ranks files");
    argList::addOption("ranksPerNode", "label", "Ranks per node (default 104)");
    argList::addOption("nParallel", "label", "Candidates decomposed concurrently (default 1)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const fileName candFile(args.get<fileName>("candidates"));
    const fileName outDir(args.get<fileName>("outDir"));
    const label ranksPerNode = args.getOrDefault<label>("ranksPerNode", 104);
    const label nParallel = max(label(1), args.getOrDefault<label>("nParallel", 1));
    mkDir(outDir);

    // Build everything the children share before forking
    mesh.faceOwner();
    mesh.faceNeighbour();
    const pointField& cc = mesh.cellCentres();
    Info<< "Mesh: cells " << mesh.nCells() << " faces " << mesh.nFaces()
        << " internal " << mesh.nInternalFaces() << nl << endl;

    DynamicList<word> names;
    DynamicList<fileName> dicts;
    PtrList<scalarField> coeffs;
    {
        IFstream cis(candFile);
        std::string line;
        while (cis.good() && cis.getLine(line))
        {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            std::string nm, dp;
            if (!(iss >> nm >> dp)) continue;
            names.push_back(word(nm));
            dicts.push_back(fileName(dp));
            double a, b, c;
            coeffs.resize(names.size());
            if (iss >> a >> b >> c)
            {
                coeffs.set(names.size() - 1, new scalarField({a, b, c}));
            }
        }
    }

    std::map<pid_t, label> running;
    label nFailed = 0;
    auto reapOne = [&]()
    {
        int status = 0;
        const pid_t pid = ::waitpid(-1, &status, 0);
        if (pid <= 0) return;
        const label i = running[pid];
        running.erase(pid);
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        {
            ++nFailed;
            Info<< "CANDIDATE_FAILED " << names[i] << " status " << status << endl;
        }
    };

    forAll(names, i)
    {
        while (label(running.size()) >= nParallel) reapOne();
        Info.stdStream().flush();
        std::cout.flush();
        const pid_t pid = ::fork();
        if (pid < 0)
        {
            FatalErrorInFunction << "fork failed" << exit(FatalError);
        }
        if (pid == 0)
        {
            scoreCandidate
            (
                mesh, cc, names[i], dicts[i], outDir, ranksPerNode,
                coeffs.set(i) ? &coeffs[i] : nullptr
            );
            std::cout.flush();
            ::_exit(0);
        }
        running[pid] = i;
    }
    while (!running.empty()) reapOne();

    Info<< "\nCandidates " << names.size() << " failed " << nFailed << nl
        << "End\n" << endl;
    return 0;
}
