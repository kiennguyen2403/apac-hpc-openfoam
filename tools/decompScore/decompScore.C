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

    Usage: decompScore -candidates <file: "name dictPath" per line>
                       -outDir <dir> [-ranksPerNode 104] [-nParallel 1]
    Output: <outDir>/<name>.ranks with columns
            rank cells procFaces bndFaces neighbours offNodeFaces
            internalFaces wallFaces wallCells
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

using namespace Foam;

// Decompose with one candidate dictionary and write its per-rank table
static void scoreCandidate
(
    const polyMesh& mesh,
    const pointField& cc,
    const word& name,
    const fileName& dictPath,
    const fileName& outDir,
    const label ranksPerNode
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
    const labelList cellToProc(method->decompose(mesh, cc));
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
        os  << "# candidate " << name << " dict " << dictPath << nl
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
    argList::addOption("candidates", "file", "List of 'name dictPath' lines");
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
    {
        IFstream cis(candFile);
        std::string line;
        while (cis.good() && cis.getLine(line))
        {
            if (line.empty() || line[0] == '#') continue;
            const auto sp = line.find_first_of(" \t");
            if (sp == std::string::npos) continue;
            names.push_back(word(line.substr(0, sp)));
            dicts.push_back(fileName(line.substr(line.find_first_not_of(" \t", sp))));
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
            scoreCandidate(mesh, cc, names[i], dicts[i], outDir, ranksPerNode);
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
