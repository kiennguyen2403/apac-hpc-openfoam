/*---------------------------------------------------------------------------*\
    gamgCensus: LD_PRELOAD instrumentation of the GAMG pressure solve, per
    rank and per level. No numerics change: every wrapper times the call and
    forwards it to the next definition (stock, or libstepReuse's smoother).

    STEP_REUSE_CENSUS=<dir>  enables it; each rank writes <dir>/rank<N>.txt
                             at exit (N from PMI_RANK)

    Inside GAMGSolver::solve, per level (key: nCells, nFaces of the matrix):
      self times (nested wrapped calls excluded) of
        GaussSeidelSmoother::smooth / scalarSmooth   (smoothing)
        lduMatrix::Amul                              (incl. GAMGSolver::scale)
        lduMatrix::residual
        GAMGSolver::scale                            (its own sums/update)
        lduMatrix::initMatrixInterfaces / updateMatrixInterfaces
      and MPI time (Allreduce-like, Wait/Test/Probe/Recv-like, posting of
      sends/receives, other collectives) while that level is the innermost
      one. "rest" is the rest of GAMGSolver::solve (V-cycle bookkeeping,
      restriction, prolongation, coarsest-level solve).
    Outside GAMGSolver::solve only the MPI totals ("outside").

    Step phases (gamgPhases section below): self wall time and self MPI time
    of the U solve, scalar solves (per field), the turbulence model, the
    Laplacian/convection/gradient/divDevReff operators (per field), the wall
    functions and each function object, plus the whole time step (between
    successive simpleControl::loop() calls).

    Everything is accumulated from the second simpleControl::loop() call on,
    i.e. over the time steps the wall-clock statistics measure (step 1 and
    start-up excluded).

    LD_PRELOAD order: this library first (it wraps symbols that
    libfusedVectorSolve and libstepReuse define).
\*---------------------------------------------------------------------------*/

#include "GAMGSolver.H"
#include "GaussSeidelSmoother.H"
#include "lduMatrix.H"
#include "processorLduInterface.H"
#include "simpleControl.H"
#include "fvMatrices.H"
#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "RASModel.H"
#include "eddyViscosity.H"
#include "kOmegaSSTBase.H"
#include "linearViscousStress.H"
#include "gaussConvectionScheme.H"
#include "laplacianScheme.H"
#include "gaussGrad.H"
#include "cellLimitedGrad.H"
#include "minmodGradientLimiter.H"
#include "nutUSpaldingWallFunctionFvPatchScalarField.H"
#include "omegaWallFunctionFvPatchScalarField.H"
#include "timeControlFunctionObject.H"

#include <mpi.h>
#include <dlfcn.h>
#include <sched.h>
#include <sys/resource.h>
#include <unistd.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

namespace
{

enum Cat
{
    GS, AMUL, RESID, SCALE, IFINIT, IFUPD, ALLRED, WAIT, POST, COLL, REST, NCAT
};
const char* catName[NCAT] =
{
    "gs", "amul", "resid", "scale", "ifInit", "ifUpd", "allreduce", "wait",
    "post", "coll", "rest"
};

struct Level
{
    long nCells = -1, nFaces = 0, nIfaces = 0, nIfaceFaces = 0;
    // addressing: neighbour - owner distance (mean, quantiles, fractions
    // above 1k/16k/64k cells) and mean jump between consecutive faces'
    // neighbours
    double dMean = 0, d50 = 0, d90 = 0, d99 = 0, f1k = 0, f16k = 0, f64k = 0;
    double jump = 0;
    std::vector<std::pair<int, long>> nbrs;   // processor interfaces

    double t[NCAT] = {};
    long n[NCAT] = {};
};

struct Phase
{
    double wall = 0, mpi = 0;
    long n = 0;
};

struct Census
{
    bool on = false;
    long loopCalls = 0;
    double mpiAcc = 0;                      // all MPI time so far
    double loopT0 = 0, loopM0 = 0, loopT = 0, loopM = 0;
    struct rusage ru0 = {}, ru1 = {};       // at the 2nd and the last loop()
    int threads = -1;                       // at the 3rd loop()
    std::map<std::string, Phase> phases;
    double phChild = 0, phChildMpi = 0;
    std::string dir;
    std::map<std::pair<long, long>, Level> levels;
    Level rest;             // rest of GAMGSolver::solve
    Level outside;          // MPI outside GAMGSolver::solve
    Level* cur = nullptr;   // innermost level inside GAMG
    int depth = 0;          // GAMGSolver::solve nesting
    long nSolves = 0;
    double tSolve = 0, child = 0;

    Census()
    {
        const char* e = std::getenv("STEP_REUSE_CENSUS");
        if (e && *e) { on = true; dir = e; }
    }

    ~Census()
    {
        if (!on) return;
        const char* r = std::getenv("PMI_RANK");
        const std::string f = dir + "/rank" + (r ? r : "x") + ".txt";
        FILE* os = std::fopen(f.c_str(), "w");
        if (!os) return;
        std::fprintf(os, "# rank %s gamgSolves %ld gamgSolveTime %.6f\n",
            r ? r : "x", nSolves, tSolve);
        char host[256] = {};
        ::gethostname(host, sizeof(host) - 1);
        std::fprintf(os, "host %s %d\n", host, ::sched_getcpu());
        std::fprintf(os, "# level nCells nFaces nIfaces nIfaceFaces");
        for (int c = 0; c < NCAT; ++c)
            std::fprintf(os, " t_%s n_%s", catName[c], catName[c]);
        std::fprintf(os, "\n");
        auto line = [&](const char* name, const Level& L)
        {
            std::fprintf(os, "%s %ld %ld %ld %ld", name,
                L.nCells, L.nFaces, L.nIfaces, L.nIfaceFaces);
            for (int c = 0; c < NCAT; ++c)
                std::fprintf(os, " %.6f %ld", L.t[c], L.n[c]);
            std::fprintf(os, "\n");
        };
        for (auto it = levels.rbegin(); it != levels.rend(); ++it)
            line("level", it->second);
        for (auto it = levels.rbegin(); it != levels.rend(); ++it)
        {
            const Level& L = it->second;
            std::fprintf(os,
                "addr %ld %ld %.3f %.0f %.0f %.0f %.6f %.6f %.6f %.3f\n",
                L.nCells, L.nFaces, L.dMean, L.d50, L.d90, L.d99,
                L.f1k, L.f16k, L.f64k, L.jump);
        }
        if (!levels.empty())
        {
            const Level& L0 = levels.rbegin()->second;
            for (const auto& nb : L0.nbrs)
            {
                std::fprintf(os, "nbr %ld %d %ld\n", L0.nCells, nb.first, nb.second);
            }
        }
        line("rest", rest);
        line("outside", outside);
        std::fprintf(os, "steps %ld %.6f %.6f\n",
            loopCalls - 2, loopT - loopT0, loopM - loopM0);
        auto tv = [](const timeval& a, const timeval& b)
            { return (b.tv_sec - a.tv_sec) + 1e-6*(b.tv_usec - a.tv_usec); };
        std::fprintf(os, "rusage %.6f %.6f %ld %ld %ld %ld threads %d\n",
            tv(ru0.ru_utime, ru1.ru_utime), tv(ru0.ru_stime, ru1.ru_stime),
            ru1.ru_minflt - ru0.ru_minflt, ru1.ru_majflt - ru0.ru_majflt,
            ru1.ru_nvcsw - ru0.ru_nvcsw, ru1.ru_nivcsw - ru0.ru_nivcsw, threads);
        for (const auto& p : phases)
        {
            std::fprintf(os, "phase %s %.6f %.6f %ld\n", p.first.c_str(),
                p.second.wall, p.second.mpi, p.second.n);
        }
        std::fclose(os);
    }
};

Census census;

int threadCount()
{
    FILE* f = std::fopen("/proc/self/status", "r");
    if (!f) return -1;
    char line[256];
    int n = -1;
    while (std::fgets(line, sizeof(line), f))
    {
        if (std::sscanf(line, "Threads: %d", &n) == 1) break;
    }
    std::fclose(f);
    return n;
}

inline double now()
{
    return std::chrono::duration<double>
    (
        std::chrono::steady_clock::now().time_since_epoch()
    ).count();
}

void* nextSym(const char* name)
{
    void* p = ::dlsym(RTLD_NEXT, name);
    if (!p)
    {
        std::fprintf(stderr, "gamgCensus: no next definition of %s\n", name);
        std::abort();
    }
    return p;
}

Level& levelOf
(
    const Foam::lduMatrix& m,
    const Foam::lduInterfaceFieldPtrsList& interfaces
)
{
    const long nCells = m.diag().size();
    const long nFaces = m.lduAddr().lowerAddr().size();
    Level& L = census.levels[{nCells, nFaces}];
    if (L.nCells < 0)
    {
        L.nCells = nCells;
        L.nFaces = nFaces;
        forAll(interfaces, i)
        {
            if (interfaces.set(i))
            {
                ++L.nIfaces;
                const Foam::lduInterface& li = interfaces[i].interface();
                L.nIfaceFaces += li.faceCells().size();
                const auto* pli =
                    dynamic_cast<const Foam::processorLduInterface*>(&li);
                if (pli)
                {
                    L.nbrs.emplace_back(pli->neighbProcNo(), li.faceCells().size());
                }
            }
        }

        const Foam::labelUList& l = m.lduAddr().lowerAddr();
        const Foam::labelUList& u = m.lduAddr().upperAddr();
        std::vector<long> d(nFaces);
        double sum = 0, jump = 0;
        long n1k = 0, n16k = 0, n64k = 0;
        for (long f = 0; f < nFaces; ++f)
        {
            d[f] = long(u[f]) - long(l[f]);
            sum += d[f];
            n1k += d[f] > 1024;
            n16k += d[f] > 16384;
            n64k += d[f] > 65536;
            if (f) jump += std::labs(long(u[f]) - long(u[f - 1]));
        }
        if (nFaces > 1)
        {
            L.dMean = sum/nFaces;
            L.f1k = double(n1k)/nFaces;
            L.f16k = double(n16k)/nFaces;
            L.f64k = double(n64k)/nFaces;
            L.jump = jump/(nFaces - 1);
            auto q = [&](double p)
            {
                const long k = std::min(nFaces - 1, long(p*nFaces));
                std::nth_element(d.begin(), d.begin() + k, d.end());
                return double(d[k]);
            };
            L.d50 = q(0.5);
            L.d90 = q(0.9);
            L.d99 = q(0.99);
        }
    }
    return L;
}

// Time f() as category c of level L (self time: nested timed calls excluded)
template<class F>
inline void timed(Cat c, Level& L, F&& f)
{
    const double saved = census.child;
    Level* savedCur = census.cur;
    census.child = 0;
    census.cur = &L;
    const double t0 = now();
    f();
    const double el = now() - t0;
    L.t[c] += el - census.child;
    ++L.n[c];
    census.child = saved + el;
    census.cur = savedCur;
}

inline bool measuring()
{
    return census.on && census.loopCalls >= 2;
}

inline bool active()
{
    return measuring() && census.depth > 0;
}

template<class F>
inline int timedMpi(Cat c, F&& f)
{
    if (!census.on)
    {
        return f();
    }
    const double t0 = now();
    const int ret = f();
    const double el = now() - t0;
    census.mpiAcc += el;
    if (measuring())
    {
        Level& L = census.depth > 0 ? *census.cur : census.outside;
        L.t[c] += el;
        ++L.n[c];
        census.child += el;
    }
    return ret;
}

// Time f() as phase "name": self wall time and self MPI time (nested phases
// excluded)
template<class F>
inline void phased(const std::string& name, F&& f)
{
    if (!measuring())
    {
        f();
        return;
    }
    const double sC = census.phChild, sM = census.phChildMpi;
    census.phChild = 0;
    census.phChildMpi = 0;
    const double t0 = now(), m0 = census.mpiAcc;
    f();
    const double el = now() - t0, em = census.mpiAcc - m0;
    Phase& P = census.phases[name];
    P.wall += el - census.phChild;
    P.mpi += em - census.phChildMpi;
    ++P.n;
    census.phChild = sC + el;
    census.phChildMpi = sM + em;
}

} // End anonymous namespace


namespace Foam
{

solverPerformance GAMGSolver::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    typedef solverPerformance (*fn)
        (const GAMGSolver*, scalarField&, const scalarField&, direction);
    static const fn stock = reinterpret_cast<fn>
        (nextSym("_ZNK4Foam10GAMGSolver5solveERNS_5FieldIdEERKS2_h"));

    if (!measuring())
    {
        return stock(this, psi, source, cmpt);
    }

    solverPerformance perf;
    ++census.depth;
    const double t0 = now();
    timed(REST, census.rest, [&]{ perf = stock(this, psi, source, cmpt); });
    census.tSolve += now() - t0;
    ++census.nSolves;
    --census.depth;
    return perf;
}


void GAMGSolver::scale
(
    solveScalarField& field,
    solveScalarField& Acf,
    const lduMatrix& A,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const solveScalarField& source,
    const direction cmpt
) const
{
    typedef void (*fn)
    (
        const GAMGSolver*, solveScalarField&, solveScalarField&,
        const lduMatrix&, const FieldField<Field, scalar>&,
        const lduInterfaceFieldPtrsList&, const solveScalarField&, direction
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZNK4Foam10GAMGSolver5scaleERNS_5FieldIdEES3_RKNS_9lduMatrixERKNS_"
            "10FieldFieldIS1_dEERKNS_8UPtrListIKNS_17lduInterfaceFieldEEERKS2_h"
        )
    );

    if (!active())
    {
        stock(this, field, Acf, A, interfaceLevelBouCoeffs, interfaceLevel,
            source, cmpt);
        return;
    }

    timed(SCALE, levelOf(A, interfaceLevel), [&]
    {
        stock(this, field, Acf, A, interfaceLevelBouCoeffs, interfaceLevel,
            source, cmpt);
    });
}


void GaussSeidelSmoother::smooth
(
    solveScalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    typedef void (*fn)
    (
        const GaussSeidelSmoother*, solveScalarField&, const scalarField&,
        direction, label
    );
    static const fn stock = reinterpret_cast<fn>
        (nextSym("_ZNK4Foam19GaussSeidelSmoother6smoothERNS_5FieldIdEERKS2_hi"));

    if (!active())
    {
        stock(this, psi, source, cmpt, nSweeps);
        return;
    }

    timed(GS, levelOf(matrix_, interfaces_), [&]
    {
        stock(this, psi, source, cmpt, nSweeps);
    });
}


void GaussSeidelSmoother::scalarSmooth
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    typedef void (*fn)
    (
        const GaussSeidelSmoother*, solveScalarField&,
        const solveScalarField&, direction, label
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        ("_ZNK4Foam19GaussSeidelSmoother12scalarSmoothERNS_5FieldIdEERKS2_hi")
    );

    if (!active())
    {
        stock(this, psi, source, cmpt, nSweeps);
        return;
    }

    timed(GS, levelOf(matrix_, interfaces_), [&]
    {
        stock(this, psi, source, cmpt, nSweeps);
    });
}


void lduMatrix::Amul
(
    solveScalarField& Apsi,
    const tmp<solveScalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    typedef void (*fn)
    (
        const lduMatrix*, solveScalarField&, const tmp<solveScalarField>&,
        const FieldField<Field, scalar>&, const lduInterfaceFieldPtrsList&,
        direction
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZNK4Foam9lduMatrix4AmulERNS_5FieldIdEERKNS_3tmpIS2_EERKNS_"
            "10FieldFieldIS1_dEERKNS_8UPtrListIKNS_17lduInterfaceFieldEEEh"
        )
    );

    if (!active())
    {
        stock(this, Apsi, tpsi, interfaceBouCoeffs, interfaces, cmpt);
        return;
    }

    timed(AMUL, levelOf(*this, interfaces), [&]
    {
        stock(this, Apsi, tpsi, interfaceBouCoeffs, interfaces, cmpt);
    });
}


void lduMatrix::residual
(
    solveScalarField& rA,
    const solveScalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    typedef void (*fn)
    (
        const lduMatrix*, solveScalarField&, const solveScalarField&,
        const scalarField&, const FieldField<Field, scalar>&,
        const lduInterfaceFieldPtrsList&, direction
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZNK4Foam9lduMatrix8residualERNS_5FieldIdEERKS2_S5_RKNS_"
            "10FieldFieldIS1_dEERKNS_8UPtrListIKNS_17lduInterfaceFieldEEEh"
        )
    );

    if (!active())
    {
        stock(this, rA, psi, source, interfaceBouCoeffs, interfaces, cmpt);
        return;
    }

    timed(RESID, levelOf(*this, interfaces), [&]
    {
        stock(this, rA, psi, source, interfaceBouCoeffs, interfaces, cmpt);
    });
}


void lduMatrix::initMatrixInterfaces
(
    const bool add,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt
) const
{
    typedef void (*fn)
    (
        const lduMatrix*, bool, const FieldField<Field, scalar>&,
        const lduInterfaceFieldPtrsList&, const solveScalarField&,
        solveScalarField&, direction
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZNK4Foam9lduMatrix20initMatrixInterfacesEbRKNS_10FieldFieldINS_"
            "5FieldEdEERKNS_8UPtrListIKNS_17lduInterfaceFieldEEERKNS2_IdEERSC_h"
        )
    );

    if (!active())
    {
        stock(this, add, interfaceCoeffs, interfaces, psiif, result, cmpt);
        return;
    }

    timed(IFINIT, levelOf(*this, interfaces), [&]
    {
        stock(this, add, interfaceCoeffs, interfaces, psiif, result, cmpt);
    });
}


void lduMatrix::updateMatrixInterfaces
(
    const bool add,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt,
    const label startRequest
) const
{
    typedef void (*fn)
    (
        const lduMatrix*, bool, const FieldField<Field, scalar>&,
        const lduInterfaceFieldPtrsList&, const solveScalarField&,
        solveScalarField&, direction, label
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZNK4Foam9lduMatrix22updateMatrixInterfacesEbRKNS_10FieldFieldINS_"
            "5FieldEdEERKNS_8UPtrListIKNS_17lduInterfaceFieldEEERKNS2_IdEERSC_hi"
        )
    );

    if (!active())
    {
        stock
        (
            this, add, interfaceCoeffs, interfaces, psiif, result, cmpt,
            startRequest
        );
        return;
    }

    timed(IFUPD, levelOf(*this, interfaces), [&]
    {
        stock
        (
            this, add, interfaceCoeffs, interfaces, psiif, result, cmpt,
            startRequest
        );
    });
}

} // End namespace Foam


// Step phases

namespace Foam
{

bool simpleControl::loop()
{
    typedef bool (*fn)(simpleControl*);
    static const fn stock =
        reinterpret_cast<fn>(nextSym("_ZN4Foam13simpleControl4loopEv"));

    if (census.on)
    {
        ++census.loopCalls;
        const double t = now();
        if (census.loopCalls == 2)
        {
            census.loopT0 = t;
            census.loopM0 = census.mpiAcc;
            ::getrusage(RUSAGE_SELF, &census.ru0);
        }
        if (census.loopCalls == 3)
        {
            census.threads = threadCount();
        }
        census.loopT = t;
        census.loopM = census.mpiAcc;
        ::getrusage(RUSAGE_SELF, &census.ru1);
    }
    return stock(this);
}


template<>
SolverPerformance<vector> fvMatrix<vector>::solveSegregated
(
    const dictionary& solverControls
)
{
    typedef SolverPerformance<vector> (*fn)(fvMatrix<vector>*, const dictionary&);
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym("_ZN4Foam8fvMatrixINS_6VectorIdEEE15solveSegregatedERKNS_10dictionaryE")
    );
    SolverPerformance<vector> perf;
    phased("solve:" + psi().name(), [&]{ perf = stock(this, solverControls); });
    return perf;
}


template<>
solverPerformance fvMatrix<scalar>::solveSegregated
(
    const dictionary& solverControls
)
{
    typedef solverPerformance (*fn)(fvMatrix<scalar>*, const dictionary&);
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym("_ZN4Foam8fvMatrixIdE15solveSegregatedERKNS_10dictionaryE")
    );
    solverPerformance perf;
    phased("solve:" + psi().name(), [&]{ perf = stock(this, solverControls); });
    return perf;
}


typedef eddyViscosity<RASModel<IncompressibleTurbulenceModel<transportModel>>>
    censusEddyViscosityRAS;

template<>
void kOmegaSSTBase<censusEddyViscosityRAS>::correct()
{
    typedef void (*fn)(kOmegaSSTBase<censusEddyViscosityRAS>*);
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZN4Foam13kOmegaSSTBaseINS_13eddyViscosityINS_8RASModelINS_"
            "29IncompressibleTurbulenceModelINS_14transportModelEEEEEEEE7correctEv"
        )
    );
    phased("turbulence", [&]{ stock(this); });
}


typedef RASModel<IncompressibleTurbulenceModel<transportModel>>
    censusIncompressibleRAS;

template<>
tmp<fvVectorMatrix>
linearViscousStress<censusIncompressibleRAS>::divDevRhoReff
(
    volVectorField& U
) const
{
    typedef tmp<fvVectorMatrix> (*fn)
        (const linearViscousStress<censusIncompressibleRAS>*, volVectorField&);
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZNK4Foam19linearViscousStressINS_8RASModelINS_"
            "29IncompressibleTurbulenceModelINS_14transportModelEEEEEE"
            "13divDevRhoReffERNS_14GeometricFieldINS_6VectorIdEENS_"
            "12fvPatchFieldENS_7volMeshEEE"
        )
    );
    tmp<fvVectorMatrix> t;
    phased("divDevReff", [&]{ t = stock(this, U); });
    return t;
}


template<>
tmp<fvMatrix<vector>> fv::gaussConvectionScheme<vector>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const volVectorField& vf
) const
{
    typedef tmp<fvMatrix<vector>> (*fn)
    (
        const fv::gaussConvectionScheme<vector>*, const surfaceScalarField&,
        const volVectorField&
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZNK4Foam2fv21gaussConvectionSchemeINS_6VectorIdEEE6fvmDivERKNS_"
            "14GeometricFieldIdNS_13fvsPatchFieldENS_11surfaceMeshEEERKNS5_IS3_"
            "NS_12fvPatchFieldENS_7volMeshEEE"
        )
    );
    tmp<fvMatrix<vector>> t;
    phased("conv:" + vf.name(), [&]{ t = stock(this, faceFlux, vf); });
    return t;
}


template<>
tmp<fvMatrix<scalar>> fv::laplacianScheme<scalar, scalar>::fvmLaplacian
(
    const volScalarField& gamma,
    const volScalarField& vf
)
{
    typedef tmp<fvMatrix<scalar>> (*fn)
    (
        fv::laplacianScheme<scalar, scalar>*, const volScalarField&,
        const volScalarField&
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZN4Foam2fv15laplacianSchemeIddE12fvmLaplacianERKNS_"
            "14GeometricFieldIdNS_12fvPatchFieldENS_7volMeshEEES8_"
        )
    );
    tmp<fvMatrix<scalar>> t;
    phased("lap:" + vf.name(), [&]{ t = stock(this, gamma, vf); });
    return t;
}


template<>
tmp<fvMatrix<vector>> fv::laplacianScheme<vector, scalar>::fvmLaplacian
(
    const volScalarField& gamma,
    const volVectorField& vf
)
{
    typedef tmp<fvMatrix<vector>> (*fn)
    (
        fv::laplacianScheme<vector, scalar>*, const volScalarField&,
        const volVectorField&
    );
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZN4Foam2fv15laplacianSchemeINS_6VectorIdEEdE12fvmLaplacianERKNS_"
            "14GeometricFieldIdNS_12fvPatchFieldENS_7volMeshEEERKNS5_IS3_S6_S7_EE"
        )
    );
    tmp<fvMatrix<vector>> t;
    phased("lap:" + vf.name(), [&]{ t = stock(this, gamma, vf); });
    return t;
}


// Gradients: returning and in-place versions, one phase per field
#define CENSUS_GRAD(Scheme, Type, mangledReturn, mangledInPlace)               \
                                                                               \
template<>                                                                     \
tmp<GeometricField<outerProduct<vector, Type>::type, fvPatchField, volMesh>>   \
Scheme::calcGrad                                                               \
(                                                                              \
    const GeometricField<Type, fvPatchField, volMesh>& vsf,                    \
    const word& name                                                           \
) const                                                                        \
{                                                                              \
    typedef GeometricField<outerProduct<vector, Type>::type, fvPatchField,     \
        volMesh> GradFieldType;                                                \
    typedef tmp<GradFieldType> (*fn)                                           \
    (                                                                          \
        const Scheme*, const GeometricField<Type, fvPatchField, volMesh>&,     \
        const word&                                                            \
    );                                                                         \
    static const fn stock = reinterpret_cast<fn>(nextSym(mangledReturn));      \
    tmp<GradFieldType> t;                                                      \
    phased("grad:" + vsf.name(), [&]{ t = stock(this, vsf, name); });          \
    return t;                                                                  \
}                                                                              \
                                                                               \
template<>                                                                     \
void Scheme::calcGrad                                                          \
(                                                                              \
    GeometricField<outerProduct<vector, Type>::type, fvPatchField, volMesh>& g,\
    const GeometricField<Type, fvPatchField, volMesh>& vsf                     \
) const                                                                        \
{                                                                              \
    typedef void (*fn)                                                         \
    (                                                                          \
        const Scheme*,                                                         \
        GeometricField<outerProduct<vector, Type>::type, fvPatchField,         \
            volMesh>&,                                                         \
        const GeometricField<Type, fvPatchField, volMesh>&                     \
    );                                                                         \
    static const fn stock = reinterpret_cast<fn>(nextSym(mangledInPlace));     \
    phased("grad:" + vsf.name(), [&]{ stock(this, g, vsf); });                 \
}

CENSUS_GRAD
(
    fv::gaussGrad<scalar>, scalar,
    "_ZNK4Foam2fv9gaussGradIdE8calcGradERKNS_14GeometricFieldIdNS_"
    "12fvPatchFieldENS_7volMeshEEERKNS_4wordE",
    "_ZNK4Foam2fv9gaussGradIdE8calcGradERNS_14GeometricFieldINS_6VectorIdEENS_"
    "12fvPatchFieldENS_7volMeshEEERKNS3_IdS6_S7_EE"
)

CENSUS_GRAD
(
    fv::gaussGrad<vector>, vector,
    "_ZNK4Foam2fv9gaussGradINS_6VectorIdEEE8calcGradERKNS_14GeometricFieldIS3_"
    "NS_12fvPatchFieldENS_7volMeshEEERKNS_4wordE",
    "_ZNK4Foam2fv9gaussGradINS_6VectorIdEEE8calcGradERNS_14GeometricFieldINS_"
    "6TensorIdEENS_12fvPatchFieldENS_7volMeshEEERKNS5_IS3_S8_S9_EE"
)

typedef fv::cellLimitedGrad<vector, fv::gradientLimiters::minmod>
    censusCellLimitedVector;

CENSUS_GRAD
(
    censusCellLimitedVector, vector,
    "_ZNK4Foam2fv15cellLimitedGradINS_6VectorIdEENS0_16gradientLimiters"
    "6minmodEE8calcGradERKNS_14GeometricFieldIS3_NS_12fvPatchFieldENS_"
    "7volMeshEEERKNS_4wordE",
    "_ZNK4Foam2fv15cellLimitedGradINS_6VectorIdEENS0_16gradientLimiters"
    "6minmodEE8calcGradERNS_14GeometricFieldINS_6TensorIdEENS_"
    "12fvPatchFieldENS_7volMeshEEERKNS7_IS3_SA_SB_EE"
)


tmp<scalarField> nutUSpaldingWallFunctionFvPatchScalarField::calcNut() const
{
    typedef tmp<scalarField> (*fn)
        (const nutUSpaldingWallFunctionFvPatchScalarField*);
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym("_ZNK4Foam42nutUSpaldingWallFunctionFvPatchScalarField7calcNutEv")
    );
    tmp<scalarField> t;
    phased("wall:nut", [&]{ t = stock(this); });
    return t;
}


void omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    typedef void (*fn)(omegaWallFunctionFvPatchScalarField*);
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym("_ZN4Foam35omegaWallFunctionFvPatchScalarField12updateCoeffsEv")
    );
    phased("wall:omegaUpdate", [&]{ stock(this); });
}


void omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    typedef void (*fn)
        (omegaWallFunctionFvPatchScalarField*, fvMatrix<scalar>&);
    static const fn stock = reinterpret_cast<fn>
    (
        nextSym
        (
            "_ZN4Foam35omegaWallFunctionFvPatchScalarField16manipulateMatrixERNS_"
            "8fvMatrixIdEE"
        )
    );
    phased("wall:omegaMatrix", [&]{ stock(this, matrix); });
}


bool functionObjects::timeControl::execute()
{
    typedef bool (*fn)(functionObjects::timeControl*);
    static const fn stock = reinterpret_cast<fn>
        (nextSym("_ZN4Foam15functionObjects11timeControl7executeEv"));
    bool ok = true;
    phased("fo:" + name(), [&]{ ok = stock(this); });
    return ok;
}


bool functionObjects::timeControl::write()
{
    typedef bool (*fn)(functionObjects::timeControl*);
    static const fn stock = reinterpret_cast<fn>
        (nextSym("_ZN4Foam15functionObjects11timeControl5writeEv"));
    bool ok = true;
    phased("fo:" + name(), [&]{ ok = stock(this); });
    return ok;
}

} // End namespace Foam


// MPI: every communication call the Pstream library makes, timed and
// forwarded to the MPI library
#define CENSUS_MPI(cat, name, params, args)                                    \
int name params                                                                \
{                                                                              \
    typedef int (*fn) params;                                                  \
    static const fn stock = reinterpret_cast<fn>(nextSym(#name));              \
    return timedMpi(cat, [&]{ return stock args; });                           \
}

extern "C"
{

CENSUS_MPI(ALLRED, MPI_Allreduce,
    (const void* s, void* r, int c, MPI_Datatype t, MPI_Op o, MPI_Comm m),
    (s, r, c, t, o, m))
CENSUS_MPI(ALLRED, MPI_Iallreduce,
    (const void* s, void* r, int c, MPI_Datatype t, MPI_Op o, MPI_Comm m,
     MPI_Request* q),
    (s, r, c, t, o, m, q))
CENSUS_MPI(ALLRED, MPI_Reduce,
    (const void* s, void* r, int c, MPI_Datatype t, MPI_Op o, int root,
     MPI_Comm m),
    (s, r, c, t, o, root, m))
CENSUS_MPI(ALLRED, MPI_Ireduce,
    (const void* s, void* r, int c, MPI_Datatype t, MPI_Op o, int root,
     MPI_Comm m, MPI_Request* q),
    (s, r, c, t, o, root, m, q))
CENSUS_MPI(ALLRED, MPI_Scan,
    (const void* s, void* r, int c, MPI_Datatype t, MPI_Op o, MPI_Comm m),
    (s, r, c, t, o, m))
CENSUS_MPI(ALLRED, MPI_Exscan,
    (const void* s, void* r, int c, MPI_Datatype t, MPI_Op o, MPI_Comm m),
    (s, r, c, t, o, m))

CENSUS_MPI(WAIT, MPI_Wait, (MPI_Request* q, MPI_Status* st), (q, st))
CENSUS_MPI(WAIT, MPI_Waitall,
    (int n, MPI_Request q[], MPI_Status st[]), (n, q, st))
CENSUS_MPI(WAIT, MPI_Waitany,
    (int n, MPI_Request q[], int* i, MPI_Status* st), (n, q, i, st))
CENSUS_MPI(WAIT, MPI_Waitsome,
    (int n, MPI_Request q[], int* o, int is[], MPI_Status* st),
    (n, q, o, is, st))
CENSUS_MPI(WAIT, MPI_Test,
    (MPI_Request* q, int* f, MPI_Status* st), (q, f, st))
CENSUS_MPI(WAIT, MPI_Testall,
    (int n, MPI_Request q[], int* f, MPI_Status* st), (n, q, f, st))
CENSUS_MPI(WAIT, MPI_Testsome,
    (int n, MPI_Request q[], int* o, int is[], MPI_Status* st),
    (n, q, o, is, st))
CENSUS_MPI(WAIT, MPI_Probe,
    (int src, int tag, MPI_Comm m, MPI_Status* st), (src, tag, m, st))
CENSUS_MPI(WAIT, MPI_Iprobe,
    (int src, int tag, MPI_Comm m, int* f, MPI_Status* st),
    (src, tag, m, f, st))
CENSUS_MPI(WAIT, MPI_Recv,
    (void* b, int c, MPI_Datatype t, int src, int tag, MPI_Comm m,
     MPI_Status* st),
    (b, c, t, src, tag, m, st))

CENSUS_MPI(POST, MPI_Isend,
    (const void* b, int c, MPI_Datatype t, int d, int tag, MPI_Comm m,
     MPI_Request* q),
    (b, c, t, d, tag, m, q))
CENSUS_MPI(POST, MPI_Issend,
    (const void* b, int c, MPI_Datatype t, int d, int tag, MPI_Comm m,
     MPI_Request* q),
    (b, c, t, d, tag, m, q))
CENSUS_MPI(POST, MPI_Irecv,
    (void* b, int c, MPI_Datatype t, int src, int tag, MPI_Comm m,
     MPI_Request* q),
    (b, c, t, src, tag, m, q))
CENSUS_MPI(POST, MPI_Send,
    (const void* b, int c, MPI_Datatype t, int d, int tag, MPI_Comm m),
    (b, c, t, d, tag, m))
CENSUS_MPI(POST, MPI_Ssend,
    (const void* b, int c, MPI_Datatype t, int d, int tag, MPI_Comm m),
    (b, c, t, d, tag, m))
CENSUS_MPI(POST, MPI_Bsend,
    (const void* b, int c, MPI_Datatype t, int d, int tag, MPI_Comm m),
    (b, c, t, d, tag, m))

CENSUS_MPI(COLL, MPI_Barrier, (MPI_Comm m), (m))
CENSUS_MPI(COLL, MPI_Ibarrier, (MPI_Comm m, MPI_Request* q), (m, q))
CENSUS_MPI(COLL, MPI_Bcast,
    (void* b, int c, MPI_Datatype t, int root, MPI_Comm m),
    (b, c, t, root, m))
CENSUS_MPI(COLL, MPI_Allgather,
    (const void* s, int sc, MPI_Datatype st, void* r, int rc,
     MPI_Datatype rt, MPI_Comm m),
    (s, sc, st, r, rc, rt, m))
CENSUS_MPI(COLL, MPI_Iallgather,
    (const void* s, int sc, MPI_Datatype st, void* r, int rc,
     MPI_Datatype rt, MPI_Comm m, MPI_Request* q),
    (s, sc, st, r, rc, rt, m, q))
CENSUS_MPI(COLL, MPI_Alltoall,
    (const void* s, int sc, MPI_Datatype st, void* r, int rc,
     MPI_Datatype rt, MPI_Comm m),
    (s, sc, st, r, rc, rt, m))
CENSUS_MPI(COLL, MPI_Gather,
    (const void* s, int sc, MPI_Datatype st, void* r, int rc,
     MPI_Datatype rt, int root, MPI_Comm m),
    (s, sc, st, r, rc, rt, root, m))
CENSUS_MPI(COLL, MPI_Igather,
    (const void* s, int sc, MPI_Datatype st, void* r, int rc,
     MPI_Datatype rt, int root, MPI_Comm m, MPI_Request* q),
    (s, sc, st, r, rc, rt, root, m, q))
CENSUS_MPI(COLL, MPI_Gatherv,
    (const void* s, int sc, MPI_Datatype st, void* r, const int rcs[],
     const int ds[], MPI_Datatype rt, int root, MPI_Comm m),
    (s, sc, st, r, rcs, ds, rt, root, m))
CENSUS_MPI(COLL, MPI_Scatter,
    (const void* s, int sc, MPI_Datatype st, void* r, int rc,
     MPI_Datatype rt, int root, MPI_Comm m),
    (s, sc, st, r, rc, rt, root, m))
CENSUS_MPI(COLL, MPI_Iscatter,
    (const void* s, int sc, MPI_Datatype st, void* r, int rc,
     MPI_Datatype rt, int root, MPI_Comm m, MPI_Request* q),
    (s, sc, st, r, rc, rt, root, m, q))
CENSUS_MPI(COLL, MPI_Scatterv,
    (const void* s, const int scs[], const int ds[], MPI_Datatype st,
     void* r, int rc, MPI_Datatype rt, int root, MPI_Comm m),
    (s, scs, ds, st, r, rc, rt, root, m))

} // extern "C"


// ************************************************************************* //
