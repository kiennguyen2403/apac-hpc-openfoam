/*---------------------------------------------------------------------------*\
    LD_PRELOAD replacement for Foam::solution::cache: answers as the stock
    function (the cache dictionary of fvSolution; the shipped case has
    cache { grad(U); }) and also yes for grad(p), grad(k) and grad(omega), as
    if they were listed there. Only gradScheme::grad asks for these names.

      grad(p)      computed several times per step while p is unchanged
                   (job 181366105: -2.30% at 4n)
      grad(k),
      grad(omega)  computed for CDkOmega, reused by the corrected laplacian
                   terms (job 181376795: ~0.6% at 4n)

    gradScheme::grad reuses a cached gradient only while it is up to date
    with its field (event number) and recomputes it otherwise, so a cached
    gradient is the field the uncached code would compute: bit-identical.
    Was "grad(p); grad(k); grad(omega);" added to the cache block of
    fvSolution; now the shipped fvSolution is used unchanged.

    STEP_REUSE_OFF=gcache: stock answers only.
\*---------------------------------------------------------------------------*/

#include "solution.H"
#include "objectRegistry.H"
#include "fixedIterSolvers.H"


bool Foam::solution::cache(const word& name) const
{
    typedef bool (*cacheFn)(const solution*, const word&);
    static const cacheFn stock = reinterpret_cast<cacheFn>
    (
        fixedIterSolvers::next("_ZNK4Foam8solution5cacheERKNS_4wordE")
    );

    if (stock(this, name))
    {
        return true;
    }

    if
    (
        (name == "grad(p)" || name == "grad(k)" || name == "grad(omega)")
     && !fixedIterSolvers::off("gcache")
     && cache_.getOrDefault<bool>("active", true)
    )
    {
        if (fixedIterSolvers::reporting())
        {
            fixedIterSolvers::report("cache: " + name + " cached in code");

            // Asked again while the field is stored: gradScheme::grad
            // reuses it, or updates it in place if its field has changed
            if (db().found(name))
            {
                fixedIterSolvers::report
                (
                    "cache: " + name + " stored when asked again"
                    " (reused while up to date)"
                );
            }
        }

        return true;
    }

    return false;
}
