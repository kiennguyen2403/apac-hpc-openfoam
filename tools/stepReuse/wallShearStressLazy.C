/*---------------------------------------------------------------------------*\
    Replacement for functionObjects::wallShearStress::execute() and write().

    Stock execute() evaluates devReff() over the whole mesh at every time step
    and stores the wall values in a registered field (NO_WRITE) that only the
    function object's own write() reads. Here execute() does nothing and
    write() runs the stock execute() first, then the stock write().

    functionObjectList calls write() right after execute() for the same time
    step (and timeControl::write() runs execute() first if it has not run
    this time step), so write() sees the same state either way and writes the
    same values, provided execute() would run at every time step
    (executeControl timeStep, the default). With an executeInterval > 1 the
    stock object would write the values of its last execute() instead.
\*---------------------------------------------------------------------------*/

#include "wallShearStress.H"
#include "stepReuse.H"

namespace
{
    typedef bool (*memberFn)(Foam::functionObjects::wallShearStress*);

    memberFn stockExecute()
    {
        static const memberFn f = reinterpret_cast<memberFn>
        (
            Foam::stepReuse::next
            (
                "_ZN4Foam15functionObjects15wallShearStress7executeEv"
            )
        );
        return f;
    }

    memberFn stockWrite()
    {
        static const memberFn f = reinterpret_cast<memberFn>
        (
            Foam::stepReuse::next
            (
                "_ZN4Foam15functionObjects15wallShearStress5writeEv"
            )
        );
        return f;
    }
}


bool Foam::functionObjects::wallShearStress::execute()
{
    if (stepReuse::off("wss"))
    {
        return stockExecute()(this);
    }

    stepReuse::report("wss", "wallShearStress: evaluated in write() only");
    return true;
}


bool Foam::functionObjects::wallShearStress::write()
{
    if (!stepReuse::off("wss"))
    {
        stockExecute()(this);
    }

    return stockWrite()(this);
}


// ************************************************************************* //
