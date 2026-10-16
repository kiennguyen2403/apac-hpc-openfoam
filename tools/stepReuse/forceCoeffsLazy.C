/*---------------------------------------------------------------------------*\
    Replacement for functionObjects::forceCoeffs::reset(), calcForceCoeffs()
    and calcMomentCoeffs().

    Stock: every execute() sets the registered volVectorFields forceCoeff and
    momentCoeff to zero over the whole mesh, then to forceScaling*force()
    and momentScaling*moment() (full-mesh operations). They are only written
    when writeFields is on; the coefficients that are logged and written
    (Cf_, Cm_) come from the summed patch forces. Here the two fields are
    only computed when writeFields is on; the coefficients are computed
    exactly as in the stock code. Output identical.

    execute() is the stock code too: in the library it has reset() inlined,
    so it is replaced as well, to call the reset() here.

    STEP_REUSE_OFF=fco restores the stock behaviour.
\*---------------------------------------------------------------------------*/

#include "forceCoeffs.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "stepReuse.H"

namespace
{
    bool coeffFieldsNeeded(const bool writeFields)
    {
        return writeFields || Foam::stepReuse::off("fco");
    }
}


void Foam::functionObjects::forceCoeffs::reset()
{
    Cf_.reset();
    Cm_.reset();

    if (coeffFieldsNeeded(writeFields_))
    {
        forceCoeff() == Zero;
        momentCoeff() == Zero;
    }
    else
    {
        stepReuse::report
        (
            "fco",
            "forceCoeffs: coefficient fields only with writeFields"
        );
    }
}


void Foam::functionObjects::forceCoeffs::calcForceCoeffs()
{
    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const dimensionedScalar forceScaling
    (
        dimless/dimForce,
        scalar(1)/(Aref_*pDyn + SMALL)
    );

    const auto& coordSys = coordSysPtr_();

    // Calculate force coefficients
    if (coeffFieldsNeeded(writeFields_))
    {
        forceCoeff() = forceScaling*force();
    }

    Cf_.reset
    (
        forceScaling.value()*coordSys.localVector(sumPatchForcesP_),
        forceScaling.value()*coordSys.localVector(sumPatchForcesV_),
        forceScaling.value()*coordSys.localVector(sumInternalForces_)
    );
}


void Foam::functionObjects::forceCoeffs::calcMomentCoeffs()
{
    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const dimensionedScalar momentScaling
    (
        dimless/(dimForce*dimLength),
        scalar(1)/(Aref_*pDyn*lRef_ + SMALL)
    );

    const auto& coordSys = coordSysPtr_();

    // Calculate moment coefficients
    if (coeffFieldsNeeded(writeFields_))
    {
        momentCoeff() = momentScaling*moment();
    }

    Cm_.reset
    (
        momentScaling.value()*coordSys.localVector(sumPatchMomentsP_),
        momentScaling.value()*coordSys.localVector(sumPatchMomentsV_),
        momentScaling.value()*coordSys.localVector(sumInternalMoments_)
    );
}


bool Foam::functionObjects::forceCoeffs::execute()
{
    forces::calcForcesMoments();

    initialise();

    reset();

    Log << type() << " " << name() << " write:" << nl
        << "    " << "Coefficient"
        << tab << "Total"
        << tab << "Pressure"
        << tab << "Viscous"
        << tab << "Internal"
        << nl;

    calcForceCoeffs();

    calcMomentCoeffs();

    auto logValues = [](const word& name, const vector& coeff, Ostream& os)
    {
        os  << "    " << name << ":"
            << tab << coeff.x() + coeff.y() + coeff.z()
            << tab << coeff.x()
            << tab << coeff.y()
            << tab << coeff.z()
            << nl;
    };

    // Always setting all results
    for (const auto& iter : coeffs_.csorted())
    {
        const word& coeffName = iter.key();
        const auto& coeff = iter.val();

        // Vectors for x:pressure, y:viscous, z:internal
        const vector coeffValue = coeff.value(Cf_, Cm_);

        if (log && coeff.active_)
        {
            logValues(coeffName, coeffValue, Info);
        }

        setResult(coeffName, coeffValue.x() + coeffValue.y() + coeffValue.z());
        setResult(coeffName & "pressure", coeffValue.x());
        setResult(coeffName & "viscous", coeffValue.y());
        setResult(coeffName & "internal", coeffValue.z());
    }

    Log  << endl;

    return true;
}


// ************************************************************************* //
