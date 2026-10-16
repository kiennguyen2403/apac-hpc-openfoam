/*---------------------------------------------------------------------------*\
    Replacement for fv::laplacianScheme<Type, scalar>::fvmLaplacian(gamma, vf)
    with a volScalarField gamma (Type scalar and vector; the call every
    fvm::laplacian(volScalarField, vf) makes: U, k, omega and p equations).

    Stock (Gauss linear corrected): interpolate gamma to the faces (a face
    field), gammaMagSf = gamma*magSf (a face field), upper = deltaCoeffs
    *gammaMagSf, negSumDiag, patch coefficients; then the non-orthogonal
    correction: dotInterpolate(nonOrthCorrectionVectors, grad(vf)) (a face
    field), gammaMagSf*correction (a face field), fvc::div of it (surface
    integral, a cell field), source -= V*div.

    Here one loop over the internal faces does the interpolation, upper, the
    diagonal (owner then neighbour, in face order, as negSumDiag) and the
    correction flux added to owner and subtracted from neighbour (in face
    order, then the patches, as surfaceIntegrate). Each value is the same
    expression as in the stock code (the interpolation lambda*(a - b) + b is
    left to the same contraction as there), so the result is bit-identical.
    The patch values use the stock functions (patchNeighbourField, lerp,
    gradient{Internal,Boundary}Coeffs). The face flux correction is stored
    for fields that need it (p) as in the stock code.

    Only for gaussLaplacianScheme with linear gamma interpolation and the
    corrected or uncorrected snGrad; otherwise, or with STEP_REUSE_OFF=lap,
    the stock body runs.
\*---------------------------------------------------------------------------*/

#include "gaussLaplacianScheme.H"
#include "gradScheme.H"
#include "fvMatrices.H"
#include "stepReuse.H"

#include <memory>

namespace Foam
{
namespace stepReuseDetail
{

template<class Type>
tmp<fvMatrix<Type>> fusedGaussLaplacian
(
    const surfaceInterpolationScheme<scalar>& gammaInterp,
    const fv::snGradScheme<Type>& snGrad,
    const volScalarField& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<surfaceScalarField> tlambdas = gammaInterp.weights(gamma);
    const surfaceScalarField& lambdas = tlambdas();
    tmp<surfaceScalarField> tdeltaCoeffs = snGrad.deltaCoeffs(vf);
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();
    const surfaceScalarField& magSf = mesh.magSf();

    const scalarField& lambda = lambdas.primitiveField();
    const scalarField& dC = deltaCoeffs.primitiveField();
    const scalarField& mSf = magSf.primitiveField();
    const scalarField& g = gamma.primitiveField();

    // Patch values of the interpolated gamma times magSf (stock dotInterpolate
    // with geometricOneField, then gamma*magSf)
    PtrList<scalarField> pGammaMagSf(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        const fvPatchScalarField& pg = gamma.boundaryField()[patchi];
        if (pg.coupled())
        {
            pGammaMagSf.set
            (
                patchi,
                new scalarField
                (
                    lerp
                    (
                        pg.patchNeighbourField(),
                        pg.patchInternalField(),
                        lambdas.boundaryField()[patchi]
                    )*magSf.boundaryField()[patchi]
                )
            );
        }
        else
        {
            pGammaMagSf.set
            (
                patchi,
                new scalarField(pg*magSf.boundaryField()[patchi])
            );
        }
    }

    const dimensionSet gammaMagSfDims(gamma.dimensions()*magSf.dimensions());

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSfDims*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalarField& upper = fvm.upper();
    scalarField& diag = fvm.diag();

    // Patch coefficients (fvmLaplacianUncorrected)
    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
        const scalarField& pGamma = pGammaMagSf[patchi];
        const fvsPatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchi];

        if (pvf.coupled())
        {
            fvm.internalCoeffs()[patchi] =
                pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            fvm.boundaryCoeffs()[patchi] =
               -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
        }
        else
        {
            fvm.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();
            fvm.boundaryCoeffs()[patchi] = -pGamma*pvf.gradientBoundaryCoeffs();
        }
    }

    const bool corrected = snGrad.corrected();

    if (!corrected)
    {
        forAll(P, facei)
        {
            const scalar gface =
                lambda[facei]*(g[P[facei]] - g[N[facei]]) + g[N[facei]];
            const scalar gms = gface*mSf[facei];
            const scalar up = dC[facei]*gms;
            upper[facei] = up;
            diag[P[facei]] -= up;
            diag[N[facei]] -= up;
        }
    }
    else
    {
        // Gradient as correctedSnGrad::fullGradCorrection gets it, after the
        // matrix is constructed (as in the stock order)
        tmp<GradFieldType> tgrad =
            fv::gradScheme<Type>::New
            (
                mesh,
                mesh.gradScheme("grad(" + vf.name() + ')')
            )().grad(vf, "grad(" + vf.name() + ')');
        const GradFieldType& gradVf = tgrad();
        const Field<GradType>& G = gradVf.primitiveField();

        const surfaceVectorField& ncv = mesh.nonOrthCorrectionVectors();
        const vectorField& ncvi = ncv.primitiveField();

        const bool storeFlux = mesh.fluxRequired(vf.name());

        std::unique_ptr<SurfaceFieldType> ffcPtr;
        if (storeFlux)
        {
            ffcPtr.reset
            (
                new SurfaceFieldType
                (
                    IOobject
                    (
                        "faceFluxCorrection(" + vf.name() + ')',
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    mesh,
                    dimensioned<Type>
                    (
                        gammaMagSfDims*ncv.dimensions()*gradVf.dimensions(),
                        Zero
                    )
                )
            );
            ffcPtr->setOriented();
        }

        Field<Type> acc(vf.size(), Zero);

        forAll(P, facei)
        {
            const scalar gface =
                lambda[facei]*(g[P[facei]] - g[N[facei]]) + g[N[facei]];
            const scalar gms = gface*mSf[facei];
            const scalar up = dC[facei]*gms;
            upper[facei] = up;
            diag[P[facei]] -= up;
            diag[N[facei]] -= up;

            const Type corr =
                ncvi[facei]
              & (lambda[facei]*(G[P[facei]] - G[N[facei]]) + G[N[facei]]);
            const Type flux = gms*corr;

            if (storeFlux)
            {
                ffcPtr->primitiveFieldRef()[facei] = flux;
            }

            acc[P[facei]] += flux;
            acc[N[facei]] -= flux;
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatchField<GradType>& pG = gradVf.boundaryField()[patchi];
            const fvsPatchVectorField& pncv = ncv.boundaryField()[patchi];

            Field<Type> pcorr;
            if (pG.coupled())
            {
                pcorr =
                    pncv
                  & lerp
                    (
                        pG.patchNeighbourField(),
                        pG.patchInternalField(),
                        lambdas.boundaryField()[patchi]
                    );
            }
            else
            {
                pcorr = pncv & pG;
            }

            const Field<Type> pflux(pGammaMagSf[patchi]*pcorr);

            if (storeFlux)
            {
                ffcPtr->boundaryFieldRef()[patchi] = pflux;
            }

            const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
            forAll(faceCells, facei)
            {
                acc[faceCells[facei]] += pflux[facei];
            }
        }

        acc /= mesh.Vsc();
        fvm.source() -= mesh.V()*acc;

        if (storeFlux)
        {
            fvm.faceFluxCorrectionPtr() = std::move(ffcPtr);
        }
    }

    return tfvm;
}


template<class Type>
bool fusedLaplacianOk
(
    const fv::laplacianScheme<Type, scalar>& scheme,
    const surfaceInterpolationScheme<scalar>& gammaInterp,
    const fv::snGradScheme<Type>& snGrad
)
{
    return
        !stepReuse::off("lap")
     && isA<fv::gaussLaplacianScheme<Type, scalar>>(scheme)
     && gammaInterp.type() == "linear"
     && (snGrad.type() == "corrected" || snGrad.type() == "uncorrected");
}

} // End namespace stepReuseDetail


template<>
tmp<fvMatrix<scalar>> fv::laplacianScheme<scalar, scalar>::fvmLaplacian
(
    const volScalarField& gamma,
    const volScalarField& vf
)
{
    if
    (
        stepReuseDetail::fusedLaplacianOk
        (
            *this, tinterpGammaScheme_(), tsnGradScheme_()
        )
    )
    {
        stepReuse::report("lap", "laplacian: fused Gauss linear assembly");
        return stepReuseDetail::fusedGaussLaplacian
        (
            tinterpGammaScheme_(), tsnGradScheme_(), gamma, vf
        );
    }

    return fvmLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<>
tmp<fvMatrix<vector>> fv::laplacianScheme<vector, scalar>::fvmLaplacian
(
    const volScalarField& gamma,
    const volVectorField& vf
)
{
    if
    (
        stepReuseDetail::fusedLaplacianOk
        (
            *this, tinterpGammaScheme_(), tsnGradScheme_()
        )
    )
    {
        stepReuse::report("lap", "laplacian: fused Gauss linear assembly");
        return stepReuseDetail::fusedGaussLaplacian
        (
            tinterpGammaScheme_(), tsnGradScheme_(), gamma, vf
        );
    }

    return fvmLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}

} // End namespace Foam


// ************************************************************************* //
