/*---------------------------------------------------------------------------*\
    Replacement for fv::gaussGrad<Type>::calcGrad (Type scalar and vector):
    the returning version (grad(U) through cellLimited, first evaluation of a
    cached gradient) and the in-place version (a cached gradient updated:
    grad(p), grad(k), grad(omega)).

    Stock: interpolate the field to the faces (a face field), then a face
    loop adding Sf*value to the owner and subtracting it from the
    neighbour, then the patches, then dividing by the cell volumes. Here the
    interpolation is done inside that face loop: the face value is the same
    expression as in the stock linear interpolation (lambda*(a - b) + b, left
    to the same contraction), the accumulation order is unchanged, and the
    patch values use the stock functions. Bit-identical.

    Only with linear interpolation; otherwise, or with STEP_REUSE_OFF=grad,
    the stock function runs (through the next definition of the symbol).
\*---------------------------------------------------------------------------*/

#include "gaussGrad.H"
#include "stepReuse.H"

namespace Foam
{
namespace stepReuseDetail
{

// igGrad += Gauss sum over the faces of Sf*linear-interpolate(vsf), then
// divided by the cell volumes (igGrad zero on entry)
template<class Type>
void fusedGaussGradSum
(
    Field<typename outerProduct<vector, Type>::type>& igGrad,
    const surfaceScalarField& lambdas,
    const GeometricField<Type, fvPatchField, volMesh>& vsf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();
    const scalarField& lambda = lambdas.primitiveField();
    const Field<Type>& v = vsf.primitiveField();

    forAll(owner, facei)
    {
        const Type ssf =
            lambda[facei]*(v[owner[facei]] - v[neighbour[facei]])
          + v[neighbour[facei]];

        const GradType Sfssf = Sf[facei]*ssf;

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const fvPatchField<Type>& pvf = vsf.boundaryField()[patchi];

        if (pvf.coupled())
        {
            const Field<Type> pssf
            (
                lerp
                (
                    pvf.patchNeighbourField(),
                    pvf.patchInternalField(),
                    lambdas.boundaryField()[patchi]
                )
            );

            forAll(mesh.boundary()[patchi], facei)
            {
                igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
            }
        }
        else
        {
            forAll(mesh.boundary()[patchi], facei)
            {
                igGrad[pFaceCells[facei]] += pSf[facei]*pvf[facei];
            }
        }
    }

    igGrad /= mesh.V();
}

} // End namespace stepReuseDetail


// Returning version (stock: gradf(interpolate(vsf), name), then
// correctBoundaryConditions(vsf, gGrad))
#define STEP_REUSE_GAUSS_GRAD(Type, mangledReturn, mangledInPlace)            \
                                                                               \
template<>                                                                     \
tmp<GeometricField<outerProduct<vector, Type>::type, fvPatchField, volMesh>>   \
fv::gaussGrad<Type>::calcGrad                                                  \
(                                                                              \
    const GeometricField<Type, fvPatchField, volMesh>& vsf,                    \
    const word& name                                                           \
) const                                                                        \
{                                                                              \
    typedef outerProduct<vector, Type>::type GradType;                         \
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;     \
                                                                               \
    if (stepReuse::off("grad") || tinterpScheme_().type() != "linear")         \
    {                                                                          \
        typedef tmp<GradFieldType> (*fn)                                       \
        (                                                                      \
            const fv::gaussGrad<Type>*,                                        \
            const GeometricField<Type, fvPatchField, volMesh>&,                \
            const word&                                                        \
        );                                                                     \
        static const fn stock =                                                \
            reinterpret_cast<fn>(stepReuse::next(mangledReturn));              \
        return stock(this, vsf, name);                                         \
    }                                                                          \
                                                                               \
    stepReuse::report("grad", "gaussGrad: interpolation fused");               \
                                                                               \
    const fvMesh& mesh = vsf.mesh();                                           \
                                                                               \
    tmp<GradFieldType> tgGrad                                                  \
    (                                                                          \
        new GradFieldType                                                      \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                name,                                                          \
                vsf.instance(),                                                \
                mesh,                                                          \
                IOobject::NO_READ,                                             \
                IOobject::NO_WRITE                                             \
            ),                                                                 \
            mesh,                                                              \
            dimensioned<GradType>(vsf.dimensions()/dimLength, Zero),           \
            fvPatchFieldBase::extrapolatedCalculatedType()                     \
        )                                                                      \
    );                                                                         \
    GradFieldType& gGrad = tgGrad.ref();                                       \
                                                                               \
    stepReuseDetail::fusedGaussGradSum                                         \
    (                                                                          \
        gGrad.primitiveFieldRef(),                                             \
        tinterpScheme_().weights(vsf)(),                                       \
        vsf                                                                    \
    );                                                                         \
                                                                               \
    gGrad.correctBoundaryConditions();                                         \
                                                                               \
    correctBoundaryConditions(vsf, gGrad);                                     \
                                                                               \
    return tgGrad;                                                             \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
void fv::gaussGrad<Type>::calcGrad                                             \
(                                                                              \
    GeometricField                                                             \
    <                                                                          \
        outerProduct<vector, Type>::type,                                      \
        fvPatchField,                                                          \
        volMesh                                                                \
    >& gGrad,                                                                  \
    const GeometricField<Type, fvPatchField, volMesh>& vsf                     \
) const                                                                        \
{                                                                              \
    typedef outerProduct<vector, Type>::type GradType;                         \
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;     \
                                                                               \
    if (stepReuse::off("grad") || tinterpScheme_().type() != "linear")         \
    {                                                                          \
        typedef void (*fn)                                                     \
        (                                                                      \
            const fv::gaussGrad<Type>*,                                        \
            GradFieldType&,                                                    \
            const GeometricField<Type, fvPatchField, volMesh>&                 \
        );                                                                     \
        static const fn stock =                                                \
            reinterpret_cast<fn>(stepReuse::next(mangledInPlace));             \
        stock(this, gGrad, vsf);                                               \
        return;                                                                \
    }                                                                          \
                                                                               \
    stepReuse::report("grad", "gaussGrad: interpolation fused");               \
                                                                               \
    gGrad = dimensioned<GradType>(vsf.dimensions()/dimLength, Zero);           \
                                                                               \
    stepReuseDetail::fusedGaussGradSum                                         \
    (                                                                          \
        gGrad.primitiveFieldRef(),                                             \
        tinterpScheme_().weights(vsf)(),                                       \
        vsf                                                                    \
    );                                                                         \
                                                                               \
    gGrad.correctBoundaryConditions();                                         \
                                                                               \
    correctBoundaryConditions(vsf, gGrad);                                     \
}


STEP_REUSE_GAUSS_GRAD
(
    scalar,
    "_ZNK4Foam2fv9gaussGradIdE8calcGradERKNS_14GeometricFieldIdNS_"
    "12fvPatchFieldENS_7volMeshEEERKNS_4wordE",
    "_ZNK4Foam2fv9gaussGradIdE8calcGradERNS_14GeometricFieldINS_6VectorIdEENS_"
    "12fvPatchFieldENS_7volMeshEEERKNS3_IdS6_S7_EE"
)

STEP_REUSE_GAUSS_GRAD
(
    vector,
    "_ZNK4Foam2fv9gaussGradINS_6VectorIdEEE8calcGradERKNS_14GeometricFieldIS3_"
    "NS_12fvPatchFieldENS_7volMeshEEERKNS_4wordE",
    "_ZNK4Foam2fv9gaussGradINS_6VectorIdEEE8calcGradERNS_14GeometricFieldINS_"
    "6TensorIdEENS_12fvPatchFieldENS_7volMeshEEERKNS5_IS3_S8_S9_EE"
)

#undef STEP_REUSE_GAUSS_GRAD

} // End namespace Foam


// ************************************************************************* //
