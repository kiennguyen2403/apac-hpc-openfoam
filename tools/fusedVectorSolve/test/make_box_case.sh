#!/bin/bash
# 3D bit-identity test case for fusedVectorSolve: 30x12x12 box with the boundary types of
# occDrivAerStaticMesh that matter for the segregated vector solve (fixedValue inlet,
# inletOutlet outlet, noSlip wall with wall functions, slip top, symmetry sides), kOmegaSST,
# fixed-iteration solver controls. usage: make_box_case.sh <case dir>
set -euo pipefail
C=$1
mkdir -p "$C"/{0,constant,system}
cd "$C"
H='FoamFile { version 2.0; format ascii; class'
cat > system/blockMeshDict <<EOF
$H dictionary; object blockMeshDict; }
scale 1;
vertices ( (0 0 0) (3 0 0) (3 1 0) (0 1 0) (0 0 1) (3 0 1) (3 1 1) (0 1 1) );
blocks ( hex (0 1 2 3 4 5 6 7) (30 12 12) simpleGrading (1 1 4) );
boundary
(
    inlet  { type patch; faces ((0 4 7 3)); }
    outlet { type patch; faces ((1 2 6 5)); }
    ground { type wall;  faces ((0 3 2 1)); }
    top    { type patch; faces ((4 5 6 7)); }
    sides  { type symmetry; faces ((0 1 5 4) (3 7 6 2)); }
);
EOF
fld() { # name class dims internal inlet outlet ground top
cat > 0/$1 <<EOF
$H $2; object $1; }
dimensions $3;
internalField uniform $4;
boundaryField
{
    inlet  { $5 }
    outlet { $6 }
    ground { $7 }
    top    { $8 }
    sides  { type symmetry; }
}
EOF
}
fld U volVectorField '[0 1 -1 0 0 0 0]' '(10 0 0)' 'type fixedValue; value uniform (10 0.5 0.2);' 'type inletOutlet; inletValue uniform (0 0 0); value uniform (10 0 0);' 'type noSlip;' 'type slip;'
fld p volScalarField '[0 2 -2 0 0 0 0]' 0 'type zeroGradient;' 'type fixedValue; value uniform 0;' 'type zeroGradient;' 'type slip;'
fld k volScalarField '[0 2 -2 0 0 0 0]' 0.24 'type fixedValue; value uniform 0.24;' 'type inletOutlet; inletValue uniform 0.24; value uniform 0.24;' 'type kqRWallFunction; value uniform 0.24;' 'type slip;'
fld omega volScalarField '[0 0 -1 0 0 0 0]' 1.78 'type fixedValue; value uniform 1.78;' 'type inletOutlet; inletValue uniform 1.78; value uniform 1.78;' 'type omegaWallFunction; value uniform 1.78;' 'type slip;'
fld nut volScalarField '[0 2 -1 0 0 0 0]' 0 'type calculated; value uniform 0;' 'type calculated; value uniform 0;' 'type nutkWallFunction; value uniform 0;' 'type slip;'
cat > constant/transportProperties <<EOF
$H dictionary; object transportProperties; }
transportModel Newtonian; nu 1.5e-05;
EOF
cat > constant/turbulenceProperties <<EOF
$H dictionary; object turbulenceProperties; }
simulationType RAS; RAS { RASModel kOmegaSST; turbulence on; printCoeffs off; }
EOF
cat > system/controlDict <<EOF
$H dictionary; object controlDict; }
application simpleFoam; startFrom startTime; startTime 0; stopAt endTime; endTime 20; deltaT 1;
writeControl timeStep; writeInterval 20; writeFormat ascii; writePrecision 17; runTimeModifiable false;
EOF
cat > system/fvSchemes <<EOF
$H dictionary; object fvSchemes; }
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwindV grad(U); div(phi,k) bounded Gauss upwind; div(phi,omega) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
EOF
cat > system/fvSolution <<EOF
$H dictionary; object fvSolution; }
solvers
{
    p { solver GAMG; smoother GaussSeidel; tolerance 0; relTol 0; maxIter 12;
        coarsestLevelCorr { solver PCG; preconditioner DIC; minIter 35; maxIter 35; } }
    "(U|k|omega)" { solver smoothSolver; smoother GaussSeidel; tolerance 0; relTol 0; maxIter 15; nSweeps 1; }
}
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; }
relaxationFactors { equations { U 0.9; ".*" 0.9; } }
EOF
printf "$H dictionary; object decomposeParDict; }\nnumberOfSubdomains 4;\nmethod scotch;\n" > system/decomposeParDict
