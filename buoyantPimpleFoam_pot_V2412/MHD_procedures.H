// *** low Re_m MHD model is solved outside the PIMPLE loop
// becuase dynamics values of U are needed in Lorentz force term in UEqn.H

surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();

// *** Poisson equation for electric potential
fvScalarMatrix PotEEqn
(
fvm::laplacian(PotE) == fvc::div(psiub)
);

PotEEqn.setReference(0, 0);

// *** Solving Poisson equation
PotEEqn.solve();

// *** Computation of current density at cell faces
surfaceScalarField jn = -(fvc::snGrad(PotE) * mesh.magSf()) + psiub;

// *** Current density at face center
surfaceVectorField jnv = jn * mesh.Cf();

// *** Interpolation of current density at cell center
volVectorField jfinal = fvc::surfaceIntegrate(jnv) - (fvc::surfaceIntegrate(jn) * mesh.C());

// *** Update current density distribution and boundary condition
jfinal.correctBoundaryConditions();

// *** Lorentz force computation
JxB = sigma* (jfinal ^ B0);
//
