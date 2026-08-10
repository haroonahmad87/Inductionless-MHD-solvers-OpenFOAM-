/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantBoussinesqPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids,
    with optional mesh motion and mesh topology changes.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

Note 
   This code does DNS of buoyant incompressible flows with Boussinesq approximation: modified by Haroon Ahmad

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"     // Commented for DNS
#include "radiationModel.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient DNS solver for buoyant, turbulent flow"
        " of incompressible fluids, with optional mesh"
        " motion and mesh topology changes.\n"
        "Uses the Boussinesq approximation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    //turbulence->validate();     // Commented for DNS

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // *** Lorentz force term initialization 
    volVectorField lorentz = sigma * (-fvc::grad(PotE) ^ B0) + sigma * ((U ^ B0) ^ B0);

    while (runTime.run())  // *** time-loop
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop  // *** outer PIMPLE loop (SIMPLE-like)
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

//                if (mesh.changing())
//                {
//                    MRF.update();
//
//                    if (correctPhi)
//                    {
//                        // Calculate absolute flux
//                        // from the mapped surface velocity
//                        phi = mesh.Sf() & Uf();
//
//                        #include "correctPhi.H"
//
//                        // Make the flux relative to the mesh motion
//                        fvc::makeRelative(phi, U);
//                    }
//
//                    if (checkMeshCourantNo)
//                    {
//                        #include "meshCourantNo.H"
//                    }
//               }
            }

            #include "UEqn.H"   // *** solves the momentum equation containing Lorentz force term
            #include "TEqn.H"   // *** solves the energy equation 

            // --- Pressure corrector loop  // *** Inner PISO loop
	    while (pimple.correct())
	    {
                #include "pEqn.H"   // *** solves the pressure-Poisson eqn. and corrects the velocity & fluxes
            }

            if (pimple.turbCorr())  // *** correction of flow properties in turbulence model  
            {
                laminarTransport.correct();
                turbulence->correct();
            }
         
	// *** low Re_m MHD model is solved inside the PIMPLE loop
	// becuase dynamics values of U are needed in Lorentz force term in UEqn.H 

	// *** Interpolating cross product u x B over mesh faces
//        surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();

        // *** Poisson equation for electric potential
//        fvScalarMatrix PotEEqn
//        (
//        fvm::laplacian(PotE) == fvc::div(psiub)
//        );

        // *** Reference potential
        //PotEEqn.setReference(PotERefCell, PotERefValue);
//        PotEEqn.setReference(pRefCell, pRefValue);

        // *** Solving Poisson equation
//        PotEEqn.solve();

        // *** Computation of current density at cell faces
//        surfaceScalarField jn = -(fvc::snGrad(PotE) * mesh.magSf()) + psiub;

        // *** Current density at face center
//        surfaceVectorField jnv = jn * mesh.Cf();

        // *** Interpolation of current density at cell center
//        volVectorField jfinal = fvc::surfaceIntegrate(jnv) - (fvc::surfaceIntegrate(jn) * mesh.C());

        // *** Update current density distribution and boundary condition
//        jfinal.correctBoundaryConditions();

        // *** Lorentz force computation
//        lorentz = sigma* (jfinal ^ B0);

        }

        // *** low Re_m MHD model is solved outside the PIMPLE loop
        // becuase dynamics values of U are needed in Lorentz force term in UEqn.H

        // *** Interpolating cross product u x B over mesh faces
        surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();

        // *** Poisson equation for electric potential
        fvScalarMatrix PotEEqn
        (
        fvm::laplacian(PotE) == fvc::div(psiub)
        );

        // *** Reference potential
        //PotEEqn.setReference(PotERefCell, PotERefValue);
        PotEEqn.setReference(pRefCell, pRefValue);

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
        lorentz = sigma* (jfinal ^ B0);	

        runTime.write();

        runTime.printExecutionTime(Info);
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
