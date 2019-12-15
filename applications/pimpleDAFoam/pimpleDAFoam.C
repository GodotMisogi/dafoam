/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Adjoint solver for pimpleFoam

\*---------------------------------------------------------------------------*/

static char help[] = "Solves a linear system in parallel with KSP in OpenFOAM.\n\n";

#include <petscksp.h>
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "AdjointIO.H"
#include "AdjointSolverRegistry.H"
#include "AdjointRASModel.H"
#include "AdjointIndexing.H"
#include "AdjointJacobianConnectivity.H"
#include "AdjointObjectiveFunction.H"
#include "AdjointDerivative.H"
#include "AdjointNewtonKrylov.H"
#include "nearWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    // #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    // Initialize the petsc solver. This needs to be called after the case
    // setup so that petsc uses the OpenFOAM MPI_COMM
    PetscInitialize(&argc,&argv,(char*)0,help);

    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "createFvOptions.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        
        adjObj.printObjFuncValues();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    adjObj.writeObjFuncValues();

    Info<< "End\n" << endl;

    // solve adjoint
    autoPtr<AdjointDerivative> adjDev
    (
        AdjointDerivative::New(mesh,adjIO,adjReg(),adjRAS(),adjIdx,adjCon(),adjObj)
    );

    adjDev->calcFlowResidualStatistics("print");

    if(adjIO.useNKSolver and !adjIO.solveAdjoint)
    {

        AdjointNewtonKrylov adjNK(mesh,adjIO,adjReg,adjRAS(),adjObj,adjDev());
        adjNK.solve();

        adjDev->calcFlowResidualStatistics("print");

        adjDev->writeStates();
        adjRAS->writeTurbStates();

        Info<< "NK Finished!\n" << endl;
        
        return 0;
    }

    if (adjIO.solveAdjoint)
    {
        adjDev->solve();
    }

    PetscFinalize();

    return 0;
}


// ************************************************************************* //
