/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

\*---------------------------------------------------------------------------*/

#include "AdjointSolverRegistryInterFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AdjointSolverRegistryInterFoam, 0);
    addToRunTimeSelectionTable(AdjointSolverRegistry, AdjointSolverRegistryInterFoam, dictionary);

    AdjointSolverRegistryInterFoam :: AdjointSolverRegistryInterFoam(const fvMesh& mesh):
        AdjointSolverRegistry(mesh)
        {
        // Register state variables
        // NOTE: do not include any turbulence state variables since they will be added 
        // in the AdjointRASModel class!
        volScalarStates.append("p");
        volScalarStates.append("alpha.water");
        volVectorStates.append("U");
        surfaceScalarStates.append("phi");

        // for debugging
        //Info<<"volScalarStates: "<<volScalarStates<<endl;
        //Info<<"volVectorStates: "<<volVectorStates<<endl;
        //Info<<"surfaceScalarStates: "<<surfaceScalarStates<<endl;
        

        // Need to call this to check if the registered states are actually created in db
        //this->validate();
        
        // setup the derived state info
        this->setDerivedInfo();
        }
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //