/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Child class for interFoam

\*---------------------------------------------------------------------------*/

#ifndef AdjointSolverRegistryInterFoam_H
#define AdjointSolverRegistryInterFoam_H

#include "AdjointSolverRegistry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class AdjointSolverRegistryInterFoam Declaration
\*---------------------------------------------------------------------------*/

class AdjointSolverRegistryInterFoam:
    public AdjointSolverRegistry
    {
        private:

        public:

        TypeName ("interDAFoam");
        // Constructors

        //- Construct from components
        AdjointSolverRegistryInterFoam
        (
            const fvMesh& mesh
        );

        //- Destructor
        virtual ~AdjointSolverRegistryInterFoam()
        {}
        
    };



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //