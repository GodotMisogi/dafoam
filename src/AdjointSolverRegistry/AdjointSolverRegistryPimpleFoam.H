/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Child class for PimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef AdjointSolverRegistryPimpleFoam_H
#define AdjointSolverRegistryPimpleFoam_H

#include "AdjointSolverRegistry.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class AdjointSolverRegistryPimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class AdjointSolverRegistryPimpleFoam
:
public AdjointSolverRegistry
{

private:


public:

    TypeName ("pimpleDAFoam");
    // Constructors

    //- Construct from components
    AdjointSolverRegistryPimpleFoam
    (
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~AdjointSolverRegistryPimpleFoam()
    {}
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
