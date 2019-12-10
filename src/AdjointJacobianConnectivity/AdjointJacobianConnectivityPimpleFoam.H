/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Child class for pimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef AdjointJacobianConnectivityPimpleFoam_H
#define AdjointJacobianConnectivityPimpleFoam_H

#include "AdjointJacobianConnectivity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class AdjointJacobianConnectivityPimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class AdjointJacobianConnectivityPimpleFoam
:
public AdjointJacobianConnectivity
{

private:


public:

    TypeName ("PimpleDAFoam");
    // Constructors

    //- Construct from components
    AdjointJacobianConnectivityPimpleFoam
    (
        const fvMesh& mesh,
        const AdjointIO& adjIO,
        const AdjointSolverRegistry& adjReg,
        AdjointRASModel& adjRAS,
        AdjointIndexing& adjIdx
    );

    //- Destructor
    virtual ~AdjointJacobianConnectivityPimpleFoam()
    {}
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
