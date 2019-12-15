/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Child class for InterFoam

\*---------------------------------------------------------------------------*/

#ifndef AdjointJacobianConnectivityInterFoam_H
#define AdjointJacobianConnectivityInterFoam_H

#include "AdjointJacobianConnectivity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class AdjointJacobianConnectivityInterFoam Declaration
\*---------------------------------------------------------------------------*/

class AdjointJacobianConnectivityInterFoam
:
public AdjointJacobianConnectivity
{

private:


public:

    TypeName ("interDAFoam");
    // Constructors

    //- Construct from components
    AdjointJacobianConnectivityInterFoam
    (
        const fvMesh& mesh,
        const AdjointIO& adjIO,
        const AdjointSolverRegistry& adjReg,
        AdjointRASModel& adjRAS,
        AdjointIndexing& adjIdx
    );

    //- Destructor
    virtual ~AdjointJacobianConnectivityInterFoam()
    {}
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
