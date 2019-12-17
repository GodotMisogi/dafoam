/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Child class for pimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef AdjointDerivativePimpleFoam_H
#define AdjointDerivativePimpleFoam_H

#include "AdjointDerivative.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class AdjointDerivativePimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class AdjointDerivativePimpleFoam
:
public AdjointDerivative
{

private:

protected:

    /// \name These are state variables, reference state variables, flow residuals, and partial derivatives
    //@{
    // volVectorField& U_;
    // volVectorField URef_;
    // volVectorField URes_;
    // volVectorField UResRef_;
    // volVectorField UResPartDeriv_;
    volVectorField& UMean_;
    volVectorField UMeanRef_;
    volVectorField UMeanRes_;
    volVectorField UMeanResRef_;
    volVectorField UMeanResPartDeriv_;

    // volScalarField& p_;
    // volScalarField pRef_;
    // volScalarField pRes_;
    // volScalarField pResRef_;
    // volScalarField pResPartDeriv_;
    volScalarField& pMean_;
    volScalarField pMeanRef_;
    volScalarField pMeanRes_;
    volScalarField pMeanResRef_;
    volScalarField pMeanResPartDeriv_;
    
    // surfaceScalarField& phi_;
    // surfaceScalarField phiRef_;  
    // surfaceScalarField phiRes_;
    // surfaceScalarField phiResRef_;
    // surfaceScalarField phiResPartDeriv_;
    surfaceScalarField& phiMean_; 
    surfaceScalarField phiMeanRef_;  
    surfaceScalarField phiMeanRes_;
    surfaceScalarField phiMeanResRef_;
    surfaceScalarField phiMeanResPartDeriv_;
    //@}

    pimpleControl pimple_;


public:

    TypeName ("pimpleDAFoam");
    // Constructors

    //- Construct from components
    AdjointDerivativePimpleFoam
    (
        fvMesh& mesh,
        const AdjointIO& adjIO,
        const AdjointSolverRegistry& adjReg,
        AdjointRASModel& adjRAS,
        AdjointIndexing& adjIdx,
        AdjointJacobianConnectivity& adjCon,
        AdjointObjectiveFunction& adjObj
    );

    //- Destructor
    virtual ~AdjointDerivativePimpleFoam()
    {}
    
    virtual void calcResiduals
    (
        const label isRef,
        const label isPC,
        const word fvMatrixName="None",
        const label updatePhi=0
    );
    
    virtual void updateIntermediateVariables();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
