/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

\*---------------------------------------------------------------------------*/

#include "AdjointDerivativePimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(AdjointDerivativePimpleFoam, 0);
addToRunTimeSelectionTable(AdjointDerivative, AdjointDerivativePimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AdjointDerivativePimpleFoam::AdjointDerivativePimpleFoam
(
    fvMesh& mesh,
    const AdjointIO& adjIO,
    const AdjointSolverRegistry& adjReg,
    AdjointRASModel& adjRAS,
    AdjointIndexing& adjIdx,
    AdjointJacobianConnectivity& adjCon,
    AdjointObjectiveFunction& adjObj
)
    :
    AdjointDerivative(mesh,adjIO,adjReg,adjRAS,adjIdx,adjCon,adjObj),
    // initialize and register state variables and their residuals, we use macros defined in macroFunctions.H
    // setResidualClassMemberVector(U,dimensionSet(0,1,-2,0,0,0,0)),
    // setResidualClassMemberScalar(p,dimensionSet(0,0,-1,0,0,0,0)),
    // setResidualClassMemberPhi(phi),
    setResidualClassMemberVector(UMean,dimensionSet(0,1,-2,0,0,0,0)),
    setResidualClassMemberScalar(pMean,dimensionSet(0,0,-1,0,0,0,0)),
    setResidualClassMemberPhi(phiMean),
    // create PimpleControl
    pimple_(mesh) 
    
{
    this->copyStates("Var2Ref"); // copy states to statesRef
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void AdjointDerivativePimpleFoam::calcResiduals
(
    const label isRef,
    const label isPC,
    const word fvMatrixName,
    const label updatePhi
)
{
    // We dont support MRF and fvOptions so all the related lines are commented 
    // out for now
    
    word divUScheme="div(phi,U)";
    if(isPC) divUScheme="div(pc)";

    // ******** U Residuals **********
    // copied and modified from UEqn.H

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(UMean_) + fvm::div(phiMean_, UMean_)
      + this->MRF_.DDt(UMean_)
      + adjRAS_.divDevReff(UMean_)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    //UEqn.relax();

    // set fvMatrix for fast PC construction in NK solver
    setFvMatrix("UMean",UEqn);

    if (!updatePhi)
    {
        if(isRef) UMeanResRef_  = (UEqn&UMean_) + fvc::grad(pMean_);
        else UMeanRes_  = (UEqn&UMean_) + fvc::grad(pMean_);
        normalizeResiduals(UMeanRes);
        scaleResiduals(UMeanRes);
    }

    // ******** p Residuals **********
    // copied and modified from pEqn.H
    // NOTE manually set pRefCell and pRefValue
    label pRefCell=0;
    scalar pRefValue=0.0;
    
    // Note: relax UEqn after the URes is calculated
    UEqn.relax();
    
    volScalarField rAU(1.0/UEqn.A());
    //volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), UMean_, pMean_));
    //***************** NOTE *******************
    // we should not use the constrainHbyA function above since it
    // will degrade the accuracy of shape derivatives. Basically, we should
    // not constrain any variable because it will create discontinuity
    volVectorField HbyA("HbyA", UMean_);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    this->MRF_.makeRelative(phiHbyA);

    if (pMean_.needReference())
    {
        fvc::makeRelative(phiHbyA, UMean_);
        adjustPhi(phiHbyA, UMean_, pMean_);
        fvc::makeAbsolute(phiHbyA, UMean_);
    }

    tmp<volScalarField> rAtU(rAU);

    if (pimple_.consistent())
    {
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
        phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(pMean_)*mesh_.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(pMean_);
    }

    if (pimple_.nCorrPISO() <= 1)
    {
        tUEqn.clear();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pMean_, UMean_, phiHbyA, rAtU(), this->MRF_);
    
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), pMean_) == fvc::div(phiHbyA)
    );
    pEqn.setReference(pRefCell, pRefValue);

    // set fvMatrix for fast PC construction in NK solver
    setFvMatrix("pMean",pEqn);

    if (!updatePhi)
    {
        if(isRef) pMeanResRef_  = pEqn&pMean_;
        else pMeanRes_  = pEqn&pMean_;
        normalizeResiduals(pMeanRes);
        scaleResiduals(pMeanRes);
    }

    if(updatePhi) phiMean_=phiHbyA - pEqn.flux();

    // ******** phi Residuals **********
    // copied and modified from pEqn.H
    if(isRef) phiMeanResRef_ = phiHbyA - pEqn.flux() - phiMean_;
    else phiMeanRes_ = phiHbyA - pEqn.flux() - phiMean_;
    // need to normalize phiRes
    normalizePhiResiduals(phiMeanRes);
    scalePhiResiduals(phiMeanRes);   
    
    return;

}

void AdjointDerivativePimpleFoam::updateIntermediateVariables()
{
    // update velocity boundary based on MRF
    this->MRF_.correctBoundaryVelocity(UMean_);
}

} // End namespace Foam

// ************************************************************************* //
