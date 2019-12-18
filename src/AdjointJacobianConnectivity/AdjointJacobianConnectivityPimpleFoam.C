/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

\*---------------------------------------------------------------------------*/

#include "AdjointJacobianConnectivityPimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(AdjointJacobianConnectivityPimpleFoam, 0);
addToRunTimeSelectionTable(AdjointJacobianConnectivity, AdjointJacobianConnectivityPimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AdjointJacobianConnectivityPimpleFoam::AdjointJacobianConnectivityPimpleFoam
(
    const fvMesh& mesh,
    const AdjointIO& adjIO,
    const AdjointSolverRegistry& adjReg,
    AdjointRASModel& adjRAS,
    AdjointIndexing& adjIdx
)
    :
    AdjointJacobianConnectivity(mesh,adjIO,adjReg,adjRAS,adjIdx)
{

    // Adj state connectivity levels for pimpleFoam, numbers denote the level of connectivity
    // N/A means this state does not connect to the corrsponding residual 
    // ********************************NOTE**********************************
    // One does not need to specify connectivity for each turbulence model, set the connectivity
    // for nut only. How is nut connected to the other turbStates will be specified in the turbulence class
    // This is done by calling correctAdjStateResidualTurbCon. For example, for SA model we just replace
    // nut with nuTilda, for SST model, we need to add extract connectivity since nut depends on grad(U), k, and omega
    // **********************************************************************
    //              U      p     nut    phi
    // URes         2      1      1      0
    // pRes         3      2      2      1
    // phiRes       2      1      1      0
    
    // adjStateResidualConInfo_.set
    // (
    //     "URes",
    //     {
    //         {"U","p","nut","phi"}, // lv0
    //         {"U","p","nut"},       // lv1
    //         {"U"}                  // lv2
    //     }
    // );
    
    // adjStateResidualConInfo_.set
    // (
    //     "pRes",
    //     {
    //         {"U","p","nut","phi"}, // lv0
    //         {"U","p","nut","phi"}, // lv1
    //         {"U","p","nut"},       // lv2
    //         {"U"}                  // lv3
    //     }
    // );
    
    // adjStateResidualConInfo_.set
    // (
    //     "phiRes",
    //     {
    //         {"U","p","nut","phi"}, // lv0
    //         {"U","p","nut"},       // lv1
    //         {"U"},                 // lv2
    //     }
    // );

    adjStateResidualConInfo_.set
    (
        "UMeanRes",
        {
            {"UMean","pMean","nut","phi"}, // lv0
            {"UMean","pMean","nut"},       // lv1
            {"UMean"}                  // lv2
        }
    );
    
    adjStateResidualConInfo_.set
    (
        "pMeanRes",
        {
            {"UMean","pMean","nut","phi"}, // lv0
            {"UMean","pMean","nut","phi"}, // lv1
            {"UMean","pMean","nut"},       // lv2
            {"UMean"}                  // lv3
        }
    );
    
    adjStateResidualConInfo_.set
    (
        "phiRes",
        {
            {"UMean","pMean","nut","phi"}, // lv0
            {"UMean","pMean","nut"},       // lv1
            {"UMean"},                 // lv2
        }
    );
    
    // need to correct turb con for each residual
    adjRAS.correctAdjStateResidualTurbCon(adjStateResidualConInfo_["UMeanRes"]);
    adjRAS.correctAdjStateResidualTurbCon(adjStateResidualConInfo_["pMeanRes"]);
    adjRAS.correctAdjStateResidualTurbCon(adjStateResidualConInfo_["phiRes"]);
    
    // add turbulence residual connectivity
    adjRAS.setAdjStateResidualTurbCon(adjStateResidualConInfo_);
    
    //Info<<adjStateResidualConInfo_<<endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
