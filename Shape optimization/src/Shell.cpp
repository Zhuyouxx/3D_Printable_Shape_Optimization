#include <Eigen/Core>
#include <iostream>
#include <map>
#include <cmath>
#include <GeometryDerivatives.h>
#include <libshell/MeshConnectivity.h>
#include <libshell/ElasticShell.h>
#include <libshell/MidedgeAngleTanFormulation.h>
#include <libshell/MidedgeAngleSinFormulation.h>
#include <libshell/MidedgeAverageFormulation.h>
#include <libshell/StVKMaterial.h>
#include <libshell/BilayerStVKMaterial.h>
#include <libshell/TensionFieldStVKMaterial.h>
#include <libshell/NeoHookeanMaterial.h>
#include <libshell/RestState.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>

//#include "../include/MeshConnectivity.h"
//#include "../include/ElasticShell.h"
//#include "../include/MidedgeAngleTanFormulation.h"
//#include "../include/MidedgeAngleSinFormulation.h"
//#include "../include/MidedgeAverageFormulation.h"
//#include "../include/StVKMaterial.h"
//#include "../include/TensionFieldStVKMaterial.h"
//#include "../include/NeoHookeanMaterial.h"
//#include "../include/RestState.h"
//#include <igl/readOBJ.h>
//#include <igl/writeOBJ.h>
//#include <igl/edges.h>



//double stretchingEnergy_t(
//    const LibShell::MeshConnectivity& mesh,
//    const Eigen::MatrixXd& curPos,
//    const LibShell::MonolayerRestState & restState,
//    int face,
//    Eigen::Matrix<double, 1, 9>* derivative, // F(face, i)
//    Eigen::Matrix<double, 9, 9>* hessian)
//{
//    using namespace Eigen;
//
//    assert(restState.type() == RestStateType::RST_MONOLAYER);
//    const LibShell::MonolayerRestState& rs = (const LibShell::MonolayerRestState&)restState;
//
//    double coeff = rs.thicknesses[face] / 4.0;
//    Matrix2d abarinv = rs.abars[face].inverse();
//    Matrix<double, 4, 9> aderiv;
//    std::vector<Matrix<double, 9, 9> > ahess;
//    //Matrix2d a = LibShell::firstFundamentalForm(mesh, curPos, face, &aderiv ,  NULL);
//
//    Eigen::Vector3d q0 = curPos.row(mesh.faceVertex(face, 0));
//    Eigen::Vector3d q1 = curPos.row(mesh.faceVertex(face, 1));
//    Eigen::Vector3d q2 = curPos.row(mesh.faceVertex(face, 2));
//    Eigen::Matrix2d a;
//    a << (q1 - q0).dot(q1 - q0), (q1 - q0).dot(q2 - q0),
//        (q2 - q0).dot(q1 - q0), (q2 - q0).dot(q2 - q0);
//    aderiv.setZero();
//    aderiv.block<1, 3>(0, 3) += 2.0 * (q1 - q0).transpose();
//    aderiv.block<1, 3>(0, 0) -= 2.0 * (q1 - q0).transpose();
//    aderiv.block<1, 3>(1, 6) += (q1 - q0).transpose();
//    aderiv.block<1, 3>(1, 3) += (q2 - q0).transpose();
//    aderiv.block<1, 3>(1, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
//    aderiv.block<1, 3>(2, 6) += (q1 - q0).transpose();
//    aderiv.block<1, 3>(2, 3) += (q2 - q0).transpose();
//    aderiv.block<1, 3>(2, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
//    aderiv.block<1, 3>(3, 6) += 2.0 * (q2 - q0).transpose();
//    aderiv.block<1, 3>(3, 0) -= 2.0 * (q2 - q0).transpose();
//
//    Matrix2d M = abarinv * (a - rs.abars[face]);
//    double dA = 0.5 * sqrt(rs.abars[face].determinant());
//    double lameAlpha = rs.lameAlpha[face];
//    double lameBeta = rs.lameBeta[face];
//
//    double StVK = 0.5 * lameAlpha * pow(M.trace(), 2) + lameBeta * (M * M).trace();
//    double result = coeff * dA * StVK;
//
//    if (derivative)
//    {
//        Matrix2d temp = lameAlpha * M.trace() * abarinv + 2 * lameBeta * M * abarinv;
//        *derivative = coeff * dA * aderiv.transpose() * Map<Vector4d>(temp.data());
//    }
//
//    if (hessian)
//    {
//        Matrix<double, 1, 9> inner = aderiv.transpose() * Map<Vector4d>(abarinv.data());
//        *hessian = lameAlpha * inner.transpose() * inner;
//
//        Matrix2d Mainv = M * abarinv;
//        for (int i = 0; i < 4; ++i) // iterate over Mainv and abarinv as if they were vectors
//            *hessian += (lameAlpha * M.trace() * abarinv(i) + 2 * lameBeta * Mainv(i)) * ahess[i];
//
//        Matrix<double, 1, 9> inner00 = abarinv(0, 0) * aderiv.row(0) + abarinv(0, 1) * aderiv.row(2);
//        Matrix<double, 1, 9> inner01 = abarinv(0, 0) * aderiv.row(1) + abarinv(0, 1) * aderiv.row(3);
//        Matrix<double, 1, 9> inner10 = abarinv(1, 0) * aderiv.row(0) + abarinv(1, 1) * aderiv.row(2);
//        Matrix<double, 1, 9> inner11 = abarinv(1, 0) * aderiv.row(1) + abarinv(1, 1) * aderiv.row(3);
//        *hessian += 2 * lameBeta * inner00.transpose() * inner00;
//        *hessian += 2 * lameBeta * (inner01.transpose() * inner10 + inner10.transpose() * inner01);
//        *hessian += 2 * lameBeta * inner11.transpose() * inner11;
//
//        *hessian *= coeff * dA;
//    }
//
//    return result;
//}


double Energy_shell(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXd& grad) {

    //Eigen::MatrixXd V;
    //Eigen::MatrixXi F;
    //igl::readOBJ("D://IDE//VisualStudio_Project//DR_manufacturing//TestOn3D//Test3D//src//bunny.obj", V, F);
    //std::cout << V.size() << std::endl;
    //igl::writeOBJ("D://IDE//VisualStudio_Project//DR_manufacturing//TestOn3D//Test3D//results//bunny.obj", V, F);
    LibShell::MeshConnectivity mesh(F);
    Eigen::VectorXd edgeDOFs;
    LibShell::MidedgeAverageFormulation::initializeExtraDOFs(edgeDOFs, mesh, V);
    LibShell::MaterialModel<LibShell::MidedgeAverageFormulation>* mat = new LibShell::StVKMaterial<LibShell::MidedgeAverageFormulation>();;
    mat = new LibShell::StVKMaterial<LibShell::MidedgeAverageFormulation>();

    double thickness = 0.1;
    double dt = 0.001;
    double possion_ratio = 0.5;
    double lameAlpha, lameBeta;
    double young = 1.0; // doesn't matter for static solves
    lameAlpha = young * possion_ratio / (1.0 - possion_ratio * possion_ratio);
    lameBeta = young / 2.0 / (1.0 + possion_ratio);

    //LibShell::MonolayerRestState* restState = new LibShell::MonolayerRestState;
    LibShell::MonolayerRestState restState;
    //std::cout << "mesh.nFaces() : " << mesh.nFaces() << std::endl;
    restState.thicknesses.resize(mesh.nFaces(), thickness);
     //initialize first fundamental forms to those of input mesh
    LibShell::ElasticShell<LibShell::MidedgeAverageFormulation>::firstFundamentalForms(mesh, V, restState.abars);
    restState.bbars.resize(mesh.nFaces());
    for (int i = 0; i < mesh.nFaces(); i++)
        restState.bbars[i].setZero();
    restState.lameAlpha.resize(mesh.nFaces(), lameAlpha);
    restState.lameBeta.resize(mesh.nFaces(), lameBeta);


    int nFaces = mesh.nFaces();
    int nVerts = (int)V.rows();
    int nEdges = mesh.nEdges();
    //Eigen::Matrix<double, 1, 9> deriv;
    double energy = 0.0;


    Eigen::VectorXd derivative;
    std::vector<Eigen::Triplet<double>> hessian;
    derivative.resize(3 * nVerts + LibShell::MidedgeAverageFormulation::numExtraDOFs * nEdges);
    //derivative.resize(3 * nVerts);
    derivative.setZero();
    energy = LibShell::ElasticShell<LibShell::MidedgeAverageFormulation>::elasticEnergy(mesh, V, edgeDOFs, *mat, restState, &derivative, NULL);

    grad = derivative;
    return energy;
    ////optimization
    //for (int i = 0; i < nVerts; i++) {
    //    V.row(i) -= derivative.segment<3>(3 * i) * dt;
    //}
    //for (int iter = 0; iter < 10000; iter++) {
    //    Eigen::Matrix<double, 1, 9> deriv;
    //    double energy = 0.0;
    //    //double energy = mat->stretchingEnergy(mesh, V, *restState, 0, &deriv, NULL);
    //    Eigen::VectorXd derivative;
    //    int n = 3 * nVerts + LibShell::MidedgeAverageFormulation::numExtraDOFs * nEdges;
    //    derivative.resize(3 * nVerts + LibShell::MidedgeAverageFormulation::numExtraDOFs * nEdges);
    //    derivative.setZero();

    //    energy = LibShell::ElasticShell<LibShell::MidedgeAverageFormulation>::elasticEnergy(mesh, V, edgeDOFs, *mat, restState, &derivative, nullptr);
    //    std::cout << "The energy is : " << energy << std::endl;
    //    for (int i = 0; i < nVerts; i++) {
    //        V.row(i) -= derivative.segment<3>(3 * i) * dt;
    //    }
    //    std::cout << "The energy is : " << energy << std::endl;
    //    if ((iter + 1) % 1000 == 0)
    //        igl::writeOBJ("D://IDE//VisualStudio_Project//DR_manufacturing//TestOn3D//Test3D//results//bunny_" + std::to_string(int(iter / 1000)) + ".obj", V, F);
    //}

}