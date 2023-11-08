/*
* GEP：使用了Eigen的库求解了GEP
* PS:但是Eigen的GEP是不是只能算Dense Matrix,好像没有给Sparse Matrix的口，所以我稀疏矩阵转成dense了再算的。
* 另外这里用的 P=HH_T，也就是说用的是p的基向量矩阵
*
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <chrono>
using namespace std;
using namespace Eigen;
//using namespace Spectra;

#include "mmg/mmg2d/libmmg2d.h"
typedef Eigen::Matrix< double, Dynamic, Dynamic> MatrixXd;
double mu = 0.3;
double E = 1e9;
//read mesh use MMG
MMG5_pMesh ReadFromMesh(char* filename, double*& points, int*& triangles, double*& normals, int*& boundary_Pid, double*& center) {
    MMG5_pMesh      mmgMesh;
    MMG5_int        k, np, nt, ne, noe;
    MMG5_pSol       mmgSol;
    mmgMesh = NULL;
    mmgSol = NULL;
    MMG2D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
        MMG5_ARG_end);
    if (MMG2D_loadMesh(mmgMesh, filename) != 1)  exit(EXIT_FAILURE);
    if (MMG2D_Get_meshSize(mmgMesh, &np, &nt, NULL, &ne) != 1)
        exit(EXIT_FAILURE);
    //if(MMG2D_Get_numberOfNonBdyEdges(mmgMesh, &noe)!=1) exit(EXIT_FAILURE);
    int Triangle[3], Edge[2];
    double  Point[3], Sol;
    double* Points = (double*)calloc(np * 2, sizeof(double));
    int* Triangles = (int*)calloc((nt) * 3, sizeof(int));
    double* Normals = (double*)calloc(np * 2, sizeof(double));
    int* Boundary_Pid = (int*)calloc(ne, sizeof(int));
    double c[2] = { 0.0,0.0 };
    for (k = 0; k < np; k++) {
        /** b) Vertex recovering */
        if (MMG2D_Get_vertex(mmgMesh, &(Point[0]), &(Point[1]), NULL, NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        Points[k * 2] = Point[0];
        Points[k * 2 + 1] = Point[1];
        c[0] += Point[0]; c[1] += Point[1];
    }
    c[0] = c[0] / np; c[1] = c[1] / np;
    center = c;
    for (k = 0; k < nt; k++) {
        if (MMG2D_Get_triangle(mmgMesh, &(Triangle[0]), &(Triangle[1]), &(Triangle[2]), NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        Triangles[k * 3] = Triangle[0] - 1;
        Triangles[k * 3 + 1] = Triangle[1] - 1;
        Triangles[k * 3 + 2] = Triangle[2] - 1;
    }
    for (k = 0; k < ne; k++) {
        if (MMG2D_Get_edge(mmgMesh, &(Edge[0]), &(Edge[1]), NULL, NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        int pi_id = Edge[0] - 1;
        int pj_id = Edge[1] - 1;
        Boundary_Pid[k] = pi_id;
        double pi_x, pi_y, pj_x, pj_y;
        pi_x = Points[pi_id * 2]; pi_y = Points[pi_id * 2 + 1];
        pj_x = Points[pj_id * 2]; pj_y = Points[pj_id * 2 + 1];
        Vector2d Normal_Of_E, Normal_Of_pi, Normal_Of_pj;
        Normal_Of_E << pj_y - pi_y, pi_x - pj_x;
        Normal_Of_E.normalize();
        Normal_Of_pi << Normals[pi_id * 2], Normals[pi_id * 2 + 1];
        Normal_Of_pj << Normals[pj_id * 2], Normals[pj_id * 2 + 1];
        Normal_Of_pi += Normal_Of_E;
        Normal_Of_pj += Normal_Of_E;
        Normal_Of_pi.normalize();
        Normal_Of_pj.normalize();
        Normals[pi_id * 2] = Normal_Of_pi[0];
        Normals[pi_id * 2 + 1] = Normal_Of_pi[1];
        Normals[pj_id * 2] = Normal_Of_pj[0];
        Normals[pj_id * 2 + 1] = Normal_Of_pj[1];
    }
    points = Points;
    triangles = Triangles;
    normals = Normals;
    boundary_Pid = Boundary_Pid;
    return mmgMesh;
}
//Build Stiffness Matrix-Sparse
Eigen::SparseMatrix<double> Build_stiffness_Matrix(int nv, double* vertices, int nt, int* triangles) {
    Eigen::SparseMatrix<double> sparse_K(nv * 2, nv * 2);
    Eigen::Matrix<double, 3, 3>D;
    Eigen::Matrix<double, 3, 6>B;
    Eigen::Matrix<double, 6, 6>Ke;
    D << 1, mu, 0,
        mu, 1, 0,
        0, 0, (1 - mu) / 2;
    D *= (E / (1 - mu * mu));
    for (int k = 0; k < nt; k++) {
        double p1_x, p1_y, p2_x, p2_y, p3_y, p3_x;
        int t1, t2, t3;
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        Matrix<double, 3, 3> Area;
        t1 = triangles[k * 3]; t2 = triangles[k * 3 + 1]; t3 = triangles[k * 3 + 2];
        p1_x = vertices[t1 * 2]; p1_y = vertices[t1 * 2 + 1];
        p2_x = vertices[t2 * 2]; p2_y = vertices[t2 * 2 + 1];
        p3_x = vertices[t3 * 2]; p3_y = vertices[t3 * 2 + 1];
        Area << 1, p1_x, p1_y,
            1, p2_x, p2_y,
            1, p3_x, p3_y;
        double A = Area.determinant() / 2.0;
        a1 = p2_x * p3_y - p3_x * p2_y;
        a2 = p3_x * p1_y - p1_x * p3_y;
        a3 = p1_x * p2_y - p2_x * p1_y;
        b1 = p2_y - p3_y;
        b2 = p3_y - p1_y;
        b3 = p1_y - p2_y;
        c1 = p3_x - p2_x;
        c2 = p1_x - p3_x;
        c3 = p2_x - p1_x;
        B << b1, 0, b2, 0, b3, 0,
            0, c1, 0, c2, 0, c3,
            c1, b1, c2, b2, c3, b3;
        B /= (2.0 * A);
        Ke = B.transpose() * D * B * A;
        int index[] = { 2 * t1,2 * t1 + 1,2 * t2,2 * t2 + 1,2 * t3,2 * t3 + 1 };
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) {
                sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
            }
    }
    return sparse_K;
}
//After solving EigenValues p/f,and get the displacement u,we calculate stress use the parameters u;
Eigen::SparseMatrix<double> Calculate_Stresses_face(int nv, double* vertices, int nt, int* triangles, VectorXd u) {
    Eigen::SparseMatrix<double> sparse_Stress(nt * 3, 1);
    Eigen::Matrix<double, 3, 3>D;
    Eigen::Matrix<double, 3, 6>B;
    Eigen::Matrix<double, 3, 1>S_tress_e;
    Eigen::Matrix<double, 6, 1>ue;
    D << 1, mu, 0,
        mu, 1, 0,
        0, 0, (1 - mu) / 2;
    D *= (E / (1 - mu * mu));
    for (int k = 0; k < nt; k++) {
        double p1_x, p1_y, p2_x, p2_y, p3_y, p3_x;
        int t1, t2, t3;
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        //double determinant = K.determinant()
        Matrix<double, 3, 3> Area;
        t1 = triangles[k * 3]; t2 = triangles[k * 3 + 1]; t3 = triangles[k * 3 + 2];
        p1_x = vertices[t1 * 2]; p1_y = vertices[t1 * 2 + 1];
        p2_x = vertices[t2 * 2]; p2_y = vertices[t2 * 2 + 1];
        p3_x = vertices[t3 * 2]; p3_y = vertices[t3 * 2 + 1];
        Area << 1, p1_x, p1_y,
            1, p2_x, p2_y,
            1, p3_x, p3_y;
        double A = Area.determinant() / 2.0;

        b1 = p2_y - p3_y;
        b2 = p3_y - p1_y;
        b3 = p1_y - p2_y;
        c1 = p3_x - p2_x;
        c2 = p1_x - p3_x;
        c3 = p2_x - p1_x;
        B << b1, 0, b2, 0, b3, 0,
            0, c1, 0, c2, 0, c3,
            c1, b1, c2, b2, c3, b3;
        B /= (2.0 * A);
        ue << u[t1 * 2], u[t1 * 2 + 1], u[t2 * 2], u[t2 * 2 + 1], u[t3 * 2], u[t3 * 2 + 1];
        S_tress_e = D * B * ue;

        //int index[] = { 2 * t1,2 * t1 + 1,2 * t2,2 * t2 + 1,2 * t3,2 * t3 + 1 };
        int index[] = { t1,t2,t3 };
        sparse_Stress.coeffRef(3 * k, 0) += S_tress_e(0, 0);
        sparse_Stress.coeffRef(3 * k + 1, 0) += S_tress_e(1, 0);
        sparse_Stress.coeffRef(3 * k + 2, 0) += S_tress_e(2, 0) * A;
    }
    return sparse_Stress;
}


SparseMatrix<double> Build_G(int nv, double* vertices, SparseMatrix<double> N, double* center) {
    Eigen::SparseMatrix<double> G(3, nv * 2);
    for (int k = 0; k < nv; k++) {
        G.insert(0, 2 * k) = 1;
        G.insert(1, 2 * k + 1) = 1;
        G.insert(2, 2 * k) = -(vertices[2 * k + 1] - center[1]);//-r_y
        G.insert(2, 2 * k + 1) = vertices[2 * k] - center[0];//r_x
    }
    G = G * N;
    return G;
}
SparseMatrix<double> Build_N(int nv, int nbp, int* Boundary_Pid, double* Normals) {
    Eigen::SparseMatrix<double> N(nv * 2, nbp);
    for (int k = 0; k < nbp; k++) {
        int id = Boundary_Pid[k];
        N.insert(2 * id, k) = -Normals[2 * id];
        N.insert(2 * id + 1, k) = -Normals[2 * id + 1];
    }
    return N;
}

int main() {
    //read mesh
    char* filename = (char*)"../IO/test_00.mesh";
    char* File_Stress = (char*)"../Results/Eworst_stress_test_00.txt";

    MMG5_pSol       mmgSol;
    MMG5_int        k, np, nt, nbe, noe, ne;
    double* points;
    int* triangles;
    double* Normals;
    int* Boundary_Pid;
    double* center;
    MMG5_pMesh mmgMesh = ReadFromMesh(filename, points, triangles, Normals, Boundary_Pid, center);
    MMG2D_Get_meshSize(mmgMesh, &np, &nt, NULL, &nbe);
    cout << "Number of Boundary:" << nbe << endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    //build Stiffness Matrix
    Eigen::SparseMatrix<double> N = Build_N(np, nbe, Boundary_Pid, Normals);
    Eigen::SparseMatrix<double> G = Build_G(np, points, N, center);
    Eigen::SparseMatrix<double> K = Build_stiffness_Matrix(np, points, nt, triangles);

    Eigen::SparseMatrix<double> I(nbe, nbe);
    for (int i = 0; i < nbe; i++) {
        I.insert(i, i) = 1.0;
    }
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::MatrixXd A = (I - G.transpose() * (sparse_GGT_inverse)*G) * N.transpose()
        * K *
        N * (I - G.transpose() * (sparse_GGT_inverse)*G);
    Eigen::MatrixXd B = (I - G.transpose() * (sparse_GGT_inverse)*G);
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver(A, B);

    solver.compute(A, B);
    Eigen::VectorXd min_eigenvector;
    if (solver.info() == Eigen::Success) {
        cout << "The (complex) numerators of the generalzied eigenvalues are: " << solver.alphas().transpose() << endl;
        cout << "The (real) denominatore of the generalzied eigenvalues are: " << solver.betas().transpose() << endl;
        cout << "The (complex) generalzied eigenvalues are (alphas./beta): " << solver.eigenvalues().real().transpose() << endl;

        //算最小特征值和对应的特征向量
        Eigen::VectorXd all_eigenvalues = solver.eigenvalues().real(); // 考虑实部，因为特征值可能是复数
        Eigen::Index min_index;  // 用于存储最小特征值的索引
        double min_eigenvalue = all_eigenvalues.minCoeff(&min_index);

        min_eigenvector = solver.eigenvectors().col(min_index).real();  // 考虑实部

        std::cout << "The minimum eigenvalue of the generalized eigenvalue problem is:\n";
        std::cout << min_eigenvalue << std::endl;

        std::cout << "The corresponding eigenvector is:\n";
        std::cout << min_eigenvector << std::endl;

    }
    else {
        std::cout << "The eigenvalue computation did not convergence." << std::endl;
    }

    VectorXd pp = min_eigenvector;
    cout << "模长(min_eigenvector)：" << pp.norm() << endl;
    Eigen::VectorXd y = (I - G.transpose() * (sparse_GGT_inverse)*G) * pp;
    cout << "y：" << y << endl;
    cout << "模长(P * min_eigenvector)：" << y.norm() << endl;

    /// ////////////////////////////////////////
    //计算time cost
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Time taken by function: " << duration << " seconds" << std::endl;

    //After calculate Eigenvalues,we'll calculate the stress;
    VectorXd u, f, p;
    y = y / y.norm();
    f = N * y;
    cout << "模长(f)：" << f.norm() << endl;
    /*Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_Kuf;
    solver_Kuf.compute(K);
    if (solver_Kuf.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return -1;
    }
    u = solver_Kuf.solve(f);*/
    LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
    lscg.compute(K);
    u = lscg.solve(f);
    //cout << "Eigenvalues" << u << endl;

    SparseMatrix<double> Stress_face = Calculate_Stresses_face(np, points, nt, triangles, u);
    VectorXd S_face_points(nt);
    for (int k = 0; k < nt; k++) {
        double sigma_k = max(abs(Stress_face.coeff(3 * k, 0)), abs(Stress_face.coeff(3 * k + 1, 0)));
        S_face_points(k) = sigma_k;
    }

    FILE* inm_stress_face;
    if (!(inm_stress_face = fopen(File_Stress, "w"))) {
        fprintf(stderr, "  ** UNABLE TO OPEN OUTPUT MESH FILE.\n");
        exit(EXIT_FAILURE);
    }
    for (int k = 0; k < nt; k++) {
        fprintf(inm_stress_face, "%lf\n", S_face_points(k, 0));
    }
    fclose(inm_stress_face);
    ////////////////////////////////////////////////////////////////////////////////////////////

}