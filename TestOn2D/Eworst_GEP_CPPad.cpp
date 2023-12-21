/*
* 这份文件最终实现了梯度的求解 mark
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <set>
#include <algorithm>
#include <random>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseRegularInverse.h>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <chrono>
# include <cppad/cppad.hpp> // the CppAD package
# include <cppad/example/cppad_eigen.hpp>
# include <cppad/example/atomic_two/eigen_mat_inv.hpp>
# include <cppad/example/atomic_two/eigen_mat_mul.hpp>
using namespace std;
using namespace Eigen;
using namespace Spectra;

#include "mmg/mmg2d/libmmg2d.h"
typedef Eigen::Matrix< double, Dynamic, Dynamic> MatrixXd;
double mu = 0.3;
double E = 1.0; //double E = 1e9;
double dt = 1e-6;
double O_mu = 1e-3;
int count_time = 0;


void sortThreeInts(int& a, int& b, int& c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
}
int write_file(char* filename, double* points, int np) {
    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }

    // 遍历数组并写入每个点到文件
    for (int i = 0; i < 2 * np; i++) {
        file << points[i] << std::endl;
    }

    // 关闭文件流
    file.close();
    return 1;
}
int write_file(char* filename, VectorXd points, int np) {
    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }

    // 遍历数组并写入每个点到文件
    for (int i = 0; i < 2 * np; i++) {
        file << points[i] << std::endl;
    }

    // 关闭文件流
    file.close();
    return 1;
}
int write_file(char* filename, vector<double> points, int np) {
    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }

    // 遍历数组并写入每个点到文件
    for (int i = 0; i < 2 * np; i++) {
        file << points[i] << std::endl;
    }

    // 关闭文件流
    file.close();
    return 1;
}
int write_file(char* filename, int* triangles, int nt) {
    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }

    // 遍历数组并写入每个点到文件
    for (int i = 0; i < nt; i++) {
        file << triangles[3 * i] << std::endl;
        file << triangles[3 * i + 1] << std::endl;
        file << triangles[3 * i + 2] << std::endl;
    }

    // 关闭文件流
    file.close();
    return 1;
}

VectorXd calculate_Ax_b(SparseMatrix<double> K, VectorXd b) {
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The inverse of K cannot be computed !");
    }
    VectorXd col = solver.solve(b);
    return col;
}

MMG5_pMesh ReadFromMesh(char* filename, double*& points, int*& triangles, int*& edges, double*& center) {
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

    int* Edges = (int*)calloc(ne * 2, sizeof(int));
    double* c = (double*)calloc(2, sizeof(double));
    //double c[2] = { 0.0,0.0 };
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
        Edges[2 * k] = pi_id;
        Edges[2 * k + 1] = pj_id;
    }
    points = Points;
    triangles = Triangles;
    edges = Edges;
    return mmgMesh;
}
MMG5_pMesh ReadFromMesh(MMG5_pSol& MmgSol, char* filename, double*& points, int*& triangles, int*& edges, double*& center) {
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

    int* Edges = (int*)calloc(ne * 2, sizeof(int));
    double* c = (double*)calloc(2, sizeof(double));
    //double c[2] = { 0.0,0.0 };
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
        Edges[2 * k] = pi_id;
        Edges[2 * k + 1] = pj_id;
    }
    points = Points;
    triangles = Triangles;
    edges = Edges;
    MmgSol = mmgSol;
    return mmgMesh;
}
int remesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, double*& points, int np, int*& triangles, int*& edges, double*& center) {

    for (int i = 0; i < np; i++) {
        mmgMesh->point[i + 1].c[0] = points[2 * i];
        mmgMesh->point[i + 1].c[1] = points[2 * i + 1];
    }
    if (MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_angle, 1) != 1) {
        fprintf(stderr, "Unable to set angle detection parameter.\n");
        return -1;
    }
    if (MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_angleDetection, 30) != 1) {
        fprintf(stderr, "Unable to set angle detection parameter.\n");
        return -1;
    }
    if (MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmin, 0.035) != 1) {
        fprintf(stderr, "Unable to set max edge length.\n");
        return EXIT_FAILURE;
    }
    if (MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmax, 0.04) != 1) {
        fprintf(stderr, "Unable to set max edge length.\n");
        return EXIT_FAILURE;
    }
    if (MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_angle, 1) != 1) {
        fprintf(stderr, "Unable to set angle detection parameter.\n");
        return -1;
    }
    if (MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_angleDetection, 30) != 1) {
        fprintf(stderr, "Unable to set angle detection parameter.\n");
        return -1;
    }
    /** Higher verbosity level */
    MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_verbose, -1);


    /** 4) (not mandatory): check if the number of given entities match with mesh size */
    if (MMG2D_Chk_meshData(mmgMesh, mmgSol) != 1) exit(EXIT_FAILURE);
    int             ier;
    ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);
    //ier = MMG2D_mmg2dmesh(mmgMesh, mmgSol);
    if (ier == MMG5_STRONGFAILURE) {
        fprintf(stdout, "BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH\n");
        return(ier);
    }
    else if (ier == MMG5_LOWFAILURE)
        fprintf(stdout, "BAD ENDING OF MMG2DLIB\n");

    MMG5_int        k, new_np, nt, ne, noe;
    if (MMG2D_Get_meshSize(mmgMesh, &new_np, &nt, NULL, &ne) != 1)
        exit(EXIT_FAILURE);

    std::string mesh_save_str = "Results/new_mesh.mesh";
    char* outname = const_cast<char*>(mesh_save_str.c_str());

    if (MMG2D_saveMesh(mmgMesh, outname) != 1)
        exit(EXIT_FAILURE);
    if (MMG2D_saveSol(mmgMesh, mmgSol, outname) != 1)  exit(EXIT_FAILURE);
    //if(MMG2D_Get_numberOfNonBdyEdges(mmgMesh, &noe)!=1) exit(EXIT_FAILURE);
    int Triangle[3], Edge[2];
    double  Point[3], Sol;
    double* Points = (double*)calloc(new_np * 2, sizeof(double));
    int* Triangles = (int*)calloc((nt) * 3, sizeof(int));

    int* Edges = (int*)calloc(ne * 2, sizeof(int));
    double* c = (double*)calloc(2, sizeof(double));
    //double c[2] = { 0.0,0.0 };
    for (k = 0; k < new_np; k++) {
        /** b) Vertex recovering */
        if (MMG2D_Get_vertex(mmgMesh, &(Point[0]), &(Point[1]), NULL, NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        Points[k * 2] = Point[0];
        Points[k * 2 + 1] = Point[1];
        c[0] += Point[0]; c[1] += Point[1];
    }
    c[0] = c[0] / new_np; c[1] = c[1] / new_np;
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
        Edges[2 * k] = pi_id;
        Edges[2 * k + 1] = pj_id;
    }
    points = Points;
    triangles = Triangles;
    edges = Edges;
    return 1;
}
Eigen::SparseMatrix<CppAD::AD<double>> Build_stiffness_Matrix(int nv, const vector<CppAD::AD<double>>& vertices, int nt, int* triangles) {
    Eigen::SparseMatrix<CppAD::AD<double>> sparse_K(nv * 2, nv * 2);
    Eigen::Matrix<CppAD::AD<double>, 3, 3>D;
    Eigen::Matrix<CppAD::AD<double>, 3, 6>B;
    Eigen::Matrix<CppAD::AD<double>, 6, 6>Ke;
    D << 1, mu, 0,
        mu, 1, 0,
        0, 0, (1 - mu) / 2;
    D *= (E / (1 - mu * mu));
    for (int k = 0; k < nt; k++) {
        CppAD::AD<double> p1_x, p1_y, p2_x, p2_y, p3_y, p3_x;
        int t1, t2, t3;
        CppAD::AD<double> a1, a2, a3, b1, b2, b3, c1, c2, c3;
        Matrix<CppAD::AD<double>, 3, 3> Area;
        t1 = triangles[k * 3]; t2 = triangles[k * 3 + 1]; t3 = triangles[k * 3 + 2];
        p1_x = vertices[t1 * 2]; p1_y = vertices[t1 * 2 + 1];
        p2_x = vertices[t2 * 2]; p2_y = vertices[t2 * 2 + 1];
        p3_x = vertices[t3 * 2]; p3_y = vertices[t3 * 2 + 1];
        Area << 1, p1_x, p1_y,
            1, p2_x, p2_y,
            1, p3_x, p3_y;
        CppAD::AD<double> A = Area.determinant() / 2.0;
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
                //sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                if (Ke(i, j) != 0) {
                    sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                }
            }

    }
    sparse_K.makeCompressed();
    return sparse_K;
}
Eigen::SparseMatrix<double> Build_stiffness_Matrix(int nv, double*& vertices, int nt, int* triangles) {
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
                //sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                if (Ke(i, j) != 0) {
                    sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                }
            }
    }
    sparse_K.makeCompressed();
    return sparse_K;
}

SparseMatrix<double> Build_G(int nv, double* vertices, SparseMatrix<double> N) {
    Eigen::SparseMatrix<double> G(3, nv * 2);
    for (int k = 0; k < nv; k++) {
        G.insert(0, 2 * k) = 1;
        G.insert(1, 2 * k + 1) = 1;
        //G.insert(2, 2 * k) = -(vertices[2 * k + 1] - center[1]);//-r_y
        //G.insert(2, 2 * k + 1) = vertices[2 * k] - center[0];//r_x
        G.insert(2, 2 * k) = -(vertices[2 * k + 1]);//-r_y
        G.insert(2, 2 * k + 1) = vertices[2 * k];//r_x
    }
    G = G * N;
    return G;
}
SparseMatrix<CppAD::AD<double>> Build_G(int nv, vector<CppAD::AD<double>> vertices, SparseMatrix<CppAD::AD<double>> N) {
    Eigen::SparseMatrix<CppAD::AD<double>> G(3, nv * 2);
    for (int k = 0; k < nv; k++) {
        G.insert(0, 2 * k) = 1;
        G.insert(1, 2 * k + 1) = 1;
        G.insert(2, 2 * k) = -(vertices[2 * k + 1]);//-r_y
        G.insert(2, 2 * k + 1) = vertices[2 * k];//r_x
    }
    G = G * N;
    return G;
}
SparseMatrix<CppAD::AD<double>> Build_N(int nv, int ne, int* Edges, const vector<CppAD::AD<double>>& Points) {
    Eigen::SparseMatrix<CppAD::AD<double>> N(nv * 2, ne);
    //CppAD::AD<double>* Normals = (CppAD::AD<double>*)calloc(nv * 2, sizeof(CppAD::AD<double>));
    vector<CppAD::AD<double>> Normals(nv * 2, CppAD::AD<double>(0));
    int* Boundary_Pid = (int*)calloc(ne, sizeof(int));
    for (int k = 0; k < ne; k++) {
        int pi_id = Edges[2 * k];
        int pj_id = Edges[2 * k + 1];
        Boundary_Pid[k] = pi_id;
        //double pi_x, pi_y, pj_x, pj_y;
        CppAD::AD<double> pi_x, pi_y, pj_x, pj_y;
        pi_x = Points[pi_id * 2]; pi_y = Points[pi_id * 2 + 1];
        pj_x = Points[pj_id * 2]; pj_y = Points[pj_id * 2 + 1];
        //Vector2d Normal_Of_E, Normal_Of_pi, Normal_Of_pj;
        Eigen::Matrix<CppAD::AD<double>, 2, 1> Normal_Of_E, Normal_Of_pi, Normal_Of_pj;;
        Normal_Of_E << pj_y - pi_y, pi_x - pj_x;
        Normal_Of_E.normalize();
        Normal_Of_pi << Normals[pi_id * 2], Normals[pi_id * 2 + 1];
        Normal_Of_pj << Normals[pj_id * 2], Normals[pj_id * 2 + 1];
        Normal_Of_pi[0] = Normal_Of_pi[0] + Normal_Of_E[0]; Normal_Of_pi[1] = Normal_Of_pi[1] + Normal_Of_E[1];
        Normal_Of_pj[0] = Normal_Of_pj[0] + Normal_Of_E[0]; Normal_Of_pj[1] = Normal_Of_pj[1] + Normal_Of_E[1];
        Normals[pi_id * 2] = Normal_Of_pi[0];
        Normals[pi_id * 2 + 1] = Normal_Of_pi[1];
        Normals[pj_id * 2] = Normal_Of_pj[0];
        Normals[pj_id * 2 + 1] = Normal_Of_pj[1];
    }
    for (int k = 0; k < ne; k++) {
        int id = Boundary_Pid[k];
        Eigen::Matrix<CppAD::AD<double>, 2, 1> Normal_Of_P;
        Normal_Of_P << Normals[2 * id], Normals[2 * id + 1];
        Normal_Of_P.normalize();
        N.insert(2 * id, k) = -Normal_Of_P[0];
        N.insert(2 * id + 1, k) = -Normal_Of_P[1];
    }
    return N;
}
SparseMatrix<double> Build_N(int nv, int ne, int* Edges, double* points, double*& Area) {
    Eigen::SparseMatrix<double> N(nv * 2, ne);
    double* Normals = (double*)calloc(nv * 2, sizeof(double));
    double* area = (double*)calloc(nv, sizeof(double));
    //vector<CppAD::AD<double>> Normals(nv * 2, CppAD::AD<double>(0));
    int* Boundary_Pid = (int*)calloc(ne, sizeof(int));
    int* Boundary_Pid_edge_num = (int*)calloc(ne, sizeof(int));
    for (int k = 0; k < ne; k++) {
        int pi_id = Edges[2 * k];
        int pj_id = Edges[2 * k + 1];
        Boundary_Pid[k] = pi_id;
        Boundary_Pid_edge_num[k] += 1;
        Boundary_Pid_edge_num[(k + 1) % ne] += 1;
        //double pi_x, pi_y, pj_x, pj_y;
        double pi_x, pi_y, pj_x, pj_y;
        pi_x = points[pi_id * 2]; pi_y = points[pi_id * 2 + 1];
        pj_x = points[pj_id * 2]; pj_y = points[pj_id * 2 + 1];
        Vector2d Normal_Of_E, Normal_Of_pi, Normal_Of_pj;
        //Eigen::Matrix<double, 2, 1> Normal_Of_E, Normal_Of_pi, Normal_Of_pj;;
        Normal_Of_E << pj_y - pi_y, pi_x - pj_x;
        double edge_norm = Normal_Of_E.norm();
        area[pi_id] += edge_norm;
        area[pj_id] += edge_norm;
        Normal_Of_E.normalize();
        Normal_Of_pi << Normals[pi_id * 2], Normals[pi_id * 2 + 1];
        Normal_Of_pj << Normals[pj_id * 2], Normals[pj_id * 2 + 1];
        Normal_Of_pi += Normal_Of_E;
        Normal_Of_pj += Normal_Of_E;
        Normals[pi_id * 2] = Normal_Of_pi[0];
        Normals[pi_id * 2 + 1] = Normal_Of_pi[1];
        Normals[pj_id * 2] = Normal_Of_pj[0];
        Normals[pj_id * 2 + 1] = Normal_Of_pj[1];
    }
    for (int k = 0; k < ne; k++) {
        int id = Boundary_Pid[k];
        area[id] /= Boundary_Pid_edge_num[k];
        Vector2d Normal_Of_P;
        Normal_Of_P << Normals[2 * id], Normals[2 * id + 1];
        Normal_Of_P.normalize();
        N.insert(2 * id, k) = -Normal_Of_P[0];
        N.insert(2 * id + 1, k) = -Normal_Of_P[1];
    }
    Area = area;
    return N;
}

SparseMatrix<double> Build_N_for_stress(int nv, int ne, int* Edges, double* points) {
    Eigen::SparseMatrix<double> N(nv * 2, ne);
    double* Normals = (double*)calloc(nv * 2, sizeof(double));
    //vector<CppAD::AD<double>> Normals(nv * 2, CppAD::AD<double>(0));
    int* Boundary_Pid = (int*)calloc(ne, sizeof(int));
    int* Boundary_Pid_edge_num = (int*)calloc(ne, sizeof(int));
    for (int k = 0; k < ne; k++) {
        int pi_id = Edges[2 * k];
        int pj_id = Edges[2 * k + 1];
        Boundary_Pid[k] = pi_id;
        Boundary_Pid_edge_num[k] += 1;
        Boundary_Pid_edge_num[(k + 1) % ne] += 1;
        //double pi_x, pi_y, pj_x, pj_y;
        double pi_x, pi_y, pj_x, pj_y;
        pi_x = points[pi_id * 2]; pi_y = points[pi_id * 2 + 1];
        pj_x = points[pj_id * 2]; pj_y = points[pj_id * 2 + 1];
        Vector2d Normal_Of_E, Normal_Of_pi, Normal_Of_pj;
        //Eigen::Matrix<double, 2, 1> Normal_Of_E, Normal_Of_pi, Normal_Of_pj;;
        Normal_Of_E << pj_y - pi_y, pi_x - pj_x;
        double edge_norm = Normal_Of_E.norm();
        Normal_Of_E.normalize();
        Normal_Of_pi << Normals[pi_id * 2], Normals[pi_id * 2 + 1];
        Normal_Of_pj << Normals[pj_id * 2], Normals[pj_id * 2 + 1];
        Normal_Of_pi += Normal_Of_E;
        Normal_Of_pj += Normal_Of_E;
        Normals[pi_id * 2] = Normal_Of_pi[0];
        Normals[pi_id * 2 + 1] = Normal_Of_pi[1];
        Normals[pj_id * 2] = Normal_Of_pj[0];
        Normals[pj_id * 2 + 1] = Normal_Of_pj[1];
    }
    for (int k = 0; k < ne; k++) {
        int id = Boundary_Pid[k];
        Vector2d Normal_Of_P;
        Normal_Of_P << Normals[2 * id], Normals[2 * id + 1];
        Normal_Of_P.normalize();
        N.insert(2 * id, k) = -Normal_Of_P[0];
        N.insert(2 * id + 1, k) = -Normal_Of_P[1];
    }
    return N;
}


Eigen::SparseMatrix<double> Calculate_Stresses(int nv, int ne, int* Edges, double* vertices, int nt, int* triangles,
    VectorXd p, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver) {
    int n = 2 * nv;
    Eigen::SparseMatrix<double> N = Build_N_for_stress(nv, ne, Edges, vertices);
    Eigen::SparseMatrix<double> G = Build_G(nv, vertices, N);
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<double> I(ne, ne);
    for (int i = 0; i < ne; i++) {
        I.insert(i, i) = 1.0;
    }
    p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    Eigen::VectorXd u(n + 3);
    Eigen::VectorXd f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 3);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0;
    u = solver.solve(f);
    // 检查解是否成功
    std::cout << "The norm of u is : " << u.norm() << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, 1> f_solve(2 * nv);
    //f_solve = K * u;
    //cout << "F_solve - f :" << (f_solve - f).norm() << endl;

    Eigen::SparseMatrix<double> sparse_Stress(nt * 3, 1);
    Eigen::Matrix<double, 3, 1>S_tress_e;
    Eigen::Matrix<double, 6, 1>ue;
    Eigen::Matrix<double, 3, 3>D;
    Eigen::Matrix<double, 3, 6>B;

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
        ue << u[t1 * 2], u[t1 * 2 + 1], u[t2 * 2], u[t2 * 2 + 1], u[t3 * 2], u[t3 * 2 + 1];
        S_tress_e = D * B * ue;
        sparse_Stress.coeffRef(3 * k, 0) += S_tress_e(0, 0);
        sparse_Stress.coeffRef(3 * k + 1, 0) += S_tress_e(1, 0);
        sparse_Stress.coeffRef(3 * k + 2, 0) += S_tress_e(2, 0);
    }
    sparse_Stress.makeCompressed();
    return sparse_Stress;
}
Eigen::SparseMatrix<double> merge_matrix(Eigen::SparseMatrix<double> K, Eigen::SparseMatrix<double> G) {
    // 行合并
    Eigen::SparseMatrix<double> Merged(K.rows() + G.rows(), K.rows() + G.rows());
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            Merged.insert(it.row(), it.col()) = it.value();
        }
    }
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
            Merged.insert(K.rows() + it.row(), it.col()) = it.value();
        }
    }
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
            Merged.insert(it.col(), K.rows() + it.row()) = it.value();
        }
    }
    Merged.makeCompressed();
    return Merged;
}

VectorXd Calculate_sigma_A(int nv, int ne, int* Edges, double* points, int nt, int* triangles, VectorXd u, double t) {
    //Eigen::SparseMatrix<double> stress(3, 2 * nv + 3);
    VectorXd stress(2 * nv + 3);
    VectorXd ue(6);
    VectorXd se(6);
    VectorXd Stress_e(3);
    VectorXd D_e(3);
    Eigen::Matrix<double, 3, 3>D;
    Eigen::Matrix<double, 3, 6>B;
    Eigen::Matrix<double, 3, 6>sigma;
    Matrix<double, 3, 3> Area;
    D << 1, mu, 0,
        mu, 1, 0,
        0, 0, (1 - mu) / 2;
    for (int k = 0; k < nt; k++) {
        double p1_x, p1_y, p2_x, p2_y, p3_y, p3_x;
        int t1, t2, t3;
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        t1 = triangles[k * 3]; t2 = triangles[k * 3 + 1]; t3 = triangles[k * 3 + 2];
        p1_x = points[t1 * 2]; p1_y = points[t1 * 2 + 1];
        p2_x = points[t2 * 2]; p2_y = points[t2 * 2 + 1];
        p3_x = points[t3 * 2]; p3_y = points[t3 * 2 + 1];
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
        sigma = D * B;
        int index[] = { 2 * t1,2 * t1 + 1,2 * t2,2 * t2 + 1,2 * t3,2 * t3 + 1 };

        ue << u[t1 * 2], u[t1 * 2 + 1], u[t2 * 2], u[t2 * 2 + 1], u[t3 * 2], u[t3 * 2 + 1];
        Stress_e = sigma * ue;
        MatrixXd tI_stress_1(2, 2), tI_stress_2(2, 2);
        MatrixXd O(2, 2);
        tI_stress_1 << t * 1.0 - Stress_e[0], -Stress_e[2], -Stress_e[2], t * 1.0 - Stress_e[1];
        tI_stress_2 << t * 1.0 + Stress_e[0], Stress_e[2], Stress_e[2], t * 1.0 + Stress_e[1];
        //cout << "-----------------------------------------------------" << endl;
        //cout << "Stress_e : " << Stress_e << endl;
        //cout << "tI_stress_1 : " << endl << tI_stress_1 << endl;
        //cout << "tI_stress_2 : " << endl << tI_stress_2 << endl;
        tI_stress_1 = tI_stress_1.inverse();
        tI_stress_2 = tI_stress_2.inverse();
        //cout << "tI_stress_1 : "<<endl<< tI_stress_1 << endl;
        //cout << "tI_stress_2 : " << endl << tI_stress_2 << endl;
        O = tI_stress_1 - tI_stress_2;
        D_e[0] = O.coeffRef(0, 0); D_e[1] = O.coeffRef(1, 1); D_e[2] = O.coeffRef(0, 1) + O.coeffRef(1, 0);
        se = D_e.transpose() * sigma;
        for (int s_i = 0; s_i < 6; s_i++) {
            stress[index[s_i]] += se[s_i];
        }
    }
    return stress;
}

double Calculate_Stresses_AD(int nv, int ne, int* Edges, double* points, int nt, int* triangles, double t, VectorXd p, VectorXd& grad_s, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver) {
    int n = 2 * nv;
    auto start_time_AD = std::chrono::high_resolution_clock::now();
    //Eigen::SparseMatrix<CppAD::AD<double>> N = Build_N_AD(nv, ne, Edges, vertices);
    //1-计算位移 u 
    Eigen::SparseMatrix<double> K = Build_stiffness_Matrix(nv, points, nt, triangles);
    Eigen::SparseMatrix<double> N = Build_N_for_stress(nv, ne, Edges, points);
    Eigen::SparseMatrix<double> G = Build_G(nv, points, N);
    Eigen::SparseMatrix<double> grad_M(2 * nv + 3, 2 * nv + 3);
    Eigen::VectorXd grad_O_s(n);

    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<double> I(ne, ne);
    for (int i = 0; i < ne; i++) {
        I.insert(i, i) = 1.0;
    }
    p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;
    VectorXd f = N * p;
    f.conservativeResize(f.size() + 3);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0;
    VectorXd u(2 * nv + 3);
    //计算位移displacement
    // 使用分解结果解决 Ax = b
    u = solver.solve(f);
    std::cout << "The norm of u is : " << u.norm() << std::endl;

    //-计算位移 u over -1

    //2-计算每个单元的stress
    //VectorXd stress_e(3);
    VectorXd ue(6);
    Eigen::VectorXd Stress_e;
    VectorXd grad_f_s(2 * nv + 3);
    Eigen::Matrix<double, 3, 3>D;
    Eigen::Matrix<double, 3, 6>B;
    Eigen::Matrix<double, 3, 6>sigma;
    vector<double> s_x(6);
    vector<CppAD::AD<double>> s(6);

    Eigen::Matrix<CppAD::AD<double>, 3, 3>D_ad;
    Eigen::Matrix<CppAD::AD<double>, 3, 6>B_ad;
    Eigen::Matrix<CppAD::AD<double>, 3, 6>sigma_s;
    Eigen::SparseMatrix<CppAD::AD<double>> N_AD;
    Eigen::SparseMatrix<CppAD::AD<double>> G_AD;
    std::vector<CppAD::AD<double>> Sigma_AD_Vector;
    std::vector<double> jac_Sigma_s;
    std::vector<double> jac_K;
    std::vector<double> jac_G;
    std::vector<size_t> row_indices;
    std::vector<size_t> col_indices;
    std::vector<size_t> G_row_indices;
    std::vector<size_t> G_col_indices;
    std::vector<CppAD::AD<double>> K_AD_Vector;
    std::vector<CppAD::AD<double>> G_AD_Vector;
    Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1> pressure(ne);
    ///这里尚未赋值，需要在最后计算p之后再进行赋值
    for (int k = 0; k < p.size(); k++) {
        pressure[k] = p[k];
    }
    D << 1, mu, 0,
        mu, 1, 0,
        0, 0, (1 - mu) / 2;
    D *= (E / (1 - mu * mu));
    D_ad << 1, mu, 0,
        mu, 1, 0,
        0, 0, (1 - mu) / 2;
    D_ad *= (E / (1 - mu * mu));
    //calculate M_Np
    auto start_time_Np = std::chrono::high_resolution_clock::now();

    VectorXd O_B = solver.solve(f);
    cout << "O_B:" << O_B.size() << endl;
    auto end_time_Np = std::chrono::high_resolution_clock::now();
    auto duration_Np = std::chrono::duration_cast<std::chrono::seconds>(end_time_Np - start_time_Np).count();
    std::cout << "||------ The cost of calculating the M-Np : " << duration_Np << " seconds ------||" << endl << endl;

    //calculate Grad_K
    vector<CppAD::AD<double>> vertices(2 * nv);
    vector<double> vertices_x(2 * nv);
    for (int i = 0; i < 2 * nv; i++) {
        vertices[i] = CppAD::AD<double>(points[i]);
        vertices_x[i] = points[i];
    }
    cout << "number of triangles: " << nt << endl;
    cout << "The number of vertices: " << n << endl;
    auto START_TIME2 = std::chrono::high_resolution_clock::now();
    CppAD::Independent(vertices);
    Eigen::SparseMatrix<CppAD::AD<double>> K_AD = Build_stiffness_Matrix(nv, vertices, nt, triangles);
    // 遍历 K 的非零元素
    auto end_time9 = std::chrono::high_resolution_clock::now();
    auto duration_Np9 = std::chrono::duration_cast<std::chrono::microseconds>(end_time9 - START_TIME2).count() / 1e6;
    std::cout << "||------ The cost of building KAD : " << duration_Np9 << " seconds ------||" << endl << endl;
    row_indices.clear(); col_indices.clear(); K_AD_Vector.clear();
    for (int k = 0; k < K_AD.outerSize(); ++k) {
        for (Eigen::SparseMatrix<CppAD::AD<double>>::InnerIterator it(K_AD, k); it; ++it) {
            row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
            col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
            K_AD_Vector.push_back(CppAD::AD<double>(it.value()));// it.value()     // 非零元素的值
        }
    }
    CppAD::ADFun<double> func(vertices, K_AD_Vector);    // 创建 ADFun 对象
    jac_K = func.Jacobian(vertices_x);
    auto end_time11 = std::chrono::high_resolution_clock::now();
    auto duration_Np11 = std::chrono::duration_cast<std::chrono::microseconds>(end_time11 - end_time9).count() / 1e6;
    std::cout << "||------ The cost of calculating the grad_K : " << duration_Np11 << " seconds ------||" << endl << endl;
    //calculate grad_G
    CppAD::Independent(vertices);
    N_AD = Build_N(nv, ne, Edges, vertices);
    G_AD = Build_G(nv, vertices, N_AD);
    G_row_indices.clear(); G_col_indices.clear(); G_AD_Vector.clear();
    for (int k = 0; k < G_AD.outerSize(); ++k) {
        for (Eigen::SparseMatrix<CppAD::AD<double>>::InnerIterator it(G_AD, k); it; ++it) {
            G_row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
            G_col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
            G_AD_Vector.push_back(CppAD::AD<double>(it.value()));// it.value()     // 非零元素的值
        }
    }
    CppAD::ADFun<double> G_func(vertices, G_AD_Vector);    // 创建 ADFun 对象
    jac_G = func.Jacobian(vertices_x);
    auto end_time12 = std::chrono::high_resolution_clock::now();
    auto duration_Np12 = std::chrono::duration_cast<std::chrono::microseconds>(end_time12 - end_time11).count() / 1e6;
    std::cout << "||------ The cost of calculating the grad_G : " << duration_Np12 << " seconds ------||" << endl << endl;


    VectorXd sigma_A = Calculate_sigma_A(nv, ne, Edges, points, nt, triangles, u, t);
    VectorXd O_A = solver.solve(sigma_A);
    for (int i = 0; i < ne; i++) {
        cout << O_A[i] << endl;
    }
    cout << "O_A:[0]: " << O_A[0] << "  " << O_A[5] << endl;

    //calculate grad_sigma(s)
    Matrix< CppAD::AD<double>, 3, 3> Area;
    for (int k = 0; k < nt; k++) {
        double p1_x, p1_y, p2_x, p2_y, p3_y, p3_x;
        int t1, t2, t3;
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        CppAD::AD<double>b1_ad, b2_ad, b3_ad, c1_ad, c2_ad, c3_ad;
        t1 = triangles[k * 3]; t2 = triangles[k * 3 + 1]; t3 = triangles[k * 3 + 2];
        s[0] = CppAD::AD<double>(points[t1 * 2]); s[1] = CppAD::AD<double>(points[t1 * 2 + 1]);
        s[2] = CppAD::AD<double>(points[t2 * 2]); s[3] = CppAD::AD<double>(points[t2 * 2 + 1]);
        s[4] = CppAD::AD<double>(points[t3 * 2]); s[5] = CppAD::AD<double>(points[t3 * 2 + 1]);
        p1_x = points[t1 * 2]; p1_y = points[t1 * 2 + 1];
        p2_x = points[t2 * 2]; p2_y = points[t2 * 2 + 1];
        p3_x = points[t3 * 2]; p3_y = points[t3 * 2 + 1];
        s_x[0] = p1_x; s_x[1] = p1_y; s_x[2] = p2_x; s_x[3] = p2_y; s_x[4] = p3_x; s_x[5] = p3_y;
        CppAD::Independent(s);
        Area << 1, s[0], s[1],
            1, s[2], s[3],
            1, s[4], s[5];
        CppAD::AD<double> A_AD = Area.determinant() / 2.0;
        double A = Value(A_AD);
        b1 = p2_y - p3_y; b1_ad = s[3] - s[5];
        b2 = p3_y - p1_y; b2_ad = s[5] - s[1];
        b3 = p1_y - p2_y; b3_ad = s[1] - s[3];
        c1 = p3_x - p2_x; c1_ad = s[4] - s[2];
        c2 = p1_x - p3_x; c2_ad = s[0] - s[4];
        c3 = p2_x - p1_x; c3_ad = s[2] - s[0];
        B << b1, 0, b2, 0, b3, 0,
            0, c1, 0, c2, 0, c3,
            c1, b1, c2, b2, c3, b3;
        B_ad << b1_ad, 0, b2_ad, 0, b3_ad, 0,
            0, c1_ad, 0, c2_ad, 0, c3_ad,
            c1_ad, b1_ad, c2_ad, b2_ad, c3_ad, b3_ad;
        B_ad /= (2.0 * A_AD);
        sigma_s = D_ad * B_ad;
        B /= (2.0 * A);
        sigma = D * B;
        Sigma_AD_Vector.clear();
        for (int sigma_i = 0; sigma_i < 3; sigma_i++)
            for (int sigma_j = 0; sigma_j < 6; sigma_j++)
                Sigma_AD_Vector.push_back(sigma_s.coeffRef(sigma_i, sigma_j));
        CppAD::ADFun<double> sigma_fun(s, Sigma_AD_Vector);
        jac_Sigma_s = sigma_fun.Jacobian(s_x);

        int index[] = { 2 * t1,2 * t1 + 1,2 * t2,2 * t2 + 1,2 * t3,2 * t3 + 1 };

        ue << u[t1 * 2], u[t1 * 2 + 1], u[t2 * 2], u[t2 * 2 + 1], u[t3 * 2], u[t3 * 2 + 1];
        Stress_e = sigma * ue;
        MatrixXd tI_stress_1(2, 2), tI_stress_2(2, 2);
        MatrixXd O_1(2, 2), O_2(2, 2);
        tI_stress_1 << t * 1.0 - Stress_e[0], -Stress_e[2], -Stress_e[2], t * 1.0 - Stress_e[1];
        tI_stress_2 << t * 1.0 + Stress_e[0], Stress_e[2], Stress_e[2], t * 1.0 + Stress_e[1];
        tI_stress_1 = tI_stress_1.inverse();
        tI_stress_2 = tI_stress_2.inverse();
        VectorXd M_Np_select(6);
        for (int t_i = 0; t_i < 6; t_i++)
            M_Np_select[t_i] = O_B[index[t_i]];
        for (int k_i = 0; k_i < 6; k_i++) {
            Eigen::MatrixXd grad_sigma_s_k(3, 6);
            for (int sigma_i = 0; sigma_i < 3; sigma_i++)
                for (int sigma_j = 0; sigma_j < 6; sigma_j++)
                    grad_sigma_s_k.coeffRef(sigma_i, sigma_j) = jac_Sigma_s[(sigma_i * 6 + sigma_j) * 6 + k_i];

            VectorXd part_1 = grad_sigma_s_k * M_Np_select;
            MatrixXd sigma_matrix(2, 2);
            MatrixXd D_part_1(2, 2);
            sigma_matrix << part_1[0], part_1[2], part_1[2], part_1[1];
            O_1 = tI_stress_1 * sigma_matrix; O_2 = tI_stress_2 * sigma_matrix;
            D_part_1 = O_1 - O_2;
            grad_O_s[index[k_i]] += O_mu * (D_part_1.coeffRef(0, 0) + D_part_1.coeffRef(1, 1));
        }
    }
    for (int s_i = 0; s_i < n; s_i++) {
        vector<CppAD::AD<double>> s_k(1);
        vector<double> s_k_x(1);
        s_k[0] = CppAD::AD<double>(points[s_i]);
        s_k_x[0] = points[s_i];
        CppAD::Independent(s_k);
        vertices[s_i] = s_k[0];
        N_AD = Build_N(nv, ne, Edges, vertices);
        Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1> F_vector_Matrix = N_AD * pressure;
        std::vector<CppAD::AD<double>> F_vector;
        for (int k_i = 0; k_i < F_vector_Matrix.outerSize(); ++k_i)
            for (int num = 0; num < F_vector_Matrix.innerSize(); num++)
                F_vector.push_back(F_vector_Matrix(num, k_i));
        CppAD::ADFun<double> N_s(s_k, F_vector);
        vector<double> jac_N(2 * nv);
        jac_N = N_s.Jacobian(s_k_x);
        for (int f_i = 0; f_i < 2 * nv; f_i++)
            grad_f_s[f_i] = jac_N[f_i];

        grad_M.setZero();
        int non_zero = K_AD_Vector.size();
        for (int it = 0; it < non_zero; it++)
            if (jac_K[it * n + s_i] != 0) {
                grad_M.insert(row_indices[it], col_indices[it]) = jac_K[it * n + s_i];
            }
        grad_M.makeCompressed();
        non_zero = G_AD_Vector.size();
        for (int it = 0; it < non_zero; it++)
            if (jac_G[it * n + s_i] != 0) {
                grad_M.insert(G_row_indices[it] + n, G_col_indices[it]) = jac_G[it * n + s_i];
                grad_M.insert(G_col_indices[it], G_row_indices[it] + n) = jac_G[it * n + s_i];
            }
        grad_M.makeCompressed();
        double part2 = O_A.transpose() * grad_f_s;
        double part3 = O_A.transpose() * grad_M * O_B;
        grad_O_s[s_i] += O_mu * (part2 - part3);
        vertices[s_i] = CppAD::AD<double>(points[s_i]);
    }

    grad_s = grad_O_s;
    auto end_time_stress = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time_stress - start_time_AD).count();
    std::cout << "|| The cost of d_stress/d_s : " << duration << " seconds" << endl << endl << endl;

    std::ostringstream stream;
    // 设置宽度为3，左侧填充0
    stream << std::setw(3) << std::setfill('0') << count_time;
    std::string grad_save_str = "Results/Grad_AD__" + stream.str() + ".txt";
    char* grad_save = const_cast<char*>(grad_save_str.c_str());
    if (!write_file(grad_save, grad_s, nv)) {
        cout << "Save FAILED!" << endl;
        return -1;
    }
    return 1.0;
}

bool isPositiveSemiDefinite(const Eigen::MatrixXd& A) {
    // 使用Eigen的自洽特征值计算方法
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(A);
    if (eigenSolver.info() != Eigen::Success) abort();

    // 如果所有特征值都是非负的，则矩阵是半正定的
    return eigenSolver.eigenvalues().minCoeff() >= 0;
}
CppAD::AD<double> Calculate_O(const CppAD::AD<double> t, Eigen::SparseMatrix<double> stress, int nt) {
    CppAD::AD<double> O = 0.0;
    for (int i = 0; i < nt; i++) {
        Eigen::Matrix<CppAD::AD<double>, 2, 2> sigma;
        Eigen::Matrix<CppAD::AD<double>, 2, 2> I;
        CppAD::AD<double> det1, det2;
        sigma << stress.coeff(3 * i, 0), stress.coeff(3 * i + 2, 0),
            stress.coeff(3 * i + 2, 0), stress.coeff(3 * i + 1, 0);
        I << t * 1.0, 0.0,
            0.0, t * 1.0;
        det1 = (I - sigma).determinant();
        det2 = (I + sigma).determinant();
        if (det1 < 1e-3 || det2 < 1e-3)
        {
            cout << t << endl;
            cout << stress.coeff(3 * i, 0) << "  " << stress.coeff(3 * i + 1, 0) << "  " << stress.coeff(3 * i + 2, 0) << endl;
            cout << "Invalid Guess s,t" << endl;
            return -1e6;
        }
        Eigen::Matrix<double, 2, 2> Sigma, tI;
        Sigma << stress.coeff(3 * i, 0), stress.coeff(3 * i + 2, 0),
            stress.coeff(3 * i + 2, 0), stress.coeff(3 * i + 1, 0);
        tI << Value(t) * 1.0, 0.0, 0.0, Value(t) * 1.0;
        Eigen::MatrixXd mat1 = tI - Sigma;
        Eigen::MatrixXd mat2 = tI + Sigma;
        if (!isPositiveSemiDefinite(mat1) || !isPositiveSemiDefinite(mat2)) {
            cout << t << endl;
            cout << stress.coeff(3 * i, 0) << "  " << stress.coeff(3 * i + 1, 0) << "  " << stress.coeff(3 * i + 2, 0) << endl;
            cout << "NOT Positive SemiDefinite." << endl;
            return -1e6;
        }
        O += -O_mu * (CppAD::log(det1) + CppAD::log(det2));
    }
    O += t;
    return O;
}

//计算最大特征值
double calculate_max_eigenvalue(SparseMatrix<double> K) {
    SparseSymMatProd<double> opK(K);
    SymEigsSolver<SparseSymMatProd<double>> eigs_K(opK, 1, 6);
    eigs_K.init();
    eigs_K.compute(SortRule::LargestAlge);
    Eigen::MatrixXd evecs;
    Eigen::VectorXd evalues;
    if (eigs_K.info() == CompInfo::Successful)
    {
        evalues = eigs_K.eigenvalues();
        evecs = eigs_K.eigenvectors();
        std::cout << "Eigenvalues of K is found:\n" << evalues << std::endl;
        //std::cout << "Eigenvalues found:\n" << evecs << std::endl;
    }
    else {
        cout << "The Eigenvalues of K were not found.Error:NotConvergence" << endl;
        return -1;
    }
    return evalues(0);
}

int main() {
    string model_name = "test_236_1";
    std::string path = "IO/";
    std::string filename_str = "IO/" + model_name + ".mesh";
    char* filename = const_cast<char*>(filename_str.c_str());
    std::string worst_stress_faces_str = "Results/Eworst_stress_face_" + model_name + ".txt";
    char* worst_stress_faces = const_cast<char*>(worst_stress_faces_str.c_str());
    MMG5_pSol       mmgSol;
    MMG5_int        k, np, nt, nbe, noe, ne;
    double* points;
    int* triangles;
    double* Normals;
    double* Area;
    int* Edges;
    int* Boundary_Pid;
    double* center;
    MMG5_pMesh mmgMesh = ReadFromMesh(mmgSol, filename, points, triangles, Edges, center);
    MMG2D_Get_meshSize(mmgMesh, &np, &nt, NULL, &nbe);
    std::cout << "Number of Boundary:" << nbe << endl;
    //Optimize
    for (int iter_time = 0; iter_time < 50; iter_time++) {
        int n = 2 * np;
        std::cout << "----------------The iteration:" << iter_time << "---------------" << endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        Eigen::SparseMatrix<double> N;
        Eigen::SparseMatrix<double> G;
        Eigen::SparseMatrix<double> K;
        //转变为double的SparseMatrix
        K = Build_stiffness_Matrix(np, points, nt, triangles);
        N = Build_N(np, nbe, Edges, points, Area);
        G = Build_G(np, points, N);        //计算G
        Eigen::SparseMatrix<double> M = merge_matrix(K, G);
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU;
        solver_LU.compute(M);
        if (solver_LU.info() != Eigen::Success) {
            std::cerr << "Decomposition failed!" << std::endl;
            throw std::runtime_error("The inverse of K cannot be computed !");
        }
        //计算K的最大特征值.
        double MaxEigenvalue = calculate_max_eigenvalue(K) * 10000.0;

        Eigen::SparseMatrix<double> I(nbe, nbe);
        for (int i = 0; i < nbe; i++) {
            I.insert(i, i) = 1.0;
        }
        Eigen::MatrixXd GGT = G * G.transpose();
        Eigen::MatrixXd GGT_inverse = GGT.inverse();
        Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();

        Eigen::SparseMatrix<double> GEP_matrix_sparse = ((I - G.transpose() * (sparse_GGT_inverse)*G).transpose() *
            N.transpose() * K * N *
            (I - G.transpose() * (sparse_GGT_inverse)*G) +
            G.transpose() * MaxEigenvalue * G
            );
        Eigen::SparseMatrix<double> A(n, n);
        for (int i = 0; i < np; i++) {
            A.insert(2 * i, 2 * i) = Area[i];
            A.insert(2 * i + 1, 2 * i + 1) = Area[i];
        }
        Eigen::SparseMatrix<double> matrix_B = N.transpose() * A * N;
        SparseSymMatProd<double> opA(GEP_matrix_sparse);
        SparseCholesky<double>  Bop(matrix_B);
        SymGEigsSolver<SparseSymMatProd<double>, SparseCholesky<double>, GEigsMode::Cholesky>
            geigs(opA, Bop, 3, 12);

        // Initialize and compute   GEP
        geigs.init();
        int nconv = geigs.compute(SortRule::SmallestAlge);

        // Retrieve results
        Eigen::VectorXd evalues;
        Eigen::MatrixXd evecs_GEP;

        if (geigs.info() == CompInfo::Successful)
        {
            evalues = geigs.eigenvalues();
            evecs_GEP = geigs.eigenvectors();
        }

        VectorXd p = evecs_GEP.col(2);
        VectorXd p_stress = evecs_GEP.col(2);

        //optimization
        bool loop = true, loop_t = true, loop_s = true;
        double* x = (double*)calloc(np * 2, sizeof(double));
        vector<CppAD::AD<double>> t(1);
        vector<double> t_x(1);

        for (size_t i = 0; i < np * 2; ++i) {
            x[i] = points[i]; // 
        }
        t[0] = 1e4; t_x[0] = 1e4;
        double threshold = 1e-5, gamma = 0.5;

        //Initial test: feasible initial guess s,t
        CppAD::Independent(t);
        Eigen::SparseMatrix<double> sparse_Stress = Calculate_Stresses(np, nbe, Edges, points, nt, triangles, p_stress, solver_LU);
        CppAD::AD<double> O = Calculate_O(t[0], sparse_Stress, nt);
        if (O < -1e5) {
            std::cout << "Invalid Initial!!" << endl;
            loop = false;
            throw std::runtime_error("The initial s,t are invalid!");
        }
        int iter = 0;
        CppAD::AD<double> last_O = O, last_O_t = O, last_O_s = O;
        std::vector<double> jac_O_t;
        double alpha = 0.001;
        while (loop) {
            iter++;
            //optimize t 
            loop_t = true;
            sparse_Stress = Calculate_Stresses(np, nbe, Edges, points, nt, triangles, p_stress, solver_LU);
            cout << "--------------------Optimize t: ---------------------------------------------------------------------------------" << endl;
            double alpha_t = 1.0;
            while (loop_t) {
                vector<CppAD::AD<double>> O_var(1);
                CppAD::Independent(t);
                O = Calculate_O(t[0], sparse_Stress, nt);
                O_var[0] = O;
                last_O_t = O;
                CppAD::ADFun<double> func(t, O_var);    // 创建 ADFun 对象
                jac_O_t = func.Jacobian(t_x);
                std::vector<double> w(1);
                w[0] = 1.0;  // 权重
                std::vector<double> hess_O_t = func.Hessian(t_x, w);
                double maxVal = 0.0;
                for (double val : jac_O_t) {
                    maxVal = std::max(maxVal, std::abs(val));
                }
                if (maxVal < threshold) { loop_t = false; break; }
                while (true) {
                    CppAD::AD<double> temp = t[0];
                    t[0] = t[0] - alpha_t * (jac_O_t[0] / hess_O_t[0]);
                    O = Calculate_O(t[0], sparse_Stress, nt);
                    if (O< last_O_t && O > -1e5) {
                        t_x[0] = Value(t[0]);
                        alpha_t /= gamma;
                        last_O_t = O;
                        break;
                    }
                    else {
                        t[0] = temp;
                        alpha_t *= gamma;
                    }
                }
                cout << "O:" << O << " The max of O_t : " << maxVal << "  t : " << t_x[0] << " dt:" << (jac_O_t[0] / hess_O_t[0]) << " |Alpha:" << alpha_t << endl;
            }

            //optimize s
            cout << "-------------------Optimize s: ----------Optimize s: ----------------- Optimize s --------------------------------------------------" << endl;
            cout << "The value of O after optomization on t : " << last_O_t << endl;
            VectorXd grad_s;
            Calculate_Stresses_AD(np, nbe, Edges, points, nt, triangles, t_x[0], p_stress, grad_s, solver_LU);
            sparse_Stress = Calculate_Stresses(np, nbe, Edges, points, nt, triangles, p_stress, solver_LU);
            O = Calculate_O(t[0], sparse_Stress, nt);
            last_O_s = O;
            VectorXd direction(n);
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> dis(-1.0f, 1.0f);
            for (int i = 0; i < np; i++) {
                double directionX = dis(gen);
                double directionY = dis(gen);
                // 计算向量长度
                double length = std::sqrt(directionX * directionX + directionY * directionY);
                // 规范化向量
                directionX /= length;
                directionY /= length;
                direction[2 * i] = directionX;
                direction[2 * i + 1] = directionY;
            }
            double eps = 1e-6;
            double O_last = Value(O);
            //做有限差分的验证
            for (int k = 0; k < 10; k++) {
                eps /= 2.0;
                for (int i = 0; i < n; i++) {
                    points[i] = (x[i] + eps * direction[i]);
                }
                //After updating points,we need to update K,N,G->M,so we need to rebuild a solver for new M.
                K = Build_stiffness_Matrix(np, points, nt, triangles);
                N = Build_N(np, nbe, Edges, points, Area);
                G = Build_G(np, points, N);        //计算G
                M = merge_matrix(K, G);
                Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU_new;
                solver_LU_new.compute(M);
                sparse_Stress = Calculate_Stresses(np, nbe, Edges, points, nt, triangles, p_stress, solver_LU_new);
                O = Calculate_O(t[0], sparse_Stress, nt);
                double O_new = Value(O);
                double dO = O_new - O_last;
                double dOdt = dO / eps;
                double grad_dot = grad_s.dot(direction);
                cout << "----------------------------------------- EPS: " << eps << "---------------------------------------" << endl;
                cout << "O_new : " << O_new << "   O : " << O_last << endl;
                cout << "dO : " << dO << endl;
                cout << "dOdt : " << dOdt << endl;
                cout << "Grad_dot : " << grad_dot << endl;
                cout << "--------------------------------------------------------------------------------" << endl;
            }

            double maxVal = 0.0;
            for (double val : grad_s) {
                maxVal = std::max(maxVal, std::abs(val));
            }
            if (maxVal < threshold) { loop = false; break; }
            if (maxVal * alpha > 0.05) {
                alpha *= gamma;
                continue;
            }
            while (true) {
                for (int i = 0; i < n; i++) {
                    points[i] = x[i] - alpha * grad_s[i];
                }
                //After updating points,we need to update K,N,G->M,so we need to rebuild a solver for new M.
                K = Build_stiffness_Matrix(np, points, nt, triangles);
                N = Build_N(np, nbe, Edges, points, Area);
                G = Build_G(np, points, N);        //计算G
                M = merge_matrix(K, G);
                Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU_new;
                solver_LU_new.compute(M);
                sparse_Stress = Calculate_Stresses(np, nbe, Edges, points, nt, triangles, p_stress, solver_LU_new);
                O = Calculate_O(t[0], sparse_Stress, nt);
                if (O< last_O_s && O > -1e5) {
                    for (int i = 0; i < n; i++) {
                        x[i] = points[i];
                    }
                    alpha /= gamma;
                    last_O_s = O;
                    break;
                }
                else {
                    cout << "The last O:" << last_O_s << "  O_new:" << O << endl;
                    for (int i = 0; i < n; i++) {
                        points[i] = x[i];
                    }
                    alpha *= gamma;
                }
            }

            cout << "The iteration : " << iter << " |O: " << last_O_s << " |max of O_s: " << maxVal << " | alpha :" << alpha << endl;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::cout << "The iteration:" << iter_time << "||Time taken by function: " << duration << " seconds" << endl << endl << endl;
    }


    return 0;

}