/*
* EP：没有使用GEP，仅仅求解EP，且使用的是p的基向量矩阵，也就是说是Pp
* 为什么解不为0，是因为：
* f的能正确求解应该是：f = NPp
* 而这份代码直接求解了：f=Np，故而f不为0，甚至模长为1
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
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseRegularInverse.h>
#include <chrono>
using namespace std;
using namespace Eigen;
using namespace Spectra;

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
//Build Sum Matrix G
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
//Build Normal Matrix N
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
    char* Stress = (char*)"../Results/Eworst_stress_test_00.txt";

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
    //start time
    auto start_time = std::chrono::high_resolution_clock::now();

    //build Stiffness Matrix
    Eigen::SparseMatrix<double> N = Build_N(np, nbe, Boundary_Pid, Normals);
    Eigen::SparseMatrix<double> G = Build_G(np, points, N, center);
    Eigen::SparseMatrix<double> K = Build_stiffness_Matrix(np, points, nt, triangles);


    // 定义矩阵-向量乘法操作
    //构造矩阵B，但是由于spectra要求B为正定矩阵，所以没用到
    class MySparseMatrixOpB {
    public:
        using Scalar = double;
    private:
        Eigen::SparseMatrix<double> N;  
        Eigen::SparseMatrix<double> G;  
        int nbe;
        int np;
        using Index = Eigen::Index;
        using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        using MapConstVec = Eigen::Map<const Vector>;
        using MapVec = Eigen::Map<Vector>;
        mutable CompInfo m_info;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> m_cg;

    public:
        using Scalar = double;  // A typedef named "Scalar" is required
        int rows() const { return nbe; }
        int cols() const { return nbe; }
        // y_out = M * x_in
        // Constructor to initialize your sparse matrices
        MySparseMatrixOpB() {};
        MySparseMatrixOpB(const Eigen::SparseMatrix<double>& inputMatrixA,
            const Eigen::SparseMatrix<double>& inputMatrixB, int NumberOfBoundary, int NumberOfPoints)
            : N(inputMatrixA), G(inputMatrixB), nbe(NumberOfBoundary), np(NumberOfPoints) {
            Eigen::MatrixXd GGT = G * G.transpose();
            Eigen::MatrixXd GGT_inverse = GGT.inverse();
            Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
            Eigen::SparseMatrix<double> I(nbe, nbe);
            for (int i = 0; i < nbe; i++) {
                I.insert(i, i) = 1.0;
            }
            m_cg.compute((I - G.transpose() * (sparse_GGT_inverse)*G));

        }

        void solve(const double* x_in, double* y_out) const
        {
            Eigen::Map<const Vector> x(x_in, N.cols());
            MapVec y(y_out, N.cols());

            y.noalias() = m_cg.solve(x);

            m_info = (m_cg.info() == Eigen::Success) ?
                CompInfo::Successful :
                CompInfo::NotConverging;
            if (m_info != CompInfo::Successful)
                throw std::runtime_error("SparseRegularInverse: CG solver does not converge");
        }
        // Define the multiplication operation
        void perform_op(const double* x_in, double* y_out)const {
            Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(x_in, N.cols());
            Eigen::SparseMatrix<double> I(nbe, nbe);
            for (int i = 0; i < nbe; i++) {
                I.insert(i, i) = 1.0;
            }
            Eigen::MatrixXd GGT = G * G.transpose();
            Eigen::MatrixXd GGT_inverse = GGT.inverse();
            Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
            Eigen::VectorXd y = (I - G.transpose() * (sparse_GGT_inverse)*G) * x;  // Example operation: first multiply with B, then with A
            std::copy(y.data(), y.data() + y.size(), y_out);
        }
    };
    //构造矩阵A
    class MySparseMatrixOpA {
    private:
        Eigen::SparseMatrix<double> N;
        Eigen::SparseMatrix<double> G;
        Eigen::SparseMatrix<double> K;
        int nbe;
        int np;

    public:
        using Scalar = double;  // A typedef named "Scalar" is required
        int rows() const { return nbe; }
        int cols() const { return nbe; }
        // y_out = M * x_in
        // Constructor to initialize your sparse matrices
        MySparseMatrixOpA() {};
        MySparseMatrixOpA(const Eigen::SparseMatrix<double>& Stiffness_Matrix,
            const Eigen::SparseMatrix<double>& inputMatrixA,
            const Eigen::SparseMatrix<double>& inputMatrixB, int NumberOfBoundary, int NumberOfPoints)
            :K(Stiffness_Matrix), N(inputMatrixA), G(inputMatrixB), nbe(NumberOfBoundary), np(NumberOfPoints) {}
        // Define the multiplication operation
        void perform_op(const double* x_in, double* y_out)const {
            Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(x_in, N.cols());
            Eigen::SparseMatrix<double> I(nbe, nbe);
            for (int i = 0; i < nbe; i++) {
                I.insert(i, i) = 1.0;
            }
            Eigen::MatrixXd GGT = G * G.transpose();
            Eigen::MatrixXd GGT_inverse = GGT.inverse();
            Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
            Eigen::VectorXd y = (I - G.transpose() * (sparse_GGT_inverse)*G) *
                N.transpose() * K * N *
                (I - G.transpose() * (sparse_GGT_inverse)*G) * x;
            std::copy(y.data(), y.data() + y.size(), y_out);
        }
    };

    MySparseMatrixOpA opA(K, N, G, nbe, np);
    //MySparseMatrixOpB opB(N,G,nbe, np);
    SymEigsSolver<MySparseMatrixOpA> eigs(opA, 1, 6);
    //SymGEigsSolver<MySparseMatrixOpA, MySparseMatrixOpB, GEigsMode::RegularInverse>
    //    eigs(opA, opB, 3, 6);
    eigs.init();
    eigs.compute(SortRule::SmallestAlge);
    Eigen::MatrixXd evecs;
    Eigen::VectorXd evalues;
    if (eigs.info() == CompInfo::Successful)
    {
        evalues = eigs.eigenvalues();
        evecs = eigs.eigenvectors();
        std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    }
    else {
        cout << "The Eigenvalues were not found.Error:NotConvergence" << endl;
        return 0;
    }
    // calculate the time cost;
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Time taken by function: " << duration << " seconds" << std::endl;

    //After calculate Eigenvalues,we'll calculate the stress;
    VectorXd u, f, p;
    f = N * evecs.col(0);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return -1;
    }
    u = solver.solve(f);
    //cout << "Eigenvalues" << u << endl;

    SparseMatrix<double> Stress_face = Calculate_Stresses_face(np, points, nt, triangles, u);

    VectorXd S_face_points(nt);
    for (int k = 0; k < nt; k++) {
        double sigma_k = max(abs(Stress_face.coeff(3 * k, 0)), abs(Stress_face.coeff(3 * k + 1, 0)));
        S_face_points(k) = sigma_k;
    }

    FILE* inm_stress_face;
    if (!(inm_stress_face = fopen(Stress, "w"))) {
        fprintf(stderr, "  ** UNABLE TO OPEN OUTPUT MESH FILE.\n");
        exit(EXIT_FAILURE);
    }
    for (int k = 0; k < nt; k++) {
        fprintf(inm_stress_face, "%lf\n", S_face_points(k, 0));
    }
    fclose(inm_stress_face);
    ////////////////////////////////////////////////////////////////////////////////////////////

}