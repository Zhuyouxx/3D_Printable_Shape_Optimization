/*
* 该文件用于计算微分，且将模型顶点按照梯度的方向进行移动。
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
using namespace std;
using namespace Eigen;
using namespace Spectra;

#include "mmg/mmg2d/libmmg2d.h"
typedef Eigen::Matrix< double, Dynamic, Dynamic> MatrixXd;
double mu = 0.3;
double E = 1e9; //double E = 1e9;
double dt = 1e5;
class MySparseMatrixOpA {
private:
    Eigen::SparseMatrix<double> N;
    Eigen::SparseMatrix<double> G;
    Eigen::SparseMatrix<double> K;
    int nbe;
    int np;
    double M;

public:
    using Scalar = double;  // A typedef named "Scalar" is required
    int rows() const { return nbe; }
    int cols() const { return nbe; }
    // y_out = M * x_in
    // Constructor to initialize your sparse matrices
    MySparseMatrixOpA() {};
    MySparseMatrixOpA(const Eigen::SparseMatrix<double>& Stiffness_Matrix,
        const Eigen::SparseMatrix<double>& inputMatrixA,
        const Eigen::SparseMatrix<double>& inputMatrixB, int NumberOfBoundary, int NumberOfPoints, double MaxEigenValueOfK)
        :K(Stiffness_Matrix), N(inputMatrixA), G(inputMatrixB), nbe(NumberOfBoundary), np(NumberOfPoints), M(MaxEigenValueOfK) {}
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
        Eigen::VectorXd y = ((I - G.transpose() * (sparse_GGT_inverse)*G) *
            N.transpose() * K * N *
            (I - G.transpose() * (sparse_GGT_inverse)*G) +
            G.transpose() * M * G
            ) * x;
        std::copy(y.data(), y.data() + y.size(), y_out);
    }
};

template<typename Scalar>
Eigen::SparseMatrix<Scalar> ConvertSparseMatrix(const Eigen::SparseMatrix<CppAD::AD<Scalar>>& input) {
    Eigen::SparseMatrix<Scalar> output(input.rows(), input.cols());
    output.reserve(input.nonZeros());

    // 遍历每个元素进行转换
    for (int k = 0; k < input.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<CppAD::AD<Scalar>>::InnerIterator it(input, k); it; ++it) {
            // 仅提取CppAD::AD<Scalar>值的实数部分
            output.insert(it.row(), it.col()) = CppAD::Value(it.value());
        }
    }

    return output;
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
        //G.insert(2, 2 * k) = -(vertices[2 * k + 1] - center[1]);//-r_y
        //G.insert(2, 2 * k + 1) = vertices[2 * k] - center[0];//r_x
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
        Normal_Of_pi.normalize();
        Normal_Of_pj.normalize();
        Normals[pi_id * 2] = Normal_Of_pi[0];
        Normals[pi_id * 2 + 1] = Normal_Of_pi[1];
        Normals[pj_id * 2] = Normal_Of_pj[0];
        Normals[pj_id * 2 + 1] = Normal_Of_pj[1];
    }
    for (int k = 0; k < ne; k++) {
        int id = Boundary_Pid[k];
        N.insert(2 * id, k) = -Normals[2 * id];
        N.insert(2 * id + 1, k) = -Normals[2 * id + 1];
    }
    return N;
}

SparseMatrix<double> Build_N_pure(int nv, int ne, int* Edges, double* points) {
    Eigen::SparseMatrix<double> N(nv * 2, ne);
    double* Normals = (double*)calloc(nv * 2, sizeof(double));
    //vector<CppAD::AD<double>> Normals(nv * 2, CppAD::AD<double>(0));
    int* Boundary_Pid = (int*)calloc(ne, sizeof(int));
    for (int k = 0; k < ne; k++) {
        int pi_id = Edges[2 * k];
        int pj_id = Edges[2 * k + 1];
        Boundary_Pid[k] = pi_id;
        //double pi_x, pi_y, pj_x, pj_y;
        double pi_x, pi_y, pj_x, pj_y;
        pi_x = points[pi_id * 2]; pi_y = points[pi_id * 2 + 1];
        pj_x = points[pj_id * 2]; pj_y = points[pj_id * 2 + 1];
        Vector2d Normal_Of_E, Normal_Of_pi, Normal_Of_pj;
        //Eigen::Matrix<double, 2, 1> Normal_Of_E, Normal_Of_pi, Normal_Of_pj;;
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
    for (int k = 0; k < ne; k++) {
        int id = Boundary_Pid[k];
        N.insert(2 * id, k) = -Normals[2 * id];
        N.insert(2 * id + 1, k) = -Normals[2 * id + 1];
    }
    return N;
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

Eigen::MatrixXd calculate_max_pressure(SparseMatrix<double> K, SparseMatrix<double> N,
    SparseMatrix<double> G, int nbe, int np, double MaxEigenvalue) {
    MySparseMatrixOpA opA(K, N, G, nbe, np, MaxEigenvalue);
    //MySparseMatrixOpB opB(N,G,nbe, np);
    SymEigsSolver<MySparseMatrixOpA> eigs(opA, 1, 6);
    eigs.init();
    eigs.compute(SortRule::SmallestAlge, 10000, 1e-6);
    Eigen::MatrixXd evecs;
    Eigen::VectorXd evalues;
    if (eigs.info() == CompInfo::Successful)
    {
        evalues = eigs.eigenvalues();
        evecs = eigs.eigenvectors();
        std::cout << "Eigenvalues found:\n" << evalues << std::endl;
        //std::cout << "Eigenvalues found:\n" << evecs << std::endl;
    }
    else {
        cout << "The Eigenvalues were not found.Error:NotConvergence" << endl;
        throw std::runtime_error("Cannot compute the result.");
    }
    return evecs;
}

//输出结果↓
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
int cal_stress_save(char* filename, double* points, int np, int nt, int nbe, int* triangles,
    SparseMatrix<double> N, SparseMatrix<double> G, SparseMatrix<double> K, VectorXd u) {
    SparseMatrix<double> Stress_face = Calculate_Stresses_face(np, points, nt, triangles, u);
    VectorXd S_face_points(nt);
    for (int k = 0; k < nt; k++) {
        double sigma_k = max(abs(Stress_face.coeff(3 * k, 0)), abs(Stress_face.coeff(3 * k + 1, 0)));
        S_face_points(k) = sigma_k;
    }
    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }

    // 遍历数组并写入每个点到文件
    for (int i = 0; i < nt; i++) {
        file << S_face_points(i, 0) << std::endl;
    }
    // 关闭文件流
    file.close();
    return 1;
}


int main() {
    string model_name = "test_20_out";
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
    int* Edges;
    int* Boundary_Pid;
    double* center;
    MMG5_pMesh mmgMesh = ReadFromMesh(filename, points, triangles, Edges, center);
    MMG2D_Get_meshSize(mmgMesh, &np, &nt, NULL, &nbe);
    cout << "Number of Boundary:" << nbe << endl;
    //build Stiffness Matrix

    int n = 2 * np;
    // 定义矩阵-向量乘法操作


    for (int iter_time = 0; iter_time < 50; iter_time++) {
        cout << "----------------The iteration:" << iter_time << "---------------" << endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        vector<CppAD::AD<double>> V(np * 2);
        vector<double> x(n);
        std::vector<size_t> row_indices;
        std::vector<size_t> col_indices;
        std::vector<CppAD::AD<double>> K_AD_Vector;
        std::vector<double> jac_K;

        Eigen::SparseMatrix<CppAD::AD<double>> N_AD;
        Eigen::SparseMatrix<CppAD::AD<double>> K_AD;
        Eigen::SparseMatrix<double> N;
        Eigen::SparseMatrix<double> G;
        Eigen::SparseMatrix<double> K;

        //points值赋给V和x
        for (size_t i = 0; i < np * 2; ++i) {
            V[i] = CppAD::AD<double>(points[i]); // 初始值
            x[i] = Value(V[i]); // 假设 V 已经被赋值
        }
        //定义变量
        CppAD::Independent(V);
        //计算V和K的关系
        K_AD = Build_stiffness_Matrix(np, V, nt, triangles);
        // 遍历 K 的非零元素
        row_indices.clear(); col_indices.clear(); K_AD_Vector.clear();
        for (int k = 0; k < K_AD.outerSize(); ++k) {
            for (Eigen::SparseMatrix<CppAD::AD<double>>::InnerIterator it(K_AD, k); it; ++it) {
                row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
                col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
                K_AD_Vector.push_back(CppAD::AD<double>(it.value()));// it.value()     // 非零元素的值
            }
        }
        CppAD::ADFun<double> func(V, K_AD_Vector);    // 创建 ADFun 对象
        jac_K = func.Jacobian(x);
        //转变为double的SparseMatrix
        //N = ConvertSparseMatrix(N_AD);
        K = ConvertSparseMatrix(K_AD);
        //计算G
        N = Build_N_pure(np, nbe, Edges, points);
        G = Build_G(np, points, N, center);

        //计算K的最大特征值.
        double MaxEigenvalue = calculate_max_eigenvalue(K);
        // 构造EP,计算最大势能的pressure
        MatrixXd evecs = calculate_max_pressure(K, N, G, nbe, np, MaxEigenvalue);
        VectorXd p = evecs.col(0);

        Eigen::SparseMatrix<double> I(nbe, nbe);
        for (int i = 0; i < nbe; i++) {
            I.insert(i, i) = 1.0;
        }
        Eigen::MatrixXd GGT = G * G.transpose();
        Eigen::MatrixXd GGT_inverse = GGT.inverse();
        Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
        Eigen::VectorXd y = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;
        //cout << "y：" << y << endl;
        cout << "模长(P * min_eigenvector)：" << y.norm() << endl;
        VectorXd u, f;
        y = y / y.norm();
        p = y;
        f = N * p;
        cout << "模长(f)：" << f.norm() << endl;
        LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
        lscg.compute(K);
        u = lscg.solve(f);
        //std::cout << "The solution is:\n" << u << std::endl;
        std::cout << "The norm is:\n" << u.norm() << std::endl;
        std::cout << "#iterations:     " << lscg.iterations() << std::endl;
        std::cout << "estimated error: " << lscg.error() << std::endl;

        CppAD::Independent(V);
        N_AD = Build_N(np, nbe, Edges, V);
        Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1> pressure(nbe);
        ///这里尚未赋值，需要在最后计算p之后再进行赋值
        for (int k = 0; k < p.size(); k++) {
            pressure[k] = p[k];
        }
        double Eworst = u.transpose() * f;
        cout << "Eworst of the " << iter_time << " iteration is :" << Eworst << endl;
        Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1> F_vector_Matrix = N_AD * pressure;
        // 创建一个 CppAD::vector 用于存储 F 的所有非零元素
        std::vector<CppAD::AD<double>> F_vector;
        for (int k = 0; k < F_vector_Matrix.outerSize(); ++k) {
            for (int num = 0; num < F_vector_Matrix.innerSize(); num++) {
                F_vector.push_back(F_vector_Matrix(num, k));
            }
        }
        CppAD::ADFun<double> N_V(V, F_vector);
        // compute derivative using operation sequence stored in f

        vector<double> jac_N(n * n); // Jacobian of f (m by n matrix)

        jac_N = N_V.Jacobian(x);      // Jacobian for operation sequence
        VectorXd KNp = u;
        VectorXd dEds_K(n);
        VectorXd dEds_N(n);
        cout << "KNp.norm():" << KNp.norm() << endl;
        //构造dNds
        Eigen::SparseMatrix<double> dNds(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                // 检查 jac_N 中对应于 (i,j) 位置的元素是否非零
                double value = jac_N[i * n + j]; // 行优先排列
                if (value != 0) { // 只处理非零元素
                    // 直接插入到稀疏矩阵中
                    //dNds.insert(j, i) = value;
                    dNds.insert(i, j) = value;
                }
            }
        }
        // 最后，调用makeCompressed来优化稀疏矩阵的存储
        dNds.makeCompressed();
        dEds_N = 2 * KNp.transpose() * dNds;
        //计算dKds
        int non_zero = K_AD_Vector.size();
        for (size_t i = 0; i < n; ++i) {
            Eigen::SparseMatrix<double> dKdsi(n, n);
            for (int it = 0; it < non_zero; it++) {
                if (jac_K[it * n + i] != 0) {
                    dKdsi.insert(row_indices[it], col_indices[it]) = jac_K[it * n + i];
                }
            }

            dKdsi.makeCompressed();
            double dEds_K_si = KNp.transpose() * dKdsi * KNp;
            dEds_K[i] = dEds_K_si;
        }
        VectorXd dEds = dEds_N - dEds_K;

        for (int i = 0; i < n; ++i) {
            if (abs(dEds[i] * dt) > 0.08) {
                cout << "check" << endl;
                cout << dEds[i] << endl;
                cout << dEds_N[i] << endl;
                cout << dEds_K[i] << endl;
            }
            points[i] -= dEds[i] * dt;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        cout << "The iteration:" << iter_time << "||Time taken by function: " << duration << " seconds" << endl << endl << endl;
    }


    return 0;

}