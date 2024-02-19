/**
 * Example of input output for the mmg3d library for multiple solutions at mesh
 * vertices
 *
 * \author Algiane Froehly (InriaSoft)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <string>
#include <iostream>
#include <vector>
#include <random>
#include <unordered_map>
#include <unordered_set>
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
#include <cppad/cppad.hpp> // the CppAD package
#include <cppad/example/cppad_eigen.hpp>
#include <cppad/example/atomic_two/eigen_mat_inv.hpp>
#include <cppad/example/atomic_two/eigen_mat_mul.hpp>
#include <boost/multiprecision/mpfr.hpp>
//#include <unsupported/Eigen/MPRealSupport>

using namespace std;
using namespace Eigen;
using namespace Spectra;
#include "mmg/mmg3d/libmmg3d.h"
namespace mp = boost::multiprecision;
typedef Eigen::Matrix< double, Dynamic, Dynamic> MatrixXd;
typedef Eigen::Matrix< CppAD::AD<double>, Dynamic, Dynamic> MatrixXd_AD;
typedef Eigen::Matrix<CppAD::AD<double>, 3, 1> Vector3d_AD;

typedef mp::mpfr_float mpfr_float;
typedef mp::mpfr_float_100 mpfr_float_100;
typedef CppAD::AD<mpfr_float> AD_mpfr_float;

typedef Eigen::Matrix< mpfr_float, Dynamic, Dynamic> MatrixX_mp;
typedef Eigen::Matrix<AD_mpfr_float, 3, 1> Vector3_mp_AD;
typedef Eigen::Matrix<mpfr_float, 3, 1> Vector3_mp;
typedef Eigen::Matrix<mpfr_float, Dynamic, 1> VectorX_mp;

typedef CppAD::AD<long double> AD_ld;
typedef Eigen::Matrix< long double , Dynamic, Dynamic > MatrixX_ld;
typedef Eigen::Matrix< long double, Dynamic, Dynamic> MatrixX_ld_AD;
typedef Eigen::Matrix<AD_ld, 3, 1> Vector3_ld_AD;
typedef Eigen::Matrix<long double, 3, 1> Vector3_ld;
typedef Eigen::Matrix<long double, Dynamic, 1> VectorX_ld;

double mu = 0.3;
double E = 1.0; //double E = 1e9;
double dt = 1e-6;
double O_mu = 1e-3;
int count_time = 0;
double edge_length = 0.05;
std::unordered_map<int, std::set<int>> adjacency_p;
std::unordered_map<int, std::set<int>> adjacency_tri;
unordered_map<int, int> Old_2_new;
int write_file(char* filename, double* points, int np) {
    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }

    // 遍历数组并写入每个点到文件
    for (int i = 0; i < 3 * np; i++) {
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
    for (int i = 0; i < 3 * np; i++) {
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
    for (int i = 0; i < 3 * np; i++) {
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
int read_mesh(MMG5_pMesh& mmgMesh_out, MMG5_pSol& mmgSol_out, char* filename, mpfr_float*& points, int*& triangles, int*& tetrahedras) {
    MMG5_pMesh      mmgMesh;
    MMG5_pSol       mmgSol, mmgMet, tmpSol;
    int             i, j, k;
    MMG5_int        np, n_tet, nprism, n_tri, nquad, na;

    mmgMesh = NULL;
    mmgSol = NULL;
    MMG3D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
        MMG5_ARG_end);
    mmgMesh_out = mmgMesh;
    mmgSol_out = mmgSol;
    if (MMG3D_loadMesh(mmgMesh, filename) != 1) {
        std::cerr << "Mesh loading failed!" << std::endl;
        return 0;
    }
    if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, &nprism, &n_tri, &nquad, &na) != 1)  exit(EXIT_FAILURE);
    cout << "number of points:" << np << endl;
    cout << "number of tetraheras:" << n_tet << endl;
    cout << "number of triagles on surface:" << n_tri << endl;
    cout << "nquad:" << nquad << endl;
    cout << "na:" << na << endl;
    int Triangle[3], Edge[2], Tetrahedra[4];
    double  Point[3], Sol;
    mpfr_float* Points = (mpfr_float*)calloc(np * 3, sizeof(mpfr_float));
    int* Triangles = (int*)calloc((n_tri) * 3, sizeof(int));
    int* Tetrahedras = (int*)calloc((n_tet) * 4, sizeof(int));

    for (k = 0; k < np; k++) {
        /** b) Vertex recovering */
        if (MMG3D_Get_vertex(mmgMesh, &(Point[0]), &(Point[1]), &(Point[2]), NULL, NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        Points[k * 3] = Point[0];
        Points[k * 3 + 1] = Point[1];
        Points[k * 3 + 2] = Point[2];
    }

    for (k = 0; k < n_tet; k++) {
        if (MMG3D_Get_tetrahedron(mmgMesh, &(Tetrahedra[0]), &(Tetrahedra[1]), &(Tetrahedra[2]), &(Tetrahedra[3]), NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        Tetrahedras[k * 4] = Tetrahedra[0] - 1;
        Tetrahedras[k * 4 + 1] = Tetrahedra[1] - 1;
        Tetrahedras[k * 4 + 2] = Tetrahedra[2] - 1;
        Tetrahedras[k * 4 + 3] = Tetrahedra[3] - 1;
    }
    for (k = 0; k < n_tri; k++) {
        if (MMG3D_Get_triangle(mmgMesh, &(Triangle[0]), &(Triangle[1]), &(Triangle[2]), NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        Triangles[k * 3] = Triangle[0] - 1;
        Triangles[k * 3 + 1] = Triangle[1] - 1;
        Triangles[k * 3 + 2] = Triangle[2] - 1;
    }
    points = Points;
    triangles = Triangles;
    tetrahedras = Tetrahedras;
    return 1;
}

void search_neighbor(int n_sp,int *index_sp,int n_tri,int *triangles) {
    for (int k = 0; k < n_tri; k++) {
        int v1 = triangles[3 * k];
        int v2 = triangles[3 * k + 1];
        int v3 = triangles[3 * k + 2];

        // 更新邻接关系
        adjacency_p[v1].insert(v2);
        adjacency_p[v1].insert(v3);
        adjacency_p[v2].insert(v1);
        adjacency_p[v2].insert(v3);
        adjacency_p[v3].insert(v1);
        adjacency_p[v3].insert(v2);

        adjacency_tri[v1].insert(k);
        adjacency_tri[v2].insert(k);
        adjacency_tri[v3].insert(k);

    }
        
}

Eigen::SparseMatrix<mpfr_float> merge_matrix(Eigen::SparseMatrix<mpfr_float> K, Eigen::SparseMatrix<mpfr_float> G) {
    // 行合并
    Eigen::SparseMatrix<mpfr_float> Merged(K.rows() + G.rows(), K.rows() + G.rows());
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<mpfr_float>::InnerIterator it(K, k); it; ++it) {
            Merged.insert(it.row(), it.col()) = it.value();
        }
    }
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<mpfr_float>::InnerIterator it(G, k); it; ++it) {
            Merged.insert(K.rows() + it.row(), it.col()) = it.value();
        }
    }
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<mpfr_float>::InnerIterator it(G, k); it; ++it) {
            Merged.insert(it.col(), K.rows() + it.row()) = it.value();
        }
    }
    Merged.makeCompressed();
    return Merged;
}
//计算行列式
template <typename T>
T determinant3x3(const Eigen::Matrix<T, 3, 3>& matrix) {
    return matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)) -
        matrix(0, 1) * (matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0)) +
        matrix(0, 2) * (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
}

template <typename T>
T determinant4x4(const Eigen::Matrix<T, 4, 4>& matrix) {
    T det = T(0);
    for (int i = 0; i < 4; i++) {
        Eigen::Matrix<T, 3, 3> subMatrix;
        for (int j = 1; j < 4; j++) {
            int columnIndex = 0;
            for (int k = 0; k < 4; k++) {
                if (k != i) {
                    subMatrix(j - 1, columnIndex) = matrix(j, k);
                    columnIndex++;
                }
            }
        }
        det += (i % 2 == 0 ? 1 : -1) * matrix(0, i) * determinant3x3(subMatrix);
    }
    return det;
}



// 计算矩阵的代数余子式
mpfr_float algebraicCofactor(Eigen::Matrix< mpfr_float, 4, 4>& matrix, int i, int j) {
    int n = matrix.rows();
    //Eigen::MatrixXd subMatrix(n - 1, n - 1);
    Eigen::Matrix<mpfr_float, 3, 3>subMatrix(3, 3);
    // 创建余子矩阵
    int rowIndex = 0;
    for (int row = 0; row < n; ++row) {
        if (row == i) continue;
        int colIndex = 0;
        for (int col = 0; col < n; ++col) {
            if (col == j) continue;
            subMatrix(rowIndex, colIndex) = matrix(row, col);
            colIndex++;
        }
        rowIndex++;
    }
    // 计算代数余子式
    mpfr_float cofactor = std::pow(-1, i + j) * determinant3x3(subMatrix);
    //mpfr_float cofactor = std::pow(-1, i + j) * subMatrix.determinant();
    //T cofactor = 1.0;
    return cofactor;
}

// 计算矩阵的代数余子式
AD_mpfr_float algebraicCofactor(Eigen::Matrix< AD_mpfr_float, 4, 4>& matrix, int i, int j) {
    int n = matrix.rows();
    //Eigen::MatrixXd subMatrix(n - 1, n - 1);
    Eigen::Matrix<AD_mpfr_float, 3, 3>subMatrix(3, 3);
    // 创建余子矩阵
    int rowIndex = 0;
    for (int row = 0; row < n; ++row) {
        if (row == i) continue;
        int colIndex = 0;
        for (int col = 0; col < n; ++col) {
            if (col == j) continue;
            subMatrix(rowIndex, colIndex) = matrix(row, col);
            colIndex++;
        }
        rowIndex++;
    }
    // 计算代数余子式
    AD_mpfr_float cofactor = std::pow(-1, i + j) * determinant3x3(subMatrix);
    //mpfr_float cofactor = std::pow(-1, i + j) * subMatrix.determinant();
    //T cofactor = 1.0;
    return cofactor;
}


Eigen::SparseMatrix<mpfr_float> Build_stiffness_Matrix(int nv, mpfr_float*& vertices, int n_tet, int* tetrahedras) {
    Eigen::SparseMatrix<mpfr_float> sparse_K(nv * 3, nv * 3);
    Eigen::Matrix<mpfr_float, 6, 6>D;
    Eigen::Matrix<mpfr_float, 6, 12>B; B.setZero();
    Eigen::Matrix<mpfr_float, 12, 12>Ke;
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        mpfr_float p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        Matrix<mpfr_float, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        mpfr_float Volume = determinant4x4(Lambda) / 6.0;
        //mpfr_float Volume = Lambda.determinant() / 6.0;
        for (int i = 0; i < 4; i++) {
            mpfr_float b_i = algebraicCofactor(Lambda, i, 1);
            mpfr_float c_i = algebraicCofactor(Lambda, i, 2);
            mpfr_float d_i = algebraicCofactor(Lambda, i, 3);
            B(0, 3 * i) = b_i; B(1, 3 * i + 1) = c_i; B(2, 3 * i + 2) = d_i;
            B(3, 3 * i) = c_i; B(3, 3 * i + 1) = b_i; B(4, 3 * i + 1) = d_i; B(4, 3 * i + 2) = c_i; B(5, 3 * i) = d_i; B(5, 3 * i + 2) = b_i;
        }
        B /= 6.0 * Volume;
        Ke = B.transpose() * D * B * Volume;
        int index[] = { 3 * t1, 3 * t1 + 1, 3 * t1 + 2,
                        3 * t2, 3 * t2 + 1, 3 * t2 + 2,
                        3 * t3, 3 * t3 + 1, 3 * t3 + 2,
                        3 * t4, 3 * t4 + 1, 3 * t4 + 2 };
        for (int i = 0; i < 12; i++)
            for (int j = 0; j < 12; j++) {
                for (int j = 0; j < 12; j++) {
                    sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                }
            }
    }
    sparse_K.makeCompressed();
    return sparse_K;
}

int Calculate_Area(int nv, int n_tri, mpfr_float* points, int* triangles, mpfr_float*& Area, int& n_sp, int*& Index_of_sp) {
    cout << nv << endl;
    std::unordered_set<int> index_sp;
    mpfr_float* area = (mpfr_float*)calloc(nv, sizeof(mpfr_float));
    for (int i = 0; i < nv; ++i) {
        new (&area[i]) mpfr_float;  // 使用placement new调用mpfr_float的构造函数
    }
    cout << area[0] << endl;
    int* nt_of_points = (int*)calloc(nv, sizeof(int));
    for (int k = 0; k < n_tri; k++) {
        int pi_id = triangles[3 * k];
        int pj_id = triangles[3 * k + 1];
        int pk_id = triangles[3 * k + 2];
        index_sp.insert(pi_id); index_sp.insert(pj_id); index_sp.insert(pk_id);
        Vector3_mp p_i(points[3 * pi_id], points[3 * pi_id + 1], points[3 * pi_id + 2]);
        Vector3_mp p_j(points[3 * pj_id], points[3 * pj_id + 1], points[3 * pj_id + 2]);
        Vector3_mp p_k(points[3 * pk_id], points[3 * pk_id + 1], points[3 * pk_id + 2]);

        Vector3_mp edge1 = p_j - p_i;
        Vector3_mp edge2 = p_k - p_i;
        // 计算叉积
        Vector3_mp crossProduct = edge1.cross(edge2);

        // 计算面积
        mpfr_float area_tri = 0.5 * crossProduct.norm();

        nt_of_points[pi_id] += 1; nt_of_points[pj_id] += 1; nt_of_points[pk_id] += 1;
        area[pi_id] += area_tri; area[pj_id] += area_tri; area[pk_id] += area_tri;
    }
    for (int k = 0; k < nv; k++) {
        if (nt_of_points[k] != 0)
            area[k] /= nt_of_points[k];
    }
    int* index_of_points = (int*)calloc(index_sp.size(), sizeof(int));
    int index = 0;
    Old_2_new.clear();
    for (const auto& element : index_sp) {
        Old_2_new[element] = index;
        index_of_points[index++] = element;
    }
    Area = area;
    n_sp = index_sp.size();
    Index_of_sp = index_of_points;

    free(nt_of_points);
    return 1;
}

SparseMatrix<mpfr_float> Build_N(int nv, int n_tri, mpfr_float* points, int* triangles, int n_sp, int* index_sp, mpfr_float** normals) {
    Eigen::SparseMatrix<mpfr_float> N(nv * 3, n_sp);
    mpfr_float* Normals = (mpfr_float*)calloc(nv * 3, sizeof(mpfr_float));
    for (int i = 0; i < 3 * nv; i++)
        new (&Normals[i]) mpfr_float;
    for (int k = 0; k < n_tri; k++) {
        int pi_id = triangles[3 * k];
        int pj_id = triangles[3 * k + 1];
        int pk_id = triangles[3 * k + 2];
        Vector3_mp p_i(points[3 * pi_id], points[3 * pi_id + 1], points[3 * pi_id + 2]);
        Vector3_mp p_j(points[3 * pj_id], points[3 * pj_id + 1], points[3 * pj_id + 2]);
        Vector3_mp p_k(points[3 * pk_id], points[3 * pk_id + 1], points[3 * pk_id + 2]);
        Vector3_mp edge1 = p_i - p_j;
        Vector3_mp edge2 = p_k - p_j;
        Vector3_mp n_of_tirangle = edge1.cross(edge2);

        Normals[3 * pi_id] += n_of_tirangle[0]; Normals[3 * pi_id + 1] += n_of_tirangle[1]; Normals[3 * pi_id + 2] += n_of_tirangle[2];
        Normals[3 * pj_id] += n_of_tirangle[0]; Normals[3 * pj_id + 1] += n_of_tirangle[1]; Normals[3 * pj_id + 2] += n_of_tirangle[2];
        Normals[3 * pk_id] += n_of_tirangle[0]; Normals[3 * pk_id + 1] += n_of_tirangle[1]; Normals[3 * pk_id + 2] += n_of_tirangle[2];
    }
    for (int k = 0; k < n_sp; k++) {
        int id = index_sp[k];
        Vector3_mp Normal_Of_P(Normals[3 * id], Normals[3 * id + 1], Normals[3 * id + 2]);
        Normal_Of_P.normalize();
        N.insert(3 * id, k) = Normal_Of_P[0];
        N.insert(3 * id + 1, k) = Normal_Of_P[1];
        N.insert(3 * id + 2, k) = Normal_Of_P[2];
        Normals[3 * id] = Normal_Of_P[0];
        Normals[3 * id + 1] = Normal_Of_P[1];
        Normals[3 * id + 2] = Normal_Of_P[2];
    }
    if (normals != nullptr)
        *normals = Normals;
    else
    {
        for (int i = 0; i < 3*nv; ++i) {
            Normals[i].~mpfr_float();  // 显式调用析构函数
        }
        free(Normals);
    }
        
    return N;
}

//SparseMatrix<mpfr_float> Build_N(int nv, int n_tri, mpfr_float* points, int* triangles, int n_sp, int* index_sp, mpfr_float** normals) {
//    Eigen::SparseMatrix<mpfr_float> N(nv * 3, n_sp);
//    mpfr_float* Normals = (mpfr_float*)calloc(nv * 3, sizeof(mpfr_float));
//
//    for (int j = 0; j < n_sp; j++) {
//        int jd = index_sp[j];
//        Vector3_mp n_jd(mpfr_float(0.0), mpfr_float(0.0), mpfr_float(0.0));
//        for (const auto& tri_id : adjacency_tri[jd]) {
//            int t_2 = -1, t_1 = -1, t_3 = -1;
//            for (int k = 0; k < 3; k++) {
//                if (triangles[3 * tri_id + k] == jd) {
//                    t_2 = jd;
//                    t_3 = triangles[3 * tri_id + (k + 1) % 3];
//                    t_1 = triangles[3 * tri_id + (k + 2) % 3];
//                    break;
//                }
//            }
//            Vector3_mp v1(points[3 * t_1] - points[3 * t_2], points[3 * t_1 + 1] - points[3 * t_2 + 1], points[3 * t_1 + 2] - points[3 * t_2 + 2]);
//            Vector3_mp v2(points[3 * t_3] - points[3 * t_2], points[3 * t_3 + 1] - points[3 * t_2 + 1], points[3 * t_3 + 2] - points[3 * t_2 + 2]);
//            Vector3_mp v1_v2_cross = v1.cross(v2);
//            n_jd += v1_v2_cross;
//        }
//        n_jd = n_jd.normalized();
//        //mpfr_float n_norm = n_jd.norm();
//        int k_i = Old_2_new[jd];
//        N.insert(3 * jd, k_i) = n_jd[0];
//        N.insert(3 * jd + 1, k_i) = n_jd[1];
//        N.insert(3 * jd + 2, k_i) = n_jd[2];
//        Normals[3 * jd] = n_jd[0];
//        Normals[3 * jd + 1] = n_jd[1];
//        Normals[3 * jd + 2] = n_jd[2];
//    }
//    if (normals != nullptr)
//        *normals = Normals;
//    else
//    {
//        //for (int i = 0; i < 3 * nv; ++i) {
//        //    Normals[i].~mpfr_float();  // 显式调用析构函数
//        //}
//        free(Normals);
//    }
//
//    return N;
//}


SparseMatrix<mpfr_float> Build_G(int nv, mpfr_float* vertices, SparseMatrix<mpfr_float> N, int n_sp) {
    Eigen::SparseMatrix<mpfr_float> G(6, nv * 3);
    Eigen::SparseMatrix<mpfr_float> G_N(6, n_sp);
    for (int k = 0; k < nv; k++) {
        G.insert(0, 3 * k) = 1;
        G.insert(1, 3 * k + 1) = 1;
        G.insert(2, 3 * k + 2) = 1;

        G.insert(3, 3 * k + 1) = -(vertices[3 * k + 2]); G.insert(3, 3 * k + 2) = (vertices[3 * k + 1]);
        G.insert(4, 3 * k) = (vertices[3 * k + 2]); G.insert(4, 3 * k + 2) = -(vertices[3 * k]);
        G.insert(5, 3 * k) = -(vertices[3 * k + 1]); G.insert(5, 3 * k + 1) = (vertices[3 * k]);
    }
    G_N = G * N;
    return G_N;
}

SparseMatrix<mpfr_float> Build_M_G(int nv, mpfr_float* vertices) {
    Eigen::SparseMatrix<mpfr_float> G(6, nv * 3);
    for (int k = 0; k < nv; k++) {
        G.insert(0, 3 * k) = 1;
        G.insert(1, 3 * k + 1) = 1;
        G.insert(2, 3 * k + 2) = 1;

        G.insert(3, 3 * k + 1) = -(vertices[3 * k + 2]); G.insert(3, 3 * k + 2) = (vertices[3 * k + 1]);
        G.insert(4, 3 * k) = (vertices[3 * k + 2]); G.insert(4, 3 * k + 2) = -(vertices[3 * k]);
        G.insert(5, 3 * k) = -(vertices[3 * k + 1]); G.insert(5, 3 * k + 1) = (vertices[3 * k]);
    }
    return G;
}

//计算最大特征值
mpfr_float calculate_max_eigenvalue(SparseMatrix<mpfr_float> K) {
    SparseSymMatProd<mpfr_float> opK(K);
    SymEigsSolver<SparseSymMatProd<mpfr_float>> eigs_K(opK, 1, 6);
    eigs_K.init();
    eigs_K.compute(SortRule::LargestAlge);
    if (eigs_K.info() != CompInfo::Successful) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The Eigenvalues of K were not found. !");
    }
    //Eigen::MatrixXd evecs;
    //Eigen::VectorXd evalues;
    MatrixX_mp evecs;
    VectorX_mp evalues;
    if (eigs_K.info() == CompInfo::Successful)
    {
        evalues = eigs_K.eigenvalues();
        evecs = eigs_K.eigenvectors();
        std::cout << "Eigenvalues of K is found:" << evalues << std::endl;
        //std::cout << "Eigenvalues found:\n" << evecs << std::endl;
    }
    else {
        cout << "The Eigenvalues of K were not found.Error:NotConvergence" << endl;
        return -1;
    }
    return evalues(0);
}

Eigen::SparseMatrix<mpfr_float> Calculate_Stresses(int nv, mpfr_float* vertices, int n_sp, int* index_sp, int n_tri, int* triangles, int n_tet, int* tetrahedras,
    VectorX_mp p, Eigen::SparseLU<Eigen::SparseMatrix<mpfr_float>>& solver, int& intersect) {
    int n = 3 * nv;
    Eigen::SparseMatrix<mpfr_float> N = Build_N(nv, n_tri, vertices, triangles, n_sp, index_sp, nullptr);
    Eigen::SparseMatrix<mpfr_float> G = Build_G(nv, vertices, N, n_sp);
    MatrixX_mp GGT = G * G.transpose();
    MatrixX_mp GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<mpfr_float> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<mpfr_float> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    //p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    VectorX_mp u(n + 6);
    VectorX_mp f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver.solve(f);
    // 检查解是否成功
    std::cout << "The norm of u is : " << u.norm() << std::endl;

    Eigen::SparseMatrix<mpfr_float> sparse_Stress(n_tet * 6, 1);
    Eigen::Matrix<mpfr_float, 6, 1>S_tress_e;
    Eigen::Matrix<mpfr_float, 12, 1>ue;
    Eigen::Matrix<mpfr_float, 6, 6>D;
    Eigen::Matrix<mpfr_float, 6, 12>B; B.setZero();
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        mpfr_float p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        Matrix<mpfr_float, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        //mpfr_float Volume = Lambda.determinant() / 6.0;
        mpfr_float Volume = determinant4x4(Lambda) / 6.0;
        if (Volume < 0.0) {
            intersect = 1;
            break;
        }
        for (int i = 0; i < 4; i++) {
            mpfr_float b_i = algebraicCofactor(Lambda, i, 1);
            mpfr_float c_i = algebraicCofactor(Lambda, i, 2);
            mpfr_float d_i = algebraicCofactor(Lambda, i, 3);
            B(0, 3 * i) = b_i; B(1, 3 * i + 1) = c_i; B(2, 3 * i + 2) = d_i;
            B(3, 3 * i) = c_i; B(3, 3 * i + 1) = b_i; B(4, 3 * i + 1) = d_i; B(4, 3 * i + 2) = c_i; B(5, 3 * i) = d_i; B(5, 3 * i + 2) = b_i;
        }
        B /= 6.0 * Volume;
        ue << u[t1 * 3], u[t1 * 3 + 1], u[t1 * 3 + 2],
            u[t2 * 3], u[t2 * 3 + 1], u[t2 * 3 + 2],
            u[t3 * 3], u[t3 * 3 + 1], u[t3 * 3 + 2],
            u[t4 * 3], u[t4 * 3 + 1], u[t4 * 3 + 2];
        S_tress_e = D * B * ue;
        for (int stress_i = 0; stress_i < 6; stress_i++) {
            sparse_Stress.coeffRef(6 * k + stress_i, 0) += S_tress_e(stress_i, 0);
        }
    }
    sparse_Stress.makeCompressed();
    return sparse_Stress;
}

Eigen::SparseMatrix<mpfr_float> Calculate_Stresses_for_test(int nv, mpfr_float* points, mpfr_float* vertices, int n_sp, int* index_sp, int n_tri, int* triangles, int n_tet, int* tetrahedras,
    VectorX_mp p, Eigen::SparseLU<Eigen::SparseMatrix<mpfr_float>>& solver, int& intersect) {
    int n = 3 * nv;
    Eigen::SparseMatrix<mpfr_float> N = Build_N(nv, n_tri, vertices, triangles, n_sp, index_sp, nullptr);
    Eigen::SparseMatrix<mpfr_float> G = Build_G(nv, vertices, N, n_sp);
    MatrixX_mp GGT = G * G.transpose();
    MatrixX_mp GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<mpfr_float> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<mpfr_float> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    //p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    VectorX_mp u(n + 6);
    VectorX_mp f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver.solve(f);
    // 检查解是否成功
    std::cout << "The norm of u is : " << u.norm() << std::endl;

    Eigen::SparseMatrix<mpfr_float> sparse_Stress(n_tet * 6, 1);
    Eigen::Matrix<mpfr_float, 6, 1>S_tress_e;
    Eigen::Matrix<mpfr_float, 12, 1>ue;
    Eigen::Matrix<mpfr_float, 6, 6>D;
    Eigen::Matrix<mpfr_float, 6, 12>B; B.setZero();
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        mpfr_float p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        Matrix<mpfr_float, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        //mpfr_float Volume = Lambda.determinant() / 6.0;
        mpfr_float Volume = determinant4x4(Lambda) / 6.0;
        if (Volume < 0.0) {
            intersect = 1;
            break;
        }
        for (int i = 0; i < 4; i++) {
            mpfr_float b_i = algebraicCofactor(Lambda, i, 1);
            mpfr_float c_i = algebraicCofactor(Lambda, i, 2);
            mpfr_float d_i = algebraicCofactor(Lambda, i, 3);
            B(0, 3 * i) = b_i; B(1, 3 * i + 1) = c_i; B(2, 3 * i + 2) = d_i;
            B(3, 3 * i) = c_i; B(3, 3 * i + 1) = b_i; B(4, 3 * i + 1) = d_i; B(4, 3 * i + 2) = c_i; B(5, 3 * i) = d_i; B(5, 3 * i + 2) = b_i;
        }
        B /= 6.0 * Volume;
        ue << u[t1 * 3], u[t1 * 3 + 1], u[t1 * 3 + 2],
            u[t2 * 3], u[t2 * 3 + 1], u[t2 * 3 + 2],
            u[t3 * 3], u[t3 * 3 + 1], u[t3 * 3 + 2],
            u[t4 * 3], u[t4 * 3 + 1], u[t4 * 3 + 2];
        S_tress_e = D * B * ue;
        for (int stress_i = 0; stress_i < 6; stress_i++) {
            sparse_Stress.coeffRef(6 * k + stress_i, 0) += S_tress_e(stress_i, 0);
        }
    }
    sparse_Stress.makeCompressed();
    return sparse_Stress;
}


bool isPositiveSemiDefinite(const Eigen::MatrixXd& A) {
    // 使用Eigen的自洽特征值计算方法
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(A);
    if (eigenSolver.info() != Eigen::Success) abort();

    // 如果所有特征值都是非负的，则矩阵是半正定的
    return eigenSolver.eigenvalues().minCoeff() >= 0;
}
bool isPositiveSemiDefinite(const MatrixX_mp& A) {
    // 使用Eigen的自洽特征值计算方法
    Eigen::SelfAdjointEigenSolver< MatrixX_mp> eigenSolver(A);
    if (eigenSolver.info() != Eigen::Success) abort();

    // 如果所有特征值都是非负的，则矩阵是半正定的
    return eigenSolver.eigenvalues().minCoeff() >= 0;
}
//gradient
//CppAD::AD<double> algebraicCofactor(Eigen::Matrix< CppAD::AD<double>, 4, 4>& matrix, int i, int j) {
//    int n = matrix.rows();
//    MatrixXd_AD subMatrix(n - 1, n - 1);
//
//    // 创建余子矩阵
//    int rowIndex = 0;
//    for (int row = 0; row < n; ++row) {
//        if (row == i) continue;
//        int colIndex = 0;
//        for (int col = 0; col < n; ++col) {
//            if (col == j) continue;
//            subMatrix(rowIndex, colIndex) = matrix(row, col);
//            colIndex++;
//        }
//        rowIndex++;
//    }
//    // 计算代数余子式
//    CppAD::AD<double> cofactor = std::pow(-1, i + j) * subMatrix.determinant();
//    return cofactor;
//}

CppAD::AD<double> Calculate_O(CppAD::AD<double> t, Eigen::SparseMatrix<mpfr_float> stress, int n_tet) {
    CppAD::AD<double> O = 0.0;
    for (int i = 0; i < n_tet; i++) {
        Eigen::Matrix<CppAD::AD<double>, 3, 3> sigma;
        Eigen::Matrix<CppAD::AD<double>, 3, 3> I;
        CppAD::AD<double> det1, det2;
        sigma << stress.coeff(6 * i, 0).convert_to<double>(), stress.coeff(6 * i + 3, 0).convert_to<double>(), stress.coeff(6 * i + 5, 0).convert_to<double>(),
            stress.coeff(6 * i + 3, 0).convert_to<double>(), stress.coeff(6 * i + 1, 0).convert_to<double>(), stress.coeff(6 * i + 4, 0).convert_to<double>(),
            stress.coeff(6 * i + 5, 0).convert_to<double>(), stress.coeff(6 * i + 4, 0).convert_to<double>(), stress.coeff(6 * i + 2, 0).convert_to<double>();
        I << t * 1.0, 0.0, 0.0,
            0.0, t * 1.0, 0.0,
            0.0, 0.0, t * 1.0;
        det1 = (I - sigma).determinant();
        det2 = (I + sigma).determinant();
        if (det1 < 1e-3 || det2 < 1e-3)
        { 
            cout << "Invalid Guess s,t" << endl;
            return -1e6;
        }
        Eigen::Matrix<double, 3, 3> Sigma, tI;
        Sigma << stress.coeff(6 * i, 0).convert_to<double>(), stress.coeff(6 * i + 3, 0).convert_to<double>(), stress.coeff(6 * i + 5, 0).convert_to<double>(),
            stress.coeff(6 * i + 3, 0).convert_to<double>(), stress.coeff(6 * i + 1, 0).convert_to<double>(), stress.coeff(6 * i + 4, 0).convert_to<double>(),
            stress.coeff(6 * i + 5, 0).convert_to<double>(), stress.coeff(6 * i + 4, 0).convert_to<double>(), stress.coeff(6 * i + 2, 0).convert_to<double>();
        tI << Value(t) * 1.0, 0.0, 0.0,
            0.0, Value(t) * 1.0, 0.0,
            0.0, 0.0, Value(t) * 1.0;
        Eigen::MatrixXd mat1 = tI - Sigma;
        Eigen::MatrixXd mat2 = tI + Sigma;
        if (!isPositiveSemiDefinite(mat1) || !isPositiveSemiDefinite(mat2)) {
            cout << "NOT Positive SemiDefinite." << endl;
            return -1e6;
        }
        O += -O_mu * (CppAD::log(det1) + CppAD::log(det2));
    }
    O += t;
    return O;
}


mpfr_float Calculate_O(mpfr_float t, Eigen::SparseMatrix<mpfr_float> stress, int n_tet) {
    mpfr_float O = 0.0;
    for (int i = 0; i < n_tet; i++) {
        Eigen::Matrix<mpfr_float, 3, 3> sigma;
        Eigen::Matrix<mpfr_float, 3, 3> I;
        mpfr_float det1, det2;
        sigma << stress.coeff(6 * i, 0), stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 5, 0),
            stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 1, 0), stress.coeff(6 * i + 4, 0),
            stress.coeff(6 * i + 5, 0), stress.coeff(6 * i + 4, 0), stress.coeff(6 * i + 2, 0);
        I << t , mpfr_float(0.0), mpfr_float(0.0),
            mpfr_float(0.0), t, mpfr_float(0.0),
            mpfr_float(0.0), mpfr_float(0.0), t ;
        det1 = (I - sigma).determinant();
        det2 = (I + sigma).determinant();
        if (det1 < mpfr_float(1e-3)|| det2 < mpfr_float(1e-3))
        {
            cout << "Invalid Guess s,t" << endl;
            return mpfr_float(-1e6);
        }
        Eigen::Matrix<mpfr_float, 3, 3> Sigma, tI;
        Sigma << stress.coeff(6 * i, 0), stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 5, 0),
            stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 1, 0), stress.coeff(6 * i + 4, 0),
            stress.coeff(6 * i + 5, 0), stress.coeff(6 * i + 4, 0), stress.coeff(6 * i + 2, 0);
        tI << t, mpfr_float(0.0), mpfr_float(0.0),
            mpfr_float(0.0), t, mpfr_float(0.0),
            mpfr_float(0.0), mpfr_float(0.0), t;
        MatrixX_mp mat1 = tI - Sigma;
        MatrixX_mp mat2 = tI + Sigma;
        if (!isPositiveSemiDefinite(mat1) || !isPositiveSemiDefinite(mat2)) {
            cout << "NOT Positive SemiDefinite." << endl;
            return mpfr_float(-1e6);
        }
        O += -mpfr_float(O_mu) * (log(det1) + log(det2));
    }
    O += t;
    return O;
}

Eigen::SparseMatrix<AD_mpfr_float> Build_stiffness_Matrix(int nv, const vector<AD_mpfr_float>& vertices, int n_tet, int* tetrahedras) {
    Eigen::SparseMatrix<AD_mpfr_float> sparse_K(nv * 3, nv * 3);
    Eigen::Matrix<AD_mpfr_float, 6, 6>D;
    Eigen::Matrix<AD_mpfr_float, 6, 12>B; B.setZero();
    Eigen::Matrix<AD_mpfr_float, 12, 12>Ke;
    D << mpfr_float(1 - mu), mpfr_float(mu), mpfr_float(mu), mpfr_float(0), mpfr_float(0), mpfr_float(0),
        mpfr_float(mu), mpfr_float(1 - mu), mpfr_float(mu), mpfr_float(0), mpfr_float(0), mpfr_float(0),
        mpfr_float(mu), mpfr_float(mu), mpfr_float(1 - mu), mpfr_float(0),  mpfr_float(0), mpfr_float(0),
        mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float((1 - 2 * mu) / 2), mpfr_float(0), mpfr_float(0),
        mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float((1 - 2 * mu) / 2), mpfr_float(0),
        mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float((1 - 2 * mu) / 2);
    D *= mpfr_float(E) / mpfr_float((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        AD_mpfr_float p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        Matrix<AD_mpfr_float, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        Lambda << mpfr_float(1), p1_x, p1_y, p1_z,
            mpfr_float(1), p2_x, p2_y, p2_z,
            mpfr_float(1), p3_x, p3_y, p3_z,
            mpfr_float(1), p4_x, p4_y, p4_z;
        AD_mpfr_float Volume = determinant4x4(Lambda) / 6.0;
        for (int i = 0; i < 4; i++) {
            AD_mpfr_float b_i = algebraicCofactor(Lambda, i, 1);
            AD_mpfr_float c_i = algebraicCofactor(Lambda, i, 2);
            AD_mpfr_float d_i = algebraicCofactor(Lambda, i, 3);
            B(0, 3 * i) = b_i; B(1, 3 * i + 1) = c_i; B(2, 3 * i + 2) = d_i;
            B(3, 3 * i) = c_i; B(3, 3 * i + 1) = b_i; B(4, 3 * i + 1) = d_i; B(4, 3 * i + 2) = c_i; B(5, 3 * i) = d_i; B(5, 3 * i + 2) = b_i;
        }
        B /= 6.0 * Volume;
        Ke = B.transpose() * D * B * Volume;
        int index[] = { 3 * t1, 3 * t1 + 1, 3 * t1 + 2,
                        3 * t2, 3 * t2 + 1, 3 * t2 + 2,
                        3 * t3, 3 * t3 + 1, 3 * t3 + 2,
                        3 * t4, 3 * t4 + 1, 3 * t4 + 2 };
        for (int i = 0; i < 12; i++)
            for (int j = 0; j < 12; j++) {
                for (int j = 0; j < 12; j++) {
                    sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                }
            }
    }
    sparse_K.makeCompressed();
    return sparse_K;
}

SparseMatrix<AD_mpfr_float> Build_M_G(int nv, const vector<AD_mpfr_float>& vertices) {
    Eigen::SparseMatrix<AD_mpfr_float> G(6, nv * 3);
    for (int k = 0; k < nv; k++) {
        G.insert(0, 3 * k) = 1;
        G.insert(1, 3 * k + 1) = 1;
        G.insert(2, 3 * k + 2) = 1;

        G.insert(3, 3 * k + 1) = -(vertices[3 * k + 2]); G.insert(3, 3 * k + 2) = (vertices[3 * k + 1]);
        G.insert(4, 3 * k) = (vertices[3 * k + 2]); G.insert(4, 3 * k + 2) = -(vertices[3 * k]);
        G.insert(5, 3 * k) = -(vertices[3 * k + 1]); G.insert(5, 3 * k + 1) = (vertices[3 * k]);
    }
    return G;
}

SparseMatrix<mpfr_float> calculate_grad_N(int p_id,int x_y_z,int nv, int n_tri, mpfr_float* points, int* triangles, int n_sp) {
    Eigen::SparseMatrix<mpfr_float> N(nv * 3, n_sp);
    for (const auto& jd : adjacency_p[p_id]) {
        mpfr_float v1x, v1y, v1z, v2x, v2y, v2z;
        Vector3_mp n_jd(mpfr_float(0.0), mpfr_float(0.0), mpfr_float(0.0));
        for (const auto& tri_id : adjacency_tri[jd]) {
            int t_2 = -1, t_1 = -1, t_3 = -1;

            for (int k = 0; k < 3; k++) {
                if (triangles[3*tri_id+k]== jd) {
                    t_2 = jd;
                    t_3 = triangles[3 * tri_id + (k + 1) % 3];
                    t_1 = triangles[3 * tri_id + (k + 2) % 3];
                    break;
                }
            }

            if (t_1 == p_id) {
                v1x = points[3 * t_3] - points[3 * t_2];
                v1y = points[3 * t_3 + 1] - points[3 * t_2 + 1];
                v1z = points[3 * t_3 + 2] - points[3 * t_2 + 2];
            }
            if (t_3 == p_id) {
                v2x = points[3 * t_1] - points[3 * t_2];
                v2y = points[3 * t_1 + 1] - points[3 * t_2 + 1];
                v2z = points[3 * t_1 + 2] - points[3 * t_2 + 2];
            }
            Vector3_mp v1(points[3 * t_1] - points[3 * t_2], points[3 * t_1 + 1] - points[3 * t_2 + 1], points[3 * t_1 + 2] - points[3 * t_2 + 2]);
            Vector3_mp v2(points[3 * t_3] - points[3 * t_2], points[3 * t_3 + 1] - points[3 * t_2 + 1], points[3 * t_3 + 2] - points[3 * t_2 + 2]);
            Vector3_mp v1_v2_cross = v1.cross(v2);
            n_jd += v1_v2_cross;
        }
        //n_jd = n_jd.normalized();
        mpfr_float n_norm = n_jd.norm();
        mpfr_float n_norm_2 = n_norm* n_norm;
        int k_i = Old_2_new[jd];
        if (x_y_z == 0) {
            mpfr_float d_norm_dx = -((v2z - v1z) * n_jd[1] + (v1y - v2y) * n_jd[2]) / (n_norm * n_norm_2);
            //cout << d_norm_dx << endl;
            mpfr_float dnx_dx = n_jd[0] * d_norm_dx;
            mpfr_float dny_dx = (v2z - v1z) / n_norm + n_jd[1] * d_norm_dx;
            mpfr_float dnz_dx = (v1y - v2y) / n_norm + n_jd[2] * d_norm_dx;

            N.insert(3 * jd, k_i) = dnx_dx;
            N.insert(3 * jd + 1, k_i) = dny_dx;
            N.insert(3 * jd + 2, k_i) = dnz_dx;
        }
        else if (x_y_z == 1) {
            mpfr_float d_norm_dy = -((v1z - v2z) * n_jd[0] + (v2x - v1x) * n_jd[2]) / (n_norm * n_norm_2);
            mpfr_float dnx_dy = (v1z - v2z) / n_norm + n_jd[0] * d_norm_dy;
            mpfr_float dny_dy = n_jd[1] * d_norm_dy;
            mpfr_float dnz_dy = (v2x - v1x) / n_norm + n_jd[2] * d_norm_dy;
            N.insert(3 * jd, k_i) = dnx_dy;
            N.insert(3 * jd + 1, k_i) = dny_dy;
            N.insert(3 * jd + 2, k_i) = dnz_dy;
        }
        else if (x_y_z == 2) {
            mpfr_float d_norm_dz = -((v2y - v1y) * n_jd[0] + (v1x - v2x) * n_jd[1]) / (n_norm * n_norm_2);
            mpfr_float dnx_dz = (v2y - v1y) / n_norm + n_jd[0] * d_norm_dz;
            mpfr_float dny_dz = (v1x - v2x) / n_norm + n_jd[1] * d_norm_dz;
            mpfr_float dnz_dz = n_jd[2] * d_norm_dz;
            N.insert(3 * jd, k_i) = dnx_dz;
            N.insert(3 * jd + 1, k_i) = dny_dz;
            N.insert(3 * jd + 2, k_i) = dnz_dz;
        }
    }
    return N;
}


VectorX_mp Calculate_sigma_A(int nv, mpfr_float* points, int n_tet, int* tetrahedras, VectorX_mp u, mpfr_float t) {
    int n = 3 * nv;
    VectorX_mp stress(n + 6); stress.setZero();
    VectorX_mp ue(12);
    VectorX_mp se(6);
    VectorX_mp Stress_e(6);
    VectorX_mp D_e(6);
    Eigen::Matrix<mpfr_float, 6, 6>D;
    Eigen::Matrix<mpfr_float, 6, 12>B; B.setZero();
    Eigen::Matrix<mpfr_float, 6, 12>sigma;
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        mpfr_float p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        Matrix<mpfr_float, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = points[t1 * 3]; p1_y = points[t1 * 3 + 1]; p1_z = points[t1 * 3 + 2];
        p2_x = points[t2 * 3]; p2_y = points[t2 * 3 + 1]; p2_z = points[t2 * 3 + 2];
        p3_x = points[t3 * 3]; p3_y = points[t3 * 3 + 1]; p3_z = points[t3 * 3 + 2];
        p4_x = points[t4 * 3]; p4_y = points[t4 * 3 + 1]; p4_z = points[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        mpfr_float Volume = determinant4x4(Lambda) / 6.0;
        for (int i = 0; i < 4; i++) {
            mpfr_float b_i = algebraicCofactor(Lambda, i, 1);
            mpfr_float c_i = algebraicCofactor(Lambda, i, 2);
            mpfr_float d_i = algebraicCofactor(Lambda, i, 3);
            B(0, 3 * i) = b_i; B(1, 3 * i + 1) = c_i; B(2, 3 * i + 2) = d_i;
            B(3, 3 * i) = c_i; B(3, 3 * i + 1) = b_i; B(4, 3 * i + 1) = d_i; B(4, 3 * i + 2) = c_i; B(5, 3 * i) = d_i; B(5, 3 * i + 2) = b_i;
        }
        B /= 6.0 * Volume;
        sigma = D * B;
        int index[] = { 3 * t1, 3 * t1 + 1, 3 * t1 + 2,
                        3 * t2, 3 * t2 + 1, 3 * t2 + 2,
                        3 * t3, 3 * t3 + 1, 3 * t3 + 2,
                        3 * t4, 3 * t4 + 1, 3 * t4 + 2 };
        for (int u_i = 0; u_i < 12; u_i++)
            ue[u_i] = u[index[u_i]];
        Stress_e = sigma * ue;
        MatrixX_mp tI_stress_1(3, 3), tI_stress_2(3, 3);
        MatrixX_mp O(3, 3);
        tI_stress_1 << t * 1.0 - Stress_e[0], -Stress_e[3], -Stress_e[5],
            -Stress_e[3], t * 1.0 - Stress_e[1], -Stress_e[4],
            -Stress_e[5], -Stress_e[4], t * 1.0 - Stress_e[2];
        tI_stress_2 << t * 1.0 + Stress_e[0], Stress_e[3], Stress_e[5],
            Stress_e[3], t * 1.0 + Stress_e[1], Stress_e[4],
            Stress_e[5], Stress_e[4], t * 1.0 + Stress_e[2];

        tI_stress_1 = tI_stress_1.inverse();
        tI_stress_2 = tI_stress_2.inverse();

        O = tI_stress_1 - tI_stress_2;
        D_e << O.coeff(0, 0), O.coeff(1, 1), O.coeff(2, 2), O.coeff(0, 1) + O.coeff(1, 0), O.coeff(1, 2) + O.coeff(2, 1), O.coeff(0, 2) + O.coeff(2, 0);
        se = D_e.transpose() * sigma;
        for (int s_i = 0; s_i < 12; s_i++) {
            stress[index[s_i]] += se[s_i];
        }
    }
    return stress;
}


double Calculate_Stresses_AD(int nv, mpfr_float* points, int n_tri, int* triangles, int n_tet, int* tetrahedras, int n_sp, int* index_sp, mpfr_float t,
    VectorX_mp p, VectorX_mp& grad_s, Eigen::SparseLU<Eigen::SparseMatrix<mpfr_float>>& solver) {
    int n = 3 * nv;
    auto start_time_AD = std::chrono::high_resolution_clock::now();
    //double* Normal;

    //1-计算位移 u 
    Eigen::SparseMatrix<mpfr_float> K = Build_stiffness_Matrix(nv, points, n_tet, tetrahedras);
    Eigen::SparseMatrix<mpfr_float> N = Build_N(nv, n_tri, points, triangles, n_sp, index_sp, nullptr);
    Eigen::SparseMatrix<mpfr_float> G = Build_G(nv, points, N, n_sp);
    Eigen::SparseMatrix<mpfr_float> grad_M(n + 6, n + 6); grad_M.setZero();
    VectorX_mp grad_O_s(n); grad_O_s.setZero();

    MatrixX_mp GGT = G * G.transpose();
    MatrixX_mp GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<mpfr_float> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<mpfr_float> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    //p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    VectorX_mp u(n + 6);
    VectorX_mp f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver.solve(f);
    // 检查解是否成功
    std::cout << "The norm of u is : " << u.norm() << std::endl;
    //-计算位移 u over -1

    //2-计算gradient
    VectorX_mp ue(12);
    VectorX_mp Stress_e;
    VectorX_mp grad_f_s(n + 6); grad_f_s.setZero();
    Eigen::Matrix<mpfr_float, 6, 12>sigma; sigma.setZero();
    vector<mpfr_float> s_x(12);
    vector<AD_mpfr_float> s(12);

    Eigen::Matrix<AD_mpfr_float, 6, 6>D_ad;
    Eigen::Matrix<AD_mpfr_float, 6, 12>B_ad; B_ad.setZero();
    Eigen::Matrix<AD_mpfr_float, 6, 12>sigma_s;
    //Eigen::SparseMatrix<CppAD::AD<double>> N_AD;
    Eigen::SparseMatrix<mpfr_float> N_AD;
    Eigen::SparseMatrix<AD_mpfr_float> G_AD;
    std::vector<AD_mpfr_float> Sigma_AD_Vector;
    std::vector<mpfr_float> jac_Sigma_s;
    std::vector<mpfr_float> jac_K;
    std::vector<mpfr_float> jac_G;
    std::vector<size_t> row_indices;
    std::vector<size_t> col_indices;
    std::vector<size_t> G_row_indices;
    std::vector<size_t> G_col_indices;
    std::vector<AD_mpfr_float> K_AD_Vector;
    std::vector<AD_mpfr_float> G_AD_Vector;
    Eigen::Matrix<mpfr_float, Eigen::Dynamic, 1> pressure(n_sp);

    for (int k = 0; k < p.size(); k++) {
        pressure[k] = p[k];
    }

    auto start_time_Np = std::chrono::high_resolution_clock::now();

    //2-1:calculate O_B
    VectorX_mp O_B = solver.solve(f);
    O_B = u;

    //2-3:calculate O_A
    VectorX_mp sigma_A = Calculate_sigma_A(nv, points, n_tet, tetrahedras, u, t);
    VectorX_mp O_A = solver.solve(sigma_A);

    //calculate Grad_K
    vector<AD_mpfr_float> vertices(n);
    vector<mpfr_float> vertices_x(n);
    for (int i = 0; i < n; i++) {
        vertices[i] = AD_mpfr_float(points[i]);
        vertices_x[i] = points[i];
    }
    cout << "number of triangles: " << n_tri << endl;
    cout << "The number of vertices: " << nv << endl;
    auto START_TIME2 = std::chrono::high_resolution_clock::now();
    CppAD::Independent(vertices);
    Eigen::SparseMatrix<AD_mpfr_float> K_AD = Build_stiffness_Matrix(nv, vertices, n_tet, tetrahedras);
    // 遍历 K 的非零元素
    row_indices.clear(); col_indices.clear(); K_AD_Vector.clear();
    for (int k = 0; k < K_AD.outerSize(); ++k) {
        for (Eigen::SparseMatrix<AD_mpfr_float>::InnerIterator it(K_AD, k); it; ++it) {
            row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
            col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
            K_AD_Vector.push_back(AD_mpfr_float(it.value()));// it.value()     // 非零元素的值
        }
    }
    CppAD::ADFun<mpfr_float> func(vertices, K_AD_Vector);    // 创建 ADFun 对象
    jac_K = func.Jacobian(vertices_x);
    auto end_time9 = std::chrono::high_resolution_clock::now();
    auto duration_Np9 = std::chrono::duration_cast<std::chrono::microseconds>(end_time9 - START_TIME2).count() / 1e6;
    std::cout << "||------ The cost of calculating grad_K : " << duration_Np9 << " seconds ------||" << endl << endl;

    //calculate grad_G
    CppAD::Independent(vertices);
    G_AD = Build_M_G(nv, vertices);
    G_row_indices.clear(); G_col_indices.clear(); G_AD_Vector.clear();
    for (int k = 0; k < G_AD.outerSize(); ++k) {
        for (Eigen::SparseMatrix<AD_mpfr_float>::InnerIterator it(G_AD, k); it; ++it) {
            G_row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
            G_col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
            G_AD_Vector.push_back(AD_mpfr_float(it.value()));// it.value()     // 非零元素的值
        }
    }
    CppAD::ADFun<mpfr_float> G_func(vertices, G_AD_Vector);    // 创建 ADFun 对象
    jac_G = G_func.Jacobian(vertices_x);
    auto end_time12 = std::chrono::high_resolution_clock::now();
    auto duration_Np12 = std::chrono::duration_cast<std::chrono::microseconds>(end_time12 - end_time9).count() / 1e6;
    std::cout << "||------ The cost of calculating the grad_G : " << duration_Np12 << " seconds ------||" << endl << endl;

    //calculate grad_sigma(s)
    D_ad << mpfr_float(1 - mu), mpfr_float(mu), mpfr_float(mu), mpfr_float(0), mpfr_float(0), mpfr_float(0),
        mpfr_float(mu), mpfr_float(1 - mu), mpfr_float(mu), mpfr_float(0), mpfr_float(0), mpfr_float(0),
        mpfr_float(mu), mpfr_float(mu), mpfr_float(1 - mu), mpfr_float(0), mpfr_float(0), mpfr_float(0),
        mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float((1 - 2 * mu) / 2), mpfr_float(0), mpfr_float(0),
        mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float((1 - 2 * mu) / 2), mpfr_float(0),
        mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float(0), mpfr_float((1 - 2 * mu) / 2);
    D_ad *= mpfr_float(E) / mpfr_float((1 + mu) * (1 - 2 * mu));
    Matrix<AD_mpfr_float, 4, 4> Lambda;
    for (int k = 0; k < n_tet; k++) {
        int t1, t2, t3, t4;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        int index[] = { 3 * t1, 3 * t1 + 1, 3 * t1 + 2,
                    3 * t2, 3 * t2 + 1, 3 * t2 + 2,
                    3 * t3, 3 * t3 + 1, 3 * t3 + 2,
                    3 * t4, 3 * t4 + 1, 3 * t4 + 2 };
        for (int index_i = 0; index_i < 12; index_i++) {
            s[index_i] = AD_mpfr_float(points[index[index_i]]);
            s_x[index_i] = points[index[index_i]];
            ue[index_i] = u[index[index_i]];
        }
        CppAD::Independent(s);
        Lambda << mpfr_float(1), s[0], s[1], s[2],
            mpfr_float(1), s[3], s[4], s[5],
            mpfr_float(1), s[6], s[7], s[8],
            mpfr_float(1), s[9], s[10], s[11];
        AD_mpfr_float Volume = determinant4x4(Lambda) / mpfr_float(6.0);
        mpfr_float Volume_value = Value(Volume);
        for (int i = 0; i < 4; i++) {
            AD_mpfr_float b_i = algebraicCofactor(Lambda, i, 1);
            AD_mpfr_float c_i = algebraicCofactor(Lambda, i, 2);
            AD_mpfr_float d_i = algebraicCofactor(Lambda, i, 3);
            B_ad(0, 3 * i) = b_i; B_ad(1, 3 * i + 1) = c_i; B_ad(2, 3 * i + 2) = d_i;
            B_ad(3, 3 * i) = c_i; B_ad(3, 3 * i + 1) = b_i; B_ad(4, 3 * i + 1) = d_i; B_ad(4, 3 * i + 2) = c_i; B_ad(5, 3 * i) = d_i; B_ad(5, 3 * i + 2) = b_i;
        }
        B_ad /= 6.0 * Volume;
        sigma_s = D_ad * B_ad;
        Sigma_AD_Vector.clear();
        for (int sigma_i = 0; sigma_i < 6; sigma_i++)
            for (int sigma_j = 0; sigma_j < 12; sigma_j++) {
                Sigma_AD_Vector.push_back(sigma_s.coeffRef(sigma_i, sigma_j));
                sigma.coeffRef(sigma_i, sigma_j) = Value(sigma_s.coeffRef(sigma_i, sigma_j));
            }
        CppAD::ADFun<mpfr_float> sigma_fun(s, Sigma_AD_Vector);
        jac_Sigma_s = sigma_fun.Jacobian(s_x);

        Stress_e = sigma * ue;
        MatrixX_mp tI_stress_1(3, 3), tI_stress_2(3, 3);
        MatrixX_mp O(3, 3);
        tI_stress_1 << t * 1.0 - Stress_e[0], -Stress_e[3], -Stress_e[5],
            -Stress_e[3], t * 1.0 - Stress_e[1], -Stress_e[4],
            -Stress_e[5], -Stress_e[4], t * 1.0 - Stress_e[2];
        tI_stress_2 << t * 1.0 + Stress_e[0], Stress_e[3], Stress_e[5],
            Stress_e[3], t * 1.0 + Stress_e[1], Stress_e[4],
            Stress_e[5], Stress_e[4], t * 1.0 + Stress_e[2];
        tI_stress_1 = tI_stress_1.inverse();
        tI_stress_2 = tI_stress_2.inverse();

        VectorX_mp M_Np_select(12);
        for (int t_i = 0; t_i < 12; t_i++)
            M_Np_select[t_i] = O_B[index[t_i]];
        for (int k_i = 0; k_i < 12; k_i++) {
            MatrixX_mp grad_sigma_s_k(6, 12);
            for (int sigma_i = 0; sigma_i < 6; sigma_i++)
                for (int sigma_j = 0; sigma_j < 12; sigma_j++)
                    grad_sigma_s_k.coeffRef(sigma_i, sigma_j) = jac_Sigma_s[(sigma_i * 12 + sigma_j) * 12 + k_i];
            VectorX_mp part_1 = grad_sigma_s_k * M_Np_select;
            MatrixX_mp sigma_matrix(3, 3);
            MatrixX_mp D_part_1(3, 3), O_1(3, 3), O_2(3, 3);;
            sigma_matrix << part_1[0], part_1[3], part_1[5],
                part_1[3], part_1[1], part_1[4],
                part_1[5], part_1[4], part_1[2];
            O_1 = tI_stress_1 * sigma_matrix; O_2 = tI_stress_2 * sigma_matrix;
            D_part_1 = O_1 - O_2;
            //cout << D_part_1.coeffRef(0, 0) + D_part_1.coeffRef(1, 1) + D_part_1.coeffRef(2, 2) << endl;
            grad_O_s[index[k_i]] += mpfr_float(O_mu) * (D_part_1.coeffRef(0, 0) + D_part_1.coeffRef(1, 1) + D_part_1.coeffRef(2, 2));
        }
    }
    //for (int i = 0; i < nv; i++) {
    //    cout << grad_O_s[3 * i] << " , " << grad_O_s[3 * i + 1] << " , " << grad_O_s[3 * i + 2] << endl;
    //}
    //vector<CppAD::AD<double>> vertices_N(n);
    //vector<double> vertices_N_x(n);
    //for (int i = 0; i < n; i++) {
    //    vertices_N[i] = CppAD::AD<double>(points[i].convert_to<double>());
    //    vertices_N_x[i] = points[i].convert_to<double>();
    //}
    for (int s_i = 0; s_i < n; s_i++) {
  
        int p_id = int(s_i / 3);
        int x_y_z = s_i % 3;
        N_AD = calculate_grad_N(p_id, x_y_z, nv, n_tri, points, triangles, n_sp);
        grad_f_s = N_AD * pressure;
        //cout << grad_f_s << endl;
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
                //if (G_row_indices[it] < 2) {
                //    cout << "index : " << G_row_indices[it] << " , " << G_col_indices[it] << endl;
                //}
                grad_M.insert(G_row_indices[it] + n, G_col_indices[it]) = jac_G[it * n + s_i];
                grad_M.insert(G_col_indices[it], G_row_indices[it] + n) = jac_G[it * n + s_i];
            }
        grad_M.makeCompressed();
        mpfr_float part2 = O_A.transpose() * grad_f_s;
        mpfr_float part3 = O_A.transpose() * grad_M * O_B;

        grad_O_s[s_i] += O_mu * (part2 - part3);
        //grad_O_s[s_i] += mpfr_float(O_mu) * ( - part3);


    }
    grad_s = grad_O_s;
    //auto end_time_stress = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time_stress - start_time_AD).count();
    //std::cout << "|| The cost of d_stress/d_s : " << duration << " seconds" << endl << endl << endl;

    //std::ostringstream stream;
    //// 设置宽度为3，左侧填充0
    //stream << std::setw(3) << std::setfill('0') << count_time;
    //std::string grad_save_str = "../Results/Grad_AD__" + stream.str() + ".txt";
    //char* grad_save = const_cast<char*>(grad_save_str.c_str());
    //if (!write_file(grad_save, grad_s, nv)) {
    //    cout << "Save FAILED!" << endl;
    //    return -1;
    //}
    return 1.0;
}


int main(int argc, char* argv[]) {
    mpfr_float::default_precision(100);
    string model_name = "ring_mmg";
    std::string path = "../IO/";
    std::string filename_str = path + model_name + ".mesh";
    char* filename = const_cast<char*>(filename_str.c_str());
    std::string worst_stress_faces_str = "../Results/Eworst_stress_face_" + model_name + ".txt";
    char* worst_stress_faces = const_cast<char*>(worst_stress_faces_str.c_str());
    MMG5_pSol       mmgSol;
    MMG5_int        k, np, n_tet, n_tri, n_sp;
    //np:number of points;  n_tet:number of tetraheras;    n_tri:number of triangles on surface;    n_sp:number of points on surface;
    mpfr_float* points;
    int* triangles;
    int* tetrahedras;
    mpfr_float* Normals;
    mpfr_float* Area;
    int* Index_of_sp;
    MMG5_pMesh mmgMesh;
    if (!read_mesh(mmgMesh, mmgSol, filename, points, triangles, tetrahedras)) {
        return 0;
    }
    if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, NULL, &n_tri, NULL, NULL) != 1)  exit(EXIT_FAILURE);

    Eigen::SparseMatrix<mpfr_float> N;
    Eigen::SparseMatrix<mpfr_float> G;
    Eigen::SparseMatrix<mpfr_float> M_G;
    Eigen::SparseMatrix<mpfr_float> K;

    int n = 3 * np;
    Calculate_Area(np, n_tri, points, triangles, Area, n_sp, Index_of_sp);
    search_neighbor(n_sp, Index_of_sp, n_tri, triangles);
    K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    N = Build_N(np, n_tri, points, triangles, n_sp, Index_of_sp, nullptr);

    G = Build_G(np, points, N, n_sp);        //计算G
    M_G = Build_M_G(np, points);        //计算G
    Eigen::SparseMatrix<mpfr_float> M = merge_matrix(K, M_G);
    Eigen::SparseLU<Eigen::SparseMatrix<mpfr_float>> solver_LU;

    //计算K的最大特征值.
    mpfr_float MaxEigenvalue = calculate_max_eigenvalue(K) * 100.0;
    cout << "n_sp:" << n_sp << endl;
    Eigen::SparseMatrix<mpfr_float> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    MatrixX_mp GGT = G * G.transpose();
    MatrixX_mp GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<mpfr_float> sparse_GGT_inverse = GGT_inverse.sparseView();
    cout << sparse_GGT_inverse << endl;
    Eigen::SparseMatrix<mpfr_float> GEP_matrix_sparse = ((I - G.transpose() * (sparse_GGT_inverse)*G).transpose() *
        N.transpose() * K * N *
        (I - G.transpose() * (sparse_GGT_inverse)*G) +
        G.transpose() * MaxEigenvalue * G
        );

    solver_LU.compute(M);
    if (solver_LU.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The inverse of K cannot be computed !");
    }
    Eigen::SparseMatrix<mpfr_float> A(n, n);
    for (int i = 0; i < np; i++) {
        A.insert(3 * i, 3 * i) = Area[i];
        A.insert(3 * i + 1, 3 * i + 1) = Area[i];
        A.insert(3 * i + 2, 3 * i + 2) = Area[i];
    }
    Eigen::SparseMatrix<mpfr_float> P = I - G.transpose() * (sparse_GGT_inverse)*G;
    Eigen::SparseMatrix<mpfr_float> matrix_B = N.transpose() * A * N;
    SparseSymMatProd<mpfr_float> opA(GEP_matrix_sparse);
    SparseCholesky<mpfr_float>  Bop(matrix_B);

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    SymGEigsSolver<SparseSymMatProd<mpfr_float>, SparseCholesky<mpfr_float>, GEigsMode::Cholesky>
        geigs(opA, Bop, 3, 18);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute(SortRule::SmallestAlge);

    // Retrieve results
    VectorX_mp evalues;
    MatrixX_mp evecs_GEP;

    if (geigs.info() == CompInfo::Successful)
    {
        evalues = geigs.eigenvalues();
        evecs_GEP = geigs.eigenvectors();
    }
    else {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("GEP cannot be computed !");
    }
    VectorX_mp p, u, f, y;
    mpfr_float max_energy = 0.0; int max_energy_p_i = 0;
    for (int p_i = 0; p_i < 3; p_i++) {
        p = evecs_GEP.col(p_i);
        VectorX_mp y = P * p;
        p = y;
        f = N * p;
        f.conservativeResize(n + 6);
        f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
        u = solver_LU.solve(f);
        mpfr_float energy = u.transpose() * f;
        if (energy > max_energy) {
            max_energy = energy;
            max_energy_p_i = p_i;
        }
    }
    p = evecs_GEP.col(max_energy_p_i);
    VectorX_mp p_stress(n_sp);
    p_stress = evecs_GEP.col(max_energy_p_i);
    y = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;
    VectorX_mp Gp = G * y;
    Gp = G * p;
    p = y;
    f = N * p;
    std::cout << "模长(f)：" << f.norm() << endl;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver_LU.solve(f);
    cout << "Energy E : " << u.transpose() * f << endl;
    //optimization
    int intersect = 0;
    bool loop = true, loop_t = true, loop_s = true;
    mpfr_float* x = (mpfr_float*)calloc(np * 3, sizeof(mpfr_float));
    vector<CppAD::AD<double>> t(1);
    vector<double> t_x(1);
    //points值赋给V和x
    for (size_t i = 0; i < n; ++i) {
        x[i] = points[i];
    }
    t[0] = CppAD::AD<double>(1e4); t_x[0] = 1e4;
    mpfr_float threshold = 1e-3, gamma = 0.5;
    double gamma_t = 0.5;
    mpfr_float threshold_s = 0.1;
    //initial
    Eigen::SparseMatrix<mpfr_float> sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
    cout << sparse_Stress.coeffRef(0, 0) << endl;
    CppAD::AD<double> O = Calculate_O(t[0], sparse_Stress, n_tet);
    while (O < -1e5) {
        std::cout << "Invalid Initial!!" << endl;
        t[0] *= 2.0; t_x[0] = Value(t[0]);
        sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
        O = Calculate_O(t[0], sparse_Stress, n_tet);
    }
    CppAD::AD<double> last_O = O, last_O_t = O, last_O_s = O;
    std::vector<double> jac_O_t;
    int loop_t_time = 0, loop_s_time = 0;
    mpfr_float alpha = 0.001;
    bool points_updated = false;
    int iter = 0;

    auto start_iter = std::chrono::high_resolution_clock::now();
    intersect = 0;
    //save content
    K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    M_G = Build_M_G(np, points);
    M = merge_matrix(K, M_G);
    solver_LU.compute(M);
    if (solver_LU.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The inverse of K cannot be computed !");
    }
    std::cout << endl << "----------------The iteration on loop:" << iter << "---------------" << endl;
    //optimize t 
    loop_t = true;
    sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
    cout << "--------------------Optimize t: ---------------------------------------------------------------------------------" << endl;
    double alpha_t = 1.0;
    while (loop_t) {
        vector<CppAD::AD<double>> O_var(1);
        CppAD::Independent(t);
        O = Calculate_O(t[0], sparse_Stress, n_tet);
        O_var[0] = O;
        last_O_t = O;
        CppAD::ADFun<double> func(t, O_var);    // 创建 ADFun 对象
        jac_O_t = func.Jacobian(t_x);
        std::vector<double> w(1);
        w[0] = 1.0;  // 权重
        std::vector<double> hess_O_t = func.Hessian(t_x, w);
        mpfr_float maxVal = 0.0;
        for (mpfr_float val : jac_O_t) {
            maxVal = max(maxVal, abs(val));
        }
        if (maxVal < threshold) {
            cout << "The max gradient on t : " << maxVal << endl; loop_t = false; break;
        }
        while (true) {
            CppAD::AD<double> temp = t[0];
            t[0] = t[0] - CppAD::AD<double>(alpha_t * (jac_O_t[0] / hess_O_t[0]));
            O = Calculate_O(t[0], sparse_Stress, n_tet);
            if (O< last_O_t && O > -1e5) {
                t_x[0] = Value(t[0]);
                alpha_t /= gamma_t;
                last_O_t = O;
                break;
            }
            else {
                t[0] = temp;
                alpha_t *= gamma_t;
            }
        }
        cout << "O:" << O << " The max of O_t : " << maxVal << "  t : " << t_x[0] << " dt:" << (jac_O_t[0] / hess_O_t[0]) << " |Alpha:" << alpha_t << endl;
    }

    //optimize s
    cout << "-------------------Optimize s: ----------Optimize s: ----------------- Optimize s --------------------------------------------------" << endl;
    cout << "The value of O after optimization on t : " << last_O_t << endl;
    VectorX_mp grad_s;
    Calculate_Stresses_AD(np, points, n_tri, triangles, n_tet, tetrahedras, n_sp, Index_of_sp, mpfr_float(t_x[0]), p_stress, grad_s, solver_LU);
    sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
    mpfr_float O_last = Calculate_O(mpfr_float(Value(t[0])), sparse_Stress, n_tet);
    Validation:
    VectorX_mp direction(n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0f, 1.0f);
    for (int i = 0; i < np; i++) {
        double directionX = dis(gen);
        double directionY = dis(gen);
        double directionZ = dis(gen);
        // 计算向量长度
        double length = std::sqrt(directionX * directionX + directionY * directionY + directionZ * directionZ);
        // 规范化向量
        directionX /= length;
        directionY /= length;
        directionZ /= length;
        direction[3 * i] = mpfr_float(directionX);
        direction[3 * i + 1] = mpfr_float(directionY);
        direction[3 * i + 2] = mpfr_float(directionZ);
    }
    mpfr_float eps = 1e-23;
    //做有限差分的验证
    for (int k = 0; k < 50; k++) {
        eps /= 10.0;
        for (int i = 0; i < n; i++) {
            points[i] = (x[i] + eps * direction[i]);
        }
        //After updating points,we need to update K,N,G->M,so we also need to rebuild a solver for new M.
        K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
        M_G = Build_M_G(np, points);        //计算G
        M = merge_matrix(K, M_G);
        Eigen::SparseLU<Eigen::SparseMatrix<mpfr_float>> solver_LU_new;
        solver_LU_new.compute(M);
        //sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU_new, intersect);
        sparse_Stress = Calculate_Stresses_for_test(np,x, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU_new, intersect);
        mpfr_float O_new = Calculate_O(mpfr_float(Value(t[0])), sparse_Stress, n_tet);
        mpfr_float dO = O_new - O_last;
        mpfr_float dOdt = dO / eps;
        mpfr_float grad_dot = grad_s.dot(direction);
        mpfr_float diffrence = dOdt - grad_dot;
        cout << "----------------------------------------- EPS: " << eps << "---------------------------------------" << endl;
        cout <<std::setprecision(20)<< "O_new : " << O_new << "   O : " << O_last << endl;
        cout << std::setprecision(20) << "dO : " << dO << endl;
        cout << std::setprecision(20) << "dOdt : " << dOdt << endl;
        cout << std::setprecision(20) << "Grad_dot : " << grad_dot << endl;
        cout << std::setprecision(20) << "diffrence : " << diffrence << endl;
        cout << "--------------------------------------------------------------------------------" << endl;
    }
    return 0;
    //{
    //    while (loop) {
    //        auto start_iter = std::chrono::high_resolution_clock::now();
    //        intersect = 0;
    //        //save content
    //        K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    //        M_G = Build_M_G(np, points);
    //        M = merge_matrix(K, M_G);
    //        solver_LU.compute(M);
    //        if (solver_LU.info() != Eigen::Success) {
    //            std::cerr << "Decomposition failed!" << std::endl;
    //            throw std::runtime_error("The inverse of K cannot be computed !");
    //        }
    //        std::cout << endl << "----------------The iteration on loop:" << iter << "---------------" << endl;
    //        //optimize t 
    //        loop_t = true;
    //        sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
    //        cout << "--------------------Optimize t: ---------------------------------------------------------------------------------" << endl;
    //        mpfr_float alpha_t = 1.0;
    //        //while (loop_t) {
    //        //    vector<AD_mpfr_float> O_var(1);
    //        //    CppAD::Independent(t);
    //        //    O = Calculate_O(t[0], sparse_Stress, n_tet);
    //        //    O_var[0] = O;
    //        //    last_O_t = O;
    //        //    CppAD::ADFun<mpfr_float> func(t, O_var);    // 创建 ADFun 对象
    //        //    jac_O_t = func.Jacobian(t_x);
    //        //    std::vector<mpfr_float> w(1);
    //        //    w[0] = 1.0;  // 权重
    //        //    std::vector<mpfr_float> hess_O_t = func.Hessian(t_x, w);
    //        //    mpfr_float maxVal = 0.0;
    //        //    for (mpfr_float val : jac_O_t) {
    //        //        maxVal = max(maxVal, abs(val));
    //        //    }
    //        //    if (maxVal < threshold) {
    //        //        cout << "The max gradient on t : " << maxVal << endl; loop_t = false; break;
    //        //    }
    //        //    while (true) {
    //        //        AD_mpfr_float temp = t[0];
    //        //        t[0] = t[0] - CppAD::AD<mpfr_float>(alpha_t * (jac_O_t[0] / hess_O_t[0]));
    //        //        O = Calculate_O(t[0], sparse_Stress, n_tet);
    //        //        if (O< last_O_t && O > -1e5) {
    //        //            t_x[0] = Value(t[0]);
    //        //            alpha_t /= gamma;
    //        //            last_O_t = O;
    //        //            break;
    //        //        }
    //        //        else {
    //        //            t[0] = temp;
    //        //            alpha_t *= gamma;
    //        //        }
    //        //    }
    //        //    cout << "O:" << O << " The max of O_t : " << maxVal << "  t : " << t_x[0] << " dt:" << (jac_O_t[0] / hess_O_t[0]) << " |Alpha:" << alpha_t << endl;
    //        //}
    //        //optimize s
    //        cout << "-------------------Optimize s: ----------Optimize s: ----------------- Optimize s --------------------------------------------------" << endl;
    //        cout << "The value of O after optimization on t : " << last_O_t << endl;
    //        //    VectorXd grad_s;
    //        //    Calculate_Stresses_AD(np, points, n_tri, triangles, n_tet, tetrahedras, n_sp, Index_of_sp, t_x[0], p_stress, grad_s, solver_LU);
    //        //    sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
    //        //    O = Calculate_O(t[0], sparse_Stress, n_tet);
    //        //    last_O_s = O;
    //        //    //Validation:
    //        //    VectorXd direction(n);
    //        //    std::random_device rd;
    //        //    std::mt19937 gen(rd());
    //        //    std::uniform_real_distribution<double> dis(-1.0f, 1.0f);
    //        //    for (int i = 0; i < np; i++) {
    //        //        double directionX = dis(gen);
    //        //        double directionY = dis(gen);
    //        //        double directionZ = dis(gen);
    //        //        // 计算向量长度
    //        //        double length = std::sqrt(directionX * directionX + directionY * directionY + directionZ * directionZ);
    //        //        // 规范化向量
    //        //        directionX /= length;
    //        //        directionY /= length;
    //        //        directionZ /= length;
    //        //        direction[3 * i] = directionX;
    //        //        direction[3 * i + 1] = directionY;
    //        //        direction[3 * i + 2] = directionZ;
    //        //    }
    //        //    double eps = 1e-1;
    //        //    double O_last = Value(O);
    //        //    //做有限差分的验证
    //        //    for (int k = 0; k < 10; k++) {
    //        //        eps /= 10.0;
    //        //        for (int i = 0; i < n; i++) {
    //        //            points[i] = (x[i] + eps * direction[i]);
    //        //        }
    //        //        //After updating points,we need to update K,N,G->M,so we also need to rebuild a solver for new M.
    //        //        K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    //        //        M_G = Build_M_G(np, points);        //计算G
    //        //        M = merge_matrix(K, M_G);
    //        //        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU_new;
    //        //        solver_LU_new.compute(M);
    //        //        sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU_new, intersect);
    //        //        O = Calculate_O(t[0], sparse_Stress, n_tet);
    //        //        double O_new = Value(O);
    //        //        double dO = O_new - O_last;
    //        //        double dOdt = dO / eps;
    //        //        double grad_dot = grad_s.dot(direction);
    //        //        cout << "----------------------------------------- EPS: " << eps << "---------------------------------------" << endl;
    //        //        cout << "O_new : " << O_new << "   O : " << O_last << endl;
    //        //        cout << "dO : " << dO << endl;
    //        //        cout << "dOdt : " << dOdt << endl;
    //        //        cout << "Grad_dot : " << grad_dot << endl;
    //        //        cout << "--------------------------------------------------------------------------------" << endl;
    //            //}
    //        return 0;
    //    }
    //}

    return 0;
}
