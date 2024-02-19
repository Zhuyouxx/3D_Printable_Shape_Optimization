
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
#include <unordered_set>
#include <unordered_map>
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
#include "mmg/mmg3d/libmmg3d.h"
#include "remesh.h"

//typedef Eigen::Matrix< double, Dynamic, Dynamic> MatrixXd;
//typedef Eigen::Matrix< CppAD::AD<double>, Dynamic, Dynamic> MatrixXd_AD;
//typedef Eigen::Matrix<CppAD::AD<double>, 3, 1> Vector3d_AD;

int write_file(char* filename, double* points, int np);
int write_file(char* filename, VectorXd points, int np);
int write_file(char* filename, vector<double> points, int np);
int write_file(char* filename, int* triangles, int nt);
int read_mesh(MMG5_pMesh& mmgMesh_out, MMG5_pSol& mmgSol_out, char* filename, double*& points, int*& triangles, int*& tetrahedras);
Eigen::SparseMatrix<double> merge_matrix(Eigen::SparseMatrix<double> K, Eigen::SparseMatrix<double> G);
double algebraicCofactor(Eigen::Matrix< double, 4, 4>& matrix, int i, int j);
Eigen::SparseMatrix<double> Build_stiffness_Matrix(int nv, double*& vertices, int n_tet, int* tetrahedras);
int Calculate_Area(int nv, int n_tri, double* points, int* triangles, double*& Area, int& n_sp, int*& Index_of_sp);
template <typename T>
T Calculate_log_A(const vector<T>& points, int n_tet, int* tetrahedras);
SparseMatrix<double> Build_N(int nv, int n_tri, double* points, int* triangles, int n_sp, int* index_sp, double** normals);
SparseMatrix<double> Build_G(int nv, double* vertices, SparseMatrix<double> N, int n_sp);
SparseMatrix<double> Build_M_G(int nv, double* vertices);
double calculate_max_eigenvalue(SparseMatrix<double> K);
Eigen::SparseMatrix<double> Calculate_Stresses(int nv, double* vertices, int n_sp, int* index_sp, int n_tri, int* triangles, int n_tet, int* tetrahedras,
    VectorXd p, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver, int& intersect);
Eigen::SparseMatrix<double> Calculate_Stresses(int nv, double* vertices, int n_sp, int* index_sp, int n_tri, int* triangles, int n_tet, int* tetrahedras,
    VectorXd p, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver, int& intersect, double& log_V);
bool isPositiveSemiDefinite(const Eigen::MatrixXd& A);
CppAD::AD<double> algebraicCofactor(Eigen::Matrix< CppAD::AD<double>, 4, 4>& matrix, int i, int j);
CppAD::AD<double> Calculate_O(const CppAD::AD<double> t, Eigen::SparseMatrix<double> stress, int n_tet);
CppAD::AD<double> Calculate_O(const CppAD::AD<double> t, Eigen::SparseMatrix<double> stress, int n_tet, double log_V);
Eigen::SparseMatrix<CppAD::AD<double>> Build_stiffness_Matrix(int nv, const vector<CppAD::AD<double>>& vertices, int n_tet, int* tetrahedras);

SparseMatrix<CppAD::AD<double>> Build_M_G(int nv, const vector<CppAD::AD<double>>& vertices);
SparseMatrix<CppAD::AD<double>> Build_N(int nv, int n_tri, const vector<CppAD::AD<double>>& points, int* triangles, int n_sp, int* index_sp);
VectorXd Calculate_sigma_A(int nv, double* points, int n_tet, int* tetrahedras, VectorXd u, double t);

double Calculate_Stresses_AD(int nv, double* points, int n_tri, int* triangles, int n_tet, int* tetrahedras, int n_sp, int* index_sp, double t,
    VectorXd p, vector<double>& grad_s, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver);

int cal_stress_save(char* filename, int np, double* vertices, int n_sp, int* index_sp, int n_tri, int* triangles, int n_tet, int* tetrahedras,
    VectorXd p, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver, int& intersect);

int extract_surface(int n_sp, int* Index_of_sp, int n_tri, int* triangles, double* points, double** new_points, int** new_triangles);

int extract_info(string filename, int n_sp, int* Index_of_sp, int n_tet, int* tetrahedras, int n_tri, int* triangles, double* points, SparseMatrix<double> sparse_stress, VectorXd f = VectorXd(), vector<double> grad = vector<double>());

int extract_info_with_old_grad(string filename, int n_sp, int* Index_of_sp, int n_tet, int* tetrahedras, int n_tri, int* triangles, double* points, SparseMatrix<double> sparse_stress,
    VectorXd f, vector<double> grad, vector<double> grad_old);

int stress_grad(string filein,string fileout, std::vector<double>& grad_3d,bool info_f_g);
int smooth_grad(int n_sp, int* index_sp, int* points, int n_tri, int* triangles, double* Area, vector<double>& grad,double *Normal, double lambda);
void Smooth_AD(double lambda, int nv, double* points, int n_tri, int* triangles, int n_sp, int* index_sp, vector<double>& grad_s);
int stress_grad();
