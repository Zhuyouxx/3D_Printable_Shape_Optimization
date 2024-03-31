#pragma once
#include <ipc/ipc.hpp>
#include <ipc/collisions/collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>


double Energy_Volume(Eigen::MatrixXd points, Eigen::MatrixXi tetrahedras, Eigen::VectorXd& grad, unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old,
    double dhat_clb);
void tetVolume_first_derivative(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C, Eigen::Vector3d D, Eigen::VectorXd& grad);
double tetrahedronVolume(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C, Eigen::Vector3d D);

double Energy_shell(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXd& grad);
double energy_self_intersection(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXi edges, Eigen::VectorXd& grad, const double dhat);
bool is_intersect_free(Eigen::MatrixXd V_init, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXi edges, const double dhat);
double energy_hausdorff(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd V_target, Eigen::MatrixXi F_target, Eigen::VectorXd& grad);