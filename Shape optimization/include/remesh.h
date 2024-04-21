
#include <string>
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
#include <Eigen/Core>
#include <unordered_set>
#include <unordered_map>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
using namespace std;
/** Include the mmg3d library hader file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"
using namespace std;


//int tet_mesh(double& ideal_edge_length ,std::string filename, const std::string& output_name = "", double rel = 0.05, bool m_f = false, bool smooth = false);
int tet_mesh(double& ideal_edge_length, string filename, const std::string& output_name, double rel, bool m_f, bool smooth);
int tet_remesh_mmg(string filename, string outputname, double edge_length);
int tet_remesh_mmg_Update_Points(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, double edge_length, double* points, int n_sp, int* index_sp);
int tet_remesh_mmg_Update_Points(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, double edge_length, Eigen::MatrixXd points, int update, double& worst);
int get_info_from_mesh(MMG5_pMesh& mmgMesh, MMG5_pSol& mmgSol, double*& points, int*& triangles, int*& tetrahedras);
int get_info_from_mesh(MMG5_pMesh& mmgMesh, MMG5_pSol& mmgSol, Eigen::MatrixXd& points, Eigen::MatrixXi& triangles, Eigen::MatrixXi& tetrahedras);
int read_mesh(MMG5_pMesh& mmgMesh_out, MMG5_pSol& mmgSol_out, char* filename, Eigen::MatrixXd& points, Eigen::MatrixXi& triangles, Eigen::MatrixXi& tetrahedras);

void Calculate_Area(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXd& Area);
void extract_surface(Eigen::MatrixXd points, Eigen::MatrixXi triangles, unordered_map<int, int>& Old_2_new,
    Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& new_2_Old);

double get_worst_quality(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol);
double get_worst_quality(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, Eigen::MatrixXd points, Eigen::MatrixXi triangles, Eigen::MatrixXi tetrahedras);
//int tetgen(const std::string& filename, const std::string& file_out, string options, Eigen::MatrixXd& TV, Eigen::MatrixXi& TT, Eigen::MatrixXi& TF);
double avg_edge_length(Eigen::MatrixXd V, Eigen::MatrixXi F);