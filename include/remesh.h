
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
using namespace std;
/** Include the mmg3d library hader file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"
using namespace std;


int tet_mesh(double& ideal_edge_length ,std::string filename, const std::string& output_name = "", double rel = 0.05, bool m_f = false, bool smooth = false);
int tet_remesh_mmg(string filename, string outputname, double edge_length);
int tet_remesh_mmg_Update_Points(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, double edge_length, double* points, int n_sp, int* index_sp);