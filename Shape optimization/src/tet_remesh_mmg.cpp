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
using namespace std;
/** Include the mmg3d library hader file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"
#include "../include/remesh.h"

int tet_remesh_mmg(string filename, string outputname, double edge_length) {
    MMG5_pMesh      mmgMesh;
    MMG5_pSol       mmgSol, mmgMet, tmpSol;
    int             i, j, opt;
    MMG5_int        np, ne, nprism, nt, nquad, na;
    /* To manually recover the mesh */
    int             nsol, typSol[MMG5_NSOLS_MAX];
    double* sols;

    /* Filenames */
    char* filein, * fileout;

    filein = const_cast<char*>(filename.c_str());
    fileout = const_cast<char*>(outputname.c_str());

    mmgMesh = NULL;
    mmgSol = NULL;
    mmgMet = NULL;

    MMG3D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
        MMG5_ARG_end);

    if (MMG3D_loadMesh(mmgMesh, filein) != 1)  exit(EXIT_FAILURE);
    if (MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_debug, 1) != 1)
        exit(EXIT_FAILURE);
    if (MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_verbose, 4) != 1)
        exit(EXIT_FAILURE);
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax, edge_length) != 1)
        exit(EXIT_FAILURE);
    /* Minimal mesh size (default 0)*/
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, edge_length * 0.01) != 1)
        exit(EXIT_FAILURE);
    //if (MMG3D_Set_dparameter(mmgMesh, NULL, MMG3D_DPARAM_hsiz, edge_length) != 1)
    //    exit(EXIT_FAILURE);
    //if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hgrad, 1.25) != 1)
    //    exit(EXIT_FAILURE);

    /* Global hausdorff value (default value = 0.01) applied on the whole boundary */
    //if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, 2.5 * 1e-3) != 1)
        //exit(EXIT_FAILURE);
    //if (MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_angle, 1) != 1)
    //    exit(EXIT_FAILURE);
    //if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_angleDetection, 60) != 1)
    //    exit(EXIT_FAILURE);

    /* Global hausdorff value (default value = 0.01) applied on the whole boundary */
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, 0.008) != 1)
        exit(EXIT_FAILURE);

    ///* Gradation control*/

    /** remesh function */
    int ier;
    ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

    if (MMG3D_Get_meshSize(mmgMesh, &np, &ne, &nprism, &nt, &nquad, &na) != 1)  exit(EXIT_FAILURE);
    cout << "np:" << np << endl;
    cout << "ne:" << ne << endl;
    cout << "nt:" << nt << endl;
    cout << "nquad:" << nquad << endl;
    cout << "na:" << na << endl;

    if (ier == MMG5_STRONGFAILURE) {
        fprintf(stdout, "BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
        return(ier);
    }
    else if (ier == MMG5_LOWFAILURE)
        fprintf(stdout, "BAD ENDING OF MMG3DLIB\n");
    //if (MMG3D_saveSol(mmgMesh, mmgSol, fileout) != 1)
    //    exit(EXIT_FAILURE);

    if (MMG3D_saveMesh(mmgMesh, fileout) != 1)
        exit(EXIT_FAILURE);

    /*s ave the solutions array */
    //if (MMG3D_saveAllSols(mmgMesh, &tmpSol, fileout) != 1)
    //    exit(EXIT_FAILURE);

    /** 3) Free the MMG3D structures */
    MMG3D_Free_allSols(mmgMesh, &mmgSol);

    MMG3D_Free_all(MMG5_ARG_start,
        MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppSols, NULL,
        MMG5_ARG_end);

    //free(filein);
    //filein = NULL;

    //free(fileout);
    //fileout = NULL;

    return 0;
}


int tet_remesh_mmg_Update_Points(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, double edge_length, double* points, int n_sp, int* index_sp) {

    for (int i = 0; i < n_sp; i++) {
        mmgMesh->point[index_sp[i] + 1].c[0] = points[3 * index_sp[i]];
        mmgMesh->point[index_sp[i] + 1].c[1] = points[3 * index_sp[i] + 1];
        mmgMesh->point[index_sp[i] + 1].c[2] = points[3 * index_sp[i] + 2];
    }
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax, edge_length) != 1)
        exit(EXIT_FAILURE);
    /* Minimal mesh size (default 0)*/
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, edge_length * 0.01) != 1)
        exit(EXIT_FAILURE);
    if (MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_verbose, 4) != 1)
        exit(EXIT_FAILURE);
    //if (MMG3D_Set_dparameter(mmgMesh, NULL, MMG3D_DPARAM_hsiz, edge_length) != 1)
    //    exit(EXIT_FAILURE);
    //if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hgrad, 1.2) != 1)
    //    exit(EXIT_FAILURE);

    /* Global hausdorff value (default value = 0.01) applied on the whole boundary */
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, 0.01) != 1)
        exit(EXIT_FAILURE);

    ///* Gradation control*/
    //if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hgrad, 2) != 1)
    //    exit(EXIT_FAILURE);

    /** remesh function */
    int ier;
    ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);


    if (ier == MMG5_STRONGFAILURE) {
        fprintf(stdout, "BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
        return(ier);
    }
    else if (ier == MMG5_LOWFAILURE)
        fprintf(stdout, "BAD ENDING OF MMG3DLIB\n");
    //if (MMG3D_saveSol(mmgMesh, mmgSol, fileout) != 1)
    //    exit(EXIT_FAILURE);

    std::string mesh_save_str = "../Results/new_mesh.mesh";
    char* outname = const_cast<char*>(mesh_save_str.c_str());
    if (MMG3D_saveMesh(mmgMesh, outname) != 1)
        exit(EXIT_FAILURE);

    return 1;
}

int tet_remesh_mmg_Update_Points(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, double edge_length, Eigen::MatrixXd points, int update, double& worst) {
    int nV = points.rows();
    if (update) {
        for (int i = 0; i < nV; i++) {
            mmgMesh->point[i + 1].c[0] = points.row(i)[0];
            mmgMesh->point[i + 1].c[1] = points.row(i)[1];
            mmgMesh->point[i + 1].c[2] = points.row(i)[2];
        }
    }
    if (MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_debug, 1) != 1)
        exit(EXIT_FAILURE);
    if (MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_nosurf, 0) != 1)
        exit(EXIT_FAILURE);
    if (MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_verbose, 4) != 1)
        exit(EXIT_FAILURE);
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax, edge_length) != 1)
        exit(EXIT_FAILURE);
    /* Minimal mesh size (default 0)*/
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, edge_length * 0.01) != 1)
        exit(EXIT_FAILURE);
    if (MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, 0.008) != 1)
        exit(EXIT_FAILURE);
    /** remesh function */
    int ier;
    ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

    int ntet;
    if (MMG3D_Get_meshSize(mmgMesh, NULL, &ntet, NULL, NULL, NULL, NULL) != 1)  exit(EXIT_FAILURE);

    double worst_quality = 1.0;
    for (int i = 0; i < ntet; i++) {
        double quality = MMG3D_Get_tetrahedronQuality(mmgMesh, mmgSol, i + 1);
        if (quality < worst_quality)
            worst_quality = quality;
    }
    worst = worst_quality;
    if (ier == MMG5_STRONGFAILURE) {
        fprintf(stdout, "BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
        return(ier);
    }
    else if (ier == MMG5_LOWFAILURE)
        fprintf(stdout, "BAD ENDING OF MMG3DLIB\n");
    //if (MMG3D_saveSol(mmgMesh, mmgSol, fileout) != 1)
    //    exit(EXIT_FAILURE);

    std::string mesh_save_str = "../Results/new_mesh.mesh";
    char* outname = const_cast<char*>(mesh_save_str.c_str());
    if (MMG3D_saveMesh(mmgMesh, outname) != 1)
        exit(EXIT_FAILURE);

    return 1;
}


int read_mesh(MMG5_pMesh& mmgMesh_out, MMG5_pSol& mmgSol_out, char* filename, Eigen::MatrixXd& points, Eigen::MatrixXi& triangles, Eigen::MatrixXi& tetrahedras) {
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
    points.resize(np, 3); points.setZero();
    triangles.resize(n_tri, 3); triangles.setZero();
    tetrahedras.resize(n_tet, 4); tetrahedras.setZero();

    for (k = 0; k < np; k++) {
        /** b) Vertex recovering */
        if (MMG3D_Get_vertex(mmgMesh, &(Point[0]), &(Point[1]), &(Point[2]), NULL, NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        points(k, 0) = Point[0];
        points(k, 1) = Point[1];
        points(k, 2) = Point[2];
    }

    for (k = 0; k < n_tet; k++) {
        if (MMG3D_Get_tetrahedron(mmgMesh, &(Tetrahedra[0]), &(Tetrahedra[1]), &(Tetrahedra[2]), &(Tetrahedra[3]), NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        tetrahedras(k, 0) = Tetrahedra[0] - 1;
        tetrahedras(k, 1) = Tetrahedra[1] - 1;
        tetrahedras(k, 2) = Tetrahedra[2] - 1;
        tetrahedras(k, 3) = Tetrahedra[3] - 1;
    }
    for (k = 0; k < n_tri; k++) {
        if (MMG3D_Get_triangle(mmgMesh, &(Triangle[0]), &(Triangle[1]), &(Triangle[2]), NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        triangles(k, 0) = Triangle[0] - 1;
        triangles(k, 1) = Triangle[1] - 1;
        triangles(k, 2) = Triangle[2] - 1;
    }
    return 1;
}
int get_info_from_mesh(MMG5_pMesh& mmgMesh, MMG5_pSol& mmgSol, double*& points, int*& triangles, int*& tetrahedras) {
    if (points != nullptr)
        free(points);
    if (triangles != nullptr)
        free(triangles);
    if (tetrahedras != nullptr)
        free(tetrahedras);
    int             i, j, k;
    MMG5_int        np, n_tet, nprism, n_tri, nquad, na;
    if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, &nprism, &n_tri, &nquad, &na) != 1)  exit(EXIT_FAILURE);
    cout << "number of points:" << np << endl;
    cout << "number of tetraheras:" << n_tet << endl;
    cout << "number of triagles on surface:" << n_tri << endl;
    cout << "nquad:" << nquad << endl;
    cout << "na:" << na << endl;
    int Triangle[3], Edge[2], Tetrahedra[4];
    double  Point[3], Sol;
    double* Points = (double*)calloc(np * 3, sizeof(double));
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


int get_info_from_mesh(MMG5_pMesh& mmgMesh, MMG5_pSol& mmgSol, Eigen::MatrixXd& points, Eigen::MatrixXi& triangles, Eigen::MatrixXi& tetrahedras) {
    int             i, j, k;
    MMG5_int        np, n_tet, nprism, n_tri, nquad, na;
    if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, &nprism, &n_tri, &nquad, &na) != 1)  exit(EXIT_FAILURE);
    cout << "number of points:" << np << endl;
    cout << "number of tetraheras:" << n_tet << endl;
    cout << "number of triagles on surface:" << n_tri << endl;
    cout << "nquad:" << nquad << endl;
    cout << "na:" << na << endl;
    int Triangle[3], Edge[2], Tetrahedra[4];
    double  Point[3], Sol;
    points.resize(np, 3); points.setZero();
    triangles.resize(n_tri, 3); triangles.setZero();
    tetrahedras.resize(n_tet, 4); tetrahedras.setZero();

    for (k = 0; k < np; k++) {
        /** b) Vertex recovering */
        if (MMG3D_Get_vertex(mmgMesh, &(Point[0]), &(Point[1]), &(Point[2]), NULL, NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        points(k, 0) = Point[0];
        points(k, 1) = Point[1];
        points(k, 2) = Point[2];
    }

    for (k = 0; k < n_tet; k++) {
        if (MMG3D_Get_tetrahedron(mmgMesh, &(Tetrahedra[0]), &(Tetrahedra[1]), &(Tetrahedra[2]), &(Tetrahedra[3]), NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        tetrahedras(k, 0) = Tetrahedra[0] - 1;
        tetrahedras(k, 1) = Tetrahedra[1] - 1;
        tetrahedras(k, 2) = Tetrahedra[2] - 1;
        tetrahedras(k, 3) = Tetrahedra[3] - 1;
    }
    for (k = 0; k < n_tri; k++) {
        if (MMG3D_Get_triangle(mmgMesh, &(Triangle[0]), &(Triangle[1]), &(Triangle[2]), NULL, NULL) != 1)
            exit(EXIT_FAILURE);
        triangles(k, 0) = Triangle[0] - 1;
        triangles(k, 1) = Triangle[1] - 1;
        triangles(k, 2) = Triangle[2] - 1;
    }
    return 1;
}