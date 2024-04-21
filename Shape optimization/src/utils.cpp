#include <Grad.h>
#include <DWCO.h>
#include <remesh.h>

//#include "../include/Grad.h"
//#include "../include/DWCO.h"
//#include "../include/remesh.h"

void extract_surface(Eigen::MatrixXd points, Eigen::MatrixXi triangles, unordered_map<int, int> & Old_2_new,
    Eigen::MatrixXd &V, Eigen::MatrixXi &F,Eigen::VectorXi & new_2_Old) {
    int n_tri = triangles.rows();
    std::unordered_set<int> index_sp;
    for (int k = 0; k < n_tri; k++) {
        int pi_id = triangles(k, 0);
        int pj_id = triangles(k, 1);
        int pk_id = triangles(k, 2);
        index_sp.insert(pi_id); index_sp.insert(pj_id); index_sp.insert(pk_id);
    }
    int nV = index_sp.size();
    int* index_of_points = (int*)calloc(index_sp.size(), sizeof(int));
    int index = 0;
    V.resize(nV, 3); V.setZero();
    F.resize(n_tri, 3); F.setZero();
    Old_2_new.clear();
    new_2_Old.resize(nV); new_2_Old.setZero();
    for (const auto& element : index_sp) {
        Old_2_new[element] = index;
        V.row(index) = points.row(element);
        new_2_Old[index++] = element;
    }
    for (int i = 0; i < n_tri; i++) {
        F(i, 0) = Old_2_new[triangles(i, 0)];
        F(i, 1) = Old_2_new[triangles(i, 1)];
        F(i, 2) = Old_2_new[triangles(i, 2)];
    }

}

double avg_edge_length(Eigen::MatrixXd V, Eigen::MatrixXi F) {

    double minEdgeLength = FLT_MAX; // Initialize with the maximum float value

    for (int i = 0; i < F.rows(); ++i) {
        // Get the length of each edge in the face
        for (int j = 0; j < F.cols(); ++j) {
            // Vertices of the edge
            int v0 = F(i, j);
            int v1 = F(i, (j + 1) % 3);

            // Compute the length of the edge
            double edgeLength = (V.row(v0) - V.row(v1)).norm();
            // Update the minimum edge length
            if (edgeLength < minEdgeLength) {
                minEdgeLength = edgeLength;
            }
        }
    }

    std::cout << "Shortest Edge Length: " << minEdgeLength << std::endl;
    return minEdgeLength;
}

void Calculate_Area(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::VectorXd &Area) {
    int nV = V.rows();
    int nF = F.rows();
    Area.resize(nV); Area.setZero();
    vector<int> nt_of_points(nV);
    for (int k = 0; k < nF; k++) {
        int pi_id = F(k, 0);
        int pj_id = F(k, 1);
        int pk_id = F(k, 2);
        Vector3d p_i = V.row(pi_id);
        Vector3d p_j = V.row(pj_id);
        Vector3d p_k = V.row(pk_id);
        //Vector3d p_i (V(pi_id, 0), V(pi_id, 1), V(pi_id, 2));
        //Vector3d p_j(V(pj_id, 0), V(pj_id, 1), V(pj_id, 2));
        //Vector3d p_k(V(pk_id, 0), V(pk_id, 1), V(pk_id, 2));
        Vector3d edge1 = p_j - p_i;
        Vector3d edge2 = p_k - p_i;
         //计算叉积
        Vector3d crossProduct = edge1.cross(edge2);

        // 计算面积
        double area_tri = 0.5 * crossProduct.norm();
        nt_of_points[pi_id] += 1; nt_of_points[pj_id] += 1; nt_of_points[pk_id] += 1;
        Area[pi_id] += area_tri; Area[pj_id] += area_tri; Area[pk_id] += area_tri;
    }
    for (int k = 0; k < nV; k++) {
        if (nt_of_points[k] != 0)
            Area[k] /= nt_of_points[k];
    }
}
