#include <Grad.h>
#include <DWCO.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>


double energy_self_intersection(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXi edges, Eigen::VectorXd& grad, const double dhat) {
    int nV = V.rows();
    grad.resize(3 * nV); grad.setZero();
    //std::cout << V.outerSize() << std::endl;
    igl::edges(F, edges);
    ipc::CollisionMesh collision_mesh(V, edges, F);
    //Eigen::MatrixXd vertices = collision_mesh.rest_positions();
    //const double dhat = 1e-3;
    ipc::Collisions collisions;
    collisions.build(collision_mesh, V, dhat);
    const ipc::BarrierPotential B(dhat);
    double barrier_potential = B(collisions, collision_mesh, V);
    //std::cout << "barrier_potential : " << barrier_potential << std::endl;
    //std::cout << "The Barrier_potential is : " << barrier_potential << std::endl;
    Eigen::VectorXd barrier_potential_grad =
        B.gradient(collisions, collision_mesh, V);

    return barrier_potential;
}

bool is_intersect_free(Eigen::MatrixXd V_init, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXi edges, const double dhat) {
    ipc::CollisionMesh collision_mesh(V_init, edges, F);
    //Eigen::MatrixXd vertices = collision_mesh.rest_positions();
    //vertices.col(1) *= 0.01; // Squash the bunny in the y-direction
    //vertices(0, 0) = 0.0;
    bool is_intersect = is_step_collision_free(collision_mesh, V, V_init);
    //std::cout << "is intersect:" << is_intersect << endl;
    return is_intersect;
}