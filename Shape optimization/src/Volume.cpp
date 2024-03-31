//#include <Grad.h>
//#include <DWCO.h>
#include "../include/Grad.h"
#include "../include/DWCO.h"

double tetrahedronVolume(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C, Eigen::Vector3d D) {
    MatrixXd V_mat(3, 3);
    V_mat.col(0) = D - A;
    V_mat.col(1) = D - B;
    V_mat.col(2) = D - C;
    double volume = V_mat.determinant()/6.0;
    return volume;
}

void tetVolume_first_derivative(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C, Eigen::Vector3d D, Eigen::VectorXd &grad) {
    grad.resize(12); grad.setZero();
    MatrixXd V_mat(3, 3);
    V_mat.col(0) = D - A;
    V_mat.col(1) = D - B;
    V_mat.col(2) = D - C;

    for (int index_i = 0; index_i < 3; index_i++) {
        for (int pos = 0; pos < 3; pos++) {
            MatrixXd V_select(3, 3);
            V_select = V_mat;
            V_select.col(index_i).setZero();
            V_select.row(pos).setZero();
            V_select(pos, index_i) = - 1.0;
            double gradient = V_select.determinant() / 6.0;
            grad[3 * index_i + pos] = gradient;
        }
    }
    for (int pos = 0; pos < 3; pos++) {
        MatrixXd V_select(3, 3);
        V_select = V_mat;
        V_select.row(pos).setOnes();
        double gradient = V_select.determinant() / 6.0;
        grad[3 * 3 + pos] = gradient;
    }
}

double Energy_Volume(Eigen::MatrixXd points, Eigen::MatrixXi tetrahedras, Eigen::VectorXd& grad , unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old,
    double dhat_clb) {
    int nV = new_2_Old.size();
    int nT = tetrahedras.rows();
    grad.resize(3 * nV);    grad.setZero();
    ipc::ClampedLogBarrier clb;  // 创建ClampedLogBarrier实例
    double energy = 0.0;
    //Area.resize(nV); Area.setZero();
    //vector<int> nt_of_points(nV);
    for (int k = 0; k < nT; k++) {
        int t1 = tetrahedras(k, 0);
        int t2 = tetrahedras(k, 1);
        int t3 = tetrahedras(k, 2);
        int t4 = tetrahedras(k, 3);
        Vector3d p1 = points.row(t1);
        Vector3d p2 = points.row(t2);
        Vector3d p3 = points.row(t3);
        Vector3d p4 = points.row(t4);
        double volume = tetrahedronVolume(p1, p2, p3, p4);

        // 使用barrier函数计算值
        energy += clb(volume, dhat_clb);
        if (isinf(energy)) {
            return energy;
        }
        //std::cout << "Barrier value: " << barrier_value << std::endl;

        // 使用first_derivative计算一阶导数
        double first_derivative_value = clb.first_derivative(volume, dhat_clb);
        //std::cout << "First derivative: " << first_derivative_value << std::endl;
        VectorXd tet_grad;
        tetVolume_first_derivative(p1, p2, p3, p4, tet_grad);
        for (int t = 0; t < 4; t++) {
            if (Old_2_new.find(tetrahedras(k, t)) == Old_2_new.end())
                continue;
            int id = Old_2_new[tetrahedras(k, t)];
            grad.segment<3>(3 * id) += first_derivative_value * tet_grad.segment<3>(3 * t);
        }
    }
    return energy;
}
