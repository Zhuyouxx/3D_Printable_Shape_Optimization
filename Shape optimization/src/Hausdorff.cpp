#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <cmath>
#include <random>
#include <vector>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(0.0, 1.0);

struct Barycentric {
    double alpha;
    double beta;
    double gamma;
};
std::vector<Barycentric> generateBarycentricCoordinates(int n) {
    std::vector<Barycentric> samples;
    std::srand(std::time(0));  // 使用当前时间作为随机种子
        // Random number generation setup

    for (int i = 0; i < n; ++i) {
        double r1 = dis(gen);
        double r2 = dis(gen);
        double alpha = 1.0 - r1;
        double beta = r1 * (1.0 - r2);
        double gamma = r1 * r2;
        samples.push_back({ alpha, beta, gamma });
    }

    return samples;
}

double energy_hausdorff(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd V_target, Eigen::MatrixXi F_target, Eigen::VectorXd &grad) {
    //igl::readOBJ("D://IDE//VisualStudio_Project//DR_manufacturing//TestOn3D//Test3D//IO//ball_01.obj", V, F);
    //igl::readOBJ("D://IDE//VisualStudio_Project//DR_manufacturing//TestOn3D//Test3D//IO//ball_02.obj", V_target, F_target);
    std::vector<Triangle> triangles;
    for (int i = 0; i < F_target.rows(); ++i) {
        Point p1(V_target(F_target(i, 0), 0), V_target(F_target(i, 0), 1), V_target(F_target(i, 0), 2));
        Point p2(V_target(F_target(i, 1), 0), V_target(F_target(i, 1), 1), V_target(F_target(i, 1), 2));
        Point p3(V_target(F_target(i, 2), 0), V_target(F_target(i, 2), 1), V_target(F_target(i, 2), 2));
        triangles.push_back(Triangle(p1, p2, p3));
    }
    // 构建AABB树
    Tree tree(triangles.begin(), triangles.end());


    int nFaces = F.rows();
    int nVerts = V.rows();

    double average_area = 0.0;
    double max_area = 0.0, min_area = 100.0;
    int id = -1;
    Eigen::VectorXd count(nVerts);
    count.setZero();
    for (int i = 0; i < nFaces; i++) {
        Eigen::Vector3d p1 = V.row(F(i, 0));
        Eigen::Vector3d p2 = V.row(F(i, 1));
        Eigen::Vector3d p3 = V.row(F(i, 2));
        count(F(i, 0)) += 1.0;
        count(F(i, 1)) += 1.0;
        count(F(i, 2)) += 1.0;
        double area = 0.5 * (p2 - p1).cross(p3 - p1).norm();
        average_area += area;
        if (area < min_area) {
            min_area = area;
            id = i;
        }
        else if (area > max_area)
            max_area = area;
    }
    average_area /= nFaces;
    //std::cout << "The max area is: " << max_area << " || The min area is : " << min_area << std::endl;
    Eigen::Vector3d p1 = V.row(F(id, 0));
    Eigen::Vector3d p2 = V.row(F(id, 1));
    Eigen::Vector3d p3 = V.row(F(id, 2));
    double area = 0.5 * (p2 - p1).cross(p3 - p1).norm();
    //std::cout << "p1: " << p1 << std::endl;
    //std::cout << "p2: " << p2 << std::endl;
    //std::cout << "p3: " << p3 << std::endl;
    //std::cout << "The min area is: " << area <<" || ID : "<<id << std::endl;

    //double rho = 1.0 / min_area;
    double rho = 1.0 / average_area;
    //std::cout << "For more sampling points in the unit, the selected rho is : " << rho << std::endl;
    int max_sample_num = (int)rho * max_area;
    //std::cout << "The max sampling points is : " << max_sample_num << std::endl;
    double energy = 0.0;
    Eigen::VectorXd derivative;
    derivative.resize(3 * nVerts);
    derivative.setZero();
    for (int i = 0; i < nFaces; i++) {
        Eigen::Vector3d p1 = V.row(F(i, 0));
        Eigen::Vector3d p2 = V.row(F(i, 1));
        Eigen::Vector3d p3 = V.row(F(i, 2));
        double area = 0.5 * (p2 - p1).cross(p3 - p1).norm();
        int sample_num = (int)rho * area;
        if (sample_num <= 0) {
            sample_num = 1;
        }
        //std::cout << sample_num << std::endl;
        auto samples = generateBarycentricCoordinates(sample_num);
        double dis_sum = 0.0;
        for (auto sample : samples) {
            Eigen::Vector3d sample_point = sample.alpha * p1 + sample.beta * p2 + sample.gamma * p3;
            //// 查询最近点
            Point query(sample_point(0), sample_point(1), sample_point(2)); // 示例查询点
            Point closest_point_pos = tree.closest_point(query);
            Eigen::Vector3d closet_point(closest_point_pos.x(), closest_point_pos.y(), closest_point_pos.z());

            double dis = (sample_point - closet_point).squaredNorm();
            dis_sum += dis;

            Eigen::Vector3d gradient = 2.0 * (sample_point - closet_point);
            derivative.segment<3>(3 * F(i, 0)) += area / sample_num * sample.alpha * gradient;
            derivative.segment<3>(3 * F(i, 1)) += area / sample_num * sample.beta * gradient;
            derivative.segment<3>(3 * F(i, 2)) += area / sample_num * sample.gamma * gradient;
        }
        energy += area / sample_num * dis_sum;
    }
    grad = derivative;



    return energy;
}