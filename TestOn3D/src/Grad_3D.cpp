
#include "Grad.h"
double mu = 0.3;
double E = 1.0; //double E = 1e9;
double dt = 1e-6;
double O_mu = 1e-3;
int count_time = 0;
double edge_length = 0.05;
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
int read_mesh(MMG5_pMesh& mmgMesh_out, MMG5_pSol& mmgSol_out, char* filename, double*& points, int*& triangles, int*& tetrahedras) {
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

Eigen::SparseMatrix<double> merge_matrix(Eigen::SparseMatrix<double> K, Eigen::SparseMatrix<double> G) {
    // 行合并
    Eigen::SparseMatrix<double> Merged(K.rows() + G.rows(), K.rows() + G.rows());
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            Merged.insert(it.row(), it.col()) = it.value();
        }
    }
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
            Merged.insert(K.rows() + it.row(), it.col()) = it.value();
        }
    }
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
            Merged.insert(it.col(), K.rows() + it.row()) = it.value();
        }
    }
    Merged.makeCompressed();
    return Merged;
}

// 计算矩阵的代数余子式
double algebraicCofactor(Eigen::Matrix< double, 4, 4>& matrix, int i, int j) {
    int n = matrix.rows();
    Eigen::MatrixXd subMatrix(n - 1, n - 1);

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
    double cofactor = std::pow(-1, i + j) * subMatrix.determinant();
    return cofactor;
}

Eigen::SparseMatrix<double> Build_stiffness_Matrix(int nv, double*& vertices, int n_tet, int* tetrahedras) {
    Eigen::SparseMatrix<double> sparse_K(nv * 3, nv * 3);
    Eigen::Matrix<double, 6, 6>D;
    Eigen::Matrix<double, 6, 12>B; B.setZero();
    Eigen::Matrix<double, 12, 12>Ke;
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        Matrix<double, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        double Volume = Lambda.determinant() / 6.0;
        for (int i = 0; i < 4; i++) {
            double b_i = algebraicCofactor(Lambda, i, 1);
            double c_i = algebraicCofactor(Lambda, i, 2);
            double d_i = algebraicCofactor(Lambda, i, 3);
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
                if (Ke(i, j) != 0) {
                    sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                }
            }
    }
    sparse_K.makeCompressed();
    return sparse_K;
}

int Calculate_Area(int nv, int n_tri, double* points, int* triangles, double*& Area, int& n_sp, int*& Index_of_sp) {
    std::unordered_set<int> index_sp;
    double* area = (double*)calloc(nv, sizeof(double));
    int* nt_of_points = (int*)calloc(nv, sizeof(double));
    for (int k = 0; k < n_tri; k++) {
        int pi_id = triangles[3 * k];
        int pj_id = triangles[3 * k + 1];
        int pk_id = triangles[3 * k + 2];
        index_sp.insert(pi_id); index_sp.insert(pj_id); index_sp.insert(pk_id);
        Vector3d p_i(points[3 * pi_id], points[3 * pi_id + 1], points[3 * pi_id + 2]);
        Vector3d p_j(points[3 * pj_id], points[3 * pj_id + 1], points[3 * pj_id + 2]);
        Vector3d p_k(points[3 * pk_id], points[3 * pk_id + 1], points[3 * pk_id + 2]);

        Vector3d edge1 = p_j - p_i;
        Vector3d edge2 = p_k - p_i;
        // 计算叉积
        Vector3d crossProduct = edge1.cross(edge2);

        // 计算面积
        double area_tri = 0.5 * crossProduct.norm();
        if (area_tri < 0.0) {
            std::cerr << "Wrong area!" << std::endl;
            throw std::runtime_error("Wrong area!");
            return -1;
        }
        nt_of_points[pi_id] += 1; nt_of_points[pj_id] += 1; nt_of_points[pk_id] += 1;
        area[pi_id] += area_tri; area[pj_id] += area_tri; area[pk_id] += area_tri;
    }
    for (int k = 0; k < nv; k++) {
        if (nt_of_points[k] != 0)
            area[k] /= nt_of_points[k];
    }
    int* index_of_points = (int*)calloc(index_sp.size(), sizeof(double));
    int index = 0;
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

SparseMatrix<double> Build_N(int nv, int n_tri, double* points, int* triangles, int n_sp, int* index_sp, double** normals) {
    Eigen::SparseMatrix<double> N(nv * 3, n_sp);
    double* Normals = (double*)calloc(nv * 3, sizeof(double));
    for (int k = 0; k < n_tri; k++) {
        int pi_id = triangles[3 * k];
        int pj_id = triangles[3 * k + 1];
        int pk_id = triangles[3 * k + 2];
        Vector3d p_i(points[3 * pi_id], points[3 * pi_id + 1], points[3 * pi_id + 2]);
        Vector3d p_j(points[3 * pj_id], points[3 * pj_id + 1], points[3 * pj_id + 2]);
        Vector3d p_k(points[3 * pk_id], points[3 * pk_id + 1], points[3 * pk_id + 2]);
        Vector3d edge1 = p_j - p_i;
        Vector3d edge2 = p_k - p_i;
        Vector3d n_of_tirangle = edge1.cross(edge2).normalized();

        Normals[3 * pi_id] += n_of_tirangle[0]; Normals[3 * pi_id + 1] += n_of_tirangle[1]; Normals[3 * pi_id + 2] += n_of_tirangle[2];
        Normals[3 * pj_id] += n_of_tirangle[0]; Normals[3 * pj_id + 1] += n_of_tirangle[1]; Normals[3 * pj_id + 2] += n_of_tirangle[2];
        Normals[3 * pk_id] += n_of_tirangle[0]; Normals[3 * pk_id + 1] += n_of_tirangle[1]; Normals[3 * pk_id + 2] += n_of_tirangle[2];
    }
    for (int k = 0; k < n_sp; k++) {
        int id = index_sp[k];
        Vector3d Normal_Of_P(Normals[3 * id], Normals[3 * id + 1], Normals[3 * id + 2]);
        Normal_Of_P.normalize();
        N.insert(3 * id, k) = -Normal_Of_P[0];
        N.insert(3 * id + 1, k) = -Normal_Of_P[1];
        N.insert(3 * id + 2, k) = -Normal_Of_P[2];
        Normals[3 * id] = -Normal_Of_P[0];
        Normals[3 * id + 1] = -Normal_Of_P[1];
        Normals[3 * id + 2] = -Normal_Of_P[2];
    }
    if (normals != nullptr)
        *normals = Normals;
    else
        free(Normals);
    return N;
}

SparseMatrix<double> Build_G(int nv, double* vertices, SparseMatrix<double> N, int n_sp) {
    Eigen::SparseMatrix<double> G(6, nv * 3);
    Eigen::SparseMatrix<double> G_N(6, n_sp);
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

SparseMatrix<double> Build_M_G(int nv, double* vertices) {
    Eigen::SparseMatrix<double> G(6, nv * 3);
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
double calculate_max_eigenvalue(SparseMatrix<double> K) {
    SparseSymMatProd<double> opK(K);
    SymEigsSolver<SparseSymMatProd<double>> eigs_K(opK, 1, 6);
    eigs_K.init();
    eigs_K.compute(SortRule::LargestAlge);
    if (eigs_K.info() != CompInfo::Successful) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The Eigenvalues of K were not found. !");
    }
    Eigen::MatrixXd evecs;
    Eigen::VectorXd evalues;
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

Eigen::SparseMatrix<double> Calculate_Stresses(int nv, double* vertices, int n_sp, int* index_sp, int n_tri, int* triangles, int n_tet, int* tetrahedras,
    VectorXd p, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver, int& intersect) {
    int n = 3 * nv;
    Eigen::SparseMatrix<double> N = Build_N(nv, n_tri, vertices, triangles, n_sp, index_sp, nullptr);
    Eigen::SparseMatrix<double> G = Build_G(nv, vertices, N, n_sp);
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<double> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    Eigen::VectorXd u(n + 6);
    Eigen::VectorXd f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver.solve(f);
    // 检查解是否成功
    std::cout << "The norm of u is : " << u.norm() << std::endl;

    Eigen::SparseMatrix<double> sparse_Stress(n_tet * 6, 1);
    Eigen::Matrix<double, 6, 1>S_tress_e;
    Eigen::Matrix<double, 12, 1>ue;
    Eigen::Matrix<double, 6, 6>D;
    Eigen::Matrix<double, 6, 12>B; B.setZero();
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        Matrix<double, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        double Volume = Lambda.determinant() / 6.0;
        if (Volume < 0.0) {
            intersect = 1;
            break;
        }
        for (int i = 0; i < 4; i++) {
            double b_i = algebraicCofactor(Lambda, i, 1);
            double c_i = algebraicCofactor(Lambda, i, 2);
            double d_i = algebraicCofactor(Lambda, i, 3);
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

//gradient
CppAD::AD<double> algebraicCofactor(Eigen::Matrix< CppAD::AD<double>, 4, 4>& matrix, int i, int j) {
    int n = matrix.rows();
    MatrixXd_AD subMatrix(n - 1, n - 1);

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
    CppAD::AD<double> cofactor = std::pow(-1, i + j) * subMatrix.determinant();
    return cofactor;
}

CppAD::AD<double> Calculate_O(const CppAD::AD<double> t, Eigen::SparseMatrix<double> stress, int n_tet) {
    CppAD::AD<double> O = 0.0;
    for (int i = 0; i < n_tet; i++) {
        Eigen::Matrix<CppAD::AD<double>, 3, 3> sigma;
        Eigen::Matrix<CppAD::AD<double>, 3, 3> I;
        CppAD::AD<double> det1, det2;
        sigma << stress.coeff(3 * i, 0), stress.coeff(3 * i + 3, 0), stress.coeff(3 * i + 5, 0),
            stress.coeff(3 * i + 3, 0), stress.coeff(3 * i + 1, 0), stress.coeff(3 * i + 4, 0),
            stress.coeff(3 * i + 5, 0), stress.coeff(3 * i + 4, 0), stress.coeff(3 * i + 2, 0);
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
        Sigma << stress.coeff(3 * i, 0), stress.coeff(3 * i + 3, 0), stress.coeff(3 * i + 5, 0),
            stress.coeff(3 * i + 3, 0), stress.coeff(3 * i + 1, 0), stress.coeff(3 * i + 4, 0),
            stress.coeff(3 * i + 5, 0), stress.coeff(3 * i + 4, 0), stress.coeff(3 * i + 2, 0);
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

Eigen::SparseMatrix<CppAD::AD<double>> Build_stiffness_Matrix(int nv, const vector<CppAD::AD<double>>& vertices, int n_tet, int* tetrahedras) {
    Eigen::SparseMatrix<CppAD::AD<double>> sparse_K(nv * 3, nv * 3);
    Eigen::Matrix<CppAD::AD<double>, 6, 6>D;
    Eigen::Matrix<CppAD::AD<double>, 6, 12>B; B.setZero();
    Eigen::Matrix<CppAD::AD<double>, 12, 12>Ke;
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        CppAD::AD<double> p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        CppAD::AD<double> a1, a2, a3, b1, b2, b3, c1, c2, c3;
        Matrix<CppAD::AD<double>, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        CppAD::AD<double> Volume = Lambda.determinant() / 6.0;
        for (int i = 0; i < 4; i++) {
            CppAD::AD<double> b_i = algebraicCofactor(Lambda, i, 1);
            CppAD::AD<double> c_i = algebraicCofactor(Lambda, i, 2);
            CppAD::AD<double> d_i = algebraicCofactor(Lambda, i, 3);
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
                if (Ke(i, j) != 0) {
                    sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                }
            }
    }
    sparse_K.makeCompressed();
    return sparse_K;
}

SparseMatrix<CppAD::AD<double>> Build_M_G(int nv, const vector<CppAD::AD<double>>& vertices) {
    Eigen::SparseMatrix<CppAD::AD<double>> G(6, nv * 3);
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

SparseMatrix<CppAD::AD<double>> Build_N(int nv, int n_tri, const vector<CppAD::AD<double>>& points, int* triangles, int n_sp, int* index_sp) {
    Eigen::SparseMatrix<CppAD::AD<double>> N(nv * 3, n_sp);
    vector<CppAD::AD<double>> Normals(nv * 3, CppAD::AD<double>(0));
    for (int k = 0; k < n_tri; k++) {
        int pi_id = triangles[3 * k];
        int pj_id = triangles[3 * k + 1];
        int pk_id = triangles[3 * k + 2];
        Vector3d_AD p_i(points[3 * pi_id], points[3 * pi_id + 1], points[3 * pi_id + 2]);
        Vector3d_AD p_j(points[3 * pj_id], points[3 * pj_id + 1], points[3 * pj_id + 2]);
        Vector3d_AD p_k(points[3 * pk_id], points[3 * pk_id + 1], points[3 * pk_id + 2]);
        Vector3d_AD edge1 = p_j - p_i;
        Vector3d_AD edge2 = p_k - p_i;
        Vector3d_AD n_of_tirangle = edge1.cross(edge2).normalized();

        Normals[3 * pi_id] += n_of_tirangle[0]; Normals[3 * pi_id + 1] += n_of_tirangle[1]; Normals[3 * pi_id + 2] += n_of_tirangle[2];
        Normals[3 * pj_id] += n_of_tirangle[0]; Normals[3 * pj_id + 1] += n_of_tirangle[1]; Normals[3 * pj_id + 2] += n_of_tirangle[2];
        Normals[3 * pk_id] += n_of_tirangle[0]; Normals[3 * pk_id + 1] += n_of_tirangle[1]; Normals[3 * pk_id + 2] += n_of_tirangle[2];
    }
    for (int k = 0; k < n_sp; k++) {
        int id = index_sp[k];
        Vector3d_AD Normal_Of_P(Normals[3 * id], Normals[3 * id + 1], Normals[3 * id + 2]);
        Normal_Of_P.normalize();
        N.insert(3 * id, k) = -Normal_Of_P[0];
        N.insert(3 * id + 1, k) = -Normal_Of_P[1];
        N.insert(3 * id + 2, k) = -Normal_Of_P[2];
    }
    return N;
}


VectorXd Calculate_sigma_A(int nv, double* points, int n_tet, int* tetrahedras, VectorXd u, double t) {
    int n = 3 * nv;
    VectorXd stress(n + 6); stress.setZero();
    VectorXd ue(12);
    VectorXd se(6);
    VectorXd Stress_e(6);
    VectorXd D_e(6);
    Eigen::Matrix<double, 6, 6>D;
    Eigen::Matrix<double, 6, 12>B; B.setZero();
    Eigen::Matrix<double, 6, 12>sigma;
    D << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D *= E / ((1 + mu) * (1 - 2 * mu));
    for (int k = 0; k < n_tet; k++) {
        double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z;
        int t1, t2, t3, t4;
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        Matrix<double, 4, 4> Lambda;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        p1_x = points[t1 * 3]; p1_y = points[t1 * 3 + 1]; p1_z = points[t1 * 3 + 2];
        p2_x = points[t2 * 3]; p2_y = points[t2 * 3 + 1]; p2_z = points[t2 * 3 + 2];
        p3_x = points[t3 * 3]; p3_y = points[t3 * 3 + 1]; p3_z = points[t3 * 3 + 2];
        p4_x = points[t4 * 3]; p4_y = points[t4 * 3 + 1]; p4_z = points[t4 * 3 + 2];
        Lambda << 1, p1_x, p1_y, p1_z,
            1, p2_x, p2_y, p2_z,
            1, p3_x, p3_y, p3_z,
            1, p4_x, p4_y, p4_z;
        double Volume = Lambda.determinant() / 6.0;
        for (int i = 0; i < 4; i++) {
            double b_i = algebraicCofactor(Lambda, i, 1);
            double c_i = algebraicCofactor(Lambda, i, 2);
            double d_i = algebraicCofactor(Lambda, i, 3);
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
        MatrixXd tI_stress_1(3, 3), tI_stress_2(3, 3);
        MatrixXd O(3, 3);
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


double Calculate_Stresses_AD(int nv, double* points, int n_tri, int* triangles, int n_tet, int* tetrahedras, int n_sp, int* index_sp, double t,
    VectorXd p, vector<double>& grad_s, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver) {
    int n = 3 * nv;
    auto start_time_AD = std::chrono::high_resolution_clock::now();
    double* Normal;

    //1-计算位移 u 
    Eigen::SparseMatrix<double> K = Build_stiffness_Matrix(nv, points, n_tet, tetrahedras);
    Eigen::SparseMatrix<double> N = Build_N(nv, n_tri, points, triangles, n_sp, index_sp, &Normal);
    Eigen::SparseMatrix<double> G = Build_G(nv, points, N, n_sp);
    Eigen::SparseMatrix<double> grad_M(n + 6, n + 6); grad_M.setZero();

    vector<double> grad_O_s(n);
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<double> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    Eigen::VectorXd u(n + 6);
    Eigen::VectorXd f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver.solve(f);
    // 检查解是否成功
    std::cout << "The norm of u is : " << u.norm() << std::endl;
    //-计算位移 u over -1

    //2-计算gradient
    VectorXd ue(12);
    Eigen::VectorXd Stress_e;
    VectorXd grad_f_s(n + 6); grad_f_s.setZero();
    Eigen::Matrix<double, 6, 12>sigma; sigma.setZero();
    vector<double> s_x(12);
    vector<CppAD::AD<double>> s(12);

    Eigen::Matrix<CppAD::AD<double>, 6, 6>D_ad;
    Eigen::Matrix<CppAD::AD<double>, 6, 12>B_ad; B_ad.setZero();
    Eigen::Matrix<CppAD::AD<double>, 6, 12>sigma_s;
    Eigen::SparseMatrix<CppAD::AD<double>> N_AD;
    Eigen::SparseMatrix<CppAD::AD<double>> G_AD;
    std::vector<CppAD::AD<double>> Sigma_AD_Vector;
    std::vector<double> jac_Sigma_s;
    std::vector<double> jac_K;
    std::vector<double> jac_G;
    std::vector<size_t> row_indices;
    std::vector<size_t> col_indices;
    std::vector<size_t> G_row_indices;
    std::vector<size_t> G_col_indices;
    std::vector<CppAD::AD<double>> K_AD_Vector;
    std::vector<CppAD::AD<double>> G_AD_Vector;
    Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1> pressure(n_sp);

    for (int k = 0; k < p.size(); k++) {
        pressure[k] = p[k];
    }

    auto start_time_Np = std::chrono::high_resolution_clock::now();

    //2-1:calculate O_B
    VectorXd O_B = solver.solve(f);
    O_B = u;
    //2-3:calculate O_A
    VectorXd sigma_A = Calculate_sigma_A(nv, points, n_tet, tetrahedras, u, t);
    VectorXd O_A = solver.solve(sigma_A);

    //calculate Grad_K
    vector<CppAD::AD<double>> vertices(n);
    vector<double> vertices_x(n);
    vector<CppAD::AD<double>> vertices_sp(n_sp * 3);
    vector<double> vertices_sp_x(n_sp * 3);
    for (int i = 0; i < n; i++) {
        vertices[i] = CppAD::AD<double>(points[i]);
        vertices_x[i] = points[i];
    }
    for (int i = 0; i < n_sp; i++) {
        int id = index_sp[i];
        vertices_sp[3 * i] = CppAD::AD<double>(points[3 * id]);
        vertices_sp[3 * i + 1] = CppAD::AD<double>(points[3 * id + 1]);
        vertices_sp[3 * i + 2] = CppAD::AD<double>(points[3 * id + 2]);
        vertices_sp_x[3 * i] = points[3 * id]; vertices_sp_x[3 * i + 1] = points[3 * id + 1]; vertices_sp_x[3 * i + 2] = points[3 * id + 2];
    }
    cout << "number of triangles: " << n_tri << endl;
    cout << "The number of vertices: " << nv << endl;
    cout << "The number of vertices on the surface: " << n_sp << endl;
    auto START_TIME2 = std::chrono::high_resolution_clock::now();

    //CppAD::Independent(vertices);
    CppAD::Independent(vertices_sp);
    for (int i = 0; i < n_sp; i++) {
        int id = index_sp[i];
        vertices[3 * id] = vertices_sp[3 * i]; vertices[3 * id + 1] = vertices_sp[3 * i + 1]; vertices[3 * id + 2] = vertices_sp[3 * i + 2];
    }
    Eigen::SparseMatrix<CppAD::AD<double>> K_AD = Build_stiffness_Matrix(nv, vertices, n_tet, tetrahedras);
    // 遍历 K 的非零元素
    row_indices.clear(); col_indices.clear(); K_AD_Vector.clear();
    for (int k = 0; k < K_AD.outerSize(); ++k) {
        for (Eigen::SparseMatrix<CppAD::AD<double>>::InnerIterator it(K_AD, k); it; ++it) {
            row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
            col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
            K_AD_Vector.push_back(CppAD::AD<double>(it.value()));// it.value()     // 非零元素的值
        }
    }
    CppAD::ADFun<double> func(vertices_sp, K_AD_Vector);    // 创建 ADFun 对象
    jac_K = func.Jacobian(vertices_sp_x);
    //CppAD::ADFun<double> func(vertices, K_AD_Vector);    // 创建 ADFun 对象
    //jac_K = func.Jacobian(vertices_x);
    auto end_time9 = std::chrono::high_resolution_clock::now();
    auto duration_Np9 = std::chrono::duration_cast<std::chrono::microseconds>(end_time9 - START_TIME2).count() / 1e6;
    std::cout << "||------ The cost of calculating grad_K : " << duration_Np9 << " seconds ------||" << endl << endl;

    //calculate grad_G
    CppAD::Independent(vertices);
    G_AD = Build_M_G(nv, vertices);
    G_row_indices.clear(); G_col_indices.clear(); G_AD_Vector.clear();
    for (int k = 0; k < G_AD.outerSize(); ++k) {
        for (Eigen::SparseMatrix<CppAD::AD<double>>::InnerIterator it(G_AD, k); it; ++it) {
            G_row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
            G_col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
            G_AD_Vector.push_back(CppAD::AD<double>(it.value()));// it.value()     // 非零元素的值
        }
    }
    CppAD::ADFun<double> G_func(vertices, G_AD_Vector);    // 创建 ADFun 对象
    jac_G = G_func.Jacobian(vertices_x);
    auto end_time12 = std::chrono::high_resolution_clock::now();
    auto duration_Np12 = std::chrono::duration_cast<std::chrono::microseconds>(end_time12 - end_time9).count() / 1e6;
    std::cout << "||------ The cost of calculating the grad_G : " << duration_Np12 << " seconds ------||" << endl << endl;

    //calculate grad_sigma(s)
    D_ad << (1 - mu), mu, mu, 0, 0, 0,
        mu, (1 - mu), mu, 0, 0, 0,
        mu, mu, (1 - mu), 0, 0, 0,
        0, 0, 0, (1 - 2 * mu) / 2, 0, 0,
        0, 0, 0, 0, (1 - 2 * mu) / 2, 0,
        0, 0, 0, 0, 0, (1 - 2 * mu) / 2;
    D_ad *= E / ((1 + mu) * (1 - 2 * mu));
    Matrix<CppAD::AD<double>, 4, 4> Lambda;
    for (int k = 0; k < n_tet; k++) {
        int t1, t2, t3, t4;
        t1 = tetrahedras[k * 4]; t2 = tetrahedras[k * 4 + 1]; t3 = tetrahedras[k * 4 + 2]; t4 = tetrahedras[k * 4 + 3];
        int index[] = { 3 * t1, 3 * t1 + 1, 3 * t1 + 2,
                    3 * t2, 3 * t2 + 1, 3 * t2 + 2,
                    3 * t3, 3 * t3 + 1, 3 * t3 + 2,
                    3 * t4, 3 * t4 + 1, 3 * t4 + 2 };
        for (int index_i = 0; index_i < 12; index_i++) {
            s[index_i] = CppAD::AD<double>(points[index[index_i]]);
            s_x[index_i] = points[index[index_i]];
            ue[index_i] = u[index[index_i]];
        }
        CppAD::Independent(s);
        Lambda << 1, s[0], s[1], s[2],
            1, s[3], s[4], s[5],
            1, s[6], s[7], s[8],
            1, s[9], s[10], s[11];
        CppAD::AD<double> Volume = Lambda.determinant() / 6.0;
        double Volume_value = Value(Volume);
        for (int i = 0; i < 4; i++) {
            CppAD::AD<double> b_i = algebraicCofactor(Lambda, i, 1);
            CppAD::AD<double> c_i = algebraicCofactor(Lambda, i, 2);
            CppAD::AD<double> d_i = algebraicCofactor(Lambda, i, 3);
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
        CppAD::ADFun<double> sigma_fun(s, Sigma_AD_Vector);
        jac_Sigma_s = sigma_fun.Jacobian(s_x);

        Stress_e = sigma * ue;
        MatrixXd tI_stress_1(3, 3), tI_stress_2(3, 3);
        MatrixXd O(3, 3);
        tI_stress_1 << t * 1.0 - Stress_e[0], -Stress_e[3], -Stress_e[5],
            -Stress_e[3], t * 1.0 - Stress_e[1], -Stress_e[4],
            -Stress_e[5], -Stress_e[4], t * 1.0 - Stress_e[2];
        tI_stress_2 << t * 1.0 + Stress_e[0], Stress_e[3], Stress_e[5],
            Stress_e[3], t * 1.0 + Stress_e[1], Stress_e[4],
            Stress_e[5], Stress_e[4], t * 1.0 + Stress_e[2];
        tI_stress_1 = tI_stress_1.inverse();
        tI_stress_2 = tI_stress_2.inverse();

        VectorXd M_Np_select(12);
        for (int t_i = 0; t_i < 12; t_i++)
            M_Np_select[t_i] = O_B[index[t_i]];
        for (int k_i = 0; k_i < 12; k_i++) {
            int id = (int)index[k_i] / 3;
            if (Old_2_new.find(id) == Old_2_new.end())
                continue;
            Eigen::MatrixXd grad_sigma_s_k(6, 12);
            for (int sigma_i = 0; sigma_i < 6; sigma_i++)
                for (int sigma_j = 0; sigma_j < 12; sigma_j++)
                    grad_sigma_s_k.coeffRef(sigma_i, sigma_j) = jac_Sigma_s[(sigma_i * 12 + sigma_j) * 12 + k_i];

            VectorXd part_1 = grad_sigma_s_k * M_Np_select;
            MatrixXd sigma_matrix(3, 3);
            MatrixXd D_part_1(3, 3), O_1(3, 3), O_2(3, 3);;
            sigma_matrix << part_1[0], part_1[3], part_1[5],
                part_1[3], part_1[1], part_1[4],
                part_1[5], part_1[4], part_1[2];
            O_1 = tI_stress_1 * sigma_matrix; O_2 = tI_stress_2 * sigma_matrix;
            D_part_1 = O_1 - O_2;
            grad_O_s[index[k_i]] += O_mu * (D_part_1.coeffRef(0, 0) + D_part_1.coeffRef(1, 1) + D_part_1.coeffRef(2, 2));
        }
    }
    for (int s_i = 0; s_i < n_sp; s_i++) {
        int s_i_old = index_sp[s_i];
        Vector3d normal; normal << Normal[3 * s_i_old], Normal[3 * s_i_old + 1], Normal[3 * s_i_old + 2];
        for (int s_i_cord = 0; s_i_cord < 3; s_i_cord++) {
            int s_i_old_i = 3 * s_i_old + s_i_cord;
            int s_i_i = 3 * s_i + s_i_cord;
            vector<CppAD::AD<double>> s_k(1);
            vector<double> s_k_x(1);
            s_k[0] = CppAD::AD<double>(points[s_i_old_i]);
            s_k_x[0] = points[s_i_old_i];
            CppAD::Independent(s_k);
            vertices[s_i_old_i] = s_k[0];
            N_AD = Build_N(nv, n_tri, vertices, triangles, n_sp, index_sp);
            Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1> F_vector_Matrix = N_AD * pressure;
            std::vector<CppAD::AD<double>> F_vector;
            for (int k_i = 0; k_i < F_vector_Matrix.outerSize(); ++k_i)
                for (int num = 0; num < F_vector_Matrix.innerSize(); num++)
                    F_vector.push_back(F_vector_Matrix(num, k_i));
            CppAD::ADFun<double> N_s(s_k, F_vector);
            vector<double> jac_N(3 * nv);
            jac_N = N_s.Jacobian(s_k_x);
            for (int f_i = 0; f_i < 3 * nv; f_i++)
                grad_f_s[f_i] = jac_N[f_i];

            grad_M.setZero();
            int non_zero = K_AD_Vector.size();
            for (int it = 0; it < non_zero; it++)
                if (jac_K[it * 3 * n_sp + s_i_i] != 0) {
                    grad_M.insert(row_indices[it], col_indices[it]) = jac_K[it * 3 * n_sp + s_i_i];
                }
            grad_M.makeCompressed();
            non_zero = G_AD_Vector.size();
            for (int it = 0; it < non_zero; it++)
                if (jac_G[it * 3 * n_sp + s_i_i] != 0) {
                    grad_M.insert(G_row_indices[it] + n, G_col_indices[it]) = jac_G[it * 3 * n_sp + s_i_i];
                    grad_M.insert(G_col_indices[it], G_row_indices[it] + n) = jac_G[it * 3 * n_sp + s_i_i];
                }
            grad_M.makeCompressed();
            double part2 = O_A.transpose() * grad_f_s;
            double part3 = O_A.transpose() * grad_M * O_B;

            grad_O_s[s_i_old_i] += O_mu * (part2 - part3);
            vertices[s_i_old_i] = CppAD::AD<double>(points[s_i_old_i]);


        }
        Vector3d grad_p; grad_p << grad_O_s[3 * s_i_old], grad_O_s[3 * s_i_old + 1], grad_O_s[3 * s_i_old + 2];
        grad_p = (grad_p.dot(normal) / normal.dot(normal)) * normal;
        grad_O_s[3 * s_i_old] = grad_p[0]; grad_O_s[3 * s_i_old + 1] = grad_p[1]; grad_O_s[3 * s_i_old + 2] = grad_p[2];
    }
    grad_s = grad_O_s;
    auto end_time_stress = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time_stress - start_time_AD).count();
    std::cout << "|| The cost of d_stress/d_s : " << duration << " seconds" << endl << endl << endl;

    std::ostringstream stream;
    // 设置宽度为3，左侧填充0
    stream << std::setw(3) << std::setfill('0') << count_time;
    std::string grad_save_str = "../Results/Grad_AD__" + stream.str() + ".txt";
    char* grad_save = const_cast<char*>(grad_save_str.c_str());
    if (!write_file(grad_save, grad_s, nv)) {
        cout << "Save FAILED!" << endl;
        return -1;
    }
    return 1.0;
}

int cal_stress_save(char* filename, int np, double* vertices, int n_sp, int* index_sp, int n_tri, int* triangles, int n_tet, int* tetrahedras,
    VectorXd p, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver, int& intersect) {

    SparseMatrix<double> Stress_face = Calculate_Stresses(np, vertices, n_sp, index_sp, n_tri, triangles, n_tet, tetrahedras, p, solver, intersect);
    VectorXd S_face_points(n_tet);
    for (int k = 0; k < n_tet; k++) {
        double sigma_k = max(max(abs(Stress_face.coeff(6 * k, 0)), abs(Stress_face.coeff(6 * k + 1, 0))), abs(Stress_face.coeff(6 * k + 2, 0)));
        S_face_points(k) = sigma_k;
    }
    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }

    // 遍历数组并写入每个点到文件
    for (int i = 0; i < n_tet; i++) {
        file << S_face_points(i, 0) << std::endl;
    }
    // 关闭文件流
    file.close();
    return 1;
}

int extract_surface(int n_sp, int* Index_of_sp, int n_tri, int* triangles, double* points, double** new_points, int** new_triangles) {
    unordered_map<int, int> vertexMap;
    int* New_triangles = (int*)calloc((n_tri) * 3, sizeof(int));
    if (new_points != nullptr) {
        double* New_points = (double*)calloc((n_sp) * 3, sizeof(double));
        for (int i = 0; i < n_sp; i++) {
            vertexMap[Index_of_sp[i]] = i;
            New_points[3 * i] = points[3 * Index_of_sp[i]];
            New_points[3 * i + 1] = points[3 * Index_of_sp[i] + 1];
            New_points[3 * i + 2] = points[3 * Index_of_sp[i] + 2];
        }
    }
    else
        for (int i = 0; i < n_sp; i++)
            vertexMap[Index_of_sp[i]] = i;
    for (int t = 0; t < n_tri; t++) {
        New_triangles[3 * t] = vertexMap[triangles[3 * t]];
        New_triangles[3 * t + 1] = vertexMap[triangles[3 * t + 1]];
        New_triangles[3 * t + 2] = vertexMap[triangles[3 * t + 2]];
    }
    if (new_points != nullptr)
        *new_triangles = New_triangles;

    std::string save_str = "../Results/test.obj";
    char* filename = const_cast<char*>(save_str.c_str());

    std::ofstream file(filename);

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }
    // 遍历数组并写入每个点到文件
    for (int i = 0; i < n_sp; i++) {
        file << "v " << points[3 * Index_of_sp[i]] << " " << points[3 * Index_of_sp[i] + 1] << " " << points[3 * Index_of_sp[i] + 2] << endl;
    }
    for (int t = 0; t < n_tri; t++) {
        file << "f " << New_triangles[3 * t] + 1 << " " << New_triangles[3 * t + 1] + 1 << " " << New_triangles[3 * t + 2] + 1 << endl;
    }
    // 关闭文件流
    file.close();
}

int extract_info(string filename, int n_sp, int* Index_of_sp, int n_tet, int* tetrahedras, int n_tri, int* triangles, double* points, SparseMatrix<double> sparse_stress, VectorXd f, vector<double> grad) {
    unordered_map<int, int> vertexMap;
    int* New_triangles = (int*)calloc((n_tri) * 3, sizeof(int));
    for (int i = 0; i < n_sp; i++)
        vertexMap[Index_of_sp[i]] = i;

    for (int t = 0; t < n_tri; t++) {
        New_triangles[3 * t] = vertexMap[triangles[3 * t]];
        New_triangles[3 * t + 1] = vertexMap[triangles[3 * t + 1]];
        New_triangles[3 * t + 2] = vertexMap[triangles[3 * t + 2]];
    }

    std::string save_str = filename + ".obj";
    char* save_file = const_cast<char*>(save_str.c_str());
    std::ofstream file(save_file);
    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }
    // 遍历数组并写入每个点到文件
    for (int i = 0; i < n_sp; i++) {
        file << "v " << points[3 * Index_of_sp[i]] << " " << points[3 * Index_of_sp[i] + 1] << " " << points[3 * Index_of_sp[i] + 2] << endl;
    }
    for (int t = 0; t < n_tri; t++) {
        file << "f " << New_triangles[3 * t] + 1 << " " << New_triangles[3 * t + 1] + 1 << " " << New_triangles[3 * t + 2] + 1 << endl;
    }
    // 关闭文件流
    file.close();
    string save_path = "../Results/";
    save_str = save_path + filename + "_points.txt";
    save_file = const_cast<char*>(save_str.c_str());
    file.open(save_file);
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }
    for (int i = 0; i < n_sp; i++)
        file << points[3 * Index_of_sp[i]] << endl << points[3 * Index_of_sp[i] + 1] << endl << points[3 * Index_of_sp[i] + 2] << endl;
    file.close();

    save_str = save_path + filename + "_triangles.txt";
    save_file = const_cast<char*>(save_str.c_str());
    file.open(save_file);
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }
    for (int t = 0; t < n_tri; t++)
        file << New_triangles[3 * t] << endl << New_triangles[3 * t + 1] << endl << New_triangles[3 * t + 2] << endl;
    file.close();

    vector<int> count(n_sp);
    vector<double> stress(n_sp);
    for (int i = 0; i < n_tet; i++) {
        for (int t = 0; t < 4; t++) {
            int id = tetrahedras[4 * i + t];
            if (vertexMap.find(id) != vertexMap.end()) {
                stress[vertexMap[id]] = max(max(abs(sparse_stress.coeffRef(6 * i, 0)), abs(sparse_stress.coeffRef(6 * i + 1, 0))), abs(sparse_stress.coeffRef(6 * i + 2, 0)));
                count[vertexMap[id]]++;
            }
        }
    }
    save_str = save_path + filename + "_stress.txt";
    save_file = const_cast<char*>(save_str.c_str());
    file.open(save_file);
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 0;
    }
    for (int t = 0; t < n_sp; t++) {
        double stress_id = stress[t] / count[t];
        file << stress_id << endl;
    }
    file.close();


    if (f.size() > 0) {
        save_str = save_path + filename + "_force.txt";
        save_file = const_cast<char*>(save_str.c_str());
        file.open(save_file);
        if (!file.is_open()) {
            std::cerr << "Failed to open file." << std::endl;
            return 0;
        }
        for (int i = 0; i < n_sp; i++)
            file << f[3 * Index_of_sp[i]] << endl << f[3 * Index_of_sp[i] + 1] << endl << f[3 * Index_of_sp[i] + 2] << endl;
        file.close();
    }

    if (grad.size() > 0) {
        save_str = save_path + filename + "_grad.txt";
        save_file = const_cast<char*>(save_str.c_str());
        file.open(save_file);
        if (!file.is_open()) {
            std::cerr << "Failed to open file." << std::endl;
            return 0;
        }
        for (int i = 0; i < n_sp; i++)
            file << grad[3 * Index_of_sp[i]] << endl << grad[3 * Index_of_sp[i] + 1] << endl << grad[3 * Index_of_sp[i] + 2] << endl;
        file.close();
    }

}

int stress_grad(string filein, string fileout, std::vector<double>& vec_grad, bool info_f_g) {
    double edge_length;
    //string model = "cone";
    string model = filein;
    string remesh_model = model + ".obj";
    string remesh_out = model + ".mesh";
    string remesh_mmg_out = model + "_mmg.mesh";
    //tet_remesh(0.01, remesh_model, remesh_out);
    auto start_iter1 = std::chrono::high_resolution_clock::now();
    tet_mesh(edge_length, remesh_model, remesh_out, 0.05, true, true);
    auto end_time_iter = std::chrono::high_resolution_clock::now();
    auto duration_iter = std::chrono::duration_cast<std::chrono::seconds>(end_time_iter - start_iter1).count();
    std::cout << "||Time taken by function1: " << duration_iter << " seconds" << endl << endl << endl;

    tet_remesh_mmg(remesh_out, remesh_mmg_out, edge_length * 0.5);
    auto end_time_iter2 = std::chrono::high_resolution_clock::now();
    duration_iter = std::chrono::duration_cast<std::chrono::seconds>(end_time_iter2 - end_time_iter).count();
    std::cout << "||Time taken by function2: " << duration_iter << " seconds" << endl << endl << endl;
    string model_name = model + "_mmg";
    std::string path = "../IO/";
    std::string filename_str = model_name + ".mesh";
    char* filename = const_cast<char*>(filename_str.c_str());
    std::string worst_stress_faces_str = "../Results/Eworst_stress_face_" + model_name + ".txt";
    char* worst_stress_faces = const_cast<char*>(worst_stress_faces_str.c_str());
    MMG5_pSol       mmgSol;
    MMG5_int        k, np, n_tet, n_tri, n_sp;
    //np:number of points;  n_tet:number of tetraheras;    n_tri:number of triangles on surface;    n_sp:number of points on surface;
    double* points;
    int* triangles;
    int* tetrahedras;
    double* Normals;
    double* Area;
    int* Index_of_sp;
    MMG5_pMesh mmgMesh;
    if (!read_mesh(mmgMesh, mmgSol, filename, points, triangles, tetrahedras)) {
        return 0;
    }
    if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, NULL, &n_tri, NULL, NULL) != 1)  exit(EXIT_FAILURE);

    Eigen::SparseMatrix<double> N;
    Eigen::SparseMatrix<double> G;
    Eigen::SparseMatrix<double> M_G;
    Eigen::SparseMatrix<double> K;

    int n = 3 * np;
    Calculate_Area(np, n_tri, points, triangles, Area, n_sp, Index_of_sp);

    double* new_points;
    int* new_triangles;
    extract_surface(n_sp, Index_of_sp, n_tri, triangles, points, &new_points, &new_triangles);
    //转变为double的SparseMatrix
    K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    N = Build_N(np, n_tri, points, triangles, n_sp, Index_of_sp, nullptr);
    G = Build_G(np, points, N, n_sp);        //计算G
    M_G = Build_M_G(np, points);        //计算G
    Eigen::SparseMatrix<double> M = merge_matrix(K, M_G);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU;

    //计算K的最大特征值.
    double MaxEigenvalue = calculate_max_eigenvalue(K);
    cout << "n_sp:" << n_sp << endl;
    Eigen::SparseMatrix<double> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();

    Eigen::SparseMatrix<double> GEP_matrix_sparse = ((I - G.transpose() * (sparse_GGT_inverse)*G).transpose() *
        N.transpose() * K * N *
        (I - G.transpose() * (sparse_GGT_inverse)*G) +
        G.transpose() * MaxEigenvalue * G
        );

    solver_LU.compute(M);
    if (solver_LU.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The inverse of K cannot be computed !");
    }
    Eigen::SparseMatrix<double> A(n, n);
    for (int i = 0; i < np; i++) {
        A.insert(3 * i, 3 * i) = Area[i];
        A.insert(3 * i + 1, 3 * i + 1) = Area[i];
        A.insert(3 * i + 2, 3 * i + 2) = Area[i];
    }
    Eigen::SparseMatrix<double> P = I - G.transpose() * (sparse_GGT_inverse)*G;
    Eigen::SparseMatrix<double> matrix_B = N.transpose() * A * N;
    SparseSymMatProd<double> opA(GEP_matrix_sparse);
    SparseCholesky<double>  Bop(matrix_B);

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    SymGEigsSolver<SparseSymMatProd<double>, SparseCholesky<double>, GEigsMode::Cholesky>
        geigs(opA, Bop, 3, 15);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute(SortRule::SmallestAlge);

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs_GEP;

    if (geigs.info() == CompInfo::Successful)
    {
        evalues = geigs.eigenvalues();
        evecs_GEP = geigs.eigenvectors();
    }
    else {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("GEP cannot be computed !");
    }
    VectorXd p, u, f, y;
    double max_energy = 0.0; int max_energy_p_i = 0;
    for (int p_i = 0; p_i < 3; p_i++) {
        p = evecs_GEP.col(p_i);
        Eigen::VectorXd y = P * p;
        p = y;
        f = N * p;
        f.conservativeResize(n + 6);
        f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
        u = solver_LU.solve(f);
        double energy = u.transpose() * f;
        if (energy > max_energy) {
            max_energy = energy;
            max_energy_p_i = p_i;
        }
    }
    p = evecs_GEP.col(max_energy_p_i);
    VectorXd p_stress(n_sp);
    p_stress = evecs_GEP.col(max_energy_p_i);
    y = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;
    VectorXd Gp = G * y;
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
    double* x = (double*)calloc(np * 3, sizeof(double));
    vector<CppAD::AD<double>> t(1);
    vector<double> t_x(1);
    //points值赋给V和x
    for (size_t i = 0; i < n; ++i) {
        x[i] = points[i]; // 
    }
    t[0] = 1e4; t_x[0] = 1e4;
    double threshold = 1e-3, gamma = 0.5;
    double threshold_s = 0.1;
    //initial
    //Eigen::SparseMatrix<double> sparse_Stress = Calculate_Stresses(np);
    Eigen::SparseMatrix<double> sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
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
    double alpha = 0.001;
    bool points_updated = false;
    int iter = 0;
    while (loop) {
        auto start_iter = std::chrono::high_resolution_clock::now();
        intersect = 0;
        //save content
        //std::ostringstream stream;
        //// 设置宽度为3，左侧填充0
        //stream << std::setw(3) << std::setfill('0') << iter;
        //std::string point_save_str = "../Results/new_Points_" + model_name + "_" + stream.str() + ".txt";
        //char* point_save = const_cast<char*>(point_save_str.c_str());
        //if (!write_file(point_save, points, np)) {
        //    std::cout << "Save FAILED!" << endl;
        //    return -1;
        //}
        //std::string T_save_str = "../Results/Triangles_" + model_name + "_" + stream.str() + ".txt";
        //char* T_save = const_cast<char*>(T_save_str.c_str());
        //if (!write_file(T_save, triangles, n_tri)) {
        //    std::cout << "Save FAILED!" << endl;
        //    return -1;
        //}
        //std::ostringstream stream_force;
        //// 设置宽度为3，左侧填充0
        //stream_force << std::setw(3) << std::setfill('0') << iter;
        //std::string force_save_str = "../Results/Force_" + model_name + "_" + stream.str() + ".txt";
        //char* force_save = const_cast<char*>(force_save_str.c_str());
        //if (!write_file(force_save, f, np)) {
        //    std::cout << "Save FAILED!" << endl;
        //    return -1;
        //}
        //std::ostringstream stream_stress;
        //stream_stress << std::setw(3) << std::setfill('0') << iter;
        //std::string stress_save_str = "../Results/Stress" + model_name + "_" + stream_stress.str() + ".txt";
        //char* stress_save = const_cast<char*>(stress_save_str.c_str());
        //if (!cal_stress_save(stress_save, np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect)) {
        //    std::cout << "Stress calculate Faild!" << endl;
        //    return -1;
        //}
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
            double maxVal = 0.0;
            for (double val : jac_O_t) {
                maxVal = std::max(maxVal, std::abs(val));
            }
            if (maxVal < threshold) {
                cout << "The max gradient on t : " << maxVal << endl; loop_t = false; break;
            }
            while (true) {
                CppAD::AD<double> temp = t[0];
                t[0] = t[0] - alpha_t * (jac_O_t[0] / hess_O_t[0]);
                O = Calculate_O(t[0], sparse_Stress, n_tet);
                if (O< last_O_t && O > -1e5) {
                    t_x[0] = Value(t[0]);
                    alpha_t /= gamma;
                    last_O_t = O;
                    break;
                }
                else {
                    t[0] = temp;
                    alpha_t *= gamma;
                }
            }
            cout << "O:" << O << " The max of O_t : " << maxVal << "  t : " << t_x[0] << " dt:" << (jac_O_t[0] / hess_O_t[0]) << " |Alpha:" << alpha_t << endl;
        }
        //optimize s
        cout << "-------------------Optimize s: ----------Optimize s: ----------------- Optimize s --------------------------------------------------" << endl;
        cout << "The value of O after optimization on t : " << last_O_t << endl;

        vector<double> grad_s;
        Calculate_Stresses_AD(np, points, n_tri, triangles, n_tet, tetrahedras, n_sp, Index_of_sp, t_x[0], p_stress, grad_s, solver_LU);
        vec_grad = grad_s;
        if (info_f_g)
            extract_info(fileout, n_sp, Index_of_sp, n_tet, tetrahedras, n_tri, triangles, points, sparse_Stress, f, grad_s);
        else
            extract_info(fileout, n_sp, Index_of_sp, n_tet, tetrahedras, n_tri, triangles, points, sparse_Stress);


        loop = false;
        double maxVal = 0.0;
        for (double val : grad_s) {
            maxVal = std::max(maxVal, std::abs(val));
        }
        cout << "||max of O_s: " << maxVal << endl;
        auto end_time_iter = std::chrono::high_resolution_clock::now();
        auto duration_iter = std::chrono::duration_cast<std::chrono::seconds>(end_time_iter - start_iter).count();
        std::cout << "||Time taken by function: " << duration_iter << " seconds" << endl << endl << endl;
        return 0;
        //{
        //    std::ostringstream stream;
        //    // 设置宽度为3，左侧填充0
        //    stream << std::setw(3) << std::setfill('0') << iter;
        //    std::string grad_save_str = "../Results/Grad_AD_" + model_name + "_" + stream.str() + ".txt";
        //    char* grad_save = const_cast<char*>(grad_save_str.c_str());
        //    if (!write_file(grad_save, grad_s, np)) {
        //        cout << "Save FAILED!" << endl;
        //        return -1;
        //    }
        //}

        //sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU, intersect);
        //O = Calculate_O(t[0], sparse_Stress, n_tet);
        //last_O_s = O;

        //double maxVal = 0.0;
        //for (double val : grad_s) {
        //    maxVal = std::max(maxVal, std::abs(val));
        //}
        //if (maxVal < threshold_s) { loop = false; break; }
        //if (maxVal * alpha < edge_length * 0.05) {
        //    if (points_updated) {
        //        cout << "|-The maximum update step size is: " << maxVal * alpha << ", which is much smaller than the edge length, recalculate the pressure p.-|" << endl;
        //        break;
        //    }
        //}
        //while (maxVal * alpha > edge_length) {
        //    alpha *= gamma;
        //    loop_s_time++;
        //}
        //while (true) {
        //    intersect = 0;
        //    for (int i = 0; i < n; i++) {
        //        points[i] = x[i] - alpha * grad_s[i];
        //    }
        //    //After updating points,we need to update K,N,G->M,so we also need to rebuild a solver for new M.
        //    K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
        //    M_G = Build_M_G(np, points);        //计算G
        //    M = merge_matrix(K, M_G);
        //    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU_new;
        //    solver_LU_new.compute(M);
        //    sparse_Stress = Calculate_Stresses(np, points, n_sp, Index_of_sp, n_tri, triangles, n_tet, tetrahedras, p_stress, solver_LU_new, intersect);
        //    O = Calculate_O(t[0], sparse_Stress, n_tet);
        //    if (O< last_O_s && O > -1e5 && !intersect) {
        //        for (int i = 0; i < n; i++) {
        //            x[i] = points[i];
        //        }
        //        points_updated = true;
        //        alpha /= gamma;
        //        last_O_s = O;
        //        break;
        //    }
        //    else {
        //        if (intersect) {
        //            cout << "Self-intersection !" << endl;
        //        }
        //        cout << "The last O:" << last_O_s << "  O_new:" << O << endl;
        //        for (int i = 0; i < n; i++) {
        //            points[i] = x[i];
        //        }
        //        alpha *= gamma;
        //    }


    }

}
