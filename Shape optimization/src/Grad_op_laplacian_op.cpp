
#include <Grad.h>

typedef Eigen::Matrix< double, Dynamic, Dynamic> MatrixXd;
typedef Eigen::Matrix< CppAD::AD<double>, Dynamic, Dynamic> MatrixXd_AD;
typedef Eigen::Matrix<CppAD::AD<double>, 3, 1> Vector3d_AD;
double mu = 0.3;
double E = 1.0; //double E = 1e9;
double dt = 1e-6;
double O_mu = 1e-3;
double A_mu = 0.0;
double smooth_grad_norm = 0.0;
double Smooth_mu = 0.0;
int count_time = 0;
double edge_length = 0.05;

Eigen::SparseMatrix<double> merge_matrix(Eigen::SparseMatrix<double> K, Eigen::SparseMatrix<double> G) {
    // 行合并
    int k_rows = K.rows();
    Eigen::SparseMatrix<double> Merged(K.rows() + G.rows(), K.rows() + G.rows());
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve((K.nonZeros() + G.nonZeros()) * 2);
    // 添加K矩阵的非零元素
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            tripletList.emplace_back(it.row(), it.col(), it.value());
        }
    }
    // 添加G矩阵的非零元素
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
            tripletList.emplace_back(k_rows + it.row(), it.col(), it.value());
            tripletList.emplace_back(it.col(), k_rows + it.row(), it.value());
        }
    }
    Merged.setFromTriplets(tripletList.begin(), tripletList.end());
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

Eigen::SparseMatrix<double> Build_stiffness_Matrix(int nv, Eigen::MatrixXd vertices, int n_tet, Eigen::MatrixXi tetrahedras) {
    Eigen::SparseMatrix<double> sparse_K(nv * 3, nv * 3);
    vector<map<int, double>> sparse_k_map(3 * nv);
    VectorXd sparse_pri(3 * nv); sparse_pri.setZero();
    //sparse_K.reserve(nv * 3 +6* n_tet);
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
        t1 = tetrahedras(k, 0); t2 = tetrahedras(k, 1); t3 = tetrahedras(k, 2); t4 = tetrahedras(k, 3);
        p1_x = vertices(t1, 0); p1_y = vertices(t1, 1); p1_z = vertices(t1, 2);
        p2_x = vertices(t2, 0); p2_y = vertices(t2, 1); p2_z = vertices(t2, 2);
        p3_x = vertices(t3, 0); p3_y = vertices(t3, 1); p3_z = vertices(t3, 2);
        p4_x = vertices(t4, 0); p4_y = vertices(t4, 1); p4_z = vertices(t4, 2);
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
        for (int j = 0; j < 12; j++)
            for (int i = j; i < 12; i++) {
                if (index[i] == index[j])
                    sparse_pri[index[i]] += Ke(i, i);
                else if (index[i] > index[j])
                    sparse_k_map[index[i]][index[j]] += Ke(i, j);
                else
                    sparse_k_map[index[j]][index[i]] += Ke(j, i);
            }
    }
    int totalElements = 0;
    for (const auto& map : sparse_k_map) {  // 遍历vector中的每个map
        totalElements += map.size();     // 累加每个map的大小
    }
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(totalElements * 2 + 3 * nv);
    for (int i = 0; i < 3 * nv; i++)
        tripletList.emplace_back(i, i, sparse_pri[i]);
    // 添加K矩阵的非零元素
    for (int k = 0; k < sparse_k_map.size(); ++k) {
        for (const auto& elem : sparse_k_map[k]) {
            int colIndex = elem.first;       // 列索引
            double value = elem.second;      // 值
            tripletList.emplace_back(k, colIndex, value);
            tripletList.emplace_back(colIndex, k, value);
        }
    }
    sparse_K.setFromTriplets(tripletList.begin(), tripletList.end());
    sparse_K.makeCompressed();
    return sparse_K;
}

SparseMatrix<double> Build_N(int nv, Eigen::MatrixXd V, Eigen::MatrixXi F, VectorXi index_sp, double** normals) {
    int n_sp = V.rows();
    int n_tri = F.rows();
    Eigen::SparseMatrix<double> N(nv * 3, n_sp);
    Eigen::MatrixXd Normal;
    Normal.resize(n_sp, 3); Normal.setZero();
    for (int k = 0; k < n_tri; k++) {
        int pi_id = F(k, 0);
        int pj_id = F(k, 1);
        int pk_id = F(k, 2);
        Vector3d p_i = V.row(pi_id);
        Vector3d p_j = V.row(pj_id);
        Vector3d p_k = V.row(pk_id);
        Vector3d edge1 = p_j - p_i;
        Vector3d edge2 = p_k - p_i;
        Vector3d n_of_tirangle = edge1.cross(edge2);

        Normal.row(pi_id) += n_of_tirangle;
        Normal.row(pj_id) += n_of_tirangle;
        Normal.row(pk_id) += n_of_tirangle;
    }
    for (int k = 0; k < n_sp; k++) {
        int id = index_sp[k];
        Vector3d Normal_Of_P = Normal.row(k);
        //std::cout << Normal_Of_P << std::endl;
        Normal_Of_P.normalize();
        Normal.row(k) = Normal_Of_P;
        N.insert(3 * id, k) = -Normal_Of_P[0];
        N.insert(3 * id + 1, k) = -Normal_Of_P[1];
        N.insert(3 * id + 2, k) = -Normal_Of_P[2];
    }
    if (normals != nullptr) {
        double* Normals = (double*)calloc(nv * 3, sizeof(double));
        for (int k = 0; k < n_sp; k++) {
            int id = index_sp[k];
            Normals[3 * id] = -Normal(k, 0);
            Normals[3 * id + 1] = -Normal(k, 1);
            Normals[3 * id + 2] = -Normal(k, 2);
        }
        *normals = Normals;
    }
    return N;
}

SparseMatrix<double> Build_G(int nv, Eigen::MatrixXd vertices, SparseMatrix<double> N, int n_sp) {
    Eigen::SparseMatrix<double> G(6, nv * 3);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nv * 3);
    // 添加K矩阵的非零元素
    for (int k = 0; k < nv; k++) {
        tripletList.emplace_back(0, 3 * k, 1.0);
        tripletList.emplace_back(1, 3 * k + 1, 1.0);
        tripletList.emplace_back(2, 3 * k + 2, 1.0);

        tripletList.emplace_back(3, 3 * k + 1, -(vertices(k, 2))); tripletList.emplace_back(3, 3 * k + 2, (vertices(k, 1)));
        tripletList.emplace_back(4, 3 * k, vertices(k, 2)); tripletList.emplace_back(4, 3 * k + 2, -(vertices(k, 0)));
        tripletList.emplace_back(5, 3 * k, -(vertices(k, 1))); tripletList.emplace_back(5, 3 * k + 1, (vertices(k, 0)));
    }
    G.setFromTriplets(tripletList.begin(), tripletList.end());
    G.makeCompressed();
    Eigen::SparseMatrix<double> G_N(6, n_sp);
    G_N = G * N;
    return G_N;
}

SparseMatrix<double> Build_M_G(int nv, Eigen::MatrixXd vertices) {
    Eigen::SparseMatrix<double> G(6, nv * 3);
    G.reserve(nv * 9);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nv * 3);
    // 添加K矩阵的非零元素
    for (int k = 0; k < nv; k++) {
        tripletList.emplace_back(0, 3 * k, 1.0);
        tripletList.emplace_back(1, 3 * k + 1, 1.0);
        tripletList.emplace_back(2, 3 * k + 2, 1.0);

        tripletList.emplace_back(3, 3 * k + 1, -(vertices(k, 2))); tripletList.emplace_back(3, 3 * k + 2, (vertices(k, 1)));
        tripletList.emplace_back(4, 3 * k, vertices(k, 2)); tripletList.emplace_back(4, 3 * k + 2, -(vertices(k, 0)));
        tripletList.emplace_back(5, 3 * k, -(vertices(k, 1))); tripletList.emplace_back(5, 3 * k + 1, (vertices(k, 0)));
    }
    G.setFromTriplets(tripletList.begin(), tripletList.end());
    G.makeCompressed();
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

Eigen::SparseMatrix<double> Calculate_Stresses(Eigen::MatrixXi tetrahedras, Eigen::MatrixXd vertices, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXi new_2_Old,
    VectorXd p, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver, int& intersect) {
    int nv = vertices.rows();
    int n_sp = V.rows();
    int n_tet = tetrahedras.rows();
    int n = 3 * nv;

    Eigen::SparseMatrix<double> N = Build_N(nv, V, F, new_2_Old, nullptr);
    Eigen::SparseMatrix<double> G = Build_G(nv, vertices, N, n_sp);
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<double> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    //p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

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
        t1 = tetrahedras(k, 0); t2 = tetrahedras(k, 1); t3 = tetrahedras(k, 2); t4 = tetrahedras(k, 3);
        p1_x = vertices(t1, 0); p1_y = vertices(t1, 1); p1_z = vertices(t1, 2);
        p2_x = vertices(t2, 0); p2_y = vertices(t2, 1); p2_z = vertices(t2, 2);
        p3_x = vertices(t3, 0); p3_y = vertices(t3, 1); p3_z = vertices(t3, 2);
        p4_x = vertices(t4, 0); p4_y = vertices(t4, 1); p4_z = vertices(t4, 2);
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

Eigen::SparseMatrix<double> Calculate_Stresses(Eigen::MatrixXi tetrahedras, Eigen::MatrixXd vertices, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXi new_2_Old,
    VectorXd p, Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper>& solver, int& intersect, VectorXd& u_back) {
    int nv = vertices.rows();
    int n_sp = V.rows();
    int n_tet = tetrahedras.rows();
    int n = 3 * nv;

    Eigen::SparseMatrix<double> N = Build_N(nv, V, F, new_2_Old, nullptr);
    Eigen::SparseMatrix<double> G = Build_G(nv, vertices, N, n_sp);
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<double> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    //p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    Eigen::VectorXd u(n + 6);
    Eigen::VectorXd f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver.solveWithGuess(f, u_back);
    // 检查解是否成功
    std::cout << "The norm of u is : " << u.norm() << std::endl;
    u_back = u;
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
        t1 = tetrahedras(k, 0); t2 = tetrahedras(k, 1); t3 = tetrahedras(k, 2); t4 = tetrahedras(k, 3);
        p1_x = vertices(t1, 0); p1_y = vertices(t1, 1); p1_z = vertices(t1, 2);
        p2_x = vertices(t2, 0); p2_y = vertices(t2, 1); p2_z = vertices(t2, 2);
        p3_x = vertices(t3, 0); p3_y = vertices(t3, 1); p3_z = vertices(t3, 2);
        p4_x = vertices(t4, 0); p4_y = vertices(t4, 1); p4_z = vertices(t4, 2);
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

double Calculate_O(int& notPositiveSemiDefinite, double t, Eigen::SparseMatrix<double> stress, int n_tet) {
    double O = 0.0;
    for (int i = 0; i < n_tet; i++) {
        double det1, det2;

        Eigen::Matrix<double, 3, 3> Sigma, tI;
        Sigma << stress.coeff(6 * i, 0), stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 5, 0),
            stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 1, 0), stress.coeff(6 * i + 4, 0),
            stress.coeff(6 * i + 5, 0), stress.coeff(6 * i + 4, 0), stress.coeff(6 * i + 2, 0);
        tI << t * 1.0, 0.0, 0.0,
            0.0, t * 1.0, 0.0,
            0.0, 0.0, t * 1.0;
        Eigen::MatrixXd mat1 = tI - Sigma;
        Eigen::MatrixXd mat2 = tI + Sigma;
        det1 = mat1.determinant();
        det2 = mat2.determinant();
        if (det1 < 1e-3 || det2 < 1e-3)
        {
            notPositiveSemiDefinite = 1;
            //cout << "Invalid Guess s,t" << endl;
            return -1e6;
        }
        if (!isPositiveSemiDefinite(mat1) || !isPositiveSemiDefinite(mat2)) {
            notPositiveSemiDefinite = 1;
            cout << "NOT Positive SemiDefinite." << endl;
            return -1e6;
        }
        O += -O_mu * (log(det1) + log(det2));
    }
    O += t;
    notPositiveSemiDefinite = 0;
    //cout << "Double || The value of O : " << O << "  t :" << t << endl;
    return O;
}


double Calculate_O_dt(double t, Eigen::SparseMatrix<double> stress, int n_tet, double& DO_DT, double& DO_DT2) {
    double dOdt = 0.0;
    double dOdt2 = 0.0;
    for (int i = 0; i < n_tet; i++) {
        double det1, det2;
        Eigen::Matrix<double, 3, 3> Sigma, tI;
        Sigma << stress.coeff(6 * i, 0), stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 5, 0),
            stress.coeff(6 * i + 3, 0), stress.coeff(6 * i + 1, 0), stress.coeff(6 * i + 4, 0),
            stress.coeff(6 * i + 5, 0), stress.coeff(6 * i + 4, 0), stress.coeff(6 * i + 2, 0);
        tI << t * 1.0, 0.0, 0.0,
            0.0, t * 1.0, 0.0,
            0.0, 0.0, t * 1.0;
        Eigen::MatrixXd mat1 = tI - Sigma;
        Eigen::MatrixXd mat2 = tI + Sigma;
        det1 = mat1.determinant();
        det2 = mat2.determinant();
        if (det1 < 1e-3 || det2 < 1e-3)
        {
            //cout << "Invalid Guess s,t" << endl;
            return -1e6;
        }
        if (!isPositiveSemiDefinite(mat1) || !isPositiveSemiDefinite(mat2)) {
            cout << "NOT Positive SemiDefinite." << endl;
            return -1e6;
        }
        Eigen::MatrixXd mat1_inv = mat1.inverse();
        Eigen::MatrixXd mat2_inv = mat2.inverse();
        Eigen::MatrixXd mat1_inv_2 = mat1_inv * mat1_inv;
        Eigen::MatrixXd mat2_inv_2 = mat2_inv * mat2_inv;
        mat1_inv = mat1_inv + mat2_inv;
        mat1_inv_2 = mat1_inv_2 + mat2_inv_2;
        dOdt += -O_mu * (mat1_inv(0, 0) + mat1_inv(1, 1) + mat1_inv(2, 2));
        dOdt2 += O_mu * (mat1_inv_2(0, 0) + mat1_inv_2(1, 1) + mat1_inv_2(2, 2));
    }
    dOdt += 1.0;
    DO_DT = dOdt;
    DO_DT2 = dOdt2;
    return 0.0;
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

Eigen::SparseMatrix<CppAD::AD<double>> Build_stiffness_Matrix(int nv, const vector<CppAD::AD<double>>& vertices, int n_tet, Eigen::MatrixXi tetrahedras) {
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
        t1 = tetrahedras(k, 0); t2 = tetrahedras(k, 1); t3 = tetrahedras(k, 2); t4 = tetrahedras(k, 3);
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
                sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
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

SparseMatrix<CppAD::AD<double>> Build_N(int nv, const vector<CppAD::AD<double>>& points, Eigen::MatrixXi F, VectorXi index_sp) {
    int n_sp = index_sp.size();
    int n_tri = F.rows();
    Eigen::SparseMatrix<CppAD::AD<double>> N(nv * 3, n_sp);
    vector<CppAD::AD<double>> Normals(n_sp * 3, CppAD::AD<double>(0));
    for (int k = 0; k < n_tri; k++) {
        int pi_id = F(k, 0);
        int pj_id = F(k, 1);
        int pk_id = F(k, 2);
        Vector3d_AD p_i(points[3 * pi_id], points[3 * pi_id + 1], points[3 * pi_id + 2]);
        Vector3d_AD p_j(points[3 * pj_id], points[3 * pj_id + 1], points[3 * pj_id + 2]);
        Vector3d_AD p_k(points[3 * pk_id], points[3 * pk_id + 1], points[3 * pk_id + 2]);
        Vector3d_AD edge1 = p_j - p_i;
        Vector3d_AD edge2 = p_k - p_i;
        Vector3d_AD n_of_tirangle = edge1.cross(edge2);

        Normals[3 * pi_id] += n_of_tirangle[0]; Normals[3 * pi_id + 1] += n_of_tirangle[1]; Normals[3 * pi_id + 2] += n_of_tirangle[2];
        Normals[3 * pj_id] += n_of_tirangle[0]; Normals[3 * pj_id + 1] += n_of_tirangle[1]; Normals[3 * pj_id + 2] += n_of_tirangle[2];
        Normals[3 * pk_id] += n_of_tirangle[0]; Normals[3 * pk_id + 1] += n_of_tirangle[1]; Normals[3 * pk_id + 2] += n_of_tirangle[2];
    }
    for (int k = 0; k < n_sp; k++) {
        int id = index_sp[k];
        Vector3d_AD Normal_Of_P(Normals[3 * k], Normals[3 * k + 1], Normals[3 * k + 2]);
        Normal_Of_P.normalize();
        N.insert(3 * id, k) = -Normal_Of_P[0];
        N.insert(3 * id + 1, k) = -Normal_Of_P[1];
        N.insert(3 * id + 2, k) = -Normal_Of_P[2];
    }
    return N;
}

int Stiffness_grad_id(int nv, const vector<CppAD::AD<double>> vertices, int n_tet, Eigen::MatrixXi tetrahedras, int id,
    vector<int>& row_indices_output, vector<int>& col_indices_output, vector<double>& grad_K) {
    // 遍历 K 的非零元素
    auto START_TIME2 = std::chrono::high_resolution_clock::now();
    vector<int>row_indices;
    vector<int>col_indices;
    vector<CppAD::AD<double>>K_AD_Vector;
    vector<CppAD::AD<double>> verteces_id(3);
    vector<double>verteces_id_s(3);

    for (int i = 0; i < 3; i++) {
        verteces_id[i] = vertices[3 * id + i];
        verteces_id_s[i] = Value(vertices[3 * id + i]);
    }
    CppAD::Independent(verteces_id);
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
        int t_id;
        CppAD::AD<double> a1, a2, a3, b1, b2, b3, c1, c2, c3;
        Matrix<CppAD::AD<double>, 4, 4> Lambda;
        t1 = tetrahedras(k, 0); t2 = tetrahedras(k, 1); t3 = tetrahedras(k, 2); t4 = tetrahedras(k, 3);
        if ((t1 != id) && (t2 != id) && (t3 != id) && (t4 != id))
            continue;
        if (t1 == id) {
            p1_x = verteces_id[0]; p1_y = verteces_id[1]; p1_z = verteces_id[2];
        }
        else {
            p1_x = vertices[t1 * 3]; p1_y = vertices[t1 * 3 + 1]; p1_z = vertices[t1 * 3 + 2];
        }
        if (t2 == id) {
            p2_x = verteces_id[0]; p2_y = verteces_id[1]; p2_z = verteces_id[2];
        }
        else {
            p2_x = vertices[t2 * 3]; p2_y = vertices[t2 * 3 + 1]; p2_z = vertices[t2 * 3 + 2];
        }
        if (t3 == id) {
            p3_x = verteces_id[0]; p3_y = verteces_id[1]; p3_z = verteces_id[2];
        }
        else {
            p3_x = vertices[t3 * 3]; p3_y = vertices[t3 * 3 + 1]; p3_z = vertices[t3 * 3 + 2];
        }
        if (t4 == id) {
            p4_x = verteces_id[0]; p4_y = verteces_id[1]; p4_z = verteces_id[2];
        }
        else {
            p4_x = vertices[t4 * 3]; p4_y = vertices[t4 * 3 + 1]; p4_z = vertices[t4 * 3 + 2];
        }
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
        for (int j = 0; j < 12; j++)
            for (int i = j; i < 12; i++) {
                //sparse_K.coeffRef(0, index[i]) += 1.0;
                if (index[i] >= index[j])
                    sparse_K.coeffRef(index[i], index[j]) += Ke(i, j);
                else
                    sparse_K.coeffRef(index[j], index[i]) += Ke(j, i);
            }
    }
    sparse_K.makeCompressed();
    // 遍历 K 的非零元素
    row_indices.clear(); col_indices.clear(); K_AD_Vector.clear();
    for (int k = 0; k < sparse_K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<CppAD::AD<double>>::InnerIterator it(sparse_K, k); it; ++it) {
            row_indices.push_back(it.row());// it.row()       // 非零元素的行索引
            col_indices.push_back(it.col());// it.col()       // 非零元素的列索引
            K_AD_Vector.push_back(CppAD::AD<double>(it.value()));// it.value()     // 非零元素的值
        }
    }
    CppAD::ADFun<double> func(verteces_id, K_AD_Vector);    // 创建 ADFun 对象
    vector<double>jac_K = func.Jacobian(verteces_id_s);
    row_indices_output = row_indices;
    col_indices_output = col_indices;
    grad_K = jac_K;
    return 1;
}

void grad_M_G_id(int id, vector<int>& row_indices_output, vector<int>& col_indices_output, vector<double>& grad_G) {
    vector<int> row_indices;
    vector<int> col_indices;
    vector<double> jac_G;
    row_indices.push_back(4); col_indices.push_back(3 * id + 2); jac_G.push_back(-1);
    row_indices.push_back(5); col_indices.push_back(3 * id + 1); jac_G.push_back(1);

    row_indices.push_back(3); col_indices.push_back(3 * id + 2); jac_G.push_back(1);
    row_indices.push_back(5); col_indices.push_back(3 * id); jac_G.push_back(-1);

    row_indices.push_back(3); col_indices.push_back(3 * id + 1); jac_G.push_back(-1);
    row_indices.push_back(4); col_indices.push_back(3 * id); jac_G.push_back(1);

    row_indices_output = row_indices;
    col_indices_output = col_indices;
    grad_G = jac_G;
}



VectorXd Calculate_sigma_A(int nv, Eigen::MatrixXd vertices, int n_tet, Eigen::MatrixXi tetrahedras, VectorXd u, double t) {
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
        t1 = tetrahedras(k, 0); t2 = tetrahedras(k, 1); t3 = tetrahedras(k, 2); t4 = tetrahedras(k, 3);
        p1_x = vertices(t1, 0); p1_y = vertices(t1, 1); p1_z = vertices(t1, 2);
        p2_x = vertices(t2, 0); p2_y = vertices(t2, 1); p2_z = vertices(t2, 2);
        p3_x = vertices(t3, 0); p3_y = vertices(t3, 1); p3_z = vertices(t3, 2);
        p4_x = vertices(t4, 0); p4_y = vertices(t4, 1); p4_z = vertices(t4, 2);
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


double Calculate_Stresses_AD(Eigen::MatrixXi tetrahedras, Eigen::MatrixXd points, Eigen::MatrixXd V, Eigen::MatrixXi F, unordered_map<int, int> Old_2_new, Eigen::VectorXi index_sp, double t,
    VectorXd p, vector<double>& grad_s, Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper>& solver, Eigen::VectorXd& u_back) {
    int nv = points.rows();
    int n_sp = V.rows();
    int n_tet = tetrahedras.rows();
    int n_tri = F.rows();
    int n = 3 * nv;
    auto start_time_AD = std::chrono::high_resolution_clock::now();
    double* Normal;

    //1-计算位移 u 
    Eigen::SparseMatrix<double> K = Build_stiffness_Matrix(nv, points, n_tet, tetrahedras);
    Eigen::SparseMatrix<double> N = Build_N(nv, V, F, index_sp, &Normal);
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
    //p = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;

    Eigen::VectorXd u(n + 6);
    Eigen::VectorXd f(n);
    //计算位移displacement
    f = N * p;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    //u = solver.solve(f);
    u = solver.solveWithGuess(f, u_back);
    u_back = u;
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
    std::vector<double> jac_V;
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
    for (int i = 0; i < nv; i++) {
        for (int j = 0; j < 3; j++) {
            vertices[3 * i + j] = CppAD::AD<double>(points(i, j));
            vertices_x[3 * i + j] = points(i, j);
        }
    }
    cout << "number of triangles: " << n_tri << endl;
    cout << "The number of vertices: " << nv << endl;
    cout << "The number of vertices on the surface: " << n_sp << endl;

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
        t1 = tetrahedras(k, 0); t2 = tetrahedras(k, 1); t3 = tetrahedras(k, 2); t4 = tetrahedras(k, 3);
        int index[] = { 3 * t1, 3 * t1 + 1, 3 * t1 + 2,
                    3 * t2, 3 * t2 + 1, 3 * t2 + 2,
                    3 * t3, 3 * t3 + 1, 3 * t3 + 2,
                    3 * t4, 3 * t4 + 1, 3 * t4 + 2 };
        int t_index[] = { t1,t2,t3,t4 };
        for (int index_i = 0; index_i < 4; index_i++) {
            for (int index_j = 0; index_j < 3; index_j++) {
                s[3 * index_i + index_j] = CppAD::AD<double>(points(t_index[index_i], index_j));
                s_x[3 * index_i + index_j] = points(t_index[index_i], index_j);
            }
        }
        for (int index_i = 0; index_i < 12; index_i++)
            ue[index_i] = u[index[index_i]];
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
            int id = (int)(index[k_i] / 3);
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
    double grad_norm = 0.0;

    vector<CppAD::AD<double>> V_ad(n_sp * 3);
    vector<double> V_ad_x(n_sp * 3);
    for (int i = 0; i < n_sp; i++) {
        for (int j = 0; j < 3; j++) {
            V_ad[3 * i + j] = CppAD::AD<double>(V(i, j));
            V_ad_x[3 * i + j] = V(i, j);
        }
    }
    for (int s_i = 0; s_i < n_sp; s_i++) {
        int s_i_old = index_sp[s_i];
        vector<int> row_indices;
        vector<int> col_indices;
        vector<double> grad_K;
        vector<int> row_indices_G;
        vector<int> col_indices_G;
        vector<double> grad_G;
        Stiffness_grad_id(nv, vertices, n_tet, tetrahedras, s_i_old, row_indices, col_indices, grad_K);
        grad_M_G_id(s_i_old, row_indices_G, col_indices_G, grad_G);
        Vector3d normal; normal << Normal[3 * s_i_old], Normal[3 * s_i_old + 1], Normal[3 * s_i_old + 2];
        for (int s_i_cord = 0; s_i_cord < 3; s_i_cord++) {
            int s_i_old_i = 3 * s_i_old + s_i_cord;
            int s_i_i = 3 * s_i + s_i_cord;
            vector<CppAD::AD<double>> s_k(1);
            vector<double> s_k_x(1);
            s_k[0] = CppAD::AD<double>(points(s_i_old, s_i_cord));
            s_k_x[0] = points(s_i_old, s_i_cord);
            CppAD::Independent(s_k);
            V_ad[s_i_i] = s_k[0];
            N_AD = Build_N(nv, V_ad, F, index_sp);
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
            int non_zero = row_indices.size();
            for (int it = 0; it < non_zero; it++)
                if (grad_K[it * 3 + s_i_cord] != 0.0) {
                    grad_M.insert(row_indices[it], col_indices[it]) = grad_K[it * 3 + s_i_cord];
                    if (row_indices[it] != col_indices[it])
                        grad_M.insert(col_indices[it], row_indices[it]) = grad_K[it * 3 + s_i_cord];
                }
            grad_M.makeCompressed();
            for (int it = 0; it < 2; it++) {
                //   cout << "id :" << row_indices_G[2 * s_i_cord + it] + n << " , " << col_indices_G[2 * s_i_cord + it] << " : " << grad_G[2 * s_i_cord + it] << endl;
                grad_M.insert(row_indices_G[2 * s_i_cord + it] + n, col_indices_G[2 * s_i_cord + it]) = grad_G[2 * s_i_cord + it];
                grad_M.insert(col_indices_G[2 * s_i_cord + it], row_indices_G[2 * s_i_cord + it] + n) = grad_G[2 * s_i_cord + it];
            }
            grad_M.makeCompressed();

            double part2 = O_A.transpose() * grad_f_s;
            double part3 = O_A.transpose() * grad_M * O_B;
            grad_O_s[s_i_old_i] += O_mu * (part2 - part3);
            //grad_O_s[s_i_old_i] += -jac_V[s_i_i];
            //grad_O_s[s_i_old_i] += O_mu * (part2);
            //vertices[s_i_old_i] = CppAD::AD<double>(points(s_i_old, s_i_cord));
        }
        Vector3d grad_p; grad_p << grad_O_s[3 * s_i_old], grad_O_s[3 * s_i_old + 1], grad_O_s[3 * s_i_old + 2];
        grad_p = (grad_p.dot(normal) / normal.dot(normal)) * normal;
        grad_O_s[3 * s_i_old] = grad_p[0]; grad_O_s[3 * s_i_old + 1] = grad_p[1]; grad_O_s[3 * s_i_old + 2] = grad_p[2];
        grad_norm += grad_p.norm();
    }
    grad_s = grad_O_s;
    auto end_time_stress = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time_stress - start_time_AD).count();
    std::cout << "|| The cost of stress grad : " << duration << " seconds" << endl << endl << endl;
    cout << "The norm of grad : " << grad_norm << endl;
    return grad_norm;
}

double Stresses_AD_Filter(Eigen::MatrixXi tetrahedras, SparseMatrix<double> sparse_stress,
    unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old, Eigen::VectorXd& grad) {
    int n_sp = new_2_Old.size();
    int n_tet = tetrahedras.rows();
    vector<int> count(n_sp);
    vector<double> stress(n_sp);
    for (int i = 0; i < n_tet; i++) {
        for (int t = 0; t < 4; t++) {
            int id = tetrahedras(i, t);
            if (Old_2_new.find(id) != Old_2_new.end()) {
                stress[Old_2_new[id]] += max(max(abs(sparse_stress.coeffRef(6 * i, 0)), abs(sparse_stress.coeffRef(6 * i + 1, 0))), abs(sparse_stress.coeffRef(6 * i + 2, 0)));
                count[Old_2_new[id]]++;
            }
        }
    }
    double max_stress = 0.0;
    double grad_norm = 0.0;
    for (int t = 0; t < n_sp; t++) {
        stress[t] = stress[t] / count[t];
        if (stress[t] > max_stress)
            max_stress = stress[t];
    }
    //cout << "max stress : " << max_stress << endl;
    for (int t = 0; t < n_sp; t++) {
        int id = new_2_Old[t];
        if (stress[t] < max_stress * 0.1) {
            grad.segment<3>(3 * id) *= 0.0;
        }
        /* Vector3d grad_t(grad_s[3 * id], grad_s[3 * id + 1], grad_s[3 * id + 2]);*/
        grad_norm += grad.segment<3>(3 * id).norm();
    }
    cout << "-- | The norm of grad after filtering : " << grad_norm << " | --" << endl;

}

double calculate_t(SparseMatrix<double> sparse_Stress, int n_tet, double& O_stress) {
    bool loop_t = true, loop_s = true;
    double t = 1e4;
    double threshold = 1e-15, gamma = 0.5;
    int notPositiveSemiDefinite = 0;
    double O = Calculate_O(notPositiveSemiDefinite, t, sparse_Stress, n_tet);
    while (notPositiveSemiDefinite) {
        std::cout << "Invalid Initial!!" << endl;
        t *= 2.0;
        O = Calculate_O(notPositiveSemiDefinite, t, sparse_Stress, n_tet);
        std::cout << O << endl;
    }
    double last_O = O, last_O_t = O, last_O_s = O;
    int loop_t_time = 0, loop_s_time = 0;
    double alpha = 0.001;
    auto start_iter = std::chrono::high_resolution_clock::now();
    //optimize t 
    loop_t = true;
    cout << "--------------------Optimize t: ---------------------------------------------------------------------------------" << endl;
    double alpha_t = 1.0;
    while (loop_t) {
        O = Calculate_O(notPositiveSemiDefinite, t, sparse_Stress, n_tet);
        last_O_t = O;
        double dodt, dodt2;
        Calculate_O_dt(t, sparse_Stress, n_tet, dodt, dodt2);
        double new_ton_delta = dodt / dodt2;
        double maxVal = abs(dodt);
        if (maxVal < threshold) {
            cout << "The max gradient on t : " << maxVal << endl; loop_t = false; break;
        }
        double t_temp = t;
        while (true) {
            if (abs(new_ton_delta * alpha_t) < 1e-15) {
                cout << "|-The maximum update step size is: " << abs(new_ton_delta * alpha_t) << ", which is much small.-|" << endl;
                loop_t = false;
                break;
            }
            t = t - alpha_t * (new_ton_delta);
            O = Calculate_O(notPositiveSemiDefinite, t, sparse_Stress, n_tet);
            if (O < last_O_t && !notPositiveSemiDefinite) {
                t_temp = t;
                alpha_t /= gamma;
                last_O_t = O;
                break;
            }
            else {
                t = t_temp;
                alpha_t *= gamma;
            }
        }
        //cout << "O:" << O << " The max of O_t : " << maxVal << "  t : " << t << " new_ton_delta:" << new_ton_delta << " |Alpha:" << alpha_t << endl;
    }
    //optimize s
    cout << "-------------------Optimize s: ----------Optimize s: ----------------- Optimize s --------------------------------------------------" << endl;
    cout << "The value of O after optimization on t : " << last_O_t << endl;
    cout << "The value of O t : " << t << endl;
    O_stress = O;
    return t;
}

int update_pressure(int n_sp, int n, SparseMatrix<double> K, SparseMatrix<double> N, SparseMatrix<double> G, Eigen::VectorXd Area, Eigen::VectorXi new_2_Old,
    VectorXd& pressure, Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper>& solver, int& successful, Eigen::VectorXd& u_back) {
    //计算K的最大特征值.
    double MaxEigenvalue = calculate_max_eigenvalue(K);
    Eigen::SparseMatrix<double> I(n_sp, n_sp);
    for (int i = 0; i < n_sp; i++) {
        I.insert(i, i) = 1.0;
    }
    Eigen::MatrixXd GGT = G * G.transpose();
    Eigen::MatrixXd GGT_inverse = GGT.inverse();
    //cout << G << endl;
    //cout << GGT_inverse << endl;
    Eigen::SparseMatrix<double> sparse_GGT_inverse = GGT_inverse.sparseView();
    Eigen::SparseMatrix<double> GEP_matrix_sparse = ((I - G.transpose() * (sparse_GGT_inverse)*G).transpose() *
        N.transpose() * K * N *
        (I - G.transpose() * (sparse_GGT_inverse)*G) +
        G.transpose() * MaxEigenvalue * G
        );

    Eigen::SparseMatrix<double> A(n, n);
    for (int i = 0; i < n_sp; i++) {
        int id = new_2_Old[i];
        A.insert(3 * id, 3 * id) = Area[i];
        A.insert(3 * id + 1, 3 * id + 1) = Area[i];
        A.insert(3 * id + 2, 3 * id + 2) = Area[i];
    }
    Eigen::SparseMatrix<double> P = I - G.transpose() * (sparse_GGT_inverse)*G;
    Eigen::SparseMatrix<double> matrix_B = N.transpose() * A * N;
    SparseSymMatProd<double> opA(GEP_matrix_sparse);
    SparseCholesky<double>  Bop(matrix_B);

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    SymGEigsSolver<SparseSymMatProd<double>, SparseCholesky<double>, GEigsMode::Cholesky>
        geigs(opA, Bop, 3, 24);
    // Initialize and compute

    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs_GEP;
    try {
        geigs.init();

        int nconv = geigs.compute(SortRule::SmallestAlge);
        //std::cout << nconv << std::endl;
    }
    catch (const std::runtime_error& e) {
        std::cerr << "TridiagEigen: eigen decomposition failed." << std::endl;
        successful = 0;
        return 0;
    }

    if (geigs.info() == CompInfo::Successful)
    {
        evalues = geigs.eigenvalues();
        evecs_GEP = geigs.eigenvectors();
    }
    else {
        std::cerr << "Decomposition failed!" << std::endl;
        successful = 0;
        return 0;
        //throw std::runtime_error("GEP cannot be computed !");
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
        u = solver.solve(f);
        double energy = u.transpose() * f;
        if (energy > max_energy) {
            max_energy = energy;
            max_energy_p_i = p_i;
        }
    }
    p = evecs_GEP.col(max_energy_p_i);
    VectorXd p_stress(n_sp);
    p_stress = evecs_GEP.col(max_energy_p_i);
    //y = (I - G.transpose() * (sparse_GGT_inverse)*G) * p;
    //VectorXd Gp = G * y;
    //Gp = G * p;
    //p = y;
    f = N * p;
    std::cout << "模长(f)：" << f.norm() << endl;
    f.conservativeResize(n + 6);
    f[n] = 0.0; f[n + 1] = 0.0; f[n + 2] = 0.0; f[n + 3] = 0.0; f[n + 4] = 0.0; f[n + 5] = 0.0;
    u = solver.solve(f);
    pressure = p;
    successful = 1;
    u_back = u;
    cout << "Energy E : " << u.transpose() * f << endl;
    return 1;
}


int calculate_pressure(Eigen::MatrixXd points, Eigen::MatrixXi tetrahedras, Eigen::MatrixXd V, Eigen::MatrixXi F,
    Eigen::VectorXi new_2_Old, int& successful, Eigen::VectorXd& p, Eigen::VectorXd& u) {
    int iter = 0;
    int last_smooth_iter = 0;
    double grad_start_norm = 0.0;

    Eigen::SparseMatrix<double> N;
    Eigen::SparseMatrix<double> G;
    Eigen::SparseMatrix<double> M_G;
    Eigen::SparseMatrix<double> K;
    Eigen::VectorXd Area;
    int np = points.rows();
    int n_tet = tetrahedras.rows();
    int n = 3 * np;
    int n_sp = V.rows();
    Calculate_Area(V, F, Area);
    cout << "n_sp:" << n_sp << endl;
    double* new_points;
    int* new_triangles;

    //转变为double的SparseMatrix
    K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    N = Build_N(np, V, F, new_2_Old, nullptr);
    G = Build_G(np, points, N, n_sp);        //计算G
    M_G = Build_M_G(np, points);        //计算G
    Eigen::SparseMatrix<double> M = merge_matrix(K, M_G);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper>cg;
    cg.compute(M);
    cg.setTolerance(1e-16);
    //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU;
    //solver_LU.compute(M);
    if (cg.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The inverse of K cannot be computed !");
    }

    Eigen::VectorXd p_stress;
    update_pressure(n_sp, n, K, N, G, Area, new_2_Old, p_stress, cg, successful, u);
    if (!successful)
        return 0;
    successful = 1;
    p = p_stress;
    return 1;
}

double update_t(Eigen::MatrixXd points, Eigen::MatrixXi tetrahedras, Eigen::MatrixXd V, Eigen::MatrixXi F,
    Eigen::VectorXi new_2_Old, double& Max_stress, int& successful, Eigen::VectorXd p, Eigen::VectorXd& u_back) {
    int iter = 0;
    int last_smooth_iter = 0;
    double grad_start_norm = 0.0;

    Eigen::SparseMatrix<double> N;
    Eigen::SparseMatrix<double> G;
    Eigen::SparseMatrix<double> M_G;
    Eigen::SparseMatrix<double> K;
    int np = points.rows();
    int n_tet = tetrahedras.rows();
    int n = 3 * np;
    int n_sp = V.rows();
    double* new_points;
    int* new_triangles;
    //转变为double的SparseMatrix
    K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    N = Build_N(np, V, F, new_2_Old, nullptr);
    G = Build_G(np, points, N, n_sp);        //计算G
    M_G = Build_M_G(np, points);        //计算G
    Eigen::SparseMatrix<double> M = merge_matrix(K, M_G);

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper>cg;
    //ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg
    cg.compute(M);
    cg.setTolerance(1e-16);
    if (cg.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The inverse of K cannot be computed !");
    }
    //calculate stress grad
    int intersect = 0;
    double smooth_E;
    Eigen::SparseMatrix<double> sparse_Stress = Calculate_Stresses(tetrahedras, points, V, F, new_2_Old, p, cg, intersect, u_back);
    std::cout << sparse_Stress.coeffRef(0, 0) << std::endl;
    if (intersect) {
        successful = 0;
        cout << "There is a intersection!" << endl;
        return 0.0;
    }
    double max_tet_stress = 0.0;
    for (int i = 0; i < n_tet; i++) {
        for (int pos = 0; pos < 3; pos++) {
            max_tet_stress = max(max_tet_stress, abs(sparse_Stress.coeffRef(6 * i + pos, 0)));
        }
    }
    Max_stress = max_tet_stress;
    successful = 1;
    double O;
    double t = calculate_t(sparse_Stress, n_tet, O);
    return t;
}


double stress_grad(Eigen::MatrixXd points, Eigen::MatrixXi triangles, Eigen::MatrixXi tetrahedras, Eigen::MatrixXd V, Eigen::MatrixXi F,
    unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old, Eigen::VectorXd& grad, double& Max_stress, int& successful, int update_p,
    Eigen::VectorXd p, Eigen::VectorXd& u) {
    int iter = 0;
    int last_smooth_iter = 0;
    double grad_start_norm = 0.0;

    Eigen::SparseMatrix<double> N;
    Eigen::SparseMatrix<double> G;
    Eigen::SparseMatrix<double> M_G;
    Eigen::SparseMatrix<double> K;
    Eigen::VectorXd Area;
    int np = points.rows();
    int n_tet = tetrahedras.rows();
    int n = 3 * np;
    int n_sp = V.rows();
    Calculate_Area(V, F, Area);
    cout << "n_sp:" << n_sp << endl;
    double* new_points;
    int* new_triangles;

    //转变为double的SparseMatrix
    K = Build_stiffness_Matrix(np, points, n_tet, tetrahedras);
    N = Build_N(np, V, F, new_2_Old, nullptr);
    G = Build_G(np, points, N, n_sp);        //计算G
    M_G = Build_M_G(np, points);        //计算G
    Eigen::SparseMatrix<double> M = merge_matrix(K, M_G);
    //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_LU;
    //solver_LU.compute(M);
    //if (solver_LU.info() != Eigen::Success) {
    //    std::cerr << "Decomposition failed!" << std::endl;
    //    throw std::runtime_error("The inverse of K cannot be computed !");
    //}
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper>cg;
    //ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg
    cg.compute(M);
    cg.setTolerance(1e-16);
    if (cg.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("The inverse of K cannot be computed !");
    }

    Eigen::VectorXd p_stress;
    p_stress = p;
    //calculate stress grad
    int intersect = 0;
    double smooth_E;
    Eigen::SparseMatrix<double> sparse_Stress = Calculate_Stresses(tetrahedras, points, V, F, new_2_Old, p_stress, cg, intersect, u);
    std::cout << sparse_Stress.coeffRef(0, 0) << std::endl;
    if (intersect) {
        successful = 0;
        cout << "There is a intersection!" << endl;
        return 0.0;
    }
    double max_tet_stress = 0.0;
    for (int i = 0; i < n_tet; i++) {
        for (int pos = 0; pos < 3; pos++) {
            max_tet_stress = max(max_tet_stress, abs(sparse_Stress.coeffRef(6 * i + pos, 0)));
        }
    }
    Max_stress = max_tet_stress;
    smooth_E = 0.0;
    std::cout << "max_tet_stress : " << max_tet_stress << std::endl;
    double O;
    double t = calculate_t(sparse_Stress, n_tet, O);

    vector<double> grad_s;

    double grad_new_norm = Calculate_Stresses_AD(tetrahedras, points, V, F, Old_2_new, new_2_Old, t, p_stress, grad_s, cg, u);

    grad.resize(n); grad.setZero();
    for (int i = 0; i < n; i++)
        grad[i] = grad_s[i];

    Stresses_AD_Filter(tetrahedras, sparse_Stress, Old_2_new, new_2_Old, grad);
    successful = 1;
    //for (int i = 0; i < n; i++)
    //    grad_s_old = grad_s;
    //Stresses_AD_Filter(np, points, n_tet, tetrahedras, n_sp, Index_of_sp, sparse_Stress, grad_s);
    //std::ostringstream stream_index;
    // 设置宽度为3，左侧填充0 
    //stream_index << std::setw(3) << std::setfill('0') << iter;
    //string file_out_index = fileout + stream_index.str();
    //extract_info_with_old_grad(file_out_index, n_sp, Index_of_sp, n_tet, tetrahedras, n_tri, triangles, x, sparse_Stress, f, grad_s, grad_s_old);
    return O;
}
