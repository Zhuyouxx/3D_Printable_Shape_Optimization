#include <Grad.h>
#include <DWCO.h>

//#include "../include/Grad.h"
//#include "../include/DWCO.h"

double mu_p = 1e-4, mu_sc = 1e4, mu_f = 1e1, mu_h = 1e6, mu_s = 1e4;
int iter_oi = 0;
double max_stress = 0.0;
double max_stress_bound = 0.0;
int iter_smooth_times = 800;
int not_smooth_count = 0;
int iter_smooth_min_bound = 400;
int iter_smooth_max_bound = 800;
int stress_opt_count = 0;
int stress_opt_bound = 30;
Eigen::VectorXd pressure;
Eigen::VectorXd u;
Eigen::MatrixXd V_before_opt;
double t_bound = 0.0;
string smooth_file_name;
int smooth_iter = 0;
MMG5_pMesh		mmgMesh;
MMG5_pSol       mmgSol;
double optimization_I(Eigen::VectorXd& grad_O_I, Eigen::MatrixXd& points, Eigen::MatrixXi tetrahedras, Eigen::MatrixXi triangles, Eigen::MatrixXd& V, Eigen::MatrixXi F, Eigen::MatrixXi edges,
	unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old, int& successful_o1, string log_file_name) {
	int nV = V.rows();
	double O_I = 0.0;
	Eigen::VectorXd grad_I;
	grad_I.resize(points.rows() * 3); grad_I.setZero();

	//calculate stress_grad O_p
	Eigen::VectorXd grad_stress;
	int successful = 1;
	//double enery_stress = stress_grad(points, triangles, tetrahedras, V, F, Old_2_new, new_2_Old, grad_stress, max_stress, successful);
	double enery_stress = stress_grad(points, triangles, tetrahedras, V, F, Old_2_new, new_2_Old, grad_stress, max_stress, successful, 0, pressure, u);
	if (!successful) {
		successful_o1 = 0;
		return 0.0;
	}
	O_I += mu_p * enery_stress;
	std::cout << "grad_stress : " << grad_stress.size() << endl;
	std::cout << "grad_I : " << grad_I.size() << endl;
	grad_I += mu_p * grad_stress;
	std::cout << "The Stress energy is : " << enery_stress << std::endl;

	//calculate intersect_grad O_sc
	Eigen::VectorXd grad_sc;
	const double dhat = 1e-3;
	double barrier_potential = energy_self_intersection(V, F, edges, grad_sc, dhat);
	std::cout << "The Barrier_potential is : " << barrier_potential << std::endl;
	grad_I += mu_sc * grad_stress;
	O_I += mu_sc * barrier_potential;
	for (int i = 0; i < nV; i++) {
		int id = new_2_Old[i];
		grad_I.segment<3>(3 * id) += mu_sc * grad_sc.segment<3>(3 * i);
	}


	//calculate Volume_grad O_f
	Eigen::VectorXd grad_volume;
	double dhat_clb = 1e-5;
	double energy_volume = Energy_Volume(points, tetrahedras, grad_volume, Old_2_new, new_2_Old, dhat_clb);
	cout << "The Volume energy : " << energy_volume << endl;
	O_I += mu_f * energy_volume;
	//grad_I += mu_f * grad_volume;
	for (int i = 0; i < nV; i++) {
		int id = new_2_Old[i];
		grad_I.segment<3>(3 * id) += mu_sc * grad_volume.segment<3>(3 * i);
	}
	successful_o1 = 1;
	grad_O_I = grad_I;
	return O_I;

}


double optimization_II(Eigen::VectorXd& grad_O_II, Eigen::MatrixXd& points, Eigen::MatrixXi tetrahedras, Eigen::MatrixXi triangles,
	Eigen::MatrixXd& V, Eigen::MatrixXi F, Eigen::MatrixXi edges, Eigen::MatrixXd V_target, Eigen::MatrixXi F_target,
	unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old, double m_h_o2, double m_s_o2) {
	int nV = V.rows();
	double O_II = 0.0;
	Eigen::VectorXd grad_II;
	grad_II.resize(points.rows() * 3); grad_II.setZero();

	if (m_h_o2 != 0.0) {
		//calculate hausdorff energy O_h
		Eigen::VectorXd grad_hausdorff;
		double energy_h = energy_hausdorff(V, F, V_target, F_target, grad_hausdorff);
		//std::cout << "The Hausdorff distance energy is : " << energy_h << std::endl;
		O_II += m_h_o2 * energy_h;
		//std::cout << "grad_hausdorff : " << mu_h * grad_hausdorff.norm() << std::endl;
		for (int i = 0; i < nV; i++) {
			int id = new_2_Old[i];
			grad_II.segment<3>(3 * id) += m_h_o2 * grad_hausdorff.segment<3>(3 * i);
		}
	}

	if (m_s_o2 != 0.0)
	{
		//calculate smooth energy O_s
		Eigen::VectorXd grad_smooth;
		double energy_shell = Energy_shell(V, F, V_before_opt, grad_smooth);
		//std::cout << "The Discrete Shell energy is : " << energy_shell << std::endl;
		//std::cout << "grad_smooth : " << mu_s * grad_smooth.norm() << std::endl;
		O_II += m_s_o2 * energy_shell;
		for (int i = 0; i < nV; i++) {
			int id = new_2_Old[i];
			grad_II.segment<3>(3 * id) += m_s_o2 * grad_smooth.segment<3>(3 * i);
		}
	}

	//calculate intersect_grad O_sc
	Eigen::VectorXd grad_sc;
	const double dhat = 1e-3;
	double barrier_potential = energy_self_intersection(V, F, edges, grad_sc, dhat);
	//std::cout << "grad_sc : " << mu_sc * grad_sc.norm() << std::endl;
	//std::cout << "The Barrier_potential is : " << barrier_potential << std::endl;
	//grad_I += mu_sc * grad_stress;
	O_II += mu_sc * barrier_potential;
	for (int i = 0; i < nV; i++) {
		int id = new_2_Old[i];
		grad_II.segment<3>(3 * id) += mu_sc * grad_sc.segment<3>(3 * i);
	}


	//calculate Volume_grad O_f
	Eigen::VectorXd grad_volume;
	double dhat_clb = 1e-5;
	double energy_volume = Energy_Volume(points, tetrahedras, grad_volume, Old_2_new, new_2_Old, dhat_clb);
	//std::cout << "grad_volume : " << mu_f * grad_volume.norm() << std::endl;
	//cout << "The Volume energy : " << energy_volume << endl;
	O_II += mu_f * energy_volume;
	//grad_I += mu_f * grad_volume;
	for (int i = 0; i < nV; i++) {
		int id = new_2_Old[i];
		grad_II.segment<3>(3 * id) += mu_f * grad_volume.segment<3>(3 * i);
	}
	grad_O_II = grad_II;
	return O_II;

}

int optimization(int opt_option, double edge_length, Eigen::MatrixXd& points, Eigen::MatrixXi tetrahedras, Eigen::MatrixXi triangles,
	Eigen::MatrixXd& V, Eigen::MatrixXi F, Eigen::MatrixXi edges, Eigen::MatrixXd V_target, Eigen::MatrixXi F_target,
	unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old, string log_file_name, string result) {
	int nV = V.rows();
	V.resize(nV, 3);
	Eigen::MatrixXd V_temp = V;
	Eigen::MatrixXd P_temp = points;
	Eigen::VectorXd grad_O_I, grad_O_II;
	double alpha_o1 = 1e-2, alpha_o2 = 1e-2;
	int successful_o1 = 1;
	if (opt_option == 1 || opt_option == 3) {
		double O_I = optimization_I(grad_O_I, points, tetrahedras, triangles, V, F, edges, Old_2_new, new_2_Old, successful_o1, log_file_name);
		if (!successful_o1) {
			std::cout << "stress error : GEP cannot be computed !" << endl;
			return 0;
		}
		if (opt_option == 3) {
			double t_init = update_t(points, tetrahedras, V, F, new_2_Old, max_stress, successful_o1, pressure, u);
			max_stress_bound = max_stress;
		}

		std::cout << "|| Value of O_I  || " << iter_oi << " ||: " << O_I << endl;
		Eigen::Vector3d max_grad; double max_grad_norm = 0.0;
		for (int i = 0; i < nV; i++) {
			int id = new_2_Old[i];
			Eigen::Vector3d delta_s = grad_O_I.segment<3>(3 * id);
			if (delta_s.norm() > max_grad_norm) {
				max_grad = delta_s;
				max_grad_norm = max_grad.norm();
			}
		}
		while (max_grad_norm * alpha_o1 > edge_length * 0.5)
			alpha_o1 *= 0.5;
		std::cout << "max update step : " << max_grad_norm * alpha_o1 << endl;
		while (true) {
			for (int i = 0; i < nV; i++) {
				int id = new_2_Old[i];
				Eigen::Vector3d delta_s = grad_O_I.segment<3>(3 * id) * alpha_o1;
				V_temp.row(i) = V.row(i) - delta_s.transpose();
				P_temp.row(id) = points.row(id) - delta_s.transpose();
			}
			double energy_temp = 0.0;

			//calculate Volume_grad O_f
			Eigen::VectorXd grad_volume;
			double dhat_clb = 1e-5;
			double energy_volume = Energy_Volume(P_temp, tetrahedras, grad_volume, Old_2_new, new_2_Old, dhat_clb);
			//std::cout << "|| Test || Volume energy : " << energy_volume << std::endl;
			energy_temp += mu_f * energy_volume;
			if (isinf(energy_temp)) {
				alpha_o1 *= 0.5;
				if (max_grad_norm * alpha_o1 < 1e-10) {
					std::cout << "the max update is " << max_grad_norm * alpha_o1 << ", which is much small." << std::endl;
					return 0;
				}
				continue;
			}
			//calculate intersect_grad O_sc
			const double dhat = 1e-1;
			bool is_collision_free = is_intersect_free(V, V_temp, F, edges, dhat);
			if (!is_collision_free) {
				alpha_o1 *= 0.5;
				if (max_grad_norm * alpha_o1 < 1e-20)
					return 0;
				continue;
			}

			double t_new = update_t(P_temp, tetrahedras, V_temp, F, new_2_Old, max_stress, successful_o1, pressure, u);
			cout << " || Test || : max stress : " << max_stress << endl;
			if (t_new > t_bound) {
				std::cout << "t_bound : " << t_bound << " || t_new : " << t_new << endl;
				alpha_o1 *= 0.5;
				if (max_grad_norm * alpha_o1 < 1e-10) {
					std::cout << "the max update is " << max_grad_norm * alpha_o1 << ", which is much small." << std::endl;
					return 0;
				}
				continue;
			}
			else {
				t_bound = t_new;
			}
			if (!isinf(energy_temp) & is_collision_free) {
				double max_update = 0.0;
				for (int i = 0; i < nV; i++) {
					int id = new_2_Old[i];
					Eigen::Vector3d delta_s = grad_O_I.segment<3>(3 * id) * alpha_o1;
					double update_step = delta_s.norm();
					if (update_step > max_update)
						max_update = update_step;
					V.row(i) = V.row(i) - grad_O_I.segment<3>(3 * id).transpose() * alpha_o1;
					points.row(id) = points.row(id) - grad_O_I.segment<3>(3 * id).transpose() * alpha_o1;
				}
				//alpha_o1 *= 2.0;
				std::cout << "max update step real : " << max_update << endl;
				break;
			}
			alpha_o1 *= 0.5;
			if (alpha_o1 < 1e-20)
				return 0;
		}
		std::ofstream file(log_file_name, std::ios::app);

		if (!file.is_open()) {
			std::cerr << "打开文件失败！" << std::endl;
			return 1;
		}
		// 向文件写入内容
		file << "|- OPT Stress -| Iter : " << iter_oi << " || Stress : " << max_stress << " || t : " << t_bound << endl;
		file.close();
		igl::writeOBJ("..//Results//" + result + "_" + std::to_string(int(iter_oi)) + ".obj", V, F);
		cout << " || OPT || " << int(iter_oi) << " || : max stress : " << max_stress << " || t : " << t_bound << endl;
		iter_oi += 1;
		stress_opt_count += 1;
		if (stress_opt_count > stress_opt_bound) {
			std::cout << "Stress optimization has stagnated, with continuous optimization exceeding " << stress_opt_count << " times." << std::endl;
			iter_oi -= stress_opt_bound;
			return 0;
		}
		if (max_stress < max_stress_bound * 0.9) {
			std::cout << "Max_stress has been reduced to 0.9 of its initial value, and we will proceed with the smoothing process." << std::endl;
			cout << " || OPT || " << int(iter_oi - 1) << " || : max stress : " << max_stress << " || t : " << t_bound << endl;
			return 2;
		}
	}
	else if (opt_option == 2) {
		//optimization ii
		//initial smooth region 
		//Eigen::VectorXd grad_hausdorff;
		//energy_hausdorff(V, F, V_target, F_target, grad_hausdorff);
		stress_opt_count = 0;
		Eigen::VectorXd hausdorff_dis(V.rows()); hausdorff_dis.setZero();
		std::vector<std::vector<int>> A;
		igl::adjacency_list(F, A);
		for (int i = 0; i < nV; i++) {
			//Eigen::Vector3d grad_unit = grad_hausdorff.segment<3>(3 * i);
			//double h_d = grad_unit.norm();
			double h_d = (V.row(i) - V_before_opt.row(i)).norm();
			if (h_d > edge_length * 1e-2) {
				hausdorff_dis[i] = 1.0;
				for (auto id : A[i]) {
					hausdorff_dis[id] = 1.0;
				}
			}
		}
		{
			std::string save_str = "..//Info//smooth_" + smooth_file_name + std::to_string(int(iter_oi)) + "_points.txt";
			char* save_file = const_cast<char*>(save_str.c_str());
			std::ofstream file(save_file);
			if (!file.is_open()) {
				std::cerr << "Failed to open file." << std::endl;
				return 0;
			}
			for (int i = 0; i < nV; i++)
			{
				file << V(i, 0) << endl << V(i, 1) << endl << V(i, 2) << endl;
			}
			file.close();

			save_str = "..//Info//smooth_" + smooth_file_name + std::to_string(int(iter_oi)) + "_triangles.txt";
			save_file = const_cast<char*>(save_str.c_str());
			file.open(save_file);
			if (!file.is_open()) {
				std::cerr << "Failed to open file." << std::endl;
				return 0;
			}
			for (int t = 0; t < F.rows(); t++)
				file << F(t, 0) << endl << F(t, 1) << endl << F(t, 2) << endl;
			file.close();

			save_str = "..//Info//smooth_" + smooth_file_name + std::to_string(int(iter_oi)) + "_haussdorf_distance.txt";
			save_file = const_cast<char*>(save_str.c_str());
			file.open(save_file);
			if (!file.is_open()) {
				std::cerr << "Failed to open file." << std::endl;
				return 0;
			}
			for (int i = 0; i < nV; i++)
			{
				file << hausdorff_dis[i] << endl;
			}
			file.close();
		}

		double O_II_old = 0.0;
		for (int iter_smooth = 0; iter_smooth < iter_smooth_times*2; iter_smooth++) {
			double m_h_o2 = mu_h;
			double m_s_o2 = mu_s;
			if ((iter_smooth + 1) < (int)(iter_smooth_times)) {
				m_h_o2 = 0.0;
			}
			else {
				m_h_o2 = mu_h;
				m_s_o2 = 0.0;
			}
			double O_II = optimization_II(grad_O_II, points, tetrahedras, triangles, V, F, edges, V_target, F_target, Old_2_new, new_2_Old, m_h_o2, m_s_o2);
			Eigen::Vector3d max_grad; double max_grad_norm = 0.0;
			for (int i = 0; i < nV; i++) {
				int id = new_2_Old[i];
				grad_O_II.segment<3>(3 * id) *= hausdorff_dis[i];
				Eigen::Vector3d delta_s = grad_O_II.segment<3>(3 * id) * hausdorff_dis[i];
				if (delta_s.norm() > max_grad_norm) {
					max_grad = delta_s;
					max_grad_norm = max_grad.norm();
				}
			}
			while (max_grad_norm * alpha_o2 > edge_length * 1e-3) {
				alpha_o2 *= 0.5;
			}

			if ((iter_smooth + 1) % 1000 == 0)
				std::cout << "|| Value of O_II || " << iter_oi << " ||: " << O_II << endl;
			while (true) {
				for (int i = 0; i < nV; i++) {
					int id = new_2_Old[i];
					V_temp.row(i) = V.row(i) - grad_O_II.segment<3>(3 * id).transpose() * alpha_o2;
					P_temp.row(id) = points.row(id) - grad_O_II.segment<3>(3 * id).transpose() * alpha_o2;
				}
				double energy_temp = 0.0;
				//calculate Volume_grad O_f
				Eigen::VectorXd grad_volume;
				double dhat_clb = 1e-5;
				double energy_volume = Energy_Volume(P_temp, tetrahedras, grad_volume, Old_2_new, new_2_Old, dhat_clb);
				//std::cout << "|| Test || Volume energy : " << energy_volume << std::endl;
				energy_temp += mu_f * energy_volume;
				if (isinf(energy_temp)) {
					alpha_o2 *= 0.5;
					if (max_grad_norm * alpha_o2 < 1e-20)
						return 0;
					continue;
				}
				//calculate intersect_grad O_sc
				const double dhat = 1e-1;
				bool is_collision_free = is_intersect_free(V, V_temp, F, edges, dhat);
				if (!is_collision_free) {
					alpha_o2 *= 0.5;
					if (max_grad_norm * alpha_o2 < 1e-20)
						return 0;
					continue;
				}
				double max_before = max_stress;
				double t_new = update_t(P_temp, tetrahedras, V_temp, F, new_2_Old, max_stress, successful_o1, pressure, u);
				cout << " || Test O_II || : max stress : " << max_stress << endl;
				if (t_new > t_bound * (1.0)) {
					std::cout << "t_bound : " << t_bound << " || t_new : " << t_new << " | discrepancy : " << t_new - t_bound << " | the max update is " << max_grad_norm * alpha_o2 << endl;
					alpha_o2 *= 0.5;
					if (max_grad_norm * alpha_o2 < 1e-15) {
						std::cout << "the max update is " << max_grad_norm * alpha_o2 << ", which is much small." << std::endl;
						if ((iter_smooth + 1) > iter_smooth_times ) {
							not_smooth_count = 0;
							iter_smooth_times = iter_smooth_max_bound;
							//igl::writeOBJ("..//Results//smooth01_" + std::to_string(int(iter_oi)) + "_" + std::to_string(int(iter_smooth + 1)) + ".obj", V, F);
							igl::writeOBJ("..//Results//" + result + "_" + std::to_string(int(iter_oi)) + ".obj", V, F);
							std::ofstream file(log_file_name, std::ios::app);
							if (!file.is_open()) {
								std::cerr << "打开文件失败！" << std::endl;
								return 1;
							}
							// 向文件写入内容
							file << "|- OPT Smooth -| hausdorff |Smooth Iter : " << iter_oi << " || Stress : " << max_before << " || t : " << t_bound << endl;
							file.close();
							smooth_iter = iter_oi;
							iter_oi++;
							return 3;
						}
						not_smooth_count += 1;
						if (not_smooth_count > 10) {
							iter_smooth_times = max(iter_smooth_min_bound, iter_smooth_times - 50);
						}
						iter_oi = smooth_iter + 1;
						return 0;
					}
					continue;
				}

				if (!isinf(energy_temp) & is_collision_free) {
					V = V_temp;
					points = P_temp;
					alpha_o2 *= 2.0;
					igl::writeOBJ("..//Results//smooth_" + smooth_file_name + "_" + std::to_string(int(iter_oi)) + "_" + std::to_string(int(iter_smooth + 1)) + ".obj", V, F);
					std::cout << " || Update " << (iter_smooth + 1) << "|| Grad norm : " << grad_O_II.norm() << " | alpha : " << alpha_o2 << std::endl;

					if (m_s_o2 != 0.0) {
						if ((iter_smooth + 1) % 100 == 0) {
							std::ofstream file(log_file_name, std::ios::app);
							if (!file.is_open()) {
								std::cerr << "打开文件失败！" << std::endl;
								return 1;
							}
							// 向文件写入内容
							file << "|- OPT Smooth -| shell |Smooth Iter : " << iter_smooth + 1 << " || Stress : " << max_stress << " || t : " << t_bound << endl;
							file.close();
						}
					}
					else {
						std::ofstream file(log_file_name, std::ios::app);
						if (!file.is_open()) {
							std::cerr << "打开文件失败！" << std::endl;
							return 1;
						}
						// 向文件写入内容
						file << "|- OPT Smooth -| hausdorff |Smooth Iter : " << iter_smooth + 1 << " || Stress : " << max_stress << " || t : " << t_bound << endl;
						file.close();
					}
					break;
				}
				if (max_grad_norm * alpha_o2 < 1e-15)
					return 0;
			}
		}
		std::ofstream file_log(log_file_name, std::ios::app);

		if (!file_log.is_open()) {
			std::cerr << "打开文件失败！" << std::endl;
			return 1;
		}
		// 向文件写入内容
		file_log << "|- OPT Smooth -| Iter : " << iter_oi << " || Stress : " << max_stress << " || t : " << t_bound << endl;
		file_log.close();
		not_smooth_count = 0;
		iter_smooth_times = iter_smooth_max_bound;
		igl::writeOBJ("..//Results//" + result + "_" + std::to_string(int(iter_oi)) + ".obj", V, F);
		igl::writeOBJ("..//Results//" + result + "_smooth_" + std::to_string(int(iter_oi)) + ".obj", V, F);
		cout << " || OPT || " << int(iter_oi) << " || : max stress : " << max_stress << " || t : " << t_bound << endl;
		iter_oi += 1;
		return 1;
	}

	return 1;
}

int main() {
	string file_in_path = "../IO/";
	string out_path = "../Results/";
	string model_name = "hanging_ball_repair2";
	string filein = file_in_path + model_name;
	string target_model_name = model_name + "_target";
	string result_name = "hanging_ball_repair2_05";
	smooth_file_name = result_name;
	double edge_length = 0.01;
	string remesh_model = file_in_path + model_name + ".obj";
	string target_model = file_in_path + target_model_name + ".obj";
	string remesh_out = file_in_path + model_name + ".mesh";
	string remesh_mmg_out = file_in_path + model_name + "_mmg.mesh";
	string log_file = out_path + result_name + ".txt";

	tet_mesh(edge_length, remesh_model, remesh_out, 0.035, true, true);
	edge_length *= 0.75;
	tet_remesh_mmg(remesh_out, remesh_mmg_out, edge_length);
	cout << "Edge length : " << edge_length << endl;

	//read tet
	char* filename = const_cast<char*>(remesh_mmg_out.c_str());
	MMG5_int        k, np, n_tet, n_tri, n_sp;
	//np:number of points;  n_tet:number of tetraheras;    n_tri:number of triangles on surface;    n_sp:number of points on surface;
	Eigen::MatrixXd points;
	Eigen::MatrixXi triangles;
	Eigen::MatrixXi tetrahedras;

	if (!read_mesh(mmgMesh, mmgSol, filename, points, triangles, tetrahedras)) {
		return 0;
	}
	if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, NULL, &n_tri, NULL, NULL) != 1)  exit(EXIT_FAILURE);
	cout << "number of points:" << np << endl;
	cout << "number of tetrahedras:" << n_tet << endl;
	cout << "number of triangles:" << n_tri << endl;
	double worst_quality = get_worst_quality(mmgMesh, mmgSol, points, triangles, tetrahedras);
	std::cout << "--|| WORST QUALITY : " << worst_quality << " ||--" << std::endl;
	std::ofstream file(log_file);
	if (!file.is_open()) {
		std::cerr << "can't open the file . " << std::endl;
		return 0;
	}
	file << "new_log : " << std::endl;
	file << "model : " << model_name << std::endl;
	file.close();

	Eigen::MatrixXd V;
	Eigen::MatrixXi F, edges;
	Eigen::VectorXi new_2_Old;
	unordered_map<int, int> Old_2_new;
	Eigen::VectorXd Area;
	std::cout << "Edge Length : " << edge_length << std::endl;
	extract_surface(points, triangles, Old_2_new, V, F, new_2_Old);
	double min_edge_length = avg_edge_length(V, F);
	edge_length = min_edge_length;

	std::cout << "Edge Length : " << edge_length << std::endl;
	igl::edges(F, edges);
	int nV = V.rows();
	int nF = F.rows();
	Area.resize(nV); Area.setZero();

	//get target model
	igl::writeOBJ(target_model, V, F);
	Eigen::MatrixXd V_target;
	Eigen::MatrixXi F_target;
	igl::readOBJ(target_model, V_target, F_target);
	V_before_opt = V_target.eval();
	Calculate_Area(V, F, Area);
	//igl::writeOBJ("..//Results//initial.obj", V, F);
	int successful = 1;

	u.resize(3 * points.rows() + 6); u.setRandom();
	calculate_pressure(points, tetrahedras, V, F, new_2_Old, successful, pressure, u);
	if (!successful) {
		cout << "can not calculate pressure!" << endl;
		return 0;
	}
	//initial t

	double t_init = update_t(points, tetrahedras, V, F, new_2_Old, max_stress, successful, pressure, u);
	std::cout << "|| Init || t_init : " << t_init << " || t_bound : " << t_bound << std::endl;
	t_bound = t_init;
	max_stress_bound = max_stress;
	cout << " || OPT || " << int(iter_oi) << " || : max stress : " << max_stress << endl;

	file.open(log_file, std::ios::app);
	if (!file.is_open()) {
		std::cerr << "打开文件失败！" << std::endl;
		return 1;
	}
	// 向文件写入内容
	file << "|- OPT Stress -| Iter : " << iter_oi << " || Stress : " << max_stress << " || t : " << t_bound << endl;
	// 关闭文件
	file.close();
	igl::writeOBJ("..//Results//" + result_name + "_" + std::to_string(int(iter_oi)) + ".obj", V, F);
	iter_oi += 1;


	int opt_option = 1;
	int write_id = -1;
	for (int iter = 0; iter < 300; iter++) {
		int is_op = optimization(opt_option, edge_length, points, tetrahedras, triangles, V, F, edges, V_target, F_target, Old_2_new, new_2_Old, log_file, result_name);
		if (is_op) {
			opt_option = is_op;
			write_id = iter;
		}
		else {
			if (write_id > 0) {
				iter_oi -= 1;
				remesh_model = "..//Results//" + result_name + "_" + std::to_string(int(iter_oi)) + ".obj";
				remesh_out = "..//IO//test_" + std::to_string(int(iter_oi)) + ".mesh";
				remesh_mmg_out = "..//IO//test_" + std::to_string(int(iter_oi)) + "_mmg.mesh";
				tet_mesh(edge_length, remesh_model, remesh_out, 0.035, true, true);
				edge_length *= 0.75;
				tet_remesh_mmg(remesh_out, remesh_mmg_out, edge_length);
				cout << "Edge length : " << edge_length << endl;
				filename = const_cast<char*>(remesh_mmg_out.c_str());
				if (!read_mesh(mmgMesh, mmgSol, filename, points, triangles, tetrahedras)) {
					return 0;
				}
			}
			else {
				remesh_model = file_in_path + model_name + ".obj";
				remesh_out = file_in_path + model_name + ".mesh";
				remesh_mmg_out = file_in_path + model_name + "_mmg.mesh";
				tet_mesh(edge_length, remesh_model, remesh_out, 0.035, true, true);
				edge_length *= 0.75;
				tet_remesh_mmg(remesh_out, remesh_mmg_out, edge_length);
				cout << "Edge length : " << edge_length << endl;
				filename = const_cast<char*>(remesh_mmg_out.c_str());
				if (!read_mesh(mmgMesh, mmgSol, filename, points, triangles, tetrahedras)) {
					return 0;
				}
			}
			//np:number of points;  n_tet:number of tetraheras;    n_tri:number of triangles on surface;    n_sp:number of points on surface;

			if (!get_info_from_mesh(mmgMesh, mmgSol, points, triangles, tetrahedras)) {
				return 0;
			}
			if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, NULL, &n_tri, NULL, NULL) != 1)  exit(EXIT_FAILURE);
			cout << "number of points:" << np << endl;
			cout << "number of tetrahedras:" << n_tet << endl;
			cout << "number of triangles:" << n_tri << endl;
			worst_quality = get_worst_quality(mmgMesh, mmgSol);
			std::cout << "--|| WORST QUALITY : " << worst_quality << " ||--" << std::endl;
			extract_surface(points, triangles, Old_2_new, V, F, new_2_Old);
			double min_edge_length = avg_edge_length(V, F);
			edge_length = min_edge_length;
			std::cout << "Edge Length : " << edge_length << std::endl;
			igl::edges(F, edges);
			nV = V.rows();
			nF = F.rows();
			Area.resize(nV); Area.setZero();
			V_before_opt = V.eval();
			Calculate_Area(V, F, Area);
			int successful = 1;

			u.resize(points.rows() * 3 + 6); u.setRandom();
			calculate_pressure(points, tetrahedras, V, F, new_2_Old, successful, pressure, u);
			if (!successful) {
				cout << "can not calculate pressure!" << endl;
				return 0;
			}

			//initial t
			double t_init = update_t(points, tetrahedras, V, F, new_2_Old, max_stress, successful, pressure, u);
			std::cout << "|| Init || t_init : " << t_init << " || t_bound : " << t_bound << std::endl;
			t_bound = t_init;
			max_stress_bound = max_stress;
			cout << " || Remesh || OPT || " << int(iter_oi) << " || : max stress : " << max_stress << endl;

			file.open(log_file, std::ios::app);
			if (!file.is_open()) {
				std::cerr << "打开文件失败！" << std::endl;
				return 1;
			}
			// 向文件写入内容
			file << "-- Remesh --" << endl;
			file << "|- OPT Stress -| Iter : " << iter_oi << " || Stress : " << max_stress << " || t : " << t_bound << endl;
			// 关闭文件
			file.close();
			igl::writeOBJ("..//Results//" + result_name + "_remesh_" + std::to_string(int(iter_oi)) + ".obj", V, F);
			iter_oi += 1;


			opt_option = 1;

			iter--;
		}
	}

}