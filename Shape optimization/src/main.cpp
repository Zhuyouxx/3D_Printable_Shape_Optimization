#include <Grad.h>
#include <DWCO.h>

//#include "../include/Grad.h"
//#include "../include/DWCO.h"

double mu_p = 1e-4, mu_sc = 1e4, mu_f = 1e2, mu_h = 1e6, mu_s = 1e3;
int iter_oi = 0;
double max_stress = 0.0;
int iter_smooth_times = 6000;
double optimization_I(Eigen::VectorXd& grad_O_I, Eigen::MatrixXd& points, Eigen::MatrixXi tetrahedras, Eigen::MatrixXi triangles, Eigen::MatrixXd& V, Eigen::MatrixXi F, Eigen::MatrixXi edges,
	unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old, int& successful_o1, string log_file_name) {
	int nV = V.rows();
	double O_I = 0.0;
	Eigen::VectorXd grad_I;
	grad_I.resize(points.rows() * 3); grad_I.setZero();

	//calculate stress_grad O_p
	Eigen::VectorXd grad_stress;
	int successful = 1;
	double enery_stress = stress_grad(points, triangles, tetrahedras, V, F, Old_2_new, new_2_Old, grad_stress, max_stress, successful);
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
		double energy_shell = Energy_shell(V, F, grad_smooth);
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

int optimization(int iter, double edge_length, Eigen::MatrixXd& points, Eigen::MatrixXi tetrahedras, Eigen::MatrixXi triangles,
	Eigen::MatrixXd& V, Eigen::MatrixXi F, Eigen::MatrixXi edges, Eigen::MatrixXd V_target, Eigen::MatrixXi F_target,
	unordered_map<int, int> Old_2_new, Eigen::VectorXi new_2_Old, string log_file_name) {

	int nV = V.rows();
	V.resize(nV, 3);
	Eigen::MatrixXd V_temp = V;
	Eigen::MatrixXd P_temp = points;
	Eigen::VectorXd grad_O_I, grad_O_II;
	double alpha_o1 = 1e-2, alpha_o2 = 0.1;

	int successful_o1 = 1;
	double O_I = optimization_I(grad_O_I, points, tetrahedras, triangles, V, F, edges, Old_2_new, new_2_Old, successful_o1, log_file_name);
	if (!successful_o1) {
		std::cout << "stress error : GEP cannot be computed !" << endl;
		return 0;
	}
	std::cout << "|| Value of O_I  || " << iter << " ||: " << O_I << endl;
	Eigen::Vector3d max_grad; double max_grad_norm = 0.0;
	for (int i = 0; i < nV; i++) {
		int id = new_2_Old[i];
		Eigen::Vector3d delta_s = grad_O_I.segment<3>(3 * id);
		if (delta_s.norm() > max_grad_norm) {
			max_grad = delta_s;
			max_grad_norm = max_grad.norm();
		}
	}
	while (max_grad_norm * alpha_o1 > edge_length * 0.25)
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
			if (max_grad_norm * alpha_o1 < 1e-10)
				return 0;
			continue;
		}
		//calculate intersect_grad O_sc
		Eigen::VectorXd grad_sc;
		const double dhat = 1e-3;
		double barrier_potential = energy_self_intersection(V_temp, F, edges, grad_sc, dhat);
		//std::cout << "|| Test || The Barrier_potential is : " << barrier_potential << std::endl;
		//grad_I += mu_sc * grad_stress;
		energy_temp += mu_sc * barrier_potential;

		if (!isinf(energy_temp)) {
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
	file << "Iter : " << iter << " || Stress : " << max_stress << endl;
	// 关闭文件
	file.close();

	igl::writeOBJ("..//Results//test_OI_" + std::to_string(int(iter_oi)) + ".obj", V, F);
	iter_oi += 1;
	//igl::writeOBJ("..//Results//test_.obj", V, F);

	//optimization ii
	double O_II_old = 0.0;
	for (int iter_smooth = 0; iter_smooth < iter_smooth_times; iter_smooth++) {
		double m_h_o2 = mu_h;
		double m_s_o2 = mu_s;
		if ((iter_smooth + 1) < (int)(iter_smooth_times * 0.5)) {
			m_h_o2 = 0.0;
		}
		else {
			m_h_o2 = mu_h;
			m_s_o2 = mu_s * 1e-5;
		}
		double O_II = optimization_II(grad_O_II, points, tetrahedras, triangles, V, F, edges, V_target, F_target, Old_2_new, new_2_Old, m_h_o2, m_s_o2);

		Eigen::Vector3d max_grad; double max_grad_norm = 0.0;
		for (int i = 0; i < nV; i++) {
			int id = new_2_Old[i];
			Eigen::Vector3d delta_s = grad_O_II.segment<3>(3 * id);
			if (delta_s.norm() > max_grad_norm) {
				max_grad = delta_s;
				max_grad_norm = max_grad.norm();
			}
		}
		while (max_grad_norm * alpha_o2 > edge_length * 0.5) {
			alpha_o2 *= 0.5;
		}
		if ((iter_smooth + 1) % 1000 == 0)
			std::cout << "|| Value of O_II || " << iter << " ||: " << O_II << endl;
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
			if (!isinf(energy_temp) & is_collision_free) {
				V = V_temp;
				points = P_temp;
				alpha_o2 *= 2.0;
				if ((iter_smooth + 1) % 1000 == 0) {
					igl::writeOBJ("..//Results//test_OI_" + std::to_string(int(iter_smooth + 1)) + ".obj", V, F);
					std::cout << " || Update " << (iter_smooth + 1) << "|| Grad norm : " << grad_O_II.norm() << " | alpha : " << alpha_o2 << std::endl;
				}

				//return 1;
				break;
			}
			if (max_grad_norm * alpha_o2 < 1e-20)
				return 0;
		}
	}
	return 1;
}

int main() {
	string file_in_path = "../IO/";
	string out_path = "../Results/";
	string model_name = "hanging_ball";
	string filein = file_in_path + model_name;
	string target_model_name = "hanging_ball";
	string result_name = "hanging_ball";
	vector<double> grad;
	unordered_map<int, int> Old_2_new;
	//read model 
	double edge_length = 0.01;
	string remesh_model = file_in_path + model_name + ".obj";
	string target_model = file_in_path + target_model_name + ".obj";
	//string test_remesh_model = file_in_path + model_name + ".obj";
	string remesh_out = file_in_path + model_name + ".mesh";
	string remesh_mmg_out = file_in_path + model_name + "_mmg.mesh";
	string log_file = out_path + model_name + ".txt";
	tet_mesh(edge_length, remesh_model, remesh_out, 0.05, true, true);
	edge_length *= 0.75;
	tet_remesh_mmg(remesh_out, remesh_mmg_out, edge_length);
	cout << "Edge length : " << edge_length << endl;

	Eigen::MatrixXd V_target;
	Eigen::MatrixXi F_target;
	igl::readOBJ(target_model, V_target, F_target);

	//read tet
	char* filename = const_cast<char*>(remesh_mmg_out.c_str());
	MMG5_pSol       mmgSol;
	MMG5_int        k, np, n_tet, n_tri, n_sp;
	//np:number of points;  n_tet:number of tetraheras;    n_tri:number of triangles on surface;    n_sp:number of points on surface;
	Eigen::MatrixXd points;
	Eigen::MatrixXi triangles;
	Eigen::MatrixXi tetrahedras;

	MMG5_pMesh mmgMesh;
	if (!read_mesh(mmgMesh, mmgSol, filename, points, triangles, tetrahedras)) {
		return 0;
	}
	if (MMG3D_Get_meshSize(mmgMesh, &np, &n_tet, NULL, &n_tri, NULL, NULL) != 1)  exit(EXIT_FAILURE);
	cout << "number of points:" << np << endl;
	cout << "number of tetrahedras:" << n_tet << endl;
	cout << "number of triangles:" << n_tri << endl;

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
	Eigen::VectorXd Area;
	extract_surface(points, triangles, Old_2_new, V, F, new_2_Old);
	igl::edges(F, edges);
	int nV = V.rows();
	int nF = F.rows();
	Area.resize(nV); Area.setZero();

	Calculate_Area(V, F, Area);
	igl::writeOBJ("..//Results//initial.obj", V, F);
	int write_id = -1;
	for (int iter = 0; iter < 100; iter++) {
		int is_op = optimization(iter, edge_length, points, tetrahedras, triangles, V, F, edges, V_target, F_target, Old_2_new, new_2_Old, log_file);
		if (is_op) {
			write_id = iter;
			//igl::writeOBJ("D://IDE//VisualStudio_Project//DR_manufacturing//TestOn3D//Test3D//results//test_" + std::to_string(int(iter)) + ".obj", V, F);
			igl::writeOBJ("..//Results//" + result_name + "_" + std::to_string(int(iter)) + ".obj", V, F);
		}
		else {
			double worst_quality = 1.0;
			tet_remesh_mmg_Update_Points(mmgMesh, mmgSol, edge_length, points, 1, worst_quality);
			cout << "Edge length : " << edge_length << endl;
			int count = 0;
			while (worst_quality < 0.1) {
				tet_remesh_mmg_Update_Points(mmgMesh, mmgSol, edge_length, points, 0, worst_quality);
				cout << "Edge length : " << edge_length << endl;
				count++;
				if (count > 5)
					break;
			}
			if (count > 5 && write_id > 0) {
				remesh_model = "..//Results//" + result_name + "_" + std::to_string(int(write_id)) + ".obj";
				remesh_out = "..//IO//test_" + std::to_string(int(write_id)) + ".mesh";
				remesh_mmg_out = "..//IO//test_" + std::to_string(int(write_id)) + "_mmg.mesh";
				tet_mesh(edge_length, remesh_model, remesh_out, 0.05, true, true);
				edge_length *= 0.75;
				tet_remesh_mmg(remesh_out, remesh_mmg_out, edge_length);
				cout << "Edge length : " << edge_length << endl;
				filename = const_cast<char*>(remesh_mmg_out.c_str());
				if (!read_mesh(mmgMesh, mmgSol, filename, points, triangles, tetrahedras)) {
					return 0;
				}
			}
			else if (count > 5) {
				remesh_model = file_in_path + model_name + ".obj";
				remesh_out = file_in_path + model_name + ".mesh";
				remesh_mmg_out = file_in_path + model_name + "_mmg.mesh";
				tet_mesh(edge_length, remesh_model, remesh_out, 0.05, true, true);
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

			extract_surface(points, triangles, Old_2_new, V, F, new_2_Old);
			igl::edges(F, edges);
			nV = V.rows();
			nF = F.rows();
			Area.resize(nV); Area.setZero();

			Calculate_Area(V, F, Area);
			iter--;
		}
	}

}