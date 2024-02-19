#include "../include/Grad.h"

int main() {
	string path = "../IO/";
	string out_path = "../Results/";
	string model_name = "Test4E_H_003";
	string filein = path + model_name;
	string fileout = out_path + "Test4E_H_003_grad_04_";
	vector<double> grad;
	//MatrixXd t(3, 3);
	//t << 687992, -12421.4, 24146.7,
	//	-12421.4, 49104.1, 193569,
	//	24146.7, 193569, 770855;
	//MatrixXd t_inv = t.inverse();
	//cout << t_inv << endl;
	//cout << t_inv * t << endl;
	//t << 687992.24468083310723, -12421.409938080788368, 24146.667623258339486,
	//	-12421.409938080788368, 49104.063962409501037, 193568.73069354675253,
	//	24146.667623258339486, 193568.73069354675253, 770855.12440018510677;
	//t_inv = t.inverse();
	//cout << t_inv << endl;
	//cout << t_inv * t << endl;
	stress_grad(filein, fileout, grad,true);
	//stress_grad();
	//cout << "Grad:" << endl;
	//for (int i = 0; i < 10; i++) {
	//	cout << grad[3 * i] << "    " << grad[3 * i + 1] << "    " << grad[3 * i + 2] << endl;
	//}
}