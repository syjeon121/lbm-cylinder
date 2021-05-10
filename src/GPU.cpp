#include "GPU.h"
#include <ctime>


void GPU() {
	LBM_GPU lbm;
	clock_t t0, t1;
	Error_Lp error_Lp;
//	float eps = 0.000001;
	float e = 1.0;

	int n = 0;
	int n_max = 1;
	int print_step = 100000;
	int print_step2 = 1000;

	cout << endl << "Uniform Flow around Circular Cylinder for LBM using GPU (WALL AND CURVED BOUNDARY)" << endl << endl;
	cout << "// ================== Physical property =================== //" << endl;
	cout << "Physical length in X direction = " << lbm.Lx << endl;
	cout << "Physical length in Y direction = " << lbm.Ly << endl;
	cout << "Physical grid step in X direction = " << lbm.del_x_p << endl;
	cout << "Physical grid step in Y direction = " << lbm.del_y_p << endl;
	cout << "Physical time step = " << lbm.del_t_p << endl;
	cout << "Physical mean velocity at Inlet boundary = " << (2.0 / 3.0)*lbm.Um_p << endl;
	cout << "Physical kinematic viscosity = " << lbm.nu_p << endl;
	cout << "Reynolds number = " << (2.0 / 3.0)*lbm.Um_p*lbm.D / lbm.nu_p << endl;
	cout << "// ======================================================== //" << endl << endl;

	cout << "// ===================== LBM property ===================== //" << endl;
	cout << "(LBM) The number of Node in X direction = " << lbm.nx << endl;
	cout << "(LBM) The number of Node in Y direction = " << lbm.ny << endl;
	cout << "(LBM) The mean velocity at Inlet boundary = " << (2.0 / 3.0)*lbm.Um << endl;
	cout << "(LBM) kinematic viscosity = " << lbm.nu << endl;
	cout << "Relaxation time = " << lbm.tau << endl;
	cout << "Reynolds number = " << (2.0 / 3.0)*lbm.Um*lbm.sny / lbm.nu << endl;
	cout << "The number of node at Solid = " << lbm.sn << endl;
	cout << "The number of node of diameter = " << lbm.sny << endl;
	cout << "// ======================================================== //" << endl << endl;
	
	cout << "Print step = " << print_step << endl;
	cout << "The number of maximum iteration = " << n_max << endl;

	t0 = clock();

	while (n < n_max){
		lbm.Streaming();
		lbm.BC_bounceback();
//		lbm.BC_extra();
		lbm.Collision();

		n++;
		if (n%print_step == 0) {
			lbm.Error();
			e = lbm.error;
			cout << "Iteration : [" << n << "], Error check : [" << e << "]" << endl;

		}

		lbm.Update();
		lbm.Momentum();

	//	if (n >= n_max - 10000)
	//		if (n%print_step2 == 0) lbm.Print();
		
	}
	t1 = clock();


	lbm.Print();
	cout << endl;
	cout << "The force in X direction : " << lbm.Fx << endl;
	cout << "The force in Y direction : " << lbm.Fy << endl;
	cout << "The coefficient of Drag : " << lbm.Cd << endl;
	cout << "The coefficient of lift : " << lbm.Cl << endl;
	cout << "Computation time : " << double(t1 - t0) / CLOCKS_PER_SEC << "[s]" << endl;
	cout << "Total Iteration : [" << n << "], Error check : [" << e << "]" << endl;
}