
// ============================================================================ //
//  SUBROUTINE
// ============================================================================ //
//#include "CPU.h"
//#include <ctime>
//
//void CPU() {
//	LBM_CPU lbm;
//	clock_t t0, t1;
//	Error_Lp error_Lp;
//	float eps = 0.000001;
//	float e = 1.0;
//	
//	int n = 1;
//	int print_step = 1000;
//
//	cout << endl << "Uniform Flow around Circular Cylinder for LBM using CPU" << endl << endl;
//	cout << "// ================== Physical property =================== //" << endl;
//	cout << "Physical length in X direction = " << lbm.Lx << endl;
//	cout << "Physical length in Y direction = " << lbm.Ly << endl;
//	cout << "Physical velocity at Inlet boundary = " << lbm.Ux0_p << endl;
//	cout << "Physical kinematic viscosity = " << lbm.nu << endl;
//	cout << "Reynolds number = " << lbm.Ux0_p*lbm.Ly / lbm.nu << endl;
//	cout << "// ======================================================== //" << endl << endl;
//
//	cout << "// ===================== LBM property ===================== //" << endl;
//	cout << "(LBM) The number of Node in X direction = " << lbm.nx << endl;
//	cout << "(LBM) The number of Node in Y direction = " << lbm.ny << endl;
//	cout << "(LBM) The velocity at Inlet boundary = " << lbm.Ux0 << endl;
//	cout << "(LBM) kinematic viscosity = " << lbm.nu << endl;
//	cout << "Relaxation time = " << lbm.tau << endl;
//	cout << "Reynolds number = " << lbm.Ux0*lbm.ny / lbm.nu << endl;
//	cout << "The number of node at Solid = " << lbm.snx*lbm.sny << endl;
//	cout << "// ======================================================== //" << endl << endl;
//
//	cout << "eps = " << eps << endl;
//	t0 = clock();
//	while (e > eps) {
//
//		
//		lbm.Streaming();
//		lbm.BC_bounceback();
//		lbm.BC_vel();
//		lbm.Collision();
//		lbm.Error();
//		e = lbm.error;
//		lbm.Update();
//		n++;
//		if (n%print_step == 0) cout << "Iteration : [" << n << "], Error check : [" << e << "]" << endl;
//	}
//	t1 = clock();
//
//	lbm.Print();
//	
//	cout << endl;
//	cout << "Computation time : " << double(t1 - t0) / CLOCKS_PER_SEC << "[s]" << endl;
//	cout << "Total Iteration : [" << n - 1 << "], Error check : [" << e << "]" << endl;
//}
// ============================================================================ //







// ============================================================================ //
//  NORMAL
// ============================================================================ //
#include "CPU.h"
#include <ctime>

ifstream fin_CPU("in_CPU.txt");

ofstream fout_CPU("out_CPU.dat");
ofstream fout_CPU_Cd("out_CPU_Cd.dat");
ofstream fout_CPU_Ux0("out_CPU_Ux0.dat");

void CPU() {

// ============================================================================ //
//  VARIABLES
// ============================================================================ //
	int nx;								//The number of node in X direction
	int ny;								//The number of node in Y direction
	int a;								//The number of direction in D2Q9 model

	float Lx;							//The length in X direction
	float Ly;							//The length in Y direction
	float del_x_p;						//The grid step in X direction (Physical property)
	float del_y_p;						//The grid step in Y direction (Physical property)
	float del_t_p;						//The time step (Physical property)

	float del_x;						//The grid step in X direction (LBM property)
	float del_y;						//The grid step in Y direction (LBM property)
	float del_t;						//The time step (LBM property)

	float* Ux;							//The macroscopic velocity in X direction
	float* Uy;							//The macroscopic velocity in Y direction
	float* U;							//The macroscopic velocity
	float* rho;							//The macroscopic density

	float* Ux_p;						//The macroscopic velocity in X direction (Physical property)
	float* Uy_p;						//The macroscopic velocity in Y direction (Physical property)
	float* U_p;							//The macroscopic velocity (Physical property)
	float* P;							//The macroscopic pressure (Physical property)
	float Um_p;							//The maximum velocity (Physical property) [m/s]

	float* UxN;							//The macroscopic velocity in X direction at n + 1 step
	float* UyN;							//The macroscopic velocity in Y direction at n + 1 step
	float* UN;							//The macroscopic velocity at n + 1 step
	float* rhoN;						//The macroscopic density
	float Um;							//The maximum velocity (LBM property)

	float* f;							//The distribution function at n step
	float* ftemp;						//The distribution function at temp step
	float* fN;							//The distribution function at n + 1 step
	float* feq;							//The equilibrium distribution funtion
	float* ex;							//The microscopic velocity in X direction
	float* ey;							//The microscopic velocity in Y direction


	float tau;							//The relaxation time (LBM property)
	float nu;							//The kinematic viscosity (LBM property)
	float nu_p;							//The kinematic viscosity (Physical property)
	float Re;							//The Reynolds number (Physical property = LBM property)


	float* Ux0_p;						//The velocity at inlet boundary (Physical property)
	float* Ux0;							//The velocity at inlet boundary (LBM property)
	float rho1;							//The density at outlet boundary (LBM property)

	int* is_boundary_node;
	int* is_solid_node;					//The solid node
	int* is_solid_near_node;
	int i, j, k, in, ip, jn, jp;

	int sIm, sIM;						//The Solid node at Minimum node and at Maximum node in X direction
	int sJm, sJM;						//The Solid node at Minimum node and at Maximum node in Y direction
	int sn;
	int snx, sny;
	float D;							//The diameter of cylinder
	float ic, jc;						//The center node of cylinder
	float dist;							//The distance between center node of cylinder and any node

	float error;
	float sum;

	float sum_Fx1, sum_Fx3, sum_Fx5, sum_Fx6, sum_Fx7, sum_Fx8;
	float sum_Fy2, sum_Fy4, sum_Fy5, sum_Fy6, sum_Fy7, sum_Fy8;
	float sum_Fx, sum_Fy;
	float Fx;							//The total force acting on solid nodes by fluid in X direction
	float Fy;							//The total force acting on solid nodes by fluid in Y direction
	float Cd;							//The coefficient of drag
	float Cl;							//The coefficient of lift

	float rho0, ru;

	char comment[60];

	float c;							//The lattice speed (LBM property)
	float c_s;							//The speed of sound (LBM property)
// ============================================================================ //


// ============================================================================ //
//  LOAD THE PARAMETERS
// ============================================================================ //
	fin_CPU >> nx;				fin_CPU >> comment;
	fin_CPU >> ny;				fin_CPU >> comment;
	fin_CPU >> Lx;				fin_CPU >> comment;
	fin_CPU >> Ly;				fin_CPU >> comment;
	fin_CPU >> a;				fin_CPU >> comment;
	fin_CPU >> rho1;			fin_CPU >> comment;
	fin_CPU >> D;				fin_CPU >> comment;
	fin_CPU >> Um_p;			fin_CPU >> comment;
	fin_CPU >> tau;				fin_CPU >> comment;
	fin_CPU >> nu_p;			fin_CPU >> comment;
// ============================================================================ //


// ============================================================================ //
//  NEW & CUDAMALLOC
// ============================================================================ //
	is_boundary_node = new int[nx*ny];	
	is_solid_node = new int[nx*ny];		
	is_solid_near_node = new int[nx*ny];

	U = new float[nx*ny];				
	Ux = new float[nx*ny];				
	Uy = new float[nx*ny];				
	rho = new float[nx*ny];				

	UN = new float[nx*ny];			
	UxN = new float[nx*ny];				
	UyN = new float[nx*ny];				
	rhoN = new float[nx*ny];			
	f = new float[nx*ny*a];				
	ftemp = new float[nx*ny*a];			
	fN = new float[nx*ny*a];			
	feq = new float[nx*ny*a];			
	ex = new float[a];					
	ey = new float[a];					
	U_p = new float[nx*ny];
	Ux_p = new float[nx*ny];
	Uy_p = new float[nx*ny];
	P = new float[nx*ny];

	Ux0_p = new float[ny];
	Ux0 = new float[ny];				
// ============================================================================ //


// ============================================================================ //
//  MICROSCOPIC VELOCITY
// ============================================================================ //
	ex[0] = 0.0,	ey[0] = 0.0;
	ex[1] = 1.0,	ey[1] = 0.0;
	ex[2] = 0.0,	ey[2] = 1.0;
	ex[3] = -1.0,	ey[3] = 0.0;
	ex[4] = 0.0,	ey[4] = -1.0;
	ex[5] = 1.0,	ey[5] = 1.0;
	ex[6] = -1.0,	ey[6] = 1.0;
	ex[7] = -1.0,	ey[7] = -1.0;
	ex[8] = 1.0,	ey[8] = -1.0;
// ============================================================================ //



// ============================================================================ //
//  SET BOUNDARY NODE
// ============================================================================ //
	sIm = nx / Lx * 0.15;
	sIM = nx / Lx * (0.15 + D) - 1;
	sJm = ny / Ly * 0.15;
	sJM = ny / Ly * (0.15 + D) - 1;

	snx = (sIM - sIm) + 1;
	sny = (sJM - sJm) + 1;

	sn = 0;

	ic = (float)sIm + ((float)sIM - (float)sIm) / 2;
	jc = (float)sJm + ((float)sJM - (float)sJm) / 2;


	//set boundary node
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) is_boundary_node[i + nx*j] = 1;
			else is_boundary_node[i + nx*j] = 0;

		}
	}



	//set solid node
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {

			dist = sqrt(pow((float)i - ic, 2) + pow((float)j - jc, 2));
			if (dist <= ((float)sIM - (float)sIm) / 2) is_solid_node[i + nx*j] = 1;
			else is_solid_node[i + nx*j] = 0;

			if (is_solid_node[i + nx*j]) sn = sn + 1;

		}
	}

	//set near solid node
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {

			is_solid_near_node[i + nx*j] = 0;

			in = i - 1;
			ip = i + 1;
			jn = j - 1;
			jp = j + 1;

			if (!is_boundary_node[i + nx*j]) {
				if (!is_solid_node[i + nx*j]) {


					if (is_solid_node[ip + nx*j]) {
						is_solid_near_node[i + nx*j] = 1;
					}
					else if (is_solid_node[i + nx*jp]) {
						is_solid_near_node[i + nx*j] = 1;
					}
					else if (is_solid_node[in + nx*j]) {
						is_solid_near_node[i + nx*j] = 1;
					}
					else if (is_solid_node[i + nx*jn]) {
						is_solid_near_node[i + nx*j] = 1;
					}
					else if (is_solid_node[ip + nx*jp]) {
						is_solid_near_node[i + nx*j] = 1;
					}
					else if (is_solid_node[in + nx*jp]) {
						is_solid_near_node[i + nx*j] = 1;
					}
					else if (is_solid_node[in + nx*jn]) {
						is_solid_near_node[i + nx*j] = 1;
					}
					else if (is_solid_node[ip + nx*jn]) {
						is_solid_near_node[i + nx*j] = 1;
					}


				}
			}
		}
	}
// ============================================================================ //




// ============================================================================ //
//  INITIAL CONDITION
// ============================================================================ //
	del_x = 1.0;
	del_y = 1.0;
	del_t = 1.0;

	c = del_y / del_t;
	c_s = (1.0 / sqrt(3.0))*c;

	del_x_p = D / (float)snx;
	del_y_p = D / (float)sny;
	//	del_t_p = pow(del_y_p, 2);
	//	del_t_p = 0.000013;



	//Uniform
	/*nu_p = 0.06 * (del_y_p / del_t_p) * D / Re;
	nu = (del_t_p / pow(del_y_p, 2))*nu_p;
	tau = (1.0 / pow(c_s, 2))*nu + (0.5*del_t);*/

	//Input Reynolds number and del_t
	//Um = Um_p * (del_t_p / del_y_p);
	//nu_p = (2.0 / 3.0) * Um_p * D / Re;
	//nu = (del_t_p / pow(del_y_p, 2))*nu_p;
	//tau = (1.0 / pow(c_s, 2))*nu + (0.5*del_t);


	//Input tau and kinematic viscosity
	del_t_p = pow(c_s, 2)*(tau - 0.5)*pow(del_y_p, 2) / nu_p;
	Re = (2.0 / 3.0) * Um_p * D / nu_p;
	Um = Um_p * del_t_p / del_y_p;
	nu = nu_p * del_t_p / pow(del_y_p, 2);




	cout << endl;
	cout << "// =================== Stability condition ================ //" << endl;
	cout << "Check 1. [tau > 0.5]" << endl;
	cout << "tau = " << tau << endl;
	cout << "Check 2. Mach number condition [Ma = Uavg/c_s << 1]" << endl;
	cout << "Ma = " << (2.0 / 3.0)*Um / c_s << endl;
	cout << "Check 3. BGK Stability. [If tau < 0.55, tau > 0.5 + 0.125*Uavg]" << endl;
	cout << "tau = " << tau << " > " << 0.5 + 0.125*(2.0 / 3.0)*Um << endl;
	cout << "// ======================================================== //" << endl;

	//intitalize variables
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {

			Ux[i + nx*j] = 0.0;
			Uy[i + nx*j] = 0.0;
			U[i + nx*j] = 0.0;
			UxN[i + nx*j] = 0.0;
			UyN[i + nx*j] = 0.0;
			UN[i + nx*j] = 0.0;
			P[i + nx*j] = 0.0;

			for (k = 0; k < a; k++) {
				ftemp[i + nx*j + nx*ny*k] = 0.0;
				feq[i + nx*j + nx*ny*k] = 0.0;
				fN[i + nx*j + nx*ny*k] = 0.0;
			}

			if (!is_solid_node[i + nx*j]) rho[i + nx*j] = 1.0;
			else rho[i + nx*j] = 1.0;

			f[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j];
			f[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j];
		}
	}


	//set velocity profile at inlet
	for (j = 0; j < ny; j++) {
		//		Ux0_p[j] = 4.0*Um_p / (pow(Ly, 2))*(((float)j + 1) - 0.5)*del_y_p*(Ly - (((float)j + 1) - 0.5)*del_y_p);
		Ux0[j] = 4.0*Um / (pow(ny, 2))*(((float)j + 1) - 0.5)*del_y*(164 - (((float)j + 1) - 0.5)*del_y);
		//		Ux0[j] = Ux0_p[j] * (del_t_p / del_y_p);

		fout_CPU_Ux0 << Ux0[j] << endl;
	}
// ============================================================================ //





	clock_t t0, t1;
	float e = 0.0;

	int n = 0;
	int n_max = 1;
	int print_step = 200000;
	int print_step2 = 10000;

	cout << endl << "Uniform Flow around Circular Cylinder for LBM using GPU (ver.2)" << endl << endl;
	cout << "// ================== Physical property =================== //" << endl;
	cout << "Physical length in X direction = " << Lx << endl;
	cout << "Physical length in Y direction = " << Ly << endl;
	cout << "Physical grid step in X direction = " << del_x_p << endl;
	cout << "Physical grid step in Y direction = " << del_y_p << endl;
	cout << "Physical time step = " << del_t_p << endl;
	cout << "Physical mean velocity at Inlet boundary = " << (2.0 / 3.0)*Um_p << endl;
	cout << "Physical kinematic viscosity = " << nu_p << endl;
	cout << "Reynolds number = " << (2.0 / 3.0)*Um_p*D / nu_p << endl;
	cout << "// ======================================================== //" << endl << endl;

	cout << "// ===================== LBM property ===================== //" << endl;
	cout << "(LBM) The number of Node in X direction = " << nx << endl;
	cout << "(LBM) The number of Node in Y direction = " << ny << endl;
	cout << "(LBM) The mean velocity at Inlet boundary = " << (2.0 / 3.0)*Um << endl;
	cout << "(LBM) kinematic viscosity = " << nu << endl;
	cout << "Relaxation time = " << tau << endl;
	cout << "Reynolds number = " << (2.0 / 3.0)*Um*sny / nu << endl;
	cout << "The number of node at Solid = " << sn << endl;
	cout << "The number of node of diameter = " << sny << endl;
	cout << "// ======================================================== //" << endl << endl;

	cout << "Print step = " << print_step << endl;
	cout << "The number of maximum iteration = " << n_max << endl;

	t0 = clock();

	while (n < n_max) {



		

		//STREAMING ROUTINE
		// ============================================================================ //
		//  STREAMING
		// ============================================================================ //
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {

				in = i - 1;
				ip = i + 1;
				jn = j - 1;
				jp = j + 1;

				if (!is_boundary_node[i + nx*j]) {
					if (!is_solid_node[i + nx*j]) {

						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];

						if (!is_solid_node[ip + nx*j]) ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
						else ftemp[i + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 1];

						if (!is_solid_node[i + nx*jp]) ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
						else ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];

						if (!is_solid_node[in + nx*j]) ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
						else ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3];

						if (!is_solid_node[i + nx*jn]) ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
						else ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];

						if (!is_solid_node[ip + nx*jp]) ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
						else ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];

						if (!is_solid_node[in + nx*jp]) ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
						else ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];

						if (!is_solid_node[in + nx*jn]) ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
						else ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];

						if (!is_solid_node[ip + nx*jn]) ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
						else ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
					}
				}
				else {
					if ((i == 0) && (j > 0 && j < ny - 1)) {				//INLET
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
						ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
						ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
						ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
						ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
					}
					else if ((i > 0 && i < nx - 1) && (j == ny - 1)) {			//TOP
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
						ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
						ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
						ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
						ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
					}
					else if ((i > 0 && i < nx - 1) && (j == 0)) {				//BOTTOM
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
						ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
						ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
						ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
						ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
					}
					else if ((i == nx - 1) && (j > 0 && j < ny - 1)) {			//OUTLET
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
						ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
						ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
						ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
						ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
					}
					else if ((i == 0) && (j == 0)) {							//BOTTOM-LEFT
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
						ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
						ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
					}
					else if ((i == 0) && (j == ny - 1)) {						//TOP-LEFT
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
						ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
						ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
					}
					else if ((i == nx - 1) && (j == ny - 1)) {					//TOP-RIGHT
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
						ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
						ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
					}
					else if ((i == nx - 1) && (j == 0)) {						//BOTTOM-RIGHT
						ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
						ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
						ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
						ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
					}
				}


			}
		}
		// ============================================================================ //





		//BOUNDARY ROUTINE
		// ============================================================================ //
		//  TOP BOUNDARY (HALF BOUNCEBACK)
		// ============================================================================ //
		j = ny - 1;
		for (i = 1; i < nx - 1; i++) {
			ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
			ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
			ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
		}
		// ============================================================================ //

		// ============================================================================ //
		//	BOTTOM BOUNDARY (HALF BOUNCEBACK)
		// ============================================================================ //
		j = 0;
		for (i = 1; i < nx - 1; i++) {
			ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
			ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
			ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
		}
		// ============================================================================ //

		// ============================================================================ //
		//	LEFT BOUNDARY (VELOCITY)
		// ============================================================================ //
		i = 0;
		for (j = 1; j < ny - 1; j++) {

			/*rho0 = ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
				+ 2.0*(ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 6]);

			ru = rho0 * Ux0[j];

			ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
			ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
			ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;*/

			rho0 = (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
				+ 2.0*(ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7])) / (1.0 - Ux0[j]);
			ru = rho0 * Ux0[j];

			ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
			ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru - (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
			ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru + (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
		}
		// ============================================================================ //

		// ============================================================================ //
		//  RIGHT BOUNDARY (VELOCITY)
		// ============================================================================ //
		i = nx - 1;
		for (j = 1; j < ny - 1; j++) {
			ftemp[i + nx*j + nx*ny * 0] = ftemp[(i - 1) + nx*j + nx*ny * 0];
			ftemp[i + nx*j + nx*ny * 1] = ftemp[(i - 1) + nx*j + nx*ny * 1];
			ftemp[i + nx*j + nx*ny * 2] = ftemp[(i - 1) + nx*j + nx*ny * 2];
			ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*j + nx*ny * 3];
			ftemp[i + nx*j + nx*ny * 4] = ftemp[(i - 1) + nx*j + nx*ny * 4];
			ftemp[i + nx*j + nx*ny * 5] = ftemp[(i - 1) + nx*j + nx*ny * 5];
			ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*j + nx*ny * 6];
			ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*j + nx*ny * 7];
			ftemp[i + nx*j + nx*ny * 8] = ftemp[(i - 1) + nx*j + nx*ny * 8];
		}
		// ============================================================================ //

		// ============================================================================ //
		//	TOP-LEFT CORNER (VELOCITY)
		// ============================================================================ //
		i = 0;
		j = ny - 1;

		/*ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho[i + nx*(j - 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];*/

		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

		// ============================================================================ //

		// ============================================================================ //
		//	BOTTOM-LEFT CORNER (VELOCITY)
		// ============================================================================ //
		i = 0;
		j = 0;

		/*ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho[i + nx*(j + 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];*/

		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

		// ============================================================================ //

		// ============================================================================ //
		//	TOP-RIGHT CORNER (VELOCITY)
		// ============================================================================ //
		i = nx - 1;
		j = ny - 1;

		/*ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho[i + nx*(j - 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];*/

		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

		// ============================================================================ //

		// ============================================================================ //
		//	BOTTOM-RIGHT CORNER (VELOCITY)
		// ============================================================================ //
		i = nx - 1;
		j = 0;

		/*ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8];
		ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho[i + nx*(j + 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];*/

		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

		// ============================================================================ //





		//COLLISION ROUTINE
		//Calculation of Macroscopic var 
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				Ux[i + nx*j] = 0.0;
				Uy[i + nx*j] = 0.0;
				rho[i + nx*j] = 0.0;

				if (!is_solid_node[i + nx*j]) {
					for (k = 0; k < a; k++) {
						rho[i + nx*j] += ftemp[i + nx*j + nx*ny*k];
						Ux[i + nx*j] += ftemp[i + nx*j + nx*ny*k] * ex[k];
						Uy[i + nx*j] += ftemp[i + nx*j + nx*ny*k] * ey[k];
						
					}
					Ux[i + nx*j] /= rho[i + nx*j];
					Uy[i + nx*j] /= rho[i + nx*j];
				}
			}
		}

		

		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (!is_solid_node[i + nx*j]) {


					feq[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + 3.0 * Ux[i + nx*j] + 4.5*pow(Ux[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + 3.0 * Uy[i + nx*j] + 4.5*pow(Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - 3.0 * Ux[i + nx*j] + 4.5*pow(Ux[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - 3.0 * Uy[i + nx*j] + 4.5*pow(Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (Ux[i + nx*j] + Uy[i + nx*j]) + 4.5*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (-Ux[i + nx*j] + Uy[i + nx*j]) + 4.5*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (-Ux[i + nx*j] - Uy[i + nx*j]) + 4.5*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
					feq[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (Ux[i + nx*j] - Uy[i + nx*j]) + 4.5*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

				}
			}
		}


		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {

				if (!is_solid_node[i + nx*j]) {
					for (k = 0; k < a; k++) {
						fN[i + nx*j + nx*ny*k] = ftemp[i + nx*j + nx*ny*k] - (ftemp[i + nx*j + nx*ny*k] - feq[i + nx*j + nx*ny*k]) / tau;
					}
				}
			}
		}



		//ERROR ROUTINE
		n++;

		if (n%print_step == 0) {

			for (i = 0; i < nx; i++) {
				for (j = 0; j < ny; j++) {
					Ux[i + nx*j] = 0.0;
					Uy[i + nx*j] = 0.0;
					rho[i + nx*j] = 0.0;

					if (!is_solid_node[i + nx*j]) {
						for (k = 0; k < a; k++) {

							rho[i + nx*j] += f[i + nx*j + nx*ny*k];
							Ux[i + nx*j] += f[i + nx*j + nx*ny*k] * ex[k];
							Uy[i + nx*j] += f[i + nx*j + nx*ny*k] * ey[k];
							
						}
						Ux[i + nx*j] /= rho[i + nx*j];
						Uy[i + nx*j] /= rho[i + nx*j];
						U[i + nx*j] = sqrt(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2));
					}
				}
			}

			for (i = 0; i < nx; i++) {
				for (j = 0; j < ny; j++) {
					UxN[i + nx*j] = 0.0;
					UyN[i + nx*j] = 0.0;
					rhoN[i + nx*j] = 0.0;

					if (!is_solid_node[i + nx*j]) {
						for (k = 0; k < a; k++) {
							

							rhoN[i + nx*j] += fN[i + nx*j + nx*ny*k];
							UxN[i + nx*j] += fN[i + nx*j + nx*ny*k] * ex[k];
							UyN[i + nx*j] += fN[i + nx*j + nx*ny*k] * ey[k];
						}
						UxN[i + nx*j] /= rhoN[i + nx*j];
						UyN[i + nx*j] /= rhoN[i + nx*j];
						UN[i + nx*j] = sqrt(pow(UxN[i + nx*j], 2) + pow(UyN[i + nx*j], 2));
					}
				}
			}


			sum = 0.0;
			for (i = 0; i < nx; i++) {
				for (j = 0; j < ny; j++) {

					if (!is_solid_node[i + nx*j]) {
						sum = sum + pow(abs(UN[i + nx*j] - U[i + nx*j]), 2);
					}
				}
			}

			error = sqrt(sum / (nx*ny - sn));
			e = error;
			cout << "Iteration : [" << n << "], Error check : [" << e << "]" << endl;
		}








		//UPDATE ROUTINE
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {

				if (!is_solid_node[i + nx*j]) {
					for (k = 0; k < a; k++) {
						f[i + nx*j + nx*ny*k] = fN[i + nx*j + nx*ny*k];
					}
				}
			}
		}



		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				Ux[i + nx*j] = 0.0;
				Uy[i + nx*j] = 0.0;
				rho[i + nx*j] = 0.0;

				for (k = 0; k < a; k++) {

					rho[i + nx*j] += f[i + nx*j + nx*ny*k];
					Ux[i + nx*j] += f[i + nx*j + nx*ny*k] * ex[k];
					Uy[i + nx*j] += f[i + nx*j + nx*ny*k] * ey[k];
				}
				Ux[i + nx*j] /= rho[i + nx*j];
				Uy[i + nx*j] /= rho[i + nx*j];
				U[i + nx*j] = sqrt(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2));
			}
		}


		//Cd, Cl ROUTINE
		sum_Fx1 = 0.0;
		sum_Fx3 = 0.0;
		sum_Fx5 = 0.0;
		sum_Fx6 = 0.0;
		sum_Fx7 = 0.0;
		sum_Fx8 = 0.0;
		sum_Fy2 = 0.0;
		sum_Fy4 = 0.0;
		sum_Fy5 = 0.0;
		sum_Fy6 = 0.0;
		sum_Fy7 = 0.0;
		sum_Fy8 = 0.0;

		sum_Fx = 0.0;
		sum_Fy = 0.0;

		Fx = 0.0;
		Fy = 0.0;

		Cd = 0.0;
		Cl = 0.0;

		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {

				sum_Fx1 = 0.0;
				sum_Fx3 = 0.0;
				sum_Fx5 = 0.0;
				sum_Fx6 = 0.0;
				sum_Fx7 = 0.0;
				sum_Fx8 = 0.0;
				sum_Fy2 = 0.0;
				sum_Fy4 = 0.0;
				sum_Fy5 = 0.0;
				sum_Fy6 = 0.0;
				sum_Fy7 = 0.0;
				sum_Fy8 = 0.0;

				sum_Fx = 0.0;
				sum_Fy = 0.0;

				in = i - 1;
				ip = i + 1;
				jn = j - 1;
				jp = j + 1;

				if (!is_boundary_node[i + nx*j]) {
					if (!is_solid_node[i + nx*j]) {

						if (is_solid_near_node[i + nx*j]) {

							if (is_solid_node[ip + nx*j]) {
								sum_Fx1 = f[i + nx*j + nx*ny * 1] * ex[1] * 2.0;
							}

							if (is_solid_node[i + nx*jp]) {
								sum_Fy2 = f[i + nx*j + nx*ny * 2] * ey[2] * 2.0;
							}

							if (is_solid_node[in + nx*j]) {
								sum_Fx3 = f[i + nx*j + nx*ny * 3] * ex[3] * 2.0;
							}

							if (is_solid_node[i + nx*jn]) {
								sum_Fy4 = f[i + nx*j + nx*ny * 4] * ey[4] * 2.0;
							}

							if (is_solid_node[ip + nx*jp]) {
								sum_Fx5 = f[i + nx*j + nx*ny * 5] * ex[5] * 2.0;
								sum_Fy5 = f[i + nx*j + nx*ny * 5] * ey[5] * 2.0;
							}

							if (is_solid_node[in + nx*jp]) {
								sum_Fx6 = f[i + nx*j + nx*ny * 6] * ex[6] * 2.0;
								sum_Fy6 = f[i + nx*j + nx*ny * 6] * ey[6] * 2.0;
							}

							if (is_solid_node[in + nx*jn]) {
								sum_Fx7 = f[i + nx*j + nx*ny * 7] * ex[7] * 2.0;
								sum_Fy7 = f[i + nx*j + nx*ny * 7] * ey[7] * 2.0;
							}

							if (is_solid_node[ip + nx*jn]) {
								sum_Fx8 = f[i + nx*j + nx*ny * 8] * ex[8] * 2.0;
								sum_Fy8 = f[i + nx*j + nx*ny * 8] * ey[8] * 2.0;
							}


							sum_Fx = sum_Fx1 + sum_Fx3 + sum_Fx5 + sum_Fx6 + sum_Fx7 + sum_Fx8;
							sum_Fy = sum_Fy2 + sum_Fy4 + sum_Fy5 + sum_Fy6 + sum_Fy7 + sum_Fy8;


							Fx = Fx + sum_Fx;
							Fy = Fy + sum_Fy;


						}

					}
				}
			}
		}

		Cd = 2.0*Fx / (rho1*pow((2.0 / 3.0)*Um, 2)*snx);
		Cl = 2.0*Fy / (rho1*pow((2.0 / 3.0)*Um, 2)*sny);

		fout_CPU_Cd << Cd << "\t" << Cl << endl;
		


	}
	t1 = clock();




	//PRINT ROUTINE
	// ============================================================================ //
	//  CHANGE LBM -> PHYSICAL
	// ============================================================================ //
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			Ux_p[i + nx*j] = Ux[i + nx*j];
			Uy_p[i + nx*j] = Uy[i + nx*j];
			U_p[i + nx*j] = U[i + nx*j];
			P[i + nx*j] = rho[i + nx*j] / (3.0);
		}
	}
	// ============================================================================ //


	fout_CPU << endl;
	fout_CPU << "variables = X Y Ux Uy U rho P" << endl;
	fout_CPU << "zone i=" << nx << " j=" << ny << endl;
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			fout_CPU << i << "\t" << j << "\t" << Ux_p[i + nx*j] << "\t" << Uy_p[i + nx*j] << "\t"
				<< U_p[i + nx*j] << "\t" << rho[i + nx*j] << "\t" << P[i + nx*j] << endl;
		}
	}
	fout_CPU << endl;






	cout << endl;
	cout << "The force in X direction : " << Fx << endl;
	cout << "The force in Y direction : " << Fy << endl;
	cout << "The coefficient of Drag : " << Cd << endl;
	cout << "The coefficient of lift : " << Cl << endl;
	cout << "Computation time : " << double(t1 - t0) / CLOCKS_PER_SEC << "[s]" << endl;
	cout << "Total Iteration : [" << n << "], Error check : [" << e << "]" << endl;

	delete[] Ux0;
	delete[] Ux0_p;
	delete[] P;
	delete[] Uy_p;
	delete[] Ux_p;
	delete[] U_p;
	delete[] ey;
	delete[] ex;
	delete[] fN;
	delete[] feq;
	delete[] ftemp;
	delete[] f;
	delete[] rhoN;
	delete[] UyN;
	delete[] UxN;
	delete[] UN;
	delete[] rho;
	delete[] Uy;
	delete[] Ux;
	delete[] U;
	delete[] is_boundary_node;
	delete[] is_solid_node;
	delete[] is_solid_near_node;
	cout << endl << "Done!" << endl;
}
// ============================================================================ //