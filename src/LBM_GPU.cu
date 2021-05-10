#include "LBM_GPU.cuh"
ifstream fin_GPU("in_GPU.txt");

ofstream fout_GPU("out_GPU.dat");
ofstream fout_GPU_Cd("out_GPU_Cd.dat");
ofstream fout_GPU_Ux0("out_GPU_Ux0.dat");
ofstream fout_GPU_Ux("out_GPU_Ux.dat");

LBM_GPU::LBM_GPU()
{
// ============================================================================ //
//  LOAD THE PARAMETERS
// ============================================================================ //
	fin_GPU >> nx;				fin_GPU >> comment;
	fin_GPU >> ny;				fin_GPU >> comment;
	fin_GPU >> Lx;				fin_GPU >> comment;
	fin_GPU >> Ly;				fin_GPU >> comment;
	fin_GPU >> a;				fin_GPU >> comment;
	fin_GPU >> rho1;			fin_GPU >> comment;
	fin_GPU >> BLOCK_SIZE_X;	fin_GPU >> comment;
	fin_GPU >> BLOCK_SIZE_Y;	fin_GPU >> comment;
	fin_GPU >> BLOCK_SIZE_Z;	fin_GPU >> comment;
	fin_GPU >> D;				fin_GPU >> comment;
	fin_GPU >> Um_p;			fin_GPU >> comment;
	fin_GPU >> tau;				fin_GPU >> comment;
	fin_GPU >> nu_p;			fin_GPU >> comment;
// ============================================================================ //

	
// ============================================================================ //
//  NEW & CUDAMALLOC
// ============================================================================ //
	is_boundary_node = new int[nx*ny];	cudaMalloc((void**)&d_is_boundary_node, nx*ny * sizeof(int));
	is_solid_node = new int[nx*ny];		cudaMalloc((void**)&d_is_solid_node, nx*ny * sizeof(int));
	is_solid_near_node = new int[nx*ny];

	U = new float[nx*ny];				cudaMalloc((void**)&d_U, nx*ny * sizeof(float));
	Ux = new float[nx*ny];				cudaMalloc((void**)&d_Ux, nx*ny * sizeof(float));
	Uy = new float[nx*ny];				cudaMalloc((void**)&d_Uy, nx*ny * sizeof(float));
	rho = new float[nx*ny];				cudaMalloc((void**)&d_rho, nx*ny * sizeof(float));

	UN = new float[nx*ny];				cudaMalloc((void**)&d_UN, nx*ny * sizeof(float));
	UxN = new float[nx*ny];				cudaMalloc((void**)&d_UxN, nx*ny * sizeof(float));
	UyN = new float[nx*ny];				cudaMalloc((void**)&d_UyN, nx*ny * sizeof(float));
	rhoN = new float[nx*ny];			cudaMalloc((void**)&d_rhoN, nx*ny * sizeof(float));
	f = new float[nx*ny*a];				cudaMalloc((void**)&d_f, nx*ny*a * sizeof(float));
	ftemp = new float[nx*ny*a];			cudaMalloc((void**)&d_ftemp, nx*ny*a * sizeof(float));
	fN = new float[nx*ny*a];			cudaMalloc((void**)&d_fN, nx*ny*a * sizeof(float));
	feq = new float[nx*ny*a];			cudaMalloc((void**)&d_feq, nx*ny*a * sizeof(float));
	ex = new float[a];					cudaMalloc((void**)&d_ex, a * sizeof(float));
	ey = new float[a];					cudaMalloc((void**)&d_ey, a * sizeof(float));
	U_p = new float[nx*ny];
	Ux_p = new float[nx*ny];
	Uy_p = new float[nx*ny];
	P = new float[nx*ny];
	
	Ux0_p = new float[ny];
	Ux0 = new float[ny];				cudaMalloc((void**)&d_Ux0, ny * sizeof(float));
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
	cudaMemcpy(d_ex, ex, a * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ey, ey, a * sizeof(float), cudaMemcpyHostToDevice);
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

	r = ((float)sIM - (float)sIm) / 2;

	cout << "sIm = " << sIm << endl;
	cout << "sIM = " << sIM << endl;
	cout << "sJm = " << sJm << endl;
	cout << "sJM = " << sJM << endl;
	cout << "snx = " << snx << endl;
	cout << "sny = " << sny << endl;
	cout << "ic = " << ic << endl;
	cout << "jc = " << jc << endl;
	cout << "r = " << r << endl;

	//set boundary node
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) is_boundary_node[i + nx*j] = 1;
			else is_boundary_node[i + nx*j] = 0;

		}
	}


	//Binary data
	/*for (i = 0; i < nx; i++) {
		for (j = ny - 1; j > -1; j--) {

			if ((i >= sIm && i <= sIM) && (j >= sJm && j <= sJM)) fin_grid_GPU >> is_solid_node[i + nx*j];
			else is_solid_node[i + nx*j] = 0;

			if (is_solid_node[i + nx*j]) sn = sn + 1;
		}
	}*/

	//set solid node
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {

			dist = sqrt(pow((float)i - ic, 2) + pow((float)j - jc, 2));
			if (dist <= r) is_solid_node[i + nx*j] = 1;
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


	cudaMemcpy(d_is_boundary_node, is_boundary_node, nx*ny * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_is_solid_node, is_solid_node, nx*ny * sizeof(int), cudaMemcpyHostToDevice);
// ============================================================================ //




// ============================================================================ //
//  SET PARAMETERS & INITIAL CONDITION
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
	cout << "Ma = " << (2.0/3.0)*Um/c_s << endl;
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

	cudaMemcpy(d_rho, rho, nx*ny * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_f, f, nx*ny*a * sizeof(float), cudaMemcpyHostToDevice);



	//set velocity profile at inlet
	for (j = 0; j < ny; j++) {
		Ux0_p[j] = 4.0*Um_p / (pow(Ly, 2))*(((float)j + 1) - 0.5)*del_y_p*(Ly - (((float)j + 1) - 0.5)*del_y_p);
//		Ux0[j] = 4.0*Um / (pow(ny, 2))*(((float)j + 1) - 0.5)*del_y*(ny - (((float)j + 1) - 0.5)*del_y);
		Ux0[j] = Ux0_p[j] * (del_t_p / del_y_p);

		fout_GPU_Ux0 << Ux0[j] << endl;
	}
	cudaMemcpy(d_Ux0, Ux0, ny * sizeof(float), cudaMemcpyHostToDevice);
// ============================================================================ //
}

__global__ 
void Kernel_Streaming(float* f, float* ftemp, int* is_boundary_node, int* is_solid_node, int nx, int ny, int a, float ic, float jc, float r) {

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.z * blockIdx.z + threadIdx.z;
	if (i >= nx || j >= ny || k >= a) return;

	int in, ip, jn, jp;

	in = i - 1;
	ip = i + 1;
	jn = j - 1;
	jp = j + 1;

	float dist = sqrt(pow((float)i - ic, 2) + pow((float)j - jc, 2)); 
	float q = dist - r;


	if (!is_boundary_node[i + nx*j]) {
		if (!is_solid_node[i + nx*j]) {
			
			ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];


			if (!is_solid_node[ip + nx*j]) ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 3] = 2.0 * q * f[i + nx*j + nx*ny * 1] + (1.0 - 2.0*q)*f[(i - 1) + nx*j + nx*ny * 1];
				else ftemp[i + nx*j + nx*ny * 3] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 1] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 3];
			}

			if (!is_solid_node[i + nx*jp]) ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 4] = 2.0 * q * f[i + nx*j + nx*ny * 2] + (1.0 - 2.0*q)*f[i + nx*(j - 1) + nx*ny * 2];
				else ftemp[i + nx*j + nx*ny * 4] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 2] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 4];
			}

			if (!is_solid_node[in + nx*j]) ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 1] = 2.0 * q * f[i + nx*j + nx*ny * 3] + (1.0 - 2.0*q)*f[(i + 1) + nx*j + nx*ny * 3];
				else ftemp[i + nx*j + nx*ny * 1] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 3] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 1];
			}

			if (!is_solid_node[i + nx*jn]) ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 2] = 2.0 * q * f[i + nx*j + nx*ny * 4] + (1.0 - 2.0*q)*f[i + nx*(j + 1) + nx*ny * 4];
				else ftemp[i + nx*j + nx*ny * 2] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 4] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 2];
			}

			if (!is_solid_node[ip + nx*jp]) ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 7] = 2.0 * q * f[i + nx*j + nx*ny * 5] + (1.0 - 2.0*q)*f[(i - 1) + nx*(j - 1) + nx*ny * 5];
				else ftemp[i + nx*j + nx*ny * 7] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 5] + (2.0*q - 1) / (2.0*q)*f[i + nx*j + nx*ny * 7];
			}

			if (!is_solid_node[in + nx*jp]) ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 8] = 2.0 * q * f[i + nx*j + nx*ny * 6] + (1.0 - 2.0*q)*f[(i + 1) + nx*(j - 1) + nx*ny * 6];
				else ftemp[i + nx*j + nx*ny * 8] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 6] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 8];
			}

			if (!is_solid_node[in + nx*jn]) ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 5] = 2.0 * q * f[i + nx*j + nx*ny * 7] + (1.0 - 2.0*q)*f[(i + 1) + nx*(j + 1) + nx*ny * 7];
				else ftemp[i + nx*j + nx*ny * 5] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 7] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 5];
			}

			if (!is_solid_node[ip + nx*jn]) ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
			else {
				if (q < 0.5) ftemp[i + nx*j + nx*ny * 6] = 2.0 * q * f[i + nx*j + nx*ny * 8] + (1.0 - 2.0*q)*f[(i - 1) + nx*(j + 1) + nx*ny * 8];
				else ftemp[i + nx*j + nx*ny * 6] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 8] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 6];
			}

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
void LBM_GPU::Streaming() {

	dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 dimGrid((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y, (a + BLOCK_SIZE_Z - 1) / BLOCK_SIZE_Z);
	Kernel_Streaming << < dimGrid, dimBlock >> > (d_f, d_ftemp, d_is_boundary_node, d_is_solid_node, nx, ny, a, ic, jc, r);
}

__global__ 
void Kernel_BC_bounceback(float* f, float* ftemp, float* rho, float* Ux, float* Uy, float* Ux0, float rho1, int nx, int ny, int a) {

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.z * blockIdx.z + threadIdx.z;
	if (i >= nx || j >= ny || k >= a) return;

	float rho0, ru, Ux1, Uy1, rho_extra, Ux_extra, Uy_extra;
// ============================================================================ //
//  TOP BOUNDARY (HALF-AWAY BOUNCEBACK)
// ============================================================================ //
	if ((i > 0 && i < nx - 1) && (j == ny - 1)){

		//Bounce-back boundary
		ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];

		//Periodic boundary
		/*ftemp[i + nx*0 + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
		ftemp[(i + 1) + nx*0 + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
		ftemp[(i - 1) + nx*0 + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
*/

		//Velocity boundary(first order)
		/*rho_extra = rho[i + nx*(j - 1)] + 0.5 * (rho[i + nx*(j - 1)] - rho[i + nx*(j - 2)]);
		Ux_extra = Ux[i + nx*(j - 1)];
		ru = rho_extra*Ux_extra;

		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] - (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;*/


		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[i + nx*(j - 1) + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*(j - 1) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*(j - 1) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*(j - 1) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*(j - 1) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*(j - 1) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*(j - 1) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*(j - 1) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*(j - 1) + nx*ny * 8];
*/

		//Extrapolation high order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[i + nx*(j - 1) + nx*ny * 0] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 0] - ftemp[i + nx*(j - 2) + nx*ny * 0]);
		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*(j - 1) + nx*ny * 1] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 1] - ftemp[i + nx*(j - 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*(j - 1) + nx*ny * 2] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 2] - ftemp[i + nx*(j - 2) + nx*ny * 2]);
		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*(j - 1) + nx*ny * 3] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 3] - ftemp[i + nx*(j - 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*(j - 1) + nx*ny * 4] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 4] - ftemp[i + nx*(j - 2) + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*(j - 1) + nx*ny * 5] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 5] - ftemp[i + nx*(j - 2) + nx*ny * 5]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*(j - 1) + nx*ny * 6] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 6] - ftemp[i + nx*(j - 2) + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*(j - 1) + nx*ny * 7] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 7] - ftemp[i + nx*(j - 2) + nx*ny * 7]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*(j - 1) + nx*ny * 8] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 8] - ftemp[i + nx*(j - 2) + nx*ny * 8]);
		*/

		//Extrapolation 2nd order
		/*ftemp[i + nx*j + nx*ny * 4] = 2.0 * ftemp[i + nx*(j - 1) + nx*ny * 4] - ftemp[i + nx*(j - 2) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 7] = 2.0 * ftemp[i + nx*(j - 1) + nx*ny * 7] - ftemp[i + nx*(j - 2) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = 2.0 * ftemp[i + nx*(j - 1) + nx*ny * 8] - ftemp[i + nx*(j - 2) + nx*ny * 8];
	*/


		//Equilibrium
		/*float c = 1;
		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
*/

		//NEBB method
		/*ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] + (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 1] - ftemp[i + nx*j + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] - (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 1] - ftemp[i + nx*j + nx*ny * 3]);
		*/

	}
// ============================================================================ //


// ============================================================================ //
//	BOTTOM BOUNDARY (HALF-AWAY BOUNCEBACK)
// ============================================================================ //
	if ((i > 0 && i < nx - 1) && (j == 0)){

		//Bounce-back boundary
		ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];


		//Periodic boundary
		/*ftemp[i + nx*(ny - 1) + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
		ftemp[(i - 1) + nx*(ny - 1) + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
		ftemp[(i + 1) + nx*(ny - 1) + nx*ny * 8] = f[i + nx*j + nx*ny * 8];*/


		//Velocity boundary(first order)
		/*rho_extra = rho[i + nx*(j + 1)] + 0.5 * (rho[i + nx*(j + 1)] - rho[i + nx*(j + 2)]);
		Ux_extra = Ux[i + nx*(j + 1)];
		ru = rho_extra*Ux_extra;

		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 6.0)*ru;*/


		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[i + nx*(j + 1) + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*(j + 1) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*(j + 1) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*(j + 1) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*(j + 1) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*(j + 1) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*(j + 1) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*(j + 1) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*(j + 1) + nx*ny * 8];
	*/

		//Extrapolation high order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[i + nx*(j + 1) + nx*ny * 0] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 0] - ftemp[i + nx*(j + 2) + nx*ny * 0]);
		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*(j + 1) + nx*ny * 1] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 1] - ftemp[i + nx*(j + 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*(j + 1) + nx*ny * 2] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 2] - ftemp[i + nx*(j + 2) + nx*ny * 2]);
		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*(j + 1) + nx*ny * 3] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 3] - ftemp[i + nx*(j + 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*(j + 1) + nx*ny * 4] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 4] - ftemp[i + nx*(j + 2) + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*(j + 1) + nx*ny * 5] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 5] - ftemp[i + nx*(j + 2) + nx*ny * 5]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*(j + 1) + nx*ny * 6] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 6] - ftemp[i + nx*(j + 2) + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*(j + 1) + nx*ny * 7] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 7] - ftemp[i + nx*(j + 2) + nx*ny * 7]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*(j + 1) + nx*ny * 8] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 8] - ftemp[i + nx*(j + 2) + nx*ny * 8]);
*/

		//Extrapolation 2nd order
		/*ftemp[i + nx*j + nx*ny * 2] = 2.0 * ftemp[i + nx*(j + 1) + nx*ny * 2] - ftemp[i + nx*(j + 2) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 5] = 2.0 * ftemp[i + nx*(j + 1) + nx*ny * 5] - ftemp[i + nx*(j + 2) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = 2.0 * ftemp[i + nx*(j + 1) + nx*ny * 6] - ftemp[i + nx*(j + 2) + nx*ny * 6];
	*/


		//Equilibrium
		/*float c = 1;
		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	*/

		//NEBB method
		/*ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 3] - ftemp[i + nx*j + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 3] - ftemp[i + nx*j + nx*ny * 1]);
*/

	}
// ============================================================================ //


// ============================================================================ //
//	LEFT BOUNDARY (VELOCITY)
// ============================================================================ //
	if ((i == 0) && (j > 0 && j < ny - 1)) {
		/*rho0 = rho[(i + 1) + nx*j] + 0.5*(rho[(i + 1) + nx*j] - rho[(i + 2) + nx*j]);
		ru = rho0 * Ux0;

		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;*/


		/*rho0 = ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
			+ 2.0*(ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 6]);
		

		ru = rho0 * Ux0[j];

		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;*/




		//Zou - He boundary
		rho0 = (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
			+ 2.0*(ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7])) / (1.0 - Ux0[j]);
		ru = rho0 * Ux0[j];

		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru - (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru + (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);


		//wet-node method
//		rho0 = (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
//			+ 2.0*(ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7])) / (1.0 - Ux0[j]);
//		ru = rho0 * Ux0[j];
////		ru = rho0 * 0.06;
//
//		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
//		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru - (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
//		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru + (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
		
		
		//Equilibrium
		/*float c = 1;
		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux0[j] + (4.5 / pow(c, 4)) * pow(Ux0[j], 2) - (1.5 / pow(c, 2)) * (pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux0[j] + (4.5 / pow(c, 4))*pow(Ux0[j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux0[j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux0[j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux0[j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux0[j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux0[j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux0[j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux0[j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux0[j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
	*/
	}
// ============================================================================ //


// ============================================================================ //
//  RIGHT BOUNDARY (EXTRAPOLATION)
// ============================================================================ //
	if ((i == nx - 1) && (j > 0 && j < ny - 1)) {

		//Extrapolation
//		ftemp[i + nx*j + nx*ny * 0] = ftemp[(i - 1) + nx*j + nx*ny * 0] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 0] - ftemp[(i - 2) + nx*j + nx*ny * 0]);
//		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i - 1) + nx*j + nx*ny * 1] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 1] - ftemp[(i - 2) + nx*j + nx*ny * 1]);
//		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i - 1) + nx*j + nx*ny * 2] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 2] - ftemp[(i - 2) + nx*j + nx*ny * 2]);
//		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*j + nx*ny * 3] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 3] - ftemp[(i - 2) + nx*j + nx*ny * 3]);
//		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i - 1) + nx*j + nx*ny * 4] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 4] - ftemp[(i - 2) + nx*j + nx*ny * 4]);
//		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i - 1) + nx*j + nx*ny * 5] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 5] - ftemp[(i - 2) + nx*j + nx*ny * 5]);
//		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*j + nx*ny * 6] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 6] - ftemp[(i - 2) + nx*j + nx*ny * 6]);
//		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*j + nx*ny * 7] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 7] - ftemp[(i - 2) + nx*j + nx*ny * 7]);
//		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i - 1) + nx*j + nx*ny * 8] + 0.5 * (ftemp[(i - 1) + nx*j + nx*ny * 8] - ftemp[(i - 2) + nx*j + nx*ny * 8]);
		

		//Extrapolation high order
		/*ftemp[i + nx*j + nx*ny * 0] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 0] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 0] + ftemp[(i - 3) + nx*j + nx*ny * 0]);
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 1] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 1] + ftemp[(i - 3) + nx*j + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 2] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 2] + ftemp[(i - 3) + nx*j + nx*ny * 2]);
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 3] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 3] + ftemp[(i - 3) + nx*j + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 4] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 4] + ftemp[(i - 3) + nx*j + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 5] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 5] + ftemp[(i - 3) + nx*j + nx*ny * 5]);
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 6] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 6] + ftemp[(i - 3) + nx*j + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 7] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 7] + ftemp[(i - 3) + nx*j + nx*ny * 7]);
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 8] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 8] + ftemp[(i - 3) + nx*j + nx*ny * 8]);
	*/

		//Extrapolation type2
		/*ftemp[i + nx*j + nx*ny * 3] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 3] - ftemp[(i - 2) + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 6] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 6] - ftemp[(i - 2) + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 7] - ftemp[(i - 2) + nx*j + nx*ny * 7];*/


		//Extrapolation first order
		ftemp[i + nx*j + nx*ny * 0] = ftemp[(i - 1) + nx*j + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i - 1) + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i - 1) + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i - 1) + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i - 1) + nx*j + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i - 1) + nx*j + nx*ny * 8];


		//Extrapolation second order
		/*ftemp[i + nx*j + nx*ny * 0] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 0] - ftemp[(i - 2) + nx*j + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 1] - ftemp[(i - 2) + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 2] - ftemp[(i - 2) + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 3] - ftemp[(i - 2) + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 4] - ftemp[(i - 2) + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 5] - ftemp[(i - 2) + nx*j + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 6] - ftemp[(i - 2) + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 7] - ftemp[(i - 2) + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 8] - ftemp[(i - 2) + nx*j + nx*ny * 8];*/



		//Velocity boundary (first order)
		/*rho_extra = rho[(i - 1) + nx*j] + 0.5*(rho[(i - 1) + nx*j] - rho[(i - 2) + nx*j]);
		Ux_extra = Ux[(i - 1) + nx*j];
		Uy_extra = Uy[(i - 1) + nx*j];

		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1] - (2.0 / 3.0)*rho_extra*Ux_extra;
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 6.0)*rho_extra*(Ux_extra - Uy_extra);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] - (1.0 / 6.0)*rho_extra*(Ux_extra + Uy_extra);*/


		//Pressure boundary
		/*Ux1 = Ux[(i - 1) + nx*j] + 0.5*(Ux[(i - 1) + nx*j] - Ux[(i - 2) + nx*j]);
		Uy1 = Uy[(i - 1) + nx*j] + 0.5*(Uy[(i - 1) + nx*j] - Uy[(i - 2) + nx*j]);

		ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 6] = -f[i + nx*j + nx*ny * 8] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 - Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 7] = -f[i + nx*j + nx*ny * 5] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 + Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
*/

		//wet-node method
		/*Ux1 = -1.0 + (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
			+ 2.0*(ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 8])) / rho1;
		ru = rho1 * Ux1;

		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1] - (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 6.0)*ru -(1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] - (1.0 / 6.0)*ru +(1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
*/
	}
// ============================================================================ //


// ============================================================================ //
//	TOP-LEFT CORNER (EQUILIBRIUM)
// ============================================================================ //
	if ((i == 0) && (j == ny - 1)) {
		

		//case 1
	//	ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
	//	ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];

	////	rho0 = rho[(i + 1) + nx*(j - 1)] + 0.5*(rho[(i + 1) + nx*(j - 1)] - rho[(i + 2) + nx*(j - 2)]);
	//	rho0 = rho[i + nx*(j - 1)];
	//	ru = rho0 * Ux0[j];

	//	ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
	//	ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
	//	ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;




		//case 2
		//ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		//ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
		//ftemp[i + nx*j + nx*ny * 5] = -ftemp[i + nx*j + nx*ny * 7];

		//rho0 = rho[(i + 1) + nx*(j - 1)] + 0.5*(rho[(i + 1) + nx*(j - 1)] - rho[(i + 2) + nx*(j - 2)]);
		//ru = rho0 * Ux0;

		//ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		//ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;

		//ftemp[i + nx*j + nx*ny * 0] = rho0 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
		//	+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);


		//case 3
		/*ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];*/

	

		//case 4
		/*ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];

		rho0 = rho[(i + 0) + nx*(j - 1)] + 0.5*(rho[(i + 0) + nx*(j - 1)] - rho[(i + 0) + nx*(j - 2)]);
		ru = rho0 * Ux0[j];

		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		*/
	
		//case 5
////		rho0 = rho[(i + 0) + nx*(j - 1)] + 0.5*(rho[(i + 0) + nx*(j - 1)] - rho[(i + 0) + nx*(j - 2)]);
//		rho0 = rho[i + nx*(j - 1)];
////		rho0 = 1.0;
//		ru = rho0 * Ux0[j];
////		ru = rho0* 0.005;
//
//		ftemp[i + nx*j + nx*ny * 7] = -(1.0 / 12.0) * ru;
//		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 12.0) * ru;
//
//		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
//		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
//		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;
//
//		ftemp[i + nx*j + nx*ny * 0] = rho0 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
//			+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);


		//Periodic + Velocity
		/*ftemp[i + nx*0 + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
		ftemp[(i + 1) + nx*0 + nx*ny * 5] = f[i + nx*j + nx*ny * 5];

		rho0 = rho[(i + 1) + nx*(j - 1)] + 0.5*(rho[(i + 1) + nx*(j - 1)] - rho[(i + 2) + nx*(j - 2)]);
		ru = rho0 * Ux0;

		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;*/


		//wet-node method
//		rho0 = rho[i + nx*(j - 1)];
////		rho0 = 1.001;
////		rho0 = rho[(i + 0) + nx*(j - 1)] + 0.5*(rho[(i + 0) + nx*(j - 1)] - rho[(i + 0) + nx*(j - 2)]);
////		rho0 = rho1;
//
//		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3];
//		ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
//		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
//		ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho0 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
//			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
//		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];


		//Equilibrium
		float c = 1;
		/*ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux0[j] + (4.5 / pow(c, 4)) * pow(Ux0[j], 2) - (1.5 / pow(c, 2)) * (pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux0[j] + (4.5 / pow(c, 4))*pow(Ux0[j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux0[j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux0[j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux0[j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux0[j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux0[j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux0[j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux0[j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux0[j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
	*/	

		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

	}
// ============================================================================ //


// ============================================================================ //
//	BOTTOM-LEFT CORNER (EQUILIBRIUM)
// ============================================================================ //
	if ((i == 0) && (j == 0)) {
		
		//case 1
	//	ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
	//	ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];

	////	rho0 = rho[(i + 1) + nx*(j + 1)] + 0.5*(rho[(i + 1) + nx*(j + 1)] - rho[(i + 2) + nx*(j + 2)]);
	//	rho0 = rho[i + nx*(j + 1)];
	//	ru = rho0 * Ux0[j];

	//	ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
	//	ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
	//	ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;


		//case 2
		//ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		//ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
		//ftemp[i + nx*j + nx*ny * 8] = -ftemp[i + nx*j + nx*ny * 6];

		//rho0 = rho[(i + 1) + nx*(j + 1)] + 0.5*(rho[(i + 1) + nx*(j + 1)] - rho[(i + 2) + nx*(j + 2)]);
		//ru = rho0 * Ux0;

		//ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		//ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		//
		//ftemp[i + nx*j + nx*ny * 0] = rho0 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
		//	+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);



		//case 3
		/*ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
*/

		//case 4
		/*ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];

		rho0 = rho[(i + 0) + nx*(j + 1)] + 0.5*(rho[(i + 0) + nx*(j + 1)] - rho[(i + 0) + nx*(j + 2)]);
		ru = rho0 * Ux0[j];

		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;
*/

		//case 5		
////		rho0 = rho[(i + 0) + nx*(j + 1)] + 0.5*(rho[(i + 0) + nx*(j + 1)] - rho[(i + 0) + nx*(j + 2)]);
//		rho0 = rho[i + nx*(j + 1)];
////		rho0 = 1.0;
//		ru = rho0 * Ux0[j];
////		ru = rho0* 0.005;
//
//		ftemp[i + nx*j + nx*ny * 6] = -(1.0 / 12.0) * ru;
//		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 12.0) * ru;
//
//		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
//		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
//		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
//		
//		ftemp[i + nx*j + nx*ny * 0] = rho0 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
//			+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);



		//Periodic + Velocity
		/*ftemp[i + nx*(ny - 1) + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
		ftemp[(i + 1) + nx*(ny - 1) + nx*ny * 8] = f[i + nx*j + nx*ny * 8];

		rho0 = rho[(i + 1) + nx*(j + 1)] + 0.5*(rho[(i + 1) + nx*(j + 1)] - rho[(i + 2) + nx*(j + 2)]);
		ru = rho0 * Ux0;

		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;
*/


		//wet-node method
//		rho0 = rho[i + nx*(j + 1)];
////		rho0 = 1.001;
////		rho0 = rho[(i + 0) + nx*(j + 1)] + 0.5*(rho[(i + 0) + nx*(j + 1)] - rho[(i + 0) + nx*(j + 2)]);
////		rho0 = rho1;
//
//		ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3];
//		ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
//		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
//		ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho0 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
//			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
//		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];


		//Equilibrium
		float c = 1;
		/*ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux0[j] + (4.5 / pow(c, 4)) * pow(Ux0[j], 2) - (1.5 / pow(c, 2)) * (pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux0[j] + (4.5 / pow(c, 4))*pow(Ux0[j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux0[j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux0[j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux0[j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux0[j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux0[j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux0[j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux0[j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux0[j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux0[j], 2) + pow(Uy[i + nx*j], 2)));
*/

		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

	}
// ============================================================================ //

	
// ============================================================================ //
//	TOP-RIGHT CORNER (EQUILIBRIUM)
// ============================================================================ //
	if ((i == nx - 1) && (j == ny - 1)) {

		//case 1
		/*ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];

		Ux1 = Ux[(i - 1) + nx*(j - 1)] + 0.5*(Ux[(i - 1) + nx*(j - 1)] - Ux[(i - 2) + nx*(j - 2)]);
		Uy1 = Uy[(i - 1) + nx*(j - 1)] + 0.5*(Uy[(i - 1) + nx*(j - 1)] - Uy[(i - 2) + nx*(j - 2)]);

		ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 6] = -f[i + nx*j + nx*ny * 8] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 - Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 7] = -f[i + nx*j + nx*ny * 5] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 + Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));*/


		

		//case 2
		//ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		//ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
		//ftemp[i + nx*j + nx*ny * 6] = -ftemp[i + nx*j + nx*ny * 8];

		//Ux1 = Ux[(i - 1) + nx*(j - 1)] + 0.5*(Ux[(i - 1) + nx*(j - 1)] - Ux[(i - 2) + nx*(j - 2)]);

		//ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);
		//ftemp[i + nx*j + nx*ny * 7] = -f[i + nx*j + nx*ny * 5] + (1.0 / 18.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);

		//ftemp[i + nx*j + nx*ny * 0] = rho1 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
		//	+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);


		//case 3
		/*ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];*/


		//case 4
		/*ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];

		Ux1 = Ux[(i - 1) + nx*(j - 1)] + 0.5*(Ux[(i - 1) + nx*(j - 1)] - Ux[(i - 2) + nx*(j - 2)]);

		ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);
		ftemp[i + nx*j + nx*ny * 6] = -f[i + nx*j + nx*ny * 8] + (1.0 / 18.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);
		*/


		//case 5
//		Ux1 = Ux[(i - 1) + nx*(j - 0)];
////		rho_extra = rho[(i - 1) + nx*(j - 1)] + 0.5*(rho[(i - 1) + nx*(j - 1)] - rho[(i - 2) + nx*(j - 2)]);
//		rho_extra = rho[(i - 1) + nx*(j - 0)];
//		ru = rho_extra * Ux1;
//		
//		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 12.0) * ru;
//		ftemp[i + nx*j + nx*ny * 6] = -(1.0 / 12.0) * ru;
//
//		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
//		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1] - (2.0 / 3.0) * ru;
//		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] - (1.0 / 6.0) * ru;
//
//		ftemp[i + nx*j + nx*ny * 0] = rho1 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
//			+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);


		//Periodic + Bounce back
		/*ftemp[i + nx*0 + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
		ftemp[(i - 1) + nx*0 + nx*ny * 6] = f[i + nx*j + nx*ny * 6];

		Ux1 = Ux[(i - 1) + nx*(j - 1)];
		Uy1 = Uy[(i - 1) + nx*(j - 1)];

		ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 6] = -f[i + nx*j + nx*ny * 8] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 - Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 7] = -f[i + nx*j + nx*ny * 5] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 + Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
*/


		//Periodic + Extrapolation
		/*ftemp[i + nx * 0 + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
		ftemp[(i - 1) + nx * 0 + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
*/


		//wet-node method
////		rho0 = rho[(i - 1) + nx*(j - 0)];
//		rho0 = rho1;
////		rho0 = rho[(i - 0) + nx*(j - 1)] + 0.5*(rho[(i - 0) + nx*(j - 1)] - rho[(i - 0) + nx*(j - 2)]);
//
//		ftemp[i + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 1];
//		ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
//		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
//		ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho0 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
//			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
//		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];


		//Equilibrium
		float c = 1;
		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

	}
// ============================================================================ //


// ============================================================================ //
//	BOTTOM-RIGHT CORNER (EQUILIBRIUM)
// ============================================================================ //
	if ((i == nx - 1) && (j == 0)) {

		//case 1
		/*ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];

		Ux1 = Ux[(i - 1) + nx*(j + 1)] + 0.5*(Ux[(i - 1) + nx*(j + 1)] - Ux[(i - 2) + nx*(j + 2)]);
		Uy1 = Uy[(i - 1) + nx*(j + 1)] + 0.5*(Uy[(i - 1) + nx*(j + 1)] - Uy[(i - 2) + nx*(j + 2)]);

		ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 6] = -f[i + nx*j + nx*ny * 8] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 - Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 7] = -f[i + nx*j + nx*ny * 5] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 + Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
*/

		

		//case 2
		//ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		//ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
		//ftemp[i + nx*j + nx*ny * 7] = -ftemp[i + nx*j + nx*ny * 5];

		//Ux1 = Ux[(i - 1) + nx*(j + 1)] + 0.5*(Ux[(i - 1) + nx*(j + 1)] - Ux[(i - 2) + nx*(j + 2)]);

		//ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);
		//ftemp[i + nx*j + nx*ny * 6] = -f[i + nx*j + nx*ny * 8] + (1.0 / 18.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);
		//
		//ftemp[i + nx*j + nx*ny * 0] = rho1 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
		//	+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);


		//case 3
		/*ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];*/


		//case 4
		/*ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];

		Ux1 = Ux[(i - 1) + nx*(j + 1)] + 0.5*(Ux[(i - 1) + nx*(j + 1)] - Ux[(i - 2) + nx*(j + 2)]);

		ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);
		ftemp[i + nx*j + nx*ny * 7] = -f[i + nx*j + nx*ny * 5] + (1.0 / 18.0) * rho1 * (1.0 + 3.0 * Ux1*Ux1);
*/

		//case 5
//		Ux1 = Ux[(i - 1) + nx*(j + 0)];
////		rho_extra = rho[(i - 1) + nx*(j + 1)] + 0.5*(rho[(i - 1) + nx*(j + 1)] - rho[(i - 2) + nx*(j + 2)]);
//		rho_extra = rho[(i - 1) + nx*(j + 0)];
//
//		ru = Ux1 * rho_extra;
//
//		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 12.0) * ru;
//		ftemp[i + nx*j + nx*ny * 7] = -(1.0 / 12.0) * ru;
//
//		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
//		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1] - (2.0 / 3.0) * ru;
//		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 6.0) * ru;
//		
//		ftemp[i + nx*j + nx*ny * 0] = rho1 - (ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
//			+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7] + ftemp[i + nx*j + nx*ny * 8]);
//		
		
		//Periodic + Bounce back
		/*ftemp[i + nx*(ny - 1) + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
		ftemp[(i - 1) + nx*(ny - 1) + nx*ny * 7] = f[i + nx*j + nx*ny * 7];

		Ux1 = Ux[(i - 1) + nx*(j + 1)];
		Uy1 = Uy[(i - 1) + nx*(j + 1)];

		ftemp[i + nx*j + nx*ny * 3] = -f[i + nx*j + nx*ny * 1] + (2.0 / 9.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 6] = -f[i + nx*j + nx*ny * 8] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 - Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));
		ftemp[i + nx*j + nx*ny * 7] = -f[i + nx*j + nx*ny * 5] + (1.0 / 18.0) * rho1 * (1.0 + (9.0 / 2.0)*pow(Ux1 + Uy1, 2) - (3.0 / 2.0)*(Ux1*Ux1 + Uy1*Uy1));*/


		//Periodic + Extrapolation
		/*ftemp[i + nx*(ny - 1) + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
		ftemp[(i - 1) + nx*(ny - 1) + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
*/

		//wet-node method
	////	rho0 = rho[(i - 1) + nx*(j + 0)];
	//	rho0 = rho1;
	////	rho0 = rho[(i - 0) + nx*(j + 1)] + 0.5*(rho[(i - 0) + nx*(j + 1)] - rho[(i - 0) + nx*(j + 2)]);

	//	ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
	//	ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
	//	ftemp[i + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 1];
	//	ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho0 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
	//		+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
	//	ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];



		//Equilibrium
		float c = 1;
		ftemp[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
		ftemp[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

	}
// ============================================================================ //
}
void LBM_GPU::BC_bounceback() {

	dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 dimGrid((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y, (a + BLOCK_SIZE_Z - 1) / BLOCK_SIZE_Z);
	Kernel_BC_bounceback << < dimGrid, dimBlock >> > (d_f, d_ftemp, d_rho, d_Ux, d_Uy, d_Ux0, rho1, nx, ny, a);
}

__global__
void Kernel_BC_extra(float* ftemp, float* Ux, float* rho, int nx, int ny, int a, float rho1) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.z * blockIdx.z + threadIdx.z;
	if (i >= nx || j >= ny || k >= a) return;

	float ru, Ux_extra, rho_extra;
// ============================================================================ //
//	TOP-LEFT CORNER (VELOCITY & PERIODIC)
// ============================================================================ //
	if ((i == 0) && (j == ny - 1)) {
		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i + 0) + nx*(j - 1) + nx*ny * 8];*/
		

		//Extrapolation high order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 0] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 0] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 0]);
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 1] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 1] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 2] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 2] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 2]);
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 3] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 3] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 4] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 4] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 5] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 5] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 5]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 6] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 6] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 7] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 7] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 7]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i + 1) + nx*(j - 1) + nx*ny * 8] + 0.5 * (ftemp[(i + 1) + nx*(j - 1) + nx*ny * 8] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 8]);
*/

		//Extrapolation 2nd order + moving wall
		/*ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*(j - 1) + nx*ny * 1] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 1] - ftemp[i + nx*(j - 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*(j - 1) + nx*ny * 5] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 5] - ftemp[i + nx*(j - 2) + nx*ny * 5]);
	
		rho_extra = rho[i + nx*(j - 1)] + 0.5 * (rho[i + nx*(j - 1)] - rho[i + nx*(j - 2)]);
		Ux_extra = Ux[i + nx*(j - 1)] + 0.5 * (Ux[i + nx*(j - 1)] - Ux[i + nx*(j - 2)]);
		ru = rho_extra*Ux_extra;

		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] - (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;
*/

		//Extrapolation 2nd order
		/*ftemp[i + nx*j + nx*ny * 1] = 2.0*ftemp[(i + 1) + nx*(j - 1) + nx*ny * 1] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 5] = 2.0*ftemp[(i + 1) + nx*(j - 1) + nx*ny * 5] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 4] = 2.0*ftemp[(i + 1) + nx*(j - 1) + nx*ny * 4] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 7] = 2.0*ftemp[(i + 1) + nx*(j - 1) + nx*ny * 7] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = 2.0*ftemp[(i + 1) + nx*(j - 1) + nx*ny * 8] - ftemp[(i + 2) + nx*(j - 2) + nx*ny * 8];
	*/

		//Zou - He boundary
		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho1 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];
	}
// ============================================================================ //


// ============================================================================ //
//	BOTTOM-LEFT CORNER (VELOCITY & PERIODIC)
// ============================================================================ //
	if ((i == 0) && (j == 0)) {
		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i + 0) + nx*(j + 1) + nx*ny * 8];*/

		//Extrapolation high order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 0] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 0] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 0]);
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 1] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 1] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 2] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 2] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 2]);
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 3] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 3] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 4] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 4] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 5] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 5] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 5]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 6] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 6] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 7] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 7] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 7]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i + 1) + nx*(j + 1) + nx*ny * 8] + 0.5 * (ftemp[(i + 1) + nx*(j + 1) + nx*ny * 8] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 8]);
*/

		//Extrapolation 2nd order + moving wall
		/*ftemp[i + nx*j + nx*ny * 1] = ftemp[i+ nx*(j + 1) + nx*ny * 1] + 0.5 * (ftemp[i+ nx*(j + 1) + nx*ny * 1] - ftemp[i+ nx*(j + 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i+ nx*(j + 1) + nx*ny * 8] + 0.5 * (ftemp[i+ nx*(j + 1) + nx*ny * 8] - ftemp[i+ nx*(j + 2) + nx*ny * 8]);

		rho_extra = rho[i + nx*(j + 1)] + 0.5 * (rho[i + nx*(j + 1)] - rho[i + nx*(j + 2)]);
		Ux_extra = Ux[i + nx*(j + 1)] + 0.5 * (Ux[i + nx*(j + 1)] - Ux[i + nx*(j + 2)]);
		ru = rho_extra*Ux_extra;

		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 6.0)*ru;
*/
		//Extrapolation 2nd order
		/*ftemp[i + nx*j + nx*ny * 1] = 2.0*ftemp[(i + 1) + nx*(j + 1) + nx*ny * 1] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 8] = 2.0*ftemp[(i + 1) + nx*(j + 1) + nx*ny * 8] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 8];
		ftemp[i + nx*j + nx*ny * 2] = 2.0*ftemp[(i + 1) + nx*(j + 1) + nx*ny * 2] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 5] = 2.0*ftemp[(i + 1) + nx*(j + 1) + nx*ny * 5] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = 2.0*ftemp[(i + 1) + nx*(j + 1) + nx*ny * 6] - ftemp[(i + 2) + nx*(j + 2) + nx*ny * 6];
*/

		//Zou - He boundary
		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho1 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];
	}
// ============================================================================ //


// ============================================================================ //
//	TOP-RIGHT CORNER (EXTRAPOLATION & PERIODIC)
// ============================================================================ //
	if ((i == nx - 1) && (j == ny - 1)) {

		//Extrapolation
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 0] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 0] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 0]);
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 1] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 1] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 2] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 2] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 2]);
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 3] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 3] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 4] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 4] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 5] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 5] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 5]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 6] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 6] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 7] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 7] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 7]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i - 1) + nx*(j - 1) + nx*ny * 8] + 0.5 * (ftemp[(i - 1) + nx*(j - 1) + nx*ny * 8] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 8]);*/

		//Extrapolation high order
		/*ftemp[i + nx*j + nx*ny * 3] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 3] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 3] + ftemp[(i - 3) + nx*j + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 6] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 6] + ftemp[(i - 3) + nx*j + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 7] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 7] + ftemp[(i - 3) + nx*j + nx*ny * 7]);
*/
		//Extrapolation 2nd order
		/*ftemp[i + nx*j + nx*ny * 3] = 2.0*ftemp[(i - 1) + nx*(j - 1) + nx*ny * 3] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 6] = 2.0*ftemp[(i - 1) + nx*(j - 1) + nx*ny * 6] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = 2.0*ftemp[(i - 1) + nx*(j - 1) + nx*ny * 7] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 4] = 2.0*ftemp[(i - 1) + nx*(j - 1) + nx*ny * 4] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 8] = 2.0*ftemp[(i - 1) + nx*(j - 1) + nx*ny * 8] - ftemp[(i - 2) + nx*(j - 2) + nx*ny * 8];
	*/

		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*j + nx*ny * 7];
*/

		//Extrapolation second order
		/*ftemp[i + nx*j + nx*ny * 3] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 3] - ftemp[(i - 2) + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 6] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 6] - ftemp[(i - 2) + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 7] - ftemp[(i - 2) + nx*j + nx*ny * 7];*/

		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i - 1) + nx*(j - 0) + nx*ny * 8];*/


		//Extrapolation 2nd order + moving wall
		/*ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*(j - 1) + nx*ny * 3] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 3] - ftemp[i + nx*(j - 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*(j - 1) + nx*ny * 6] + 0.5 * (ftemp[i + nx*(j - 1) + nx*ny * 6] - ftemp[i + nx*(j - 2) + nx*ny * 6]);

		rho_extra = rho[i + nx*(j - 1)] + 0.5 * (rho[i + nx*(j - 1)] - rho[i + nx*(j - 2)]);
		Ux_extra = Ux[i + nx*(j - 1)] + 0.5 * (Ux[i + nx*(j - 1)] - Ux[i + nx*(j - 2)]);
		ru = rho_extra*Ux_extra;

		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] - (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru;
*/

		//Zou - He boundary
		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho1 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];
	}
// ============================================================================ //


// ============================================================================ //
//	BOTTOM-RIGHT CORNER (EXTRAPOLATION & PERIODIC)
// ============================================================================ //
	if ((i == nx - 1) && (j == 0)) {

		//Extrapolation
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 0] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 0] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 0]);
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 1] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 1] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 1]);
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 2] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 2] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 2]);
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 3] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 3] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 4] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 4] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 4]);
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 5] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 5] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 5]);
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 6] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 6] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 7] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 7] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 7]);
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i - 1) + nx*(j + 1) + nx*ny * 8] + 0.5 * (ftemp[(i - 1) + nx*(j + 1) + nx*ny * 8] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 8]);
*/

		//Extrapolation high order
	/*	ftemp[i + nx*j + nx*ny * 3] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 3] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 3] + ftemp[(i - 3) + nx*j + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 6] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 6] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 6] + ftemp[(i - 3) + nx*j + nx*ny * 6]);
		ftemp[i + nx*j + nx*ny * 7] = (1.0 / 3.0) * (7.0*ftemp[(i - 1) + nx*j + nx*ny * 7] - 5.0*ftemp[(i - 2) + nx*j + nx*ny * 7] + ftemp[(i - 3) + nx*j + nx*ny * 7]);
*/

		//Extrapolation 2nd order
		/*ftemp[i + nx*j + nx*ny * 3] = 2.0*ftemp[(i - 1) + nx*(j + 1) + nx*ny * 3] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 6] = 2.0*ftemp[(i - 1) + nx*(j + 1) + nx*ny * 6] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = 2.0*ftemp[(i - 1) + nx*(j + 1) + nx*ny * 7] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 2] = 2.0*ftemp[(i - 1) + nx*(j + 1) + nx*ny * 2] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 5] = 2.0*ftemp[(i - 1) + nx*(j + 1) + nx*ny * 5] - ftemp[(i - 2) + nx*(j + 2) + nx*ny * 5];*/



		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*j + nx*ny * 7];*/

		//Extrapolation second order
		/*ftemp[i + nx*j + nx*ny * 3] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 3] - ftemp[(i - 2) + nx*j + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 6] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 6] - ftemp[(i - 2) + nx*j + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = 2.0*ftemp[(i - 1) + nx*j + nx*ny * 7] - ftemp[(i - 2) + nx*j + nx*ny * 7];
*/

		//Extrapolation first order
		/*ftemp[i + nx*j + nx*ny * 0] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 0];
		ftemp[i + nx*j + nx*ny * 1] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 2] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 2];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 3];
		ftemp[i + nx*j + nx*ny * 4] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 5];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 6];
		ftemp[i + nx*j + nx*ny * 7] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 7];
		ftemp[i + nx*j + nx*ny * 8] = ftemp[(i - 1) + nx*(j + 0) + nx*ny * 8];
*/

		//Extrapolation 2nd order + moving wall
		/*ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*(j + 1) + nx*ny * 3] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 3] - ftemp[i + nx*(j + 2) + nx*ny * 3]);
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*(j + 1) + nx*ny * 7] + 0.5 * (ftemp[i + nx*(j + 1) + nx*ny * 7] - ftemp[i + nx*(j + 2) + nx*ny * 7]);

		rho_extra = rho[i + nx*(j + 1)] + 0.5 * (rho[i + nx*(j + 1)] - rho[i + nx*(j + 2)]);
		Ux_extra = Ux[i + nx*(j + 1)] + 0.5 * (Ux[i + nx*(j + 1)] - Ux[i + nx*(j + 2)]);
		ru = rho_extra*Ux_extra;

		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru;
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 6.0)*ru;
*/

		//Zou - He boundary
		ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8];
		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1];
		ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho1 - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
			+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];
	}
// ============================================================================ //
}
void LBM_GPU::BC_extra() {

	dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 dimGrid((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y, (a + BLOCK_SIZE_Z - 1) / BLOCK_SIZE_Z);
	Kernel_BC_extra << < dimGrid, dimBlock >> > (d_ftemp, d_Ux, d_rho, nx, ny, a, rho1);
}

__global__ 
void Kernel_Eq(float* ftemp, float* feq, float* Ux, float* Uy, float* rho, float* ex, float* ey, int nx, int ny, int a, int* is_solid_node, float c) {

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.z * blockIdx.z + threadIdx.z;
	if (i >= nx || j >= ny || k >= a) return;


	//Calculation of Macroscopic var 
	if (!is_solid_node[i + nx*j]){
	rho[i + nx*j] = ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1]
		+ ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4]
		+ ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7]
		+ ftemp[i + nx*j + nx*ny * 8];

	Ux[i + nx*j] = ftemp[i + nx*j + nx*ny * 1] * ex[1] + ftemp[i + nx*j + nx*ny * 3] * ex[3]
		+ ftemp[i + nx*j + nx*ny * 5] * ex[5] + ftemp[i + nx*j + nx*ny * 6] * ex[6] + ftemp[i + nx*j + nx*ny * 7] * ex[7]
		+ ftemp[i + nx*j + nx*ny * 8] * ex[8];

	Uy[i + nx*j] = ftemp[i + nx*j + nx*ny * 2] * ey[2] + ftemp[i + nx*j + nx*ny * 4] * ey[4]
		+ ftemp[i + nx*j + nx*ny * 5] * ey[5] + ftemp[i + nx*j + nx*ny * 6] * ey[6] + ftemp[i + nx*j + nx*ny * 7] * ey[7]
		+ ftemp[i + nx*j + nx*ny * 8] * ey[8];

	Ux[i + nx*j] /= rho[i + nx*j];
	Uy[i + nx*j] /= rho[i + nx*j];



	feq[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4)) * pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2)) * (pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Ux[i + nx*j] + (4.5 / pow(c, 4))*pow(Ux[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - (3.0 / pow(c, 2)) * Uy[i + nx*j] + (4.5 / pow(c, 4))*pow(Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] + Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (-Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
	feq[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + (3.0 / pow(c, 2)) * (Ux[i + nx*j] - Uy[i + nx*j]) + (4.5 / pow(c, 4))*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - (1.5 / pow(c, 2))*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));

	}
}
__global__
void Kernel_Collision(float* fN, float* ftemp, float* feq, int nx, int ny, int a, float tau, int* is_solid_node) {

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.z * blockIdx.z + threadIdx.z;
	if (i >= nx || j >= ny || k >= a) return;

	if (!is_solid_node[i + nx*j]) {
		fN[i + nx*j + nx*ny*k] = ftemp[i + nx*j + nx*ny*k] - (ftemp[i + nx*j + nx*ny*k] - feq[i + nx*j + nx*ny*k]) / tau;
	}
}
void LBM_GPU::Collision() {

	dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 dimGrid((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y, (a + BLOCK_SIZE_Z - 1) / BLOCK_SIZE_Z);
	Kernel_Eq << < dimGrid, dimBlock >> > (d_ftemp, d_feq, d_Ux, d_Uy, d_rho, d_ex, d_ey, nx, ny, a, d_is_solid_node, c);
	Kernel_Collision << < dimGrid, dimBlock >> > (d_fN, d_ftemp, d_feq, nx, ny, a, tau, d_is_solid_node);
}

__global__ 
void Kernel_Error(float* f, float* Ux, float* Uy, float* U, float* rho, float* fN, float* UxN, float* UyN, float* UN, float* rhoN, float* ex, float* ey, int nx, int ny, int a, int* is_solid_node) {

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.z * blockIdx.z + threadIdx.z;
	if (i >= nx || j >= ny || k >= a) return;

	if (!is_solid_node[i + nx*j]) {
		rho[i + nx*j] = f[i + nx*j + nx*ny * 0] + f[i + nx*j + nx*ny * 1]
			+ f[i + nx*j + nx*ny * 2] + f[i + nx*j + nx*ny * 3] + f[i + nx*j + nx*ny * 4]
			+ f[i + nx*j + nx*ny * 5] + f[i + nx*j + nx*ny * 6] + f[i + nx*j + nx*ny * 7]
			+ f[i + nx*j + nx*ny * 8];

		Ux[i + nx*j] = f[i + nx*j + nx*ny * 1] * ex[1] + f[i + nx*j + nx*ny * 3] * ex[3]
			+ f[i + nx*j + nx*ny * 5] * ex[5] + f[i + nx*j + nx*ny * 6] * ex[6] + f[i + nx*j + nx*ny * 7] * ex[7]
			+ f[i + nx*j + nx*ny * 8] * ex[8];

		Uy[i + nx*j] = f[i + nx*j + nx*ny * 2] * ey[2] + f[i + nx*j + nx*ny * 4] * ey[4]
			+ f[i + nx*j + nx*ny * 5] * ey[5] + f[i + nx*j + nx*ny * 6] * ey[6] + f[i + nx*j + nx*ny * 7] * ey[7]
			+ f[i + nx*j + nx*ny * 8] * ey[8];

		Ux[i + nx*j] /= rho[i + nx*j];
		Uy[i + nx*j] /= rho[i + nx*j];
		U[i + nx*j] = sqrt(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2));




		rhoN[i + nx*j] = fN[i + nx*j + nx*ny * 0] + fN[i + nx*j + nx*ny * 1]
			+ fN[i + nx*j + nx*ny * 2] + fN[i + nx*j + nx*ny * 3] + fN[i + nx*j + nx*ny * 4]
			+ fN[i + nx*j + nx*ny * 5] + fN[i + nx*j + nx*ny * 6] + fN[i + nx*j + nx*ny * 7]
			+ fN[i + nx*j + nx*ny * 8];

		UxN[i + nx*j] = fN[i + nx*j + nx*ny * 1] * ex[1] + fN[i + nx*j + nx*ny * 3] * ex[3]
			+ fN[i + nx*j + nx*ny * 5] * ex[5] + fN[i + nx*j + nx*ny * 6] * ex[6] + fN[i + nx*j + nx*ny * 7] * ex[7]
			+ fN[i + nx*j + nx*ny * 8] * ex[8];

		UyN[i + nx*j] = fN[i + nx*j + nx*ny * 2] * ey[2] + fN[i + nx*j + nx*ny * 4] * ey[4]
			+ fN[i + nx*j + nx*ny * 5] * ey[5] + fN[i + nx*j + nx*ny * 6] * ey[6] + fN[i + nx*j + nx*ny * 7] * ey[7]
			+ fN[i + nx*j + nx*ny * 8] * ey[8];

		UxN[i + nx*j] /= rhoN[i + nx*j];
		UyN[i + nx*j] /= rhoN[i + nx*j];
		UN[i + nx*j] = sqrt(pow(UxN[i + nx*j], 2) + pow(UyN[i + nx*j], 2));
	}
}
void LBM_GPU::Error() {

	dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 dimGrid((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y, (a + BLOCK_SIZE_Z - 1) / BLOCK_SIZE_Z);
	Kernel_Error << < dimGrid, dimBlock >> > (d_f, d_Ux, d_Uy, d_U, d_rho, d_fN, d_UxN, d_UyN, d_UN, d_rhoN, d_ex, d_ey, nx, ny, a, d_is_solid_node);

	cudaMemcpy(U, d_U, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(UN, d_UN, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);

	sum = 0.0;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {

			if (!is_solid_node[i + nx*j]) {
				sum = sum + pow(abs(UN[i + nx*j] - U[i + nx*j]), 2);
			}
		}
	}
	error = sqrt(sum / (nx*ny - sn));

}

__global__ 
void Kernel_Update(float* fN, float* f, float* Ux, float* Uy, float* U, float* rho, float* ex, float* ey, int nx, int ny, int a, int* is_solid_node) {

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.z * blockIdx.z + threadIdx.z;
	if (i >= nx || j >= ny || k >= a) return;

	if(!is_solid_node[i + nx*j]) f[i + nx*j + nx*ny*k] = fN[i + nx*j + nx*ny*k];

	rho[i + nx*j] = f[i + nx*j + nx*ny * 0] + f[i + nx*j + nx*ny * 1]
		+ f[i + nx*j + nx*ny * 2] + f[i + nx*j + nx*ny * 3] + f[i + nx*j + nx*ny * 4]
		+ f[i + nx*j + nx*ny * 5] + f[i + nx*j + nx*ny * 6] + f[i + nx*j + nx*ny * 7]
		+ f[i + nx*j + nx*ny * 8];

	Ux[i + nx*j] = f[i + nx*j + nx*ny * 1] * ex[1] + f[i + nx*j + nx*ny * 3] * ex[3]
		+ f[i + nx*j + nx*ny * 5] * ex[5] + f[i + nx*j + nx*ny * 6] * ex[6] + f[i + nx*j + nx*ny * 7] * ex[7]
		+ f[i + nx*j + nx*ny * 8] * ex[8];

	Uy[i + nx*j] = f[i + nx*j + nx*ny * 2] * ey[2] + f[i + nx*j + nx*ny * 4] * ey[4]
		+ f[i + nx*j + nx*ny * 5] * ey[5] + f[i + nx*j + nx*ny * 6] * ey[6] + f[i + nx*j + nx*ny * 7] * ey[7]
		+ f[i + nx*j + nx*ny * 8] * ey[8];

	Ux[i + nx*j] /= rho[i + nx*j];
	Uy[i + nx*j] /= rho[i + nx*j];
	U[i + nx*j] = sqrt(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2));
}
void LBM_GPU::Update() {


	dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 dimGrid((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y, (a + BLOCK_SIZE_Z - 1) / BLOCK_SIZE_Z);
	Kernel_Update << < dimGrid, dimBlock >> > (d_fN, d_f, d_Ux, d_Uy, d_U, d_rho, d_ex, d_ey, nx, ny, a, d_is_solid_node);

}

void LBM_GPU::Momentum() {
	cudaMemcpy(f, d_f, nx*ny*a * sizeof(float), cudaMemcpyDeviceToHost);

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

			dist = sqrt(pow((float)i - ic, 2) + pow((float)j - jc, 2));
			q = dist - r;

			if (!is_boundary_node[i + nx*j]) {
				if (!is_solid_node[i + nx*j]) {
					if (is_solid_near_node[i + nx*j]) {


						if (is_solid_node[ip + nx*j]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 3] = 2.0 * q * f[i + nx*j + nx*ny * 1] + (1.0 - 2.0*q)*f[(i - 1) + nx*j + nx*ny * 1];
							else ftemp[i + nx*j + nx*ny * 3] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 1] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 3];

							sum_Fx1 = ex[1] * (ftemp[i + nx*j + nx*ny * 3] + f[i + nx*j + nx*ny * 1]);
						}

						if (is_solid_node[i + nx*jp]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 4] = 2.0 * q * f[i + nx*j + nx*ny * 2] + (1.0 - 2.0*q)*f[i + nx*(j - 1) + nx*ny * 2];
							else ftemp[i + nx*j + nx*ny * 4] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 2] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 4];

							sum_Fy2 = ey[2] * (ftemp[i + nx*j + nx*ny * 4] + f[i + nx*j + nx*ny * 2]);
						}

						if (is_solid_node[in + nx*j]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 1] = 2.0 * q * f[i + nx*j + nx*ny * 3] + (1.0 - 2.0*q)*f[(i + 1) + nx*j + nx*ny * 3];
							else ftemp[i + nx*j + nx*ny * 1] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 3] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 1];

							sum_Fx3 = ex[3] * (ftemp[i + nx*j + nx*ny * 1] + f[i + nx*j + nx*ny * 3]);
						}

						if (is_solid_node[i + nx*jn]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 2] = 2.0 * q * f[i + nx*j + nx*ny * 4] + (1.0 - 2.0*q)*f[i + nx*(j + 1) + nx*ny * 4];
							else ftemp[i + nx*j + nx*ny * 2] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 4] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 2];

							sum_Fy4 = ey[4] * (ftemp[i + nx*j + nx*ny * 2] + f[i + nx*j + nx*ny * 4]);
						}

						if (is_solid_node[ip + nx*jp]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 7] = 2.0 * q * f[i + nx*j + nx*ny * 5] + (1.0 - 2.0*q)*f[(i - 1) + nx*(j - 1) + nx*ny * 5];
							else ftemp[i + nx*j + nx*ny * 7] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 5] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 7];

							sum_Fx5 = ex[5] * (ftemp[i + nx*j + nx*ny * 7] + f[i + nx*j + nx*ny * 5]);
							sum_Fy5 = ey[5] * (ftemp[i + nx*j + nx*ny * 7] + f[i + nx*j + nx*ny * 5]);
						}

						if (is_solid_node[in + nx*jp]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 8] = 2.0 * q * f[i + nx*j + nx*ny * 6] + (1.0 - 2.0*q)*f[(i + 1) + nx*(j - 1) + nx*ny * 6];
							else ftemp[i + nx*j + nx*ny * 8] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 6] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 8];

							sum_Fx6 = ex[6] * (ftemp[i + nx*j + nx*ny * 8] + f[i + nx*j + nx*ny * 6]);
							sum_Fy6 = ey[6] * (ftemp[i + nx*j + nx*ny * 8] + f[i + nx*j + nx*ny * 6]);
						}

						if (is_solid_node[in + nx*jn]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 5] = 2.0 * q * f[i + nx*j + nx*ny * 7] + (1.0 - 2.0*q)*f[(i + 1) + nx*(j + 1) + nx*ny * 7];
							else ftemp[i + nx*j + nx*ny * 5] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 7] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 5];

							sum_Fx7 = ex[7] * (ftemp[i + nx*j + nx*ny * 5] + f[i + nx*j + nx*ny * 7]);
							sum_Fy7 = ey[7] * (ftemp[i + nx*j + nx*ny * 5] + f[i + nx*j + nx*ny * 7]);
						}

						if (is_solid_node[ip + nx*jn]) {

							if (q < 0.5) ftemp[i + nx*j + nx*ny * 6] = 2.0 * q * f[i + nx*j + nx*ny * 8] + (1.0 - 2.0*q)*f[(i - 1) + nx*(j + 1) + nx*ny * 8];
							else ftemp[i + nx*j + nx*ny * 6] = (1.0 / (2.0*q))*f[i + nx*j + nx*ny * 8] + (2.0*q - 1.0) / (2.0*q)*f[i + nx*j + nx*ny * 6];

							sum_Fx8 = ex[8] * (ftemp[i + nx*j + nx*ny * 6] + f[i + nx*j + nx*ny * 8]);
							sum_Fy8 = ey[8] * (ftemp[i + nx*j + nx*ny * 6] + f[i + nx*j + nx*ny * 8]);
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

	fout_GPU_Cd << Cd << "\t" << Cl << endl;
}

void LBM_GPU::Print() {

	cudaMemcpy(Ux, d_Ux, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(Uy, d_Uy, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(U, d_U, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(rho, d_rho, nx*ny * sizeof(float), cudaMemcpyDeviceToHost);


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



	fout_GPU << endl;
	fout_GPU << "variables = X Y Ux Uy U rho P" << endl;
	fout_GPU << "zone i=" << nx << " j=" << ny << endl;
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			fout_GPU << i << "\t" << j << "\t" << Ux[i + nx*j] << "\t" << Uy[i + nx*j] << "\t"
				<< U[i + nx*j] << "\t" << rho[i + nx*j] << "\t" << P[i + nx*j] << endl;
		}
	}
	fout_GPU << endl;


	i = 0;
	fout_GPU_Ux << "variables = X Y Ux" << endl;
	fout_GPU_Ux << "zone i=" << nx << " j=" << ny << endl;
	for (j = 0; j < ny; j++) {

		fout_GPU_Ux << i << "\t" << j << "\t" << Ux[i + nx*j] << endl;

		
	}
	fout_GPU_Ux << endl;
}

LBM_GPU::~LBM_GPU()
{
	cudaFree(d_Ux0);
	cudaFree(d_is_boundary_node);
	cudaFree(d_is_solid_node);
	cudaFree(d_f);
	cudaFree(d_fN);
	cudaFree(d_ftemp);
	cudaFree(d_feq);
	cudaFree(d_Ux);
	cudaFree(d_Uy);
	cudaFree(d_rho);
	cudaFree(d_ex);
	cudaFree(d_ey);
	cudaFree(d_U);
	cudaFree(d_UN);
	cudaFree(d_UxN);
	cudaFree(d_UyN);
	cudaFree(rhoN);

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
