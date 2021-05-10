#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class LBM_GPU
{
public:
// ============================================================================ //
//  VARIABLES
// ============================================================================ //
	int nx;								//The number of node in X direction
	int ny;								//The number of node in Y direction
	int a;								//The number of direction in D2Q9 model

	float Lx;							//The length in X direction (Physical property)
	float Ly;							//The length in Y direction (Physical property)
	float del_x_p;						//The grid step in X direction (Physical property)
	float del_y_p;						//The grid step in Y direction (Physical property)
	float del_t_p;						//The time step (Physical property)

	float del_x;						//The grid step in X direction (LBM property)
	float del_y;						//The grid step in Y direction (LBM property)
	float del_t;						//The time step (LBM property)


	float* Ux;							//The macroscopic velocity in X direction (LBM property)
	float* Uy;							//The macroscopic velocity in Y direction (LBM property)
	float* U;							//The macroscopic velocity (LBM property)
	float* rho;							//The macroscopic density (LBM property)

	float* Ux_p;						//The macroscopic velocity in X direction (Physical property)
	float* Uy_p;						//The macroscopic velocity in Y direction (Physical property)
	float* U_p;							//The macroscopic velocity (Physical property)
	float* P;							//The macroscopic pressure (Physical property)
	float Um_p;							//The maximum velocity (Physical property) [m/s]
	

	float* UxN;							//The macroscopic velocity in X direction at n + 1 step (LBM property)
	float* UyN;							//The macroscopic velocity in Y direction at n + 1 step (LBM property)
	float* UN;							//The macroscopic velocity at n + 1 step (LBM property)
	float* rhoN;						//The macroscopic density (LBM property)
	float Um;							//The maximum velocity (LBM property)

	float* f;							//The distribution function at n step
	float* ftemp;						//The distribution function at temp step
	float* feq;							//The equilibrium distribution funtion
	float* fN;							//The distribution function at n + 1 step
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
	

	int BLOCK_SIZE_X;					//The Size of block in X direction in GPU
	int BLOCK_SIZE_Y;					//The Size of block in Y direction in GPU
	int BLOCK_SIZE_Z;					//The Size of block in Z direction in GPU
	int* d_is_boundary_node;
	int* d_is_solid_node;

	float* d_f;							//The distribution function at n step used in GPU
	float* d_ftemp;						//The distribution function at temp step used in GPU
	float* d_feq;
	float* d_fN;
	float* d_Ux;
	float* d_Uy;
	float* d_U;		
	float* d_UN;	
	float* d_rho;
	float* d_ex;
	float* d_ey;
	float* d_UxN;
	float* d_UyN;
	float* d_rhoN;
	float* d_Ux0;

	char comment[60];

	float c;							//The lattice speed (LBM property)
	float c_s;							//The speed of sound (LBM property)
	float r;							//The radius of cylinder
	float q;
// ============================================================================ //




// ============================================================================ //
//  FUNCTIONS
// ============================================================================ //
	LBM_GPU();
	~LBM_GPU();

	void Streaming();
	void BC_bounceback();
	void BC_extra();
	void BC_vel();
	void Collision();
	void Error();
	void Update();
	void Print();
	void Momentum();
// =========================================================================== //
};