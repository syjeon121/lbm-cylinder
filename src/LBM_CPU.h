//#pragma once
//#include <iostream>
//#include <fstream>
//#include <cmath>
//using namespace std;
//
//class LBM_CPU
//{
//public:
//// ============================================================================ //
////  VARIABLES
//// ============================================================================ //
//	int nx;								//The number of node in X direction
//	int ny;								//The number of node in Y direction
//	int a;								//The number of direction in D2Q9 model
//
//	float Lx;							//The length in X direction
//	float Ly;							//The length in Y direction
//	float del_x;						//The grid step in X direction
//	float del_y;						//The grid step in Y direction
//	float del_t;						//The time step
//
//	float* Ux;							//The macroscopic velocity in X direction
//	float* Uy;							//The macroscopic velocity in Y direction
//	float* U;							//The macroscopic velocity
//	float* rho;							//The macroscopic density
//	float* Ux_p;						//The macroscopic velocity in X direction (Physical property)
//	float* Uy_p;						//The macroscopic velocity in Y direction (Physical property)
//	float* U_p;							//The macroscopic velocity (Physical property)
//	float* P;							//The macroscopic pressure
//
//	float* UxN;							//The macroscopic velocity in X direction at n + 1 step
//	float* UyN;							//The macroscopic velocity in Y direction at n + 1 step
//	float* UN;							//The macroscopic velocity at n + 1 step
//	float* rhoN;						//The macroscopic density
//
//	float* f;							//The distribution function at n step
//	float* ftemp;						//The distribution function at temp step
//	float* fN;							//The distribution function at n + 1 step
//	float* feq;							//The equilibrium distribution funtion
//	float* ex;							//The microscopic velocity in X direction
//	float* ey;							//The microscopic velocity in Y direction
//
//	int* is_boundary_node;				//The boundary node
//	int* is_solid_node;					//The solid node
//
//	float tau;							//The relaxation time
//	float nu;							//The kinematic viscosity
//	float Re;							//The Reynolds number
//
//
//	float Ux0_p;						//The velocity at inlet boundary (Physical property)
//	float Ux0;							//The velocity at inlet boundary (LBM property)
//	float rho0;							//The density at inlet boundary
//	float rho1;							//The density at outlet boundary
//	float ru;							//rho*vel
//
//	int sIm, sIM;						//The Solid node at Minimum node and at Maximum node in X direction
//	int sJm, sJM;						//The Solid node at Minimum node and at Maximum node in Y direction
//	int snx, sny;						//The number of node of Solid
//				
//
//	float error;
//	float sum;
//
//	int i, j, k, in, ip, jn, jp;
//	char comment[60];
//// ============================================================================ //
//
//
//
//
//// ============================================================================ //
////  FUNCTIONS
//// ============================================================================ //
//	LBM_CPU();
//	~LBM_CPU();
//
//	void Streaming();
//	void BC_bounceback();
//	void BC_vel();
//	void Collision();
//	void Error();
//	void Update();
//	void Print();
//// ============================================================================ //
//};
//
