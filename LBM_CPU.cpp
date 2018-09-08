//#include "LBM_CPU.h"
//ofstream fout_CPU("out_CPU.dat");
//ofstream fout_CPU_Ux("out_CPU_Ux.dat");
//ifstream fin_CPU("in_CPU.txt");
//ifstream fin_grid_CPU("grid.dat");
//
//LBM_CPU::LBM_CPU()
//{
//// ============================================================================ //
////  LOAD THE PARAMETERS
//// ============================================================================ //
//	fin_CPU >> nx;		fin_CPU >> comment;
//	fin_CPU >> ny;		fin_CPU >> comment;
//	fin_CPU >> Lx;		fin_CPU >> comment;
//	fin_CPU >> Ly;		fin_CPU >> comment;
//	fin_CPU >> a;		fin_CPU >> comment;
//	fin_CPU >> Re;		fin_CPU >> comment;
//	fin_CPU >> Ux0;		fin_CPU >> comment;
//// ============================================================================ //
//	
//
//// ============================================================================ //
////  NEW & CUDAMALLOC
//// ============================================================================ //
//	is_boundary_node = new int[nx*ny];
//	is_solid_node = new int[nx*ny];
//	U = new float[nx*ny];
//	Ux = new float[nx*ny];
//	Uy = new float[nx*ny];
//	rho = new float[nx*ny];
//	UN = new float[nx*ny];
//	UxN = new float[nx*ny];
//	UyN = new float[nx*ny];
//	rhoN = new float[nx*ny];
//	f = new float[nx*ny*a];
//	ftemp = new float[nx*ny*a];
//	fN = new float[nx*ny*a];
//	feq = new float[nx*ny*a];
//	ex = new float[a];
//	ey = new float[a];
//	U_p = new float[nx*ny];
//	Ux_p = new float[nx*ny];
//	Uy_p = new float[nx*ny];
//	P = new float[nx*ny];
//// ============================================================================ //
//
//
//// ============================================================================ //
////  Microscopic velocity
//// ============================================================================ //
//	ex[0] = 0.0,	ey[0] = 0.0;
//	ex[1] = 1.0,	ey[1] = 0.0;
//	ex[2] = 0.0,	ey[2] = 1.0;
//	ex[3] = -1.0,	ey[3] = 0.0;
//	ex[4] = 0.0,	ey[4] = -1.0;
//	ex[5] = 1.0,	ey[5] = 1.0;
//	ex[6] = -1.0,	ey[6] = 1.0;
//	ex[7] = -1.0,	ey[7] = -1.0;
//	ex[8] = 1.0,	ey[8] = -1.0;
//// ============================================================================ //
//
//
//
//// ============================================================================ //
////  SET BOUNDARY NODE & SOLID NODE
//// ============================================================================ //
//	sIm = nx / Lx * 2;
//	sIM = nx / Lx * 3;
//	sJm = ny / Ly * 3;
//	sJM = ny / Ly * 4;
//	
//	snx = (sIM - sIm) + 1;
//	sny = (sJM - sJm) + 1;
//
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//			if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) is_boundary_node[i + nx*j] = 1;
//			else is_boundary_node[i + nx*j] = 0;
//
//			if ((i >= sIm && i <= sIM) && (j >= sJm && j <= sJM)) fin_grid_CPU >> is_solid_node[i + nx*j];
//			else is_solid_node[i + nx*j] = 0;
//		}
//	}
//// ============================================================================ //
//
//
//
//
//// ============================================================================ //
////  INITIAL CONDITION
//// ============================================================================ //
//	del_x = Lx / (float)nx;
//	del_y = Ly / (float)ny;
//	del_t = pow(del_y, 2);
//
//	Ux0_p = Ux0 * (del_y / del_t);
//	tau = 3.0*(del_t / pow(del_y, 2))*(Ux0_p * Ly / Re) + 0.5;
//
//	nu = (1.0 / 3.0)*(tau - 0.5);
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//
//			Ux[i + nx*j] = 0.0;
//			Uy[i + nx*j] = 0.0;
//			U[i + nx*j] = 0.0;
//			UxN[i + nx*j] = 0.0;
//			UyN[i + nx*j] = 0.0;
//			UN[i + nx*j] = 0.0;
//			P[i + nx*j] = 0.0;
//
//			for (k = 0; k < a; k++) {
//				ftemp[i + nx*j + nx*ny*k] = 0.0;
//				feq[i + nx*j + nx*ny*k] = 0.0;
//				fN[i + nx*j + nx*ny*k] = 0.0;
//			}
//
//
//			if (!is_solid_node[i + nx*j]) {
//				rho[i + nx*j] = 1.0;
//			}
//			else rho[i + nx*j] = 10.0;
//
//			f[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j];
//			f[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j];
//
//		}
//	}
//// ============================================================================ //
//}
//
//void LBM_CPU::Streaming() {
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//
//			in = i - 1;
//			ip = i + 1;
//			jn = j - 1;
//			jp = j + 1;
//
//			if (!is_boundary_node[i + nx*j]) {
//				if (!is_solid_node[i + nx*j]) {
//
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//
//					if (!is_solid_node[ip + nx*j]) ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
//					else ftemp[i + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 1];
//
//					if (!is_solid_node[i + nx*jp]) ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
//					else ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
//
//					if (!is_solid_node[in + nx*j]) ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
//					else ftemp[i + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 3];
//
//					if (!is_solid_node[i + nx*jn]) ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
//					else ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
//
//					if (!is_solid_node[ip + nx*jp]) ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
//					else ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
//
//					if (!is_solid_node[in + nx*jp]) ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
//					else ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
//
//					if (!is_solid_node[in + nx*jn]) ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
//					else ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
//
//					if (!is_solid_node[ip + nx*jn]) ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
//					else ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
//				}
//			}
//			else {
//				if ((i == 0) && (j > 0 && j < ny - 1)) {				//INLET
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
//					ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
//					ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
//					ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
//					ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
//				}
//				else if ((i > 0 && i < nx - 1) && (j == ny - 1)) {			//TOP
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
//					ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
//					ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
//					ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
//					ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
//				}
//				else if ((i > 0 && i < nx - 1) && (j == 0)) {				//BOTTOM
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
//					ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
//					ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
//					ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
//					ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
//				}
//				else if ((i == nx - 1) && (j > 0 && j < ny - 1)) {			//OUTLET
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
//					ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
//					ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
//					ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
//					ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
//				}
//				else if ((i == 0) && (j == 0)) {							//BOTTOM-LEFT
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
//					ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
//					ftemp[ip + nx*jp + nx*ny * 5] = f[i + nx*j + nx*ny * 5];
//				}
//				else if ((i == 0) && (j == ny - 1)) {						//TOP-LEFT
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[ip + nx*j + nx*ny * 1] = f[i + nx*j + nx*ny * 1];
//					ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
//					ftemp[ip + nx*jn + nx*ny * 8] = f[i + nx*j + nx*ny * 8];
//				}
//				else if ((i == nx - 1) && (j == ny - 1)) {					//TOP-RIGHT
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
//					ftemp[i + nx*jn + nx*ny * 4] = f[i + nx*j + nx*ny * 4];
//					ftemp[in + nx*jn + nx*ny * 7] = f[i + nx*j + nx*ny * 7];
//				}
//				else if ((i == nx - 1) && (j == 0)) {						//BOTTOM-RIGHT
//					ftemp[i + nx*j + nx*ny * 0] = f[i + nx*j + nx*ny * 0];
//					ftemp[i + nx*jp + nx*ny * 2] = f[i + nx*j + nx*ny * 2];
//					ftemp[in + nx*j + nx*ny * 3] = f[i + nx*j + nx*ny * 3];
//					ftemp[in + nx*jp + nx*ny * 6] = f[i + nx*j + nx*ny * 6];
//				}
//			}
//			
//
//		}
//	}
//
//}
//
//void LBM_CPU::BC_bounceback() {
//// ============================================================================ //
////  TOP BOUNDARY (HALF BOUNCEBACK)
//// ============================================================================ //
//	j = ny - 1;
//	for (i = 1; i < nx - 1; i++) {
//		ftemp[i + nx*j + nx*ny * 4] = f[i + nx*j + nx*ny * 2];
//		ftemp[i + nx*j + nx*ny * 7] = f[i + nx*j + nx*ny * 5];
//		ftemp[i + nx*j + nx*ny * 8] = f[i + nx*j + nx*ny * 6];
//	}
//// ============================================================================ //
//
//// ============================================================================ //
////	BOTTOM BOUNDARY (HALF BOUNCEBACK)
//// ============================================================================ //
//	j = 0;
//	for (i = 1; i < nx - 1; i++) {
//		ftemp[i + nx*j + nx*ny * 2] = f[i + nx*j + nx*ny * 4];
//		ftemp[i + nx*j + nx*ny * 5] = f[i + nx*j + nx*ny * 7];
//		ftemp[i + nx*j + nx*ny * 6] = f[i + nx*j + nx*ny * 8];
//	}
//// ============================================================================ //
//}
//
//void LBM_CPU::BC_vel() {
//// ============================================================================ //
////	LEFT BOUNDARY (VELOCITY)
//// ============================================================================ //
//	i = 0;
//	for (j = 1; j < ny - 1; j++) {
//
//		rho0 = (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
//			+ 2.0*(ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 7])) / (1.0 - Ux0);
//		ru = rho0 * Ux0;
//
//		ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3] + (2.0 / 3.0)*ru;
//		ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7] + (1.0 / 6.0)*ru - (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
//		ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6] + (1.0 / 6.0)*ru + (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
//	}
//// ============================================================================ //
//
//
//// ============================================================================ //
////  RIGHT BOUNDARY (VELOCITY)
//// ============================================================================ //
//	i = nx - 1;
//	for (j = 1; j < ny - 1; j++) {
//		rho1 = (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 2] + ftemp[i + nx*j + nx*ny * 4]
//			+ 2.0*(ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 8])) / (1.0 + Ux0);
//		ru = rho1 * Ux0;
//
//		ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1] - (2.0 / 3.0)*ru;
//		ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8] - (1.0 / 6.0)*ru - (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
//		ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5] - (1.0 / 6.0)*ru + (1.0 / 2.0)*(ftemp[i + nx*j + nx*ny * 2] - ftemp[i + nx*j + nx*ny * 4]);
//	}
//// ============================================================================ //
//
//
//// ============================================================================ //
////	TOP-LEFT CORNER (VELOCITY)
//// ============================================================================ //
//	i = 0;
//	j = ny - 1;
//
//	ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3];
//	ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
//	ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];
//	ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho[i + nx*(j - 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
//		+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
//	ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];
//// ============================================================================ //
//
//// ============================================================================ //
////	BOTTOM-LEFT CORNER (VELOCITY)
//// ============================================================================ //
//	i = 0;
//	j = 0;
//
//	ftemp[i + nx*j + nx*ny * 1] = ftemp[i + nx*j + nx*ny * 3];
//	ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
//	ftemp[i + nx*j + nx*ny * 5] = ftemp[i + nx*j + nx*ny * 7];
//	ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho[i + nx*(j + 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
//		+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
//	ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];
//// ============================================================================ //
//
//
//// ============================================================================ //
////	TOP-RIGHT CORNER (VELOCITY)
//// ============================================================================ //
//	i = nx - 1;
//	j = ny - 1;
//
//	ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1];
//	ftemp[i + nx*j + nx*ny * 4] = ftemp[i + nx*j + nx*ny * 2];
//	ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];
//	ftemp[i + nx*j + nx*ny * 6] = 0.5 * (rho[i + nx*(j - 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
//		+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 5] + ftemp[i + nx*j + nx*ny * 7]));
//	ftemp[i + nx*j + nx*ny * 8] = ftemp[i + nx*j + nx*ny * 6];
//// ============================================================================ //
//
//// ============================================================================ //
////	BOTTOM-RIGHT CORNER (VELOCITY)
//// ============================================================================ //
//	i = nx - 1;
//	j = 0;
//
//	ftemp[i + nx*j + nx*ny * 3] = ftemp[i + nx*j + nx*ny * 1];
//	ftemp[i + nx*j + nx*ny * 2] = ftemp[i + nx*j + nx*ny * 4];
//	ftemp[i + nx*j + nx*ny * 6] = ftemp[i + nx*j + nx*ny * 8];
//	ftemp[i + nx*j + nx*ny * 5] = 0.5 * (rho[i + nx*(j + 1)] - (ftemp[i + nx*j + nx*ny * 0] + ftemp[i + nx*j + nx*ny * 1] + ftemp[i + nx*j + nx*ny * 2]
//		+ ftemp[i + nx*j + nx*ny * 3] + ftemp[i + nx*j + nx*ny * 4] + ftemp[i + nx*j + nx*ny * 6] + ftemp[i + nx*j + nx*ny * 8]));
//	ftemp[i + nx*j + nx*ny * 7] = ftemp[i + nx*j + nx*ny * 5];
//// ============================================================================ //
//}
//
//void LBM_CPU::Collision() {
//
//	//Calculation of Macroscopic var 
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//			Ux[i + nx*j] = 0.0;
//			Uy[i + nx*j] = 0.0;
//			rho[i + nx*j] = 0.0;
//
//			if (!is_solid_node[i + nx*j]) {
//				for (k = 0; k < a; k++) {
//					rho[i + nx*j] += ftemp[i + nx*j + nx*ny*k];
//					Ux[i + nx*j] += ftemp[i + nx*j + nx*ny*k] * ex[k];
//					Uy[i + nx*j] += ftemp[i + nx*j + nx*ny*k] * ey[k];
//				}
//				Ux[i + nx*j] /= rho[i + nx*j];
//				Uy[i + nx*j] /= rho[i + nx*j];
//			}
//		}
//	}
//
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//			if (!is_solid_node[i + nx*j]) {
//				feq[i + nx*j + nx*ny * 0] = (4.0 / 9.0) * rho[i + nx*j] * (1.0 - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 1] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + 3.0 * Ux[i + nx*j] + 4.5*pow(Ux[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 2] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 + 3.0 * Uy[i + nx*j] + 4.5*pow(Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 3] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - 3.0 * Ux[i + nx*j] + 4.5*pow(Ux[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 4] = (1.0 / 9.0) * rho[i + nx*j] * (1.0 - 3.0 * Uy[i + nx*j] + 4.5*pow(Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 5] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (Ux[i + nx*j] + Uy[i + nx*j]) + 4.5*pow(Ux[i + nx*j] + Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 6] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (-Ux[i + nx*j] + Uy[i + nx*j]) + 4.5*pow(-Ux[i + nx*j] + Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 7] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (-Ux[i + nx*j] - Uy[i + nx*j]) + 4.5*pow(-Ux[i + nx*j] - Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//				feq[i + nx*j + nx*ny * 8] = (1.0 / 36.0) * rho[i + nx*j] * (1.0 + 3.0 * (Ux[i + nx*j] - Uy[i + nx*j]) + 4.5*pow(Ux[i + nx*j] - Uy[i + nx*j], 2) - 1.5*(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2)));
//			}
//		}
//	}
//}
//
//void LBM_CPU::Error() {
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//			Ux[i + nx*j] = 0.0;
//			Uy[i + nx*j] = 0.0;
//			rho[i + nx*j] = 0.0;
//
//			if (!is_solid_node[i + nx*j]) {
//				for (k = 0; k < a; k++) {
//
//					rho[i + nx*j] += f[i + nx*j + nx*ny*k];
//					Ux[i + nx*j] += f[i + nx*j + nx*ny*k] * ex[k];
//					Uy[i + nx*j] += f[i + nx*j + nx*ny*k] * ey[k];
//				}
//				Ux[i + nx*j] /= rho[i + nx*j];
//				Uy[i + nx*j] /= rho[i + nx*j];
//				U[i + nx*j] = sqrt(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2));
//			}
//		}
//	}
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//			UxN[i + nx*j] = 0.0;
//			UyN[i + nx*j] = 0.0;
//			rhoN[i + nx*j] = 0.0;
//
//			if (!is_solid_node[i + nx*j]) {
//				for (k = 0; k < a; k++) {
//					fN[i + nx*j + nx*ny*k] = ftemp[i + nx*j + nx*ny*k] - (ftemp[i + nx*j + nx*ny*k] - feq[i + nx*j + nx*ny*k]) / tau;
//
//					rhoN[i + nx*j] += fN[i + nx*j + nx*ny*k];
//					UxN[i + nx*j] += fN[i + nx*j + nx*ny*k] * ex[k];
//					UyN[i + nx*j] += fN[i + nx*j + nx*ny*k] * ey[k];
//				}
//				UxN[i + nx*j] /= rhoN[i + nx*j];
//				UyN[i + nx*j] /= rhoN[i + nx*j];
//				UN[i + nx*j] = sqrt(pow(UxN[i + nx*j], 2) + pow(UyN[i + nx*j], 2));
//			}
//		}
//	}
//
//
//	sum = 0.0;
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//
//			if (!is_solid_node[i + nx*j]) {
//				sum = sum + pow(abs(UN[i + nx*j] - U[i + nx*j]), 2);
//			}
//		}
//	}
//
//	error = sqrt(sum / (nx*ny - snx*sny));
//
//
//}
//
//void LBM_CPU::Update() {
//
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//			if (!is_solid_node[i + nx*j]) {
//				for (k = 0; k < a; k++) {
//					f[i + nx*j + nx*ny*k] = fN[i + nx*j + nx*ny*k];
//				}
//			}
//		}
//	}
//
//
//
//	for (i = 0; i < nx; i++) {
//		for (j = 0; j < ny; j++) {
//			Ux[i + nx*j] = 0.0;
//			Uy[i + nx*j] = 0.0;
//			rho[i + nx*j] = 0.0;
//
//			for (k = 0; k < a; k++) {
//
//				rho[i + nx*j] += f[i + nx*j + nx*ny*k];
//				Ux[i + nx*j] += f[i + nx*j + nx*ny*k] * ex[k];
//				Uy[i + nx*j] += f[i + nx*j + nx*ny*k] * ey[k];
//			}
//			Ux[i + nx*j] /= rho[i + nx*j];
//			Uy[i + nx*j] /= rho[i + nx*j];
//			U[i + nx*j] = sqrt(pow(Ux[i + nx*j], 2) + pow(Uy[i + nx*j], 2));
//		}
//	}
//}
//
//void LBM_CPU::Print() {
//
//// ============================================================================ //
////  CHANGE LBM -> PHYSICAL
//// ============================================================================ //
//	for (j = 0; j < ny; j++) {
//		for (i = 0; i < nx; i++) {
//			Ux_p[i + nx*j] = Ux[i + nx*j];
//			Uy_p[i + nx*j] = Uy[i + nx*j];
//			U_p[i + nx*j] = U[i + nx*j];
//			P[i + nx*j] = rho[i + nx*j] / (3.0);
//		}
//	}
//// ============================================================================ //
//
//
//
//	fout_CPU << endl;
//	fout_CPU << "variables = X Y Ux Uy U rho P" << endl;
//	fout_CPU << "zone i=" << nx << " j=" << ny << endl;
//	for (j = 0; j < ny; j++) {
//		for (i = 0; i < nx; i++) {
//			fout_CPU << i << "\t" << j << "\t" << Ux_p[i + nx*j] << "\t" << Uy_p[i + nx*j] << "\t"
//				<< U_p[i + nx*j] << "\t" << rho[i + nx*j] << "\t" << P[i + nx*j] << endl;
//		}
//	}
//
//	fout_CPU_Ux << "variables = X Y Ux " << endl;
//	i = nx / 2;
//	for (j = 0; j < ny; j++) {
//		fout_CPU_Ux << i << "\t" << j << "\t" << Ux_p[i + nx*j] << endl;
//	}
//}
//
//LBM_CPU::~LBM_CPU()
//{
//	delete[] Uy_p;
//	delete[] Ux_p;
//	delete[] U_p;
//	delete[] ey;
//	delete[] ex;
//	delete[] fN;
//	delete[] feq;
//	delete[] ftemp;
//	delete[] f;
//	delete[] P;
//	delete[] rhoN;
//	delete[] UyN;
//	delete[] UxN;
//	delete[] UN;
//	delete[] rho;
//	delete[] Uy;
//	delete[] Ux;
//	delete[] U;
//	delete[] is_boundary_node;
//	delete[] is_solid_node;
//	cout << endl << "Done!" << endl;
//}
//
