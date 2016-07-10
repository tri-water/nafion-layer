#include "utilities.h"
#include <iostream>

void initialize(mesh& membrane, mesh& solution, nernst_equation& therm, const double& Eq, const double& Cm_bulk)
{
#pragma omp parallel for
	for (int m = 0; m < membrane.m; ++m){
		for (int n = 0; n < membrane.n; ++n){
			double K = therm.ratio_ox2red(Eq);
			membrane.Coxn[m][n] = Cm_bulk*K / (K + 1);
			membrane.Coxo[m][n] = membrane.Coxn[m][n];
			membrane.Crdn[m][n] = Cm_bulk / (K + 1);
			membrane.Crdo[m][n] = membrane.Crdn[m][n];
		}
	}
#pragma omp parallel for
	for (int m = 0; m < solution.m; ++m){
		for (int n = 0; n < solution.n; ++n){
			solution.Coxn[m][n] = 0;
			solution.Coxo[m][n] = 0;
			solution.Crdn[m][n] = 0;
			solution.Crdo[m][n] = 0;
		}
	}
}

void membrane_matrixA(mesh& membrane, const potential_signal& Esignal, const double D_Fcm, const double D_Fcpm, const double K_dis)
{
	//diffusion equation: A0*Cij + A1*Ci-1j + A2*Ci+1j + A3*Aij-1 + A4*Cij+1 = -Cij(0)
	// membrane
	//coefficient matrix of the Fc diffusion
#pragma omp parallel for
	for (int i = 1; i < membrane.m + 1; ++i){
		for (int j = 1; j < membrane.n + 1; ++j){
			if (i > 1 && i < membrane.m  && j > 1 && j < membrane.n){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
			}
			else if (i > 1 && i < membrane.m && j == membrane.n){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*K_dis*membrane.matrixA[i][j][0];
			}
			else if (i > 1 && i < membrane.m && j == 1){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.beta1[j - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
			}
			else if (i == 1 && j > 1 && j < membrane.n){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa1[i - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
			}
			else if (i == 1 && j == membrane.n) {
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa1[i - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*K_dis*membrane.matrixA[i][j][0];
			}
			else if (i == membrane.m && j > 1 && j < membrane.n){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa2[i - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
			}
			else if (i == membrane.m && j == membrane.n){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa2[i - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*K_dis*membrane.matrixA[i][j][0];
			}
			else if (i == 1 && j == 1){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa1[i - 1] + membrane.beta1[j - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
			}
			else if (i == membrane.m && j == 1){
				membrane.matrixA[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa2[i - 1] + membrane.beta1[j - 1])*Esignal.dt*D_Fcm);
				membrane.matrixA[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
				membrane.matrixA[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcm*membrane.matrixA[i][j][0];
			}
		}
	}

	//coefficient matrix of the Fc+ diffusion
#pragma omp parallel for
	for (int i = 1; i < membrane.m + 1; ++i){
		for (int j = 1; j < membrane.n + 1; ++j){
			if (i > 1 && i < membrane.m  && j > 1 && j < membrane.n){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
			}
			else if (i > 1 && i < membrane.m && j == membrane.n){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*K_dis*membrane.matrixA2[i][j][0];
			}
			else if (i > 1 && i < membrane.m && j == 1){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.beta1[j - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
			}
			else if (i == 1 && j > 1 && j < membrane.n){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa1[i - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
			}
			else if (i == 1 && j == membrane.n) {
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa1[i - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*K_dis*membrane.matrixA2[i][j][0];
			}
			else if (i == membrane.m && j > 1 && j < membrane.n){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa2[i - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
			}
			else if (i == membrane.m && j == membrane.n){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa2[i - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][3] = membrane.beta1[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*K_dis*membrane.matrixA2[i][j][0];
			}
			else if (i == 1 && j == 1){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa1[i - 1] + membrane.beta1[j - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][2] = membrane.alfa2[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
			}
			else if (i == membrane.m && j == 1){
				membrane.matrixA2[i][j][0] = 1 / (1 - (membrane.alfa3[i - 1] + membrane.beta3[j - 1] + membrane.alfa2[i - 1] + membrane.beta1[j - 1])*Esignal.dt*D_Fcpm);
				membrane.matrixA2[i][j][1] = membrane.alfa1[i - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
				membrane.matrixA2[i][j][4] = membrane.beta2[j - 1] * Esignal.dt*D_Fcpm*membrane.matrixA2[i][j][0];
			}
		}
	}
}

void solution_matrixA(mesh& solution, const potential_signal& Esignal, const mesh& membrane, const double D_Fcw, const double D_Fcpw, const double K_dis)
{
	//diffusion equation: A0*Cij + A1*Ci-1j + A2*Ci+1j + A3*Aij-1 + A4*Cij+1 = Cij(0)
	//Cij = 1/A0*Cij(0) + (-A1/A0)Ci-1j + (-A2/A0)Ci+1j + (-A3/A0)Cij-1 + (-A4/A0)Cij+1
	//coefficient matrix of Fc diffusion
#pragma omp parallel for
	for (int i = 1; i < solution.m + 1; ++i){
		for (int j = 1; j < solution.n + 1; ++j){
			
			if (i > 1 && i < solution.m + 1 && j > 1 && j < solution.n + 1){
				solution.matrixA[i][j][0] = 1/(1 - (solution.alfa3[i - 1] + solution.beta3[j - 1])*Esignal.dt*D_Fcw);
				solution.matrixA[i][j][1] = solution.alfa1[i - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
			}
			else if (i == 1 && j > 1 && j < solution.n + 1){
				solution.matrixA[i][j][0] = 1/(1 - (solution.alfa3[i - 1] + solution.beta3[j - 1] + solution.alfa1[i - 1])*Esignal.dt*D_Fcw);
				solution.matrixA[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
			}
			else if (i > 1 && i < membrane.m + 1 && j == 1){
				solution.matrixA[i][j][0] = 1/(1 - (solution.alfa3[i - 1] + solution.beta3[j - 1])*Esignal.dt*D_Fcw);
				solution.matrixA[i][j][1] = solution.alfa1[i - 1] * Esignal.dt*D_Fcw *solution.matrixA[i][j][0];
				solution.matrixA[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcw / K_dis * solution.matrixA[i][j][0];
				solution.matrixA[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
			}
			else if (i > membrane.m && i < solution.m + 1 && j == 1){
				solution.matrixA[i][j][0] = 1 / (1 - (solution.alfa3[i - 1] + solution.beta3[j - 1] + solution.beta1[j - 1])*Esignal.dt*D_Fcw);
				solution.matrixA[i][j][1] = solution.alfa1[i - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
			}
			else if (i == 1 && j == 1){
				solution.matrixA[i][j][0] = 1 / (1-(solution.alfa3[i - 1] + solution.beta3[j - 1] + solution.alfa1[i - 1])*Esignal.dt*D_Fcw);
				solution.matrixA[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
				solution.matrixA[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcw / K_dis * solution.matrixA[i][j][0];
				solution.matrixA[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcw * solution.matrixA[i][j][0];
			}
		}
	}
	//coefficient matrix of Fc+ diffusion
#pragma omp parallel for 
	for (int i = 1; i < solution.m + 1; ++i){
		for (int j = 1; j < solution.n + 1; ++j){

			if (i > 1 && i < solution.m + 1 && j > 1 && j < solution.n + 1){
				solution.matrixA2[i][j][0] = 1 / (1 - (solution.alfa3[i - 1] + solution.beta3[j - 1])*Esignal.dt*D_Fcpw);
				solution.matrixA2[i][j][1] = solution.alfa1[i - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
			}
			else if (i == 1 && j > 1 && j < solution.n + 1){
				solution.matrixA2[i][j][0] = 1 / (1 - (solution.alfa3[i - 1] + solution.beta3[j - 1] + solution.alfa1[i - 1])*Esignal.dt*D_Fcpw);
				solution.matrixA2[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
			}
			else if (i > 1 && i < membrane.m + 1 && j == 1){
				solution.matrixA2[i][j][0] = 1 / (1 - (solution.alfa3[i - 1] + solution.beta3[j - 1])*Esignal.dt*D_Fcpw);
				solution.matrixA2[i][j][1] = solution.alfa1[i - 1] * Esignal.dt*D_Fcpw *solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcpw / K_dis * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
			}
			else if (i > membrane.m && i < solution.m + 1 && j == 1){
				solution.matrixA2[i][j][0] = 1 / (1 - (solution.alfa3[i - 1] + solution.beta3[j - 1] + solution.beta1[j - 1])*Esignal.dt*D_Fcpw);
				solution.matrixA2[i][j][1] = solution.alfa1[i - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
			}
			else if (i == 1 && j == 1){
				solution.matrixA2[i][j][0] = 1 / (1 - (solution.alfa3[i - 1] + solution.beta3[j - 1] + solution.alfa1[i - 1])*Esignal.dt*D_Fcpw);
				solution.matrixA2[i][j][2] = solution.alfa2[i - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][3] = solution.beta1[j - 1] * Esignal.dt*D_Fcpw / K_dis * solution.matrixA2[i][j][0];
				solution.matrixA2[i][j][4] = solution.beta2[j - 1] * Esignal.dt*D_Fcpw * solution.matrixA2[i][j][0];
			}
		}
	}
}

void membrane_arrayb(mesh& membrane)
{
#pragma omp parallel for
	for (int i = 1; i < membrane.m + 1; ++i){
		for (int j = 1; j < membrane.n + 1; ++j){
			membrane.arrayb[i][j]= membrane.Crdo[i - 1][j - 1];
			membrane.arrayb2[i][j] = membrane.Coxo[i - 1][j - 1];
		}
	}
}

void solution_arrayb(mesh& solution)
{
#pragma omp parallel for
	for (int i = 1; i < solution.m + 1; ++i){
		for (int j = 1; j < solution.n + 1; ++j){
			solution.arrayb[i][j] = solution.Crdo[i - 1][j - 1];
			solution.arrayb2[i][j] = solution.Coxo[i - 1][j - 1];
		}
	}
}

double solver(const int ti, const double& Eq, mesh& membrane, mesh& solution, nernst_equation& therm, const potential_signal& Esignal)
{
	double K = therm.ratio_ox2red(Eq);
	double tol_I = 1e-15; // use I as the threshold
	double I, Io(0), dI(10);
	double** mx(membrane.x), ** mxk(membrane.xk), ** mxo(membrane.xo), ** mxp(membrane.xp), ** mxkp(membrane.xkp), ** mxpo(membrane.xpo);
	double** sx(solution.x), ** sxk(solution.xk), ** sxo(solution.xo), ** sxp(solution.xp), ** sxkp(solution.xkp), ** sxpo(solution.xpo);
	// use the concentration at last time node
	// initialize mx
#pragma omp parallel for
	for (int i = 0; i < membrane.m; ++i){
		for (int j = 0; j < membrane.n; ++j){
			mxk[i + 1][j + 1] = membrane.Crdn[i][j];
			mxkp[i + 1][j + 1] = membrane.Coxn[i][j];
			mxo[i + 1][j + 1] = membrane.Crdo[i][j];
			mxpo[i + 1][j + 1] = membrane.Coxo[i][j];
		}
	}
	//initialize sx
#pragma omp parallel for
	for (int i = 0; i < solution.m; ++i){
		for (int j = 0; j < solution.n; ++j){
			sxk[i + 1][j + 1] = solution.Crdn[i][j];
			sxkp[i + 1][j + 1] = solution.Coxn[i][j];
			sxo[i + 1][j + 1] = solution.Crdo[i][j];
			sxpo[i + 1][j + 1] = solution.Coxo[i][j];
		}
	}
	
	///////////////////////////////iteration////////////////////////////////
	while (dI > tol_I){
		// membrane phase
#pragma omp parallel for
		for (int i = 1; i < membrane.m + 1; ++i) {
			for (int j = 1; j < membrane.n + 1; ++j) {
				if (j == 1) {
					mx[i][j] = membrane.matrixA[i][j][0] * mxo[i][j] + membrane.matrixA[i][j][1] * mxk[i - 1][j] + membrane.matrixA[i][j][2] * mxk[i + 1][j]
						+ membrane.matrixA[i][j][3] * mxk[i][j - 1] + membrane.matrixA[i][j][4] * mxk[i][j + 1];
					mxp[i][j] = membrane.matrixA2[i][j][0] * mxpo[i][j] + membrane.matrixA2[i][j][1] * mxkp[i - 1][j] + membrane.matrixA2[i][j][2] * mxkp[i + 1][j]
						+ membrane.matrixA2[i][j][3] * mxkp[i][j - 1] + membrane.matrixA2[i][j][4] * mxkp[i][j + 1];
					mx[i][j] = 1 / (1 + K)*(mx[i][j] + mxp[i][j]);
					mxp[i][j] = K*mx[i][j];
				}
				else if (j != membrane.n) {
					mx[i][j] = membrane.matrixA[i][j][0] * mxo[i][j] + membrane.matrixA[i][j][1] * mxk[i - 1][j] + membrane.matrixA[i][j][2] * mxk[i + 1][j]
						+ membrane.matrixA[i][j][3] * mxk[i][j - 1] + membrane.matrixA[i][j][4] * mxk[i][j + 1];
					mxp[i][j] = membrane.matrixA2[i][j][0] * mxpo[i][j] + membrane.matrixA2[i][j][1] * mxkp[i - 1][j] + membrane.matrixA2[i][j][2] * mxkp[i + 1][j]
						+ membrane.matrixA2[i][j][3] * mxkp[i][j - 1] + membrane.matrixA2[i][j][4] * mxkp[i][j + 1];
				}
				else if (j == membrane.n) {
					mx[i][j] = membrane.matrixA[i][j][0] * mxo[i][j] + membrane.matrixA[i][j][1] * mxk[i - 1][j] + membrane.matrixA[i][j][2] * mxk[i + 1][j]
						+ membrane.matrixA[i][j][3] * mxk[i][j - 1] + membrane.matrixA[i][j][4] * sxk[i][1];
					mxp[i][j] = membrane.matrixA2[i][j][0] * mxpo[i][j] + membrane.matrixA2[i][j][1] * mxkp[i - 1][j] + membrane.matrixA2[i][j][2] * mxkp[i + 1][j]
						+ membrane.matrixA2[i][j][3] * mxkp[i][j - 1] + membrane.matrixA2[i][j][4] * sxkp[i][1];
				}
			}
		}
		// solution phase
#pragma omp parallel for
		for (int i = 1; i < solution.m + 1; ++i) {
			for (int j = 1; j < solution.n + 1; ++j) {
				if (j == 1 && i < membrane.m + 1) {
					sx[i][j] = solution.matrixA[i][j][0] * sxo[i][j] + solution.matrixA[i][j][1] * sxk[i - 1][j] + solution.matrixA[i][j][2] * sxk[i + 1][j]
						+ solution.matrixA[i][j][3] * mxk[i][membrane.n] + solution.matrixA[i][j][4] * sxk[i][j + 1];
					sxp[i][j] = solution.matrixA2[i][j][0] * sxpo[i][j] + solution.matrixA2[i][j][1] * sxkp[i - 1][j] + solution.matrixA2[i][j][2] * sxkp[i + 1][j]
						+ solution.matrixA2[i][j][3] * mxkp[i][membrane.n] + solution.matrixA2[i][j][4] * sxkp[i][j + 1];
				}
				else {
					sx[i][j] = solution.matrixA[i][j][0] * sxo[i][j] + solution.matrixA[i][j][1] * sxk[i - 1][j] + solution.matrixA[i][j][2] * sxk[i + 1][j]
						+ solution.matrixA[i][j][3] * sxk[i][j - 1] + solution.matrixA[i][j][4] * sxk[i][j + 1];
					sxp[i][j] = solution.matrixA2[i][j][0] * sxpo[i][j] + solution.matrixA2[i][j][1] * sxkp[i - 1][j] + solution.matrixA2[i][j][2] * sxkp[i + 1][j]
						+ solution.matrixA2[i][j][3] * sxkp[i][j - 1] + solution.matrixA2[i][j][4] * sxkp[i][j + 1];
				}
			}
		}

		// calculate current
		I = 0;
#pragma omp parallel for reduction(+ : I)
		for (int i = 1; i < membrane.m + 1; ++i) {
			I += membrane.RR[i - 1][0] *(mxp[i][1] - mxp[i][2]);
		}
		dI = I - Io;
		Io = I;
		//pass the result to medium container
#pragma omp parallel for
		for (int i = 1; i < membrane.m + 1; ++i) {
			for (int j = 1; j < membrane.n + 1; ++j) {
				mxk[i][j] = mx[i][j];
				mxkp[i][j] = mxp[i][j];
			}
		}
		for (int i = 1; i < solution.m + 1; ++i) {
			for (int j = 1; j < solution.n + 1; ++j) {
				sxk[i][j] = sx[i][j];
				sxkp[i][j] = sxp[i][j];
			}
		}
	}

	I *= 6.28*membrane.dr / (membrane.dz + membrane.dz0) * 2 * therm.n*therm.F/1000; // diffusion coefficient is multiplied in main() function

#pragma omp parallel for
	for (int i = 0; i < membrane.m; ++i) {
		for (int j = 0; j < membrane.n; ++j) {
			membrane.Crdn[i][j] = mx[i + 1][j + 1];
			membrane.Coxn[i][j] = mxp[i + 1][j + 1];
		}
	}
	for (int i = 0; i < solution.m; ++i) {
		for (int j = 0; j < solution.n; ++j) {
			solution.Crdn[i][j] = sx[i + 1][j + 1];
			solution.Coxn[i][j] = sxp[i + 1][j + 1];
		}
	}
	
	return I;
}