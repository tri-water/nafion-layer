#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <omp.h>
#include <memory>
#include "thermodynamics.h"
#include "mesh.h"
#include "potential_signal.h"
#include "utilities.h"
#include <windows.h>

using namespace std;

const double D_Fcw = 1e-5; // diffusion coefficient of fc in water, cm2s-1 (Mandal, 1993)
const double D_Fcpw = 1e-5; // diffusion coefficient of fc+ ion in water, cm2s-1
const double D_Fcm = 1e-9; // apparent diffusion coefficient of fc in nafion membrane, cm2s-1 (Dong, 1990)
const double D_Fcpm = 1e-9; // apparent diffusion coefficient of fc+ ion in nafion membrane, cm2s-1
const double K_dis = 2000; // distribution coefficient between the membrane and the solution

int main()
{
	//////////////////////////////parameters//////////////////////////
	const double Cm_bulk = 1; // the bulk concentration of the reduced active species in the membrane, mol
	const double Cs_bulk = 0; // the bulk concentration of the reduced active species in the bulk of the solution phase
	potential_signal Esignal(-0.5, 0.3, 0.002, 0.0001, 25, 0.02);

	//////////////////////////////recorder////////////////////////////
	double Eq; // applied electrode potential at time q, V

	//////////////////////////////define meshgrid/////////////////////
	mesh membrane(50, 10e-4, 10e-4, 200, 0.05e-4, 0.05e-4); // the unit of dr and dz is cm
	mesh solution(100, 10e-4, 10e-4, 250, 4e-4, 4e-4); // the unit of dr and dz is cm
	membrane_matrixA(membrane, Esignal, D_Fcm, D_Fcpm, K_dis);
	solution_matrixA(solution, Esignal, membrane, D_Fcw, D_Fcpw, K_dis);
	membrane.print_mesh("membrane");
	solution.print_mesh("solution");

	//////////////////////////////set thermodynamics//////////////////
	nernst_equation therm(0, 1); // set the formal potential at 0 and the number of transfered electrons as 1

	//////////////////////////////start the time loop/////////////////
	for (int i = 0; i < Esignal.q; ++i){

	//////////////////////////////set electrode potential/////////////
		Eq = Esignal.applied_potential(i);
		// initialise
		if (i == 0){
			initialize(membrane, solution, therm, Eq, Cm_bulk);
		}
		// start
		else{
			double Iq = solver(i, Eq, membrane, solution, therm, Esignal)*D_Fcpm;

			//////////////////////record/////////////////////////////
			if (Esignal.wave_judge == 1 && Esignal.wave_judge1 == -1) {
				Esignal.Iqr[0] = Iq;
			}
			else if (Esignal.wave_judge == -1 && Esignal.wave_judge1 == 1) {
				Esignal.Iqr[1] = Iq;
				Esignal.record_data();
				Esignal.save_peak_concentration(membrane, solution);
				++Esignal.period_counter;
			}

			membrane.Ca2Cb("Crdn", "Crdo");
			membrane.Ca2Cb("Coxn", "Coxo");
			solution.Ca2Cb("Crdn", "Crdo");
			solution.Ca2Cb("Coxn", "Coxo");
		}
	}
	Esignal.save_current();
	cout << "Completed. Enter any key to close the window\n";
	cin.get();
	return 0;
}