#include "potential_signal.h"
#include <cmath>
#include <iostream>
#include <fstream>

potential_signal::potential_signal(double fE0, double fEend, double fdE, double fdt, int fswf, double fswamp) :
E0(fE0), Eend(fEend), dE(fdE), dt(fdt), swf(fswf), swamp(fswamp), q(int(ceil((fEend - fE0) / (fswf*fdE*fdt)))), period_counter(0)
{
	Eqm = E0;
	int recd_size = int(ceil(q*dt*swf) + 5);
	Er = new double[recd_size]; // recorded electrode potential container
	tr = new double[recd_size]; // recorded time container
	Is = new double[recd_size]; // recorded square wave current container
	Io = new double[recd_size]; // recorded oxidization current container
	last_Is = new double[recd_size]; // recorded square wave current container of last period
	last_Io = new double[recd_size]; // recorded oxidization current container of last period
}

potential_signal::~potential_signal()
{
	delete Er;
	delete tr;
	delete Is;
	delete Io;
	delete last_Is;
	delete last_Io;
}

double potential_signal::applied_potential(int i)
{
	double tq = (i - 1)*dt; // current time
	double tq1 = i*dt; // next time
	double Eq; // applied potential
		
	// define the potential vs time
	// square wave voltammetry
	if (i == 0){
		Eqm = E0;
		Eq = E0;
	}
	else{
		if (wave_judge == -1 && wave_judge1 == 1){
			Eqm = Eqm + dE; // provide base potential for the current time node
		}
		
		// reset wave_judge for the next time node
		// calculate Eq
		if (tq*swf - floor(tq*swf) < 0.5){
			wave_judge = 1;
			Eq = Eqm + swamp;
		}
		else{
			wave_judge = -1;
			Eq = Eqm - swamp;
		}
	}

	// reset wave_judge1 for the next time node
	if (tq1*swf - floor(tq1*swf) < 0.5){
		wave_judge1 = 1;
	}
	else{
		wave_judge1 = -1;
	}

	return Eq;
}

void potential_signal::record_data()
{
	Is[period_counter] = Iqr[0] - Iqr[1];
	Io[period_counter] = Iqr[0];
	Er[period_counter] = Eqm;
	cout << Er[period_counter] << "V " << Is[period_counter] << "A\n ";
}

void potential_signal::save_peak_concentration(mesh& membrane, mesh& solution)
{
	if (Is[period_counter - 1] > Is[period_counter] && Is[period_counter - 1] > Is[period_counter - 2]) {
		membrane.print_concentration("Crdo", "membrane_reduced_speices_concentration@peak");
		membrane.print_concentration("Coxo", "membrane_oxidised_speices_concentration@peak");
		solution.print_concentration("Crdo", "solution_reduced_speices_concentration@peak");
		solution.print_concentration("Coxo", "solution_oxidised_speices_concentration@peak");
	}
}

void potential_signal::save_current()
{
	ofstream fcout("current.txt");

	if (fcout.is_open()) {
		fcout << "parameters: "
			<< "\nstarting potential, V: " << E0
			<< "\nend potential, V: " << Eend
			<< "\nstep potential, V: " << dE
			<< "\nfrequency, Hz: " << swf
			<< "\namplitude, V: " << swamp;

		fcout << "\npotential/V " << "SWV_Current/A " << "Ox_Current/A\n";

		for (int i = 0; i < period_counter + 1; ++i) {
			fcout << Er[i] << " " << Is[i] << " " << Io[i] << "\n";
		}
	}
	else {
		cerr << "Cannot open the file to save current";
		exit(1);
	}
}