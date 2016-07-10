#ifndef POTENTIAL_SIGNAL_H_INCLUDED
#define POTENTIAL_SIGNAL_H_INCLUDED
#include <memory>
#include "mesh.h"

using namespace std;

class potential_signal{
public:
	const double E0; // starting potential, V
	const double Eend; // end potential, V
	const double dE; // step potential, V
	const double dt; // delta time between each time node, s
	const int swf; // the frequency of square wave potential, Hz
	const double swamp; // the amplitude of square wave potential, V
	const int q; // number of time nodes
	double Eqm; // the base potential
	double* Er; // recorded electrode potential container
	double* tr; // recorded time container
	double* Is; // recorded square wave current container
	double* Io; // recorded oxidization current container
	double* last_Is; // recorded square wave current container of last period
	double* last_Io; // recorded oxidization current container of last period

	int period_counter; // count the number of period
	int wave_judge; // determine which half wave the past time is in
	int wave_judge1; // determine which half wave the comming time node is in
	double Iqr[2] = { 0, 0}; // record the current at the  end of each half period

	potential_signal(double fE0, double fEend, double fdE, double fdt, int fswf, double fswamp);
	~potential_signal();
	double applied_potential(int i);
	void record_data();
	void save_peak_concentration(mesh& membrane, mesh& solution);
	void save_current();

};
#endif