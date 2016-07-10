#ifndef THERMODYNAMICS_H_INCLUDED
#define THERMODYNAMICS_H_INCLUDED

class nernst_equation{
public:
	const double E_formal; // standard potential, V
	const double R; // gas constant, J K^-1 mol^-1
	const double T; // temperature, K
	const double F; // Faraday constant, C mol^-1
	const int n; // the number of electron transfered
	double phi;
	double alfa; // nF/RT, V^-1

	// constructor
	nernst_equation(double fE_formal, int fn);
	// destructor
	~nernst_equation(){};

	// calculate the Red concentration under thermodynamic equilibarium
	// E: the potential subjected to the reaction
	// C_total: the total concentration of the redox pair
	double ratio_ox2red(double E); 
};

#endif