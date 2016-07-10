#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#include <string>

using namespace std;

class mesh{
public:
	// define x and y interval
	const int m; // r axis node number
	const double dr0; // , cm
	const double dr; // , cm
	const int n; // z axis node number
	const double dz; // , cm
	const double dz0; // , cm

	// mesh define
	double **Crdn;
	double **Crdo;
	double **Coxn;
	double **Coxo;
	double **R; // r node for plot graph
	double **Z; // z node for plot graph
	double **RR; // the r position of the centre of each box
	double **ZZ; // the z position of the centre of each box
	double ***matrixA; // Ax = b for Fc, each element inside is initialised as 0
	double ***matrixA2; // Ax = b for Fc+, each element inside is initialised as o
	double **arrayb; // Ax = b for Fc; each element inside is initialised as 0
	double **arrayb2; // Ax = b for Fc+, each element inside is initialised as 0
	double dCrd; // the change of reduced species concentration at one point

	double *alfa1; // geometric coefficients, r direction
	double *alfa2;
	double *alfa3;
	double *beta1; // geometric coefficients, z direction
	double *beta2;
	double *beta3;

	double **x, **xk, **xo, **xp, **xkp, **xpo; //medium container for iteration

	// Constructor
	mesh(int fm, double fdr0, double fdr, int fn, double fdz0, double fdz);
	~mesh();

	// Create a mesh grid
	void mesh_grid();
	// Pass Crdn to Crdo
	void Ca2Cb(string Ca_name, string Cb_name);
	// Export the mesh grid
	void print_mesh(string ex_file_name);
	// Export the concentrations
	void print_concentration(string target_concentration_name, string ex_file_name);
};

#endif