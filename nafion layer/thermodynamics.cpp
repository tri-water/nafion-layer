#include "thermodynamics.h"
#include <cmath>

nernst_equation::nernst_equation(double fE_formal, int fn):
E_formal(fE_formal), R(8.3144621), T(293), F(9.64853399 * 10000), n(fn)
{
	alfa = n*F / R / T;
}

double nernst_equation::ratio_ox2red(double E){
	double k = exp(alfa*(E - E_formal));
	// return the Cox/Cred ratio
	return k;
	
}