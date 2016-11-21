//**************************************************************
// this program is for calculating the solution to the first 
// order differntial equation of the form y'=f(y,t)
// in this project we solve y'=y, y0=0
// the analytical solusion is y=e^t
// WRITTEN BY FARZAN MEMARIAN, FALL 2016
//************************************************************** 

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

using namespace std;

std::vector<double>  myEuler(vector<double> & y, double h, double y0, double nEuler);

void analyticalEuler(double h, double y0, double nEuler);

int main() {
	double y0=1, h=0.01, nEuler = 100; //initialize variables
	std::vector<double> y(nEuler);
	myEuler(y, h, y0, nEuler);
	cout << "these are values of y obtained from myEuler" << endl;
	for (int i=0; i<nEuler; i++){
		cout << y[i] << endl;
	}
	cout << "these are values of y obtained from analytical solution" << endl;
	analyticalEuler(h, y0, nEuler);
	// odeSover;




	return 0;

}


std::vector<double>  myEuler(vector<double> & y, double h, double y0, double nEuler)
{
	y[0]=y0;
	for (int t=0; t<nEuler; t++){
		y[t+1] = y[t] + h * y[t];
	}
	return y;
}


void analyticalEuler(double h, double y0, double nEuler)
{
	// this functin calculateds the analytical solution of the
	// first order ODE y'=y(t)
	double y=0; // the dependent variable
	double t=0; // the independent variable
	for (int i=0; i<nEuler; i++){
		t=i*h;
		y = exp(t);
		cout << y << endl;
	}
}

int func (double t, const double y[], double f[], void *params)
{
	(void)(t); 
	f[0] = y[1];
	return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t);
	// jacobian matrix, df/dy 
	*(dfdy) = 1.0;
	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array (dfdy, 0, 0);
	gsl_matrix_set (m, 0, 0, 0.0);
	dfdt[0] = 1.0;
	return GSL_SUCCESS;
}

void odeSolver(){
sage: T=ode_solver()
}

