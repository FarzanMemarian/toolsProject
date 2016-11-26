//headers.h
#ifndef ADD_H_INCLUDED
#define ADD_H_INCLUDED

std::vector<double>  myEuler(vector<double> & y, double h, double y0, double nEuler);

void analyticalEuler(double h, double y0, double nEuler);

int func (double t, const double y[], double f[], void *params);

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

void odeSolver();


int func2 (double t, const double y[], double f[], void *params);

int jac2 (double t, const double y[], double *dfdy, double dfdt[], void *params);

void odeSolver2();

#endif 
