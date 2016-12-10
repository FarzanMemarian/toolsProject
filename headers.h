//headers.h
#ifndef ADD_H_INCLUDED
#define ADD_H_INCLUDED

void  myEuler(double h, double y0 , int maxTime );

void analyticalEuler(double h, double y0 , int maxTime);

int func (double t, const double y[], double f[], void *params);

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

void gslSolver(double, int);


int func2 (double t, const double y[], double f[], void *params);

int jac2 (double t, const double y[], double *dfdy, double dfdt[], void *params);

void gslSolver_rk4(double, int, std::string&);

void gslSolver_rk2(double, int, std::string&);

void gslSolver_rkf45(double, int, std::string&);


GRVY::GRVY_Timer_Class gt;

#endif
