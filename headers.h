//headers.h
#ifndef ADD_H_INCLUDED
#define ADD_H_INCLUDED

double  myEuler(std::vector<double>&, double h, double y0 , int maxTime );

double analyticalEuler(std::vector<double>& , double h, double y0 , int maxTime);

int func (double t, const double y[], double f[], void *params);

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

double gslSolver(double, int);


int func2 (double t, const double y[], double f[], void *params);

int jac2 (double t, const double y[], double *dfdy, double dfdt[], void *params);

std::vector<double> gslSolver_rk4(  double, int, std::string&);

std::vector<double> gslSolver_rk2(  double, int, std::string&);

std::vector<double> gslSolver_rkf45( double, int, std::string&);

void verificationFunc1( double, double, double,  double , int );
void verificationFunc2 (std::vector<double> , double , int , std::string&);


GRVY::GRVY_Timer_Class gt;

#endif
