//**************************************************************
//
// this program is for calculating the solution to the first 
// order differntial equation of the form y'=f(y,t)
// in this project we solve y'=y, y0=0
// the analytical solusion is y=e^t
// WRITTEN BY FARZAN MEMARIAN, FALL 2016
//
//************************************************************** 

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <grvy.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#define FUNC_BEGIN_TIMER gt.BeginTimer(__func__);
#define FUNC_END_TIMER   gt.EndTimer  (__func__);

using namespace std;
#include "headers.h"


int numberPrint = 1;
int debug =0;
int main(int argc, char *argv[]) 
{
	// in main, the user is asked to decide which problem to choose. 
	// then for each problem, the corresponding functions are called

	// reading input using GRVY
	using namespace GRVY;

	// Initialize timing library, the global library is initialized 
	// with this call
	gt.Init("GRVY Timing");

	GRVY_Input_Class iparse;     // input parsing object
	int problem, verification;
	double h, maxTime;
	string odeMethod, inputName;

	// Initi9alize/read the file

	if (argc > 1){
		inputName = argv[1];
	}
	// Initialize/read the file 
	if(! iparse.Open(inputName.c_str()))
		exit(1);


	// Read specific variables 

	if( iparse.Read_Var("problem",&problem) );
	if( iparse.Read_Var("debug",&debug) );
	if( iparse.Read_Var("verification",&verification) );
	if( iparse.Read_Var("h",&h) );
	if( iparse.Read_Var("maxTime",&maxTime) );
	if( iparse.Read_Var("odeMethod",&odeMethod) );
	if( iparse.Read_Var("numberPrint",&numberPrint) );

	if (debug == 2)
	{
		printf("--> %-11s = %i\n","problem",problem);
		printf("--> %-11s = %i\n","debug",debug);
		printf("--> %-11s = %i\n","verification",verification);
		printf("--> %-11s = %f\n","h",h);
		printf("--> %-11s = %i\n","maxTime",maxTime);
		printf("--> %-11s = %s\n","odeMethod",odeMethod.c_str());
		printf("--> %-11s = %i\n","numberPrint",numberPrint);
	}

	if (problem == 1) 
	{
		// here we choose some of the parameters for the first problem
		double y0=1; // value of y at t_0=0

		// Define the beginning of the portion being timed 
		gt.BeginTimer("myEuler");

		int length = maxTime / (h*numberPrint);
		std::vector<double> y(length);		
		double lastEuler, lastAnalytical, lastGSL; // this is the last value for y
		lastEuler = myEuler(y ,h, y0, maxTime); // This is the function I wrote for forward euler
		// Define the beginning of the portion being timed 
		gt.EndTimer("myEuler");
		gt.Finalize();  // Finalize the myEuler Timer
		if (debug >= 1) { gt.Summarize(); } // Print performance summary to stdout
		gt.Reset();     // Reset timers for next iteration
		
		// gslSolver makes use of gsl ode solver to solve the y'=y
		// Define the beginning of the portion being timed 
		gt.BeginTimer("gslSolverProblem1");
		lastGSL = gslSolver(h,  maxTime); // this function uses gsl solver to solve problem 1 
		if (verification == 1) {
			lastAnalytical = analyticalEuler(y,h,y0,maxTime);
			verificationFunc1(lastEuler, lastAnalytical, lastGSL, h, maxTime);
		}
		// Define the beginning of the portion being timed 
		gt.EndTimer("gslSolverProblem1");
		gt.Finalize();  // Finalize the myEuler Timer
		if (debug >= 1) { gt.Summarize(); } // Print performance summary to stdout
		gt.Reset();     // Reset timers for next iteration
	}
	else if(problem == 2)
	{
		// this is for the second problem

		std::vector<double> last2(3);
		gt.BeginTimer("rk4 or rk2 or rkf45");
		if (odeMethod == "rk4") 
		{ 
			last2 = gslSolver_rk4( h, maxTime, odeMethod); 
		}	
		else if (odeMethod == "rk2") 
		{ 
			last2 = gslSolver_rk2( h, maxTime, odeMethod); 
		}	
		else if (odeMethod == "rkf45") 
		{ 
			last2 =gslSolver_rkf45(h, maxTime, odeMethod);
		}	
		else 
		{ 
			cout << "wrong method name, please enter either rk2 or rk4 or rkf45 in the input file " << endl;
		}
		// Define the beginning of the portion being timed 
		gt.EndTimer("rk4 or rk2 or rkf45");
		gt.Finalize();  // Finalize the myEuler Timer
		if (debug >= 1) { gt.Summarize(); } // Print performance summary to stdout
		gt.Reset();     // Reset timers for next iteration
		
		if (verification == 1) 
		{
			verificationFunc2(last2, h, maxTime, odeMethod);
		}
	}
	else
	{
		cout << "enter either 1 or 2 for problem in the input file" << endl;
	}

	return 0;
}

void verificationFunc2 (std::vector<double> last2, double h, int maxTime, std::string& odeMethod)
{
	double SmallStepSize = 0.00001;
	std::vector<double> lastExact = gslSolver_rk4(SmallStepSize, maxTime, odeMethod);	
	double err = sqrt( pow(last2[0]-lastExact[0] ,2) + pow(last2[1]-lastExact[1] ,2) + pow(last2[1]-lastExact[1],2) );
	int n = maxTime/h;
	printf ("%10d  %.14e \n", n, err);
}


void verificationFunc1 (double lastEuler, double lastAnalytical, double lastGSL, double h, int maxTime)
{
	double errEuler = abs(lastEuler - lastAnalytical);
	double errGSL = abs(lastGSL - lastAnalytical);
	int n = maxTime/h;	
	printf ( "%14d %.14e %.14e \n ", n , errEuler, errGSL );
}




double myEuler(std::vector<double>& y, double h, double y0, int maxTime)
{
	if (debug == 2) {cout << "the program is in function myEuler now (Problem 1)" << endl;}
	int nEuler = maxTime / h;
	y[0]=y0;
	for (int t=0; t<nEuler; t++){
		y[t+1] = y[t] + h * y[t];
	}

	// write debug output to screen
	if (debug == 2){
		cout << "these are values of y obtained from myEuler" << endl;
		for (int i=0; i<=nEuler; i++){
			printf ("%.14e %.14e \n", i*h , y[i]);
		}
	}

	// write output to a file
	FILE * pFile;
	pFile = fopen ("prob1_MyEuler.dat","w");

	for(int i=0; i<=nEuler; i++){
		if(i % numberPrint == 0) {
			fprintf (pFile, "%.14e     %.14e     %.14e    \n",i*h , y0+i*h,y[i]);
		}
	}
	fclose(pFile);
	double last;
	last = y[nEuler];
	return last;	
}

double analyticalEuler(std::vector<double>& y1, double h, double y0, int maxTime)
{
	// this functin calculateds the analytical solution of the
	// first order ODE y'=y(t)
	if (debug == 2) {cout << "the program is in function analyticalEuler now (Problem 1)" << endl;}
	double y=0; // the dependent variable
	double t=0; // the independent variable

	int nEuler = maxTime / h;
	FILE * pFile;
	pFile = fopen("prob1_AlyticalEuler.dat","w");

	if (debug == 2){
		cout << "these are values of y obtained from analytical solution, first column is time" << endl;
	}
	for (int i=0; i<=nEuler; i++){
		t=i*h;
		y = exp(t);
		y1[i] = y;
		if (debug == 2){
			printf ("%.14e %.14e \n", t, y);
		}
		if (i % numberPrint == 0) {
			fprintf (pFile, "%.14e     %.14e     %.14e   \n",i*h, t, y);
		}
	}
	fclose(pFile);
	double last;
	last = y1[nEuler];
	return last;
}

// THIS IS THE GSL VERSION
int func (double t, const double y[], double f[], void *params)
{
	(void)(t); 
	f[0] = y[0];
	return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t);
	// Jacobian matrix, df/dy 
	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array (dfdy, 0, 0);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set (m, 0, 0, 1.0);
	dfdt[0] = 1.0;
	return GSL_SUCCESS;
}




double gslSolver(double  h, int maxTime)
{
	if (debug == 2) {cout << "the program is in function gslSolver now (Problem 1)" << endl;}
	// this function solves the ODE using gsl library
	gsl_odeiv2_system sys = {func, jac, 1};

	gsl_odeiv2_driver * d = 
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
				h,10.0 , 10.0);

	int i,status;
	double t = 0.0; 
	double y[1] = { 1.0 };
	int nIter = maxTime/h;
	double last;

	FILE * pFile;
	pFile = fopen ("prob1_GSLSolver.dat","w");
	fprintf (pFile, "%.14e     %.14e     \n",t , 1.000);
	for (i = 1; i <= nIter; i++)
	{
		status =  gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y);
		last = y[0];
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		if (debug == 2){
			if (i % numberPrint == 0){ printf ("%.14e %.14e \n", t, y[0]);}
		}
		if (i % numberPrint == 0) {
			fprintf (pFile, "%.14e     %.14e     \n",t , y[0]);
		}
	}
	fclose(pFile);
	gsl_odeiv2_driver_free (d);
	return last;
}



// THIS WAS AN ATTEMPT TO VERIFY BASED ON RMS
//void verificationFunc (std::vector<double>& y, std::vector<double>& y1, int n, int m) 
//{
//	if (m == 1)
//	{
//		double sum = 0;
//		double tol = 1e-6;
//		for (int i = 0; i < n; i++) 
//		{
//			sum = sum + abs(y1[i]-y[i]);	
//		}	
//		cout << "the error between my euler function and the analytical solution is:" << endl;
//		cout << sum << endl;
//
//	} 
//	else if(m == 2)
//	{
//		std::vector<double> diff[n]
//		double sum2 = 0;
//		for (int j = 0; j < n; j++) 
//		{
//			diff[j] = y[j] - y1[j]		
//		}
//		for (int jj = 0, jj < (; jj++)
//		{
//			sum2 = sum2 + 
//		}
//	
//	}
//}


int func2 (double t, const double y[], double f[],
		void *params)
{
	(void)(t); /* avoid unused parameter warning */
	int w = 5;
	int T = 5;
	f[0] = y[3];
	f[1] = y[4];
	f[2] = y[5];
	f[3] = w*y[4]-y[3]/T;
	f[4] = -w*y[3] - y[4]/T;
	f[5] = -y[5]/T;
	return GSL_SUCCESS;
}

int jac2 (double t, const double y[], double *dfdy, 
		double dfdt[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array (dfdy, 6, 6);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 0.0);
	gsl_matrix_set (m, 0, 2, 0.0);
	gsl_matrix_set (m, 0, 3, 1.0);
	gsl_matrix_set (m, 0, 4, 0.0);
	gsl_matrix_set (m, 0, 5, 0.0);
	gsl_matrix_set (m, 1, 0, 0.0);
	gsl_matrix_set (m, 1, 1, 0.0);
	gsl_matrix_set (m, 1, 2, 0.0);
	gsl_matrix_set (m, 1, 3, 0.0);
	gsl_matrix_set (m, 1, 4, 1.0);
	gsl_matrix_set (m, 1, 5, 0.0);
	gsl_matrix_set (m, 2, 0, 0.0);
	gsl_matrix_set (m, 2, 1, 0.0);
	gsl_matrix_set (m, 2, 2, 0.0);
	gsl_matrix_set (m, 2, 3, 0.0);
	gsl_matrix_set (m, 2, 4, 0.0);
	gsl_matrix_set (m, 2, 5, 1.0);
	gsl_matrix_set (m, 3, 0, 0.0);
	gsl_matrix_set (m, 3, 1, 0.0);
	gsl_matrix_set (m, 3, 2, 0.0);
	gsl_matrix_set (m, 3, 3, -0.2);
	gsl_matrix_set (m, 3, 4, 5.0);
	gsl_matrix_set (m, 3, 5, 0.0);
	gsl_matrix_set (m, 4, 0, 0.0);
	gsl_matrix_set (m, 4, 1, 0.0);
	gsl_matrix_set (m, 4, 2, 0.0);
	gsl_matrix_set (m, 4, 3, -5.0);
	gsl_matrix_set (m, 4, 4, -0.2);
	gsl_matrix_set (m, 4, 5, 0.0);
	gsl_matrix_set (m, 5, 0, 0.0);
	gsl_matrix_set (m, 5, 1, 0.0);
	gsl_matrix_set (m, 5, 2, 0.0);
	gsl_matrix_set (m, 5, 3, 0.0);
	gsl_matrix_set (m, 5, 4, 0.0);
	gsl_matrix_set (m, 5, 5, -0.2);
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;
	dfdt[4] = 0.0;
	dfdt[5] = 0.0;

	return GSL_SUCCESS;
}


std::vector<double> gslSolver_rk2( double h, int maxTime, std::string& odeMethod)
{
	if (debug == 2) {
		cout << "the program is in function gslSolver_rk2 now " << endl;
		cout << "these are the results: time, x, y, z" << endl;
	}
	gsl_odeiv2_system sys = {func2, jac2, 6};

	gsl_odeiv2_driver * d = 
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk2,
				h, 100.0, 100.0);

	int i, status;
	double t = 0.0;
	int nIter = maxTime/h;
	double y[6] = {0.0, 0.0, 0.0, 20.0, 0.0, 2.0};

	std::vector<double> last(3);
	
	FILE * pFile;
	pFile = fopen ("prob2_rk2.dat","w");
	fprintf (pFile, "%.14e     %.14e     %.14e     %.14e     %.14e     %.14e     %.14e\n",t , 0.00, 0.00, 0.00, 20.0, 0.0, 1.0);	
	for (i = 1; i <= nIter; i++)
	{
		int status = gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y);
		last[0] = y[0];
		last[1] = y[1];
		last[2] = y[2];
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		if (debug == 2)
		{
			printf ("%.14e     %.14e     %.14e     %.14e\n", t, y[0], y[1], y[2]);
		}		

		//myfile  << t << y[0] << y[1] << y[2] << y[3] << y[4] << y[5] << endl;
		fprintf (pFile, "%.14e     %.14e     %.14e     %.14e     %.14e     %.14e     %.14e\n", t, y[0], y[1], y[2], y[3], y[4], y[5]);
		//printf ("%.14e %.14e %.14e %.14e\n", t, y[0], y[1], y[2]);
	}
	
	fclose(pFile);
	gsl_odeiv2_driver_free (d);
	return last;
}



std::vector<double> gslSolver_rk4( double h, int maxTime, std::string& odeMethod)
{
	if (debug == 2) {
		cout << "the program is in function gslSolver_rk4 now " << endl;
		cout << "these are the results: time, x, y, z" << endl;
	}
	gsl_odeiv2_system sys = {func2, jac2, 6};

	gsl_odeiv2_driver * d = 
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
				h, 100.0, 100.0);

	int i, status;
	double t = 0.0;
	int nIter = maxTime/h;
	double y[6] = {0.0, 0.0, 0.0, 20.0, 0.0, 2.0};

	std::vector<double> last(3);
	
	FILE * pFile;
	pFile = fopen ("prob2_rk4.dat","w");
	fprintf (pFile, "%.14e     %.14e     %.14e     %.14e     %.14e     %.14e     %.14e\n",t , 0.00, 0.00, 0.00, 20.0, 0.0, 1.0);	
	for (i = 1; i <= nIter; i++)
	{
		int status = gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y);
		last[0] = y[0];
		last[1] = y[1];
		last[2] = y[2];
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		if (debug == 2)
		{
			printf ("%.14e     %.14e     %.14e     %.14e\n", t, y[0], y[1], y[2]);
		}		

		//myfile  << t << y[0] << y[1] << y[2] << y[3] << y[4] << y[5] << endl;
		fprintf (pFile, "%.14e     %.14e     %.14e     %.14e     %.14e     %.14e     %.14e\n", t, y[0], y[1], y[2], y[3], y[4], y[5]);
		//printf ("%.14e %.14e %.14e %.14e\n", t, y[0], y[1], y[2]);
	}
	
	fclose(pFile);
	gsl_odeiv2_driver_free (d);
	return last;
}



std::vector<double> gslSolver_rkf45( double h, int maxTime, std::string& odeMethod)
{
	if (debug == 2) {
		cout << "the program is in function gslSolver_rkf45 now " << endl;
		cout << "these are the results: time, x, y, z" << endl;
	}
	gsl_odeiv2_system sys = {func2, jac2, 6};

	gsl_odeiv2_driver * d = 
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				h, 100.0, 100.0);

	int i, status;
	double t = 0.0;
	int nIter = maxTime/h;
	double y[6] = {0.0, 0.0, 0.0, 20.0, 0.0, 2.0};

	std::vector<double> last(3);
	
	FILE * pFile;
	pFile = fopen ("prob2_rkf45.dat","w");
	fprintf (pFile, "%.14e     %.14e     %.14e     %.14e     %.14e     %.14e     %.14e\n",t , 0.00, 0.00, 0.00, 20.0, 0.0, 1.0);	
	for (i = 1; i <= nIter; i++)
	{
		int status = gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y);
		last[0] = y[0];
		last[1] = y[1];
		last[2] = y[2];
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		if (debug == 2)
		{
			printf ("%.14e     %.14e     %.14e     %.14e\n", t, y[0], y[1], y[2]);
		}		

		//myfile  << t << y[0] << y[1] << y[2] << y[3] << y[4] << y[5] << endl;
		fprintf (pFile, "%.14e     %.14e     %.14e     %.14e     %.14e     %.14e     %.14e\n", t, y[0], y[1], y[2], y[3], y[4], y[5]);
		//printf ("%.14e %.14e %.14e %.14e\n", t, y[0], y[1], y[2]);
	}
	
	fclose(pFile);
	gsl_odeiv2_driver_free (d);
	return last;
}




















