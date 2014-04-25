/*
  nrsag.h

  Custom version of nr.h for use with SAG

  In this adaptation of the Numerical Recipes code, many programs
  have been converted to double precision.
*/

double dqromb(double (*func)(double), double a, double b);
void dpolint(double xa[], double ya[], int n, double x, double *y, 
	     double *dy);
double dtrapzd(double (*func)(double), double a, double b, int n);

float expint(int n, float x);
