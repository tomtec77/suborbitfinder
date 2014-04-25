#include <stdio.h>
#include <math.h>

#include "nrsag.h"
#include "nrutil.h"

#define EPS 1.0e-12
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

double dqromb(double (*func)(double), double a, double b)
{
  double ss,dss;
  double s[JMAXP],h[JMAXP+1];
  int j;

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=dtrapzd(func,a,b,j);
    if (j >= K) {
      dpolint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss))
	{
	  //printf("fabs(dss) <= EPS*fabs(ss):%g %g\n",                            fabs(dss),EPS*fabs(ss));
	  
	  return ss;
	}
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  nrerror("Too many steps in routine qromb");
  return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
