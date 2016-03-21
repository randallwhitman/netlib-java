/*
 *  Produced by f2java.  f2java is part of the Fortran-
 *  -to-Java project at the University of Tennessee Netlib
 *  numerical software repository.
 *
 *  FITPACK source is found in SciPy
 *  https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 *
 *  Fortran input file: fpback.f
 *  f2java version: 0.8.1
 *
 */

package com.github.fommil.netlib;
import java.lang.*;
import org.netlib.util.*;



public class Fpback {

// c  subroutine fpback calculates the solution of the system of
// c  equations a*c = z with a a n x n upper triangular matrix
// c  of bandwidth k.
// c  ..
// c  ..scalar arguments..
// c  ..array arguments..
// c  ..local scalars..
// c  ..

public static void fpback (double [] a, int _a_offset,
double [] z, int _z_offset,
int n,
int k,
double [] c, int _c_offset,
int nest)  {

double store= 0.0d;
int i= 0;
int i1= 0;
int j= 0;
int k1= 0;
int l= 0;
int m= 0;
k1 = (k-1);
c[(n-(1))+ _c_offset] = (z[(n-(1))+ _z_offset]/a[(n-(1))+(1-(1)) * (nest)+ _a_offset]);
i = (n-1);
if ((i == 0)) {
    return;  // Dummy.go_to("com/github/fommil/netlib/Fpback",30);
}
    {
for (j = 2; j <= n; j++) {
store = z[(i-(1))+ _z_offset];
i1 = k1;
if ((j <= k1)) {
    i1 = (j-1);
}
    m = i;
{
for (l = 1; l <= i1; l++) {
m = (m+1);
store = (store-(c[(m-(1))+ _c_offset]*a[(i-(1))+((l+1)-(1)) * (nest)+ _a_offset]));
//// Dummy.label("com/github/fommil/netlib/Fpback",10);
}              //  Close for() loop. 
}
c[(i-(1))+ _c_offset] = (store/a[(i-(1))+(1-(1)) * (nest)+ _a_offset]);
i = (i-1);
//// Dummy.label("com/github/fommil/netlib/Fpback",20);
}              //  Close for() loop. 
}
	/****
label30:
   Dummy.label("com/github/fommil/netlib/Fpback",30);
Dummy.go_to("com/github/fommil/netlib/Fpback",999999);
Dummy.label("com/github/fommil/netlib/Fpback",999999);
return;
	****/
   }
} // End class.
