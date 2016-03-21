/*
 *  Produced by f2java.  f2java is part of the Fortran-
 *  -to-Java project at the University of Tennessee Netlib
 *  numerical software repository.
 *
 *  FITPACK source is found in SciPy
 *  https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 *
 *  Fortran input file: fpbspl.f
 *  f2java version: 0.8.1
 *
 */

package com.github.fommil.netlib;
import java.lang.*;
import org.netlib.util.*;



public class Fpbspl {

// c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
// c  degree k at t(l) <= x < t(l+1) using the stable recurrence
// c  relation of de boor and cox.
// c  Travis Oliphant  2007
// c    changed so that weighting of 0 is used when knots with
// c      multiplicity are present.
// c    Also, notice that l+k <= n and 1 <= l+1-k
// c      or else the routine will be accessing memory outside t
// c      Thus it is imperative that that k <= l <= n-k but this
// c      is not checked.
// c  ..
// c  ..scalar arguments..
// c  ..array arguments..
// c  ..local scalars..
// c  ..local arrays..
// c  ..

public static void fpbspl (double [] t, int _t_offset,
int n,
int k,
double x,
int l,
double [] h, int _h_offset)  {

double f= 0.0d;
double one= 0.0d;
int i= 0;
int j= 0;
int li= 0;
int lj= 0;
double [] hh= new double[(19)];
one = 0.1e+01;
h[(1-(1))+ _h_offset] = one;
{
for (j = 1; j <= k; j++) {
{
for (i = 1; i <= j; i++) {
hh[(i-(1))] = h[(i-(1))+ _h_offset];
////Dummy.label("com/github/fommil/netlib/Fpbspl",10);
}              //  Close for() loop. 
}
h[(1-(1))+ _h_offset] = 0.0e0;
{
for (i = 1; i <= j; i++) {
li = (l+i);
lj = (li-j);
if ((t[(li-(1))+ _t_offset] == t[(lj-(1))+ _t_offset])) {
	////    Dummy.go_to("com/github/fommil/netlib/Fpbspl",15);
	//// }
    h[((i+1)-(1))+ _h_offset] = 0.0e0;
	//// Dummy.go_to("com/github/fommil/netlib/Fpbspl",20);
	//// label15:
} else {
	//// Dummy.label("com/github/fommil/netlib/Fpbspl",15);
f = (hh[(i-(1))]/((t[(li-(1))+ _t_offset]-t[(lj-(1))+ _t_offset])));
h[(i-(1))+ _h_offset] = (h[(i-(1))+ _h_offset]+(f*((t[(li-(1))+ _t_offset]-x))));
h[((i+1)-(1))+ _h_offset] = (f*((x-t[(lj-(1))+ _t_offset])));
} //// Dummy.label("com/github/fommil/netlib/Fpbspl",20);
}              //  Close for() loop. 
}
//// Dummy.label("com/github/fommil/netlib/Fpbspl",20);
}              //  Close for() loop. 
}
//// Dummy.go_to("com/github/fommil/netlib/Fpbspl",999999);
//// Dummy.label("com/github/fommil/netlib/Fpbspl",999999);
return;
   }
} // End class.
