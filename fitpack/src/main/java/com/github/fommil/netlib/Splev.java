/*
 *  Produced by f2java.  f2java is part of the Fortran-
 *  -to-Java project at the University of Tennessee Netlib
 *  numerical software repository.
 *
 *  FITPACK source is found in SciPy
 *  https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 *
 *  Fortran input file: splev.f
 *  f2java version: 0.8.1
 *
 */

package com.github.fommil.netlib;
import java.lang.*;
import org.netlib.util.*;



public class Splev {

// c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
// c  a spline s(x) of degree k, given in its b-spline representation.
// c
// c  calling sequence:
// c     call splev(t,n,c,k,x,y,m,e,ier)
// c
// c  input parameters:
// c    t    : array,length n, which contains the position of the knots.
// c    n    : integer, giving the total number of knots of s(x).
// c    c    : array,length n, which contains the b-spline coefficients.
// c    k    : integer, giving the degree of s(x).
// c    x    : array,length m, which contains the points where s(x) must
// c           be evaluated.
// c    m    : integer, giving the number of points where s(x) must be
// c           evaluated.
// c    e    : integer, if 0 the spline is extrapolated from the end
// c           spans for points not in the support, if 1 the spline
// c           evaluates to zero for those points, if 2 ier is set to
// c           1 and the subroutine returns, and if 3 the spline evaluates

// c           to the value of the nearest boundary point.
// c
// c  output parameter:
// c    y    : array,length m, giving the value of s(x) at the different
// c           points.
// c    ier  : error flag
// c      ier = 0 : normal return
// c      ier = 1 : argument out of bounds and e == 2
// c      ier =10 : invalid input data (see restrictions)
// c
// c  restrictions:
// c    m >= 1
// c--    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
// c
// c  other subroutines required: fpbspl.
// c
// c  references :
// c    de boor c  : on calculating with b-splines, j. approximation theory
// c                 6 (1972) 50-62.
// c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths

// c                 applics 10 (1972) 134-149.
// c    dierckx p. : curve and surface fitting with splines, monographs on

// c                 numerical analysis, oxford university press, 1993.
// c
// c  author :
// c    p.dierckx
// c    dept. computer science, k.u.leuven
// c    celestijnenlaan 200a, b-3001 heverlee, belgium.
// c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
// c
// c  latest update : march 1987
// c
// c++ pearu: 11 aug 2003
// c++   - disabled cliping x values to interval [min(t),max(t)]
// c++   - removed the restriction of the orderness of x values
// c++   - fixed initialization of sp to double precision value
// c
// c  ..scalar arguments..
// c  ..array arguments..
// c  ..local scalars..
// c++..
// c..++
// c  ..local array..
// c  ..
// c  before starting computations a data check is made. if the input data

// c  are invalid control is immediately repassed to the calling program.

public static void splev (double [] t, int _t_offset,
int n,
double [] c, int _c_offset,
int k,
double [] x, int _x_offset,
double [] y, int _y_offset,
int m,
int e,
intW ier)  {

int i= 0;
int j= 0;
int k1= 0;
int l= 0;
int ll= 0;
int l1= 0;
int nk1= 0;
int k2= 0;
double arg= 0.0d;
double sp= 0.0d;
double tb= 0.0d;
double te= 0.0d;
double [] h= new double[(20)];
ier.val = 10;
// c--      if(m-1) 100,30,10
// c++..
if ((m < 1)) {
    return;  //     Dummy.go_to("com/github/fommil/netlib/Splev",100);
}
    // c..++
// c--  10  do 20 i=2,m
// c--        if(x(i).lt.x(i-1)) go to 100
// c--  20  continue
ier.val = 0;
// c  fetch tb and te, the boundaries of the approximation interval.
k1 = (k+1);
// c++..
k2 = (k1+1);
// c..++
nk1 = (n-k1);
tb = t[(k1-(1))+ _t_offset];
te = t[((nk1+1)-(1))+ _t_offset];
l = k1;
l1 = (l+1);
// c  main loop for the different points.
{
for (i = 1; i <= m; i++) {
// c  fetch a new x-value arg.
arg = x[(i-(1))+ _x_offset];
// c  check if arg is in the support
if (((arg < tb) || (arg > te)))  {
    if ((e == 0))  {
		// Dummy.go_to("com/github/fommil/netlib/Splev",35);
	if (true) throw new RuntimeException("not-supported");
}
else if ((e == 1))  {
    y[(i-(1))+ _y_offset] = (double)(0);
	// Dummy.go_to("com/github/fommil/netlib/Splev",80);
	if (true) throw new RuntimeException("not-supported");
}              // Close else if()
else if ((e == 2))  {
    ier.val = 1;
	return;  // Dummy.go_to("com/github/fommil/netlib/Splev",100);
}              // Close else if()
else if ((e == 3))  {
    if ((arg < tb))  {
    arg = tb;
}
else  {
  arg = te;
}              //  Close else.
}              // Close else if()
}
// c  search for knot interval t(l) <= arg < t(l+1)
// c++..
/****
label35:
   Dummy.label("com/github/fommil/netlib/Splev",35);
if (((arg >= t[(l-(1))+ _t_offset]) || (l1 == k2))) {
    Dummy.go_to("com/github/fommil/netlib/Splev",40);
}
    l1 = l;
l = (l-1);
Dummy.go_to("com/github/fommil/netlib/Splev",35);
// c..++
label40:
   Dummy.label("com/github/fommil/netlib/Splev",40);
****/
while (((arg >= t[(l1-(1))+ _t_offset]) && (l != nk1))) {
    //// Dummy.go_to("com/github/fommil/netlib/Splev",50);
	//// }
    l = l1;
l1 = (l+1);
}  //// Dummy.go_to("com/github/fommil/netlib/Splev",40);
// c  evaluate the non-zero b-splines at arg.
////label50:
////   Dummy.label("com/github/fommil/netlib/Splev",50);
com.github.fommil.netlib.Fpbspl.fpbspl(t,_t_offset,n,k,arg,l,h,0);
// c  find the value of s(x) at x=arg.
sp = 0.0e0;
ll = (l-k1);
{
for (j = 1; j <= k1; j++) {
ll = (ll+1);
sp = (sp+(c[(ll-(1))+ _c_offset]*h[(j-(1))]));
////Dummy.label("com/github/fommil/netlib/Splev",60);
}              //  Close for() loop. 
}
y[(i-(1))+ _y_offset] = sp;
////Dummy.label("com/github/fommil/netlib/Splev",80);
}              //  Close for() loop. 
}
/****
label100:
   Dummy.label("com/github/fommil/netlib/Splev",100);
Dummy.go_to("com/github/fommil/netlib/Splev",999999);
Dummy.label("com/github/fommil/netlib/Splev",999999);
return;
****/
   }
} // End class.
