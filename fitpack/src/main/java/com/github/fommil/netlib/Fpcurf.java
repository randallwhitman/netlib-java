/*
 *  Produced by f2java.  f2java is part of the Fortran-
 *  -to-Java project at the University of Tennessee Netlib
 *  numerical software repository.
 *
 *  FITPACK source is found in SciPy
 *  https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 *
 *  Fortran input file: fpcurf.f
 *  f2java version: 0.8.1
 *
 */

package com.github.fommil.netlib;
import java.lang.*;
import org.netlib.util.*;



public class Fpcurf {

// c  ..
// c  ..scalar arguments..
// c  ..array arguments..
// c  ..local scalars..
// c  ..local arrays..
// c  ..function references
// c  ..subroutine references..
// c    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
// c  ..
// c  set constants

public static void fpcurf (int iopt,
double [] x, int _x_offset,
double [] y, int _y_offset,
double [] w, int _w_offset,
int m,
double xb,
double xe,
int k,
double s,
int nest,
double tol,
int maxit,
int k1,
int k2,
intW n,
double [] t, int _t_offset,
double [] c, int _c_offset,
doubleW fp,
double [] fpint, int _fpint_offset,
double [] z, int _z_offset,
double [] a, int _a_offset,
double [] b, int _b_offset,
double [] g, int _g_offset,
double [] q, int _q_offset,
int [] nrdata, int _nrdata_offset,
intW ier)  {

double acc= 0.0d;
double con1= 0.0d;
double con4= 0.0d;
double con9= 0.0d;
double cos= 0.0d;
double half= 0.0d;
double fpart= 0.0d;
double fpms= 0.0d;
double fpold= 0.0d;
double fp0= 0.0d;
double f1= 0.0d;
double f2= 0.0d;
double f3= 0.0d;
double one= 0.0d;
double p= 0.0d;
double pinv= 0.0d;
double piv= 0.0d;
double p1= 0.0d;
double p2= 0.0d;
double p3= 0.0d;
double rn= 0.0d;
double sin= 0.0d;
double store= 0.0d;
double term= 0.0d;
double wi= 0.0d;
double xi= 0.0d;
double yi= 0.0d;
int i= 0;
int ich1= 0;
int ich3= 0;
int it= 0;
int iter= 0;
int i1= 0;
int i2= 0;
int i3= 0;
int j= 0;
int k3= 0;
int l= 0;
int l0= 0;
int mk1= 0;
int New= 0;
int nk1= 0;
int nmax= 0;
int nmin= 0;
int nplus= 0;
int npl1= 0;
int nrint= 0;
int n8= 0;
double [] h= new double[(7)];
double fprati= 0.0d;
one = 0.1e+01;
con1 = 0.1e0;
con9 = 0.9e0;
con4 = 0.4e-01;
half = 0.5e0;
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// c  part 1: determination of the number of knots and their position     c
// c  **************************************************************      c
// c  given a set of knots we compute the least-squares spline sinf(x),   c
// c  and the corresponding sum of squared residuals fp=f(p=inf).         c
// c  if iopt=-1 sinf(x) is the requested approximation.                  c
// c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
// c    if fp <=s we will continue with the current set of knots.         c
// c    if fp > s we will increase the number of knots and compute the    c
// c       corresponding least-squares spline until finally fp<=s.        c
// c    the initial choice of knots depends on the value of s and iopt.   c
// c    if s=0 we have spline interpolation; in that case the number of   c
// c    knots equals nmax = m+k+1.                                        c
// c    if s > 0 and                                                      c
// c      iopt=0 we first compute the least-squares polynomial of         c
// c      degree k; n = nmin = 2*k+2                                      c
// c      iopt=1 we start with the set of knots found at the last         c
// c      call of the routine, except for the case that s > fp0; then     c
// c      we compute directly the least-squares polynomial of degree k.   c
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// c  determine nmin, the number of knots for polynomial approximation.
nmin = (2*k1);
if ((iopt >= 0)) {
    //// Dummy.go_to("com/github/fommil/netlib/Fpcurf",60);
	//// }
    // c  calculation of acc, the absolute tolerance for the root of f(p)=s.
acc = (tol*s);
// c  determine nmax, the number of knots for spline interpolation.
nmax = (m+k1);
if ((s <= 0.0e0)) {
    //// Dummy.go_to("com/github/fommil/netlib/Fpcurf",45);
	//// }
    // c  if s=0, s(x) is an interpolating spline.
// c  test whether the required storage space exceeds the available one.
n.val = nmax;
if ((nmax > nest)) {
ier.val = 1;
return;  //// Dummy.go_to("com/github/fommil/netlib/Fpcurf",420);
}
    // c  find the position of the interior knots in case of interpolation.
//// label10:
////    Dummy.label("com/github/fommil/netlib/Fpcurf",10);
mk1 = (m-k1);
if ((mk1 != 0)) {
    //// Dummy.go_to("com/github/fommil/netlib/Fpcurf",60);
	//// }
    k3 = (k/2);
i = k2;
j = (k3+2);
if (((k3*2) != k)) {  // odd
    //// Dummy.go_to("com/github/fommil/netlib/Fpcurf",30);
	//// }
    {
for (l = 1; l <= mk1; l++) {
t[(i-(1))+ _t_offset] = x[(j-(1))+ _x_offset];
i = (i+1);
j = (j+1);
////Dummy.label("com/github/fommil/netlib/Fpcurf",20);
}              //  Close for() loop. 
}
	////Dummy.go_to("com/github/fommil/netlib/Fpcurf",60);
} else {  // even  //// label30:
	////   Dummy.label("com/github/fommil/netlib/Fpcurf",30);
{
for (l = 1; l <= mk1; l++) {
t[(i-(1))+ _t_offset] = (((x[(j-(1))+ _x_offset]+x[((j-1)-(1))+ _x_offset]))*half);
i = (i+1);
j = (j+1);
////Dummy.label("com/github/fommil/netlib/Fpcurf",40);
}              //  Close for() loop. 
}
} //// Dummy.go_to("com/github/fommil/netlib/Fpcurf",60);
// c  if s>0 our initial choice of knots depends on the value of iopt.
// c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares

// c  polynomial of degree k which is a spline without interior knots.
// c  if iopt=1 and fp0>s we start computing the least squares spline
// c  according to the set of knots found at the last call of the routine.
} else {  // s>0  //// label45:   Dummy.label("com/github/fommil/netlib/Fpcurf",45);
/****if ((iopt == 0)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",50);
}
    if ((n.val == nmin)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",50);
	} ****/
if ((iopt != 0) && (n.val != nmin)) {
    fp0 = fpint[(n.val-(1))+ _fpint_offset];
fpold = fpint[((n.val-1)-(1))+ _fpint_offset];
nplus = nrdata[(n.val-(1))+ _nrdata_offset];
}
/****if ((fp0 > s)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",60);
}
    label50:
	Dummy.label("com/github/fommil/netlib/Fpcurf",50);****/
if (fp0 <= s) {
n.val = nmin;
fpold = 0.0e0;
nplus = 0;
nrdata[(1-(1))+ _nrdata_offset] = (m-2);
}
}}  // 60
// c  main loop for the different sets of knots. m is a save upper bound
// c  for the number of trials.
//// label60:
////   Dummy.label("com/github/fommil/netlib/Fpcurf",60);
{
for (iter = 1; iter <= m; iter++) {
if ((n.val == nmin)) {
    ier.val = -2;
}
    // c  find nrint, tne number of knot intervals.
nrint = ((n.val-nmin)+1);
// c  find the position of the additional knots which are needed for
// c  the b-spline representation of s(x).
nk1 = (n.val-k1);
i = n.val;
{
for (j = 1; j <= k1; j++) {
t[(j-(1))+ _t_offset] = xb;
t[(i-(1))+ _t_offset] = xe;
i = (i-1);
////Dummy.label("com/github/fommil/netlib/Fpcurf",70);
}              //  Close for() loop. 
}
// c  compute the b-spline coefficients of the least-squares spline
// c  sinf(x). the observation matrix a is built up row by row and
// c  reduced to upper triangular form by givens transformations.
// c  at the same time fp=f(p=inf) is computed.
fp.val = 0.0e0;
// c  initialize the observation matrix a.
{
for (i = 1; i <= nk1; i++) {
z[(i-(1))+ _z_offset] = 0.0e0;
{
for (j = 1; j <= k1; j++) {
a[(i-(1))+(j-(1)) * (nest)+ _a_offset] = 0.0e0;
//// Dummy.label("com/github/fommil/netlib/Fpcurf",80);
}              //  Close for() loop. 
}
//// Dummy.label("com/github/fommil/netlib/Fpcurf",80);
}              //  Close for() loop. 
}
l = k1;
{
for (it = 1; it <= m; it++) {
// c  fetch the current data point x(it),y(it).
xi = x[(it-(1))+ _x_offset];
wi = w[(it-(1))+ _w_offset];
yi = (y[(it-(1))+ _y_offset]*wi);
// c  search for knot interval t(l) <= xi < t(l+1).
/**** label85:
   Dummy.label("com/github/fommil/netlib/Fpcurf",85);
if (((xi < t[((l+1)-(1))+ _t_offset]) || (l == nk1))) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",90);
}
    l = (l+1);
	Dummy.go_to("com/github/fommil/netlib/Fpcurf",85); ****/
// //boolean flag = true;
// // while (flag) {
while (! (((xi < t[((l+1)-(1))+ _t_offset]) || (l == nk1))) ) {
	// //    flag = false;
	// // break;
    l = (l+1);
	}
// c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
//// label90:
////    Dummy.label("com/github/fommil/netlib/Fpcurf",90);
com.github.fommil.netlib.Fpbspl.fpbspl(t,_t_offset,n.val,k,xi,l,h,0);
{
for (i = 1; i <= k1; i++) {
q[(it-(1))+(i-(1)) * (m)+ _q_offset] = h[(i-(1))];
h[(i-(1))] = (h[(i-(1))]*wi);
//// Dummy.label("com/github/fommil/netlib/Fpcurf",95);
}              //  Close for() loop. 
}
// c  rotate the new row of the observation matrix into triangle.
j = (l-k1);
{
for (i = 1; i <= k1; i++) {
j = (j+1);
piv = h[(i-(1))];
////if ((piv == 0.0e0)) {
////    Dummy.go_to("com/github/fommil/netlib/Fpcurf",110);
////}
if (piv != 0.0) {
    // c  calculate the parameters of the givens transformation.
//// com.github.fommil.netlib.Fpgivs.fpgivs(piv,a,(j-(1))+(1-(1)) * (nest)+ _a_offset,cos,sin);
doubleW wrapCos = new doubleW(cos);
doubleW wrapSin = new doubleW(sin);
doubleW wrapWw = new doubleW(a[(j-(1))+(1-(1)) * (nest)+ _a_offset]);
com.github.fommil.netlib.Fpgivs.fpgivs(piv,wrapWw,wrapCos,wrapSin);
cos = wrapCos.val;
sin = wrapSin.val;
a[(j-(1))+(1-(1)) * (nest)+ _a_offset] = wrapWw.val;
// c  transformations to right hand side.
//// com.github.fommil.netlib.Fprota.fprota(cos,sin,yi,z,(j-(1))+ _z_offset);
doubleW wrapA = new doubleW(yi);
doubleW wrapB = new doubleW(z[(j-(1))+ _z_offset]);
com.github.fommil.netlib.Fprota.fprota(cos,sin,wrapA,wrapB);
yi = wrapA.val;
z[(j-(1))+ _z_offset] = wrapB.val;
if ((i != k1)) {
	////    Dummy.go_to("com/github/fommil/netlib/Fpcurf",120);
	////}
    i2 = 1;
i3 = (i+1);
{
for (i1 = i3; i1 <= k1; i1++) {
i2 = (i2+1);
// c  transformations to left hand side.
//// com.github.fommil.netlib.Fprota.fprota(cos,sin,h,(i1-(1)),a,(j-(1))+(i2-(1)) * (nest)+ _a_offset);
wrapA.val = h[(i1-(1))];
wrapB.val = a[(j-(1))+(i2-(1)) * (nest)+ _a_offset];
com.github.fommil.netlib.Fprota.fprota(cos,sin,wrapA,wrapB);
h[(i1-(1))] = wrapA.val;
a[(j-(1))+(i2-(1)) * (nest)+ _a_offset] = wrapB.val;
////Dummy.label("com/github/fommil/netlib/Fpcurf",100);
}              //  Close for() loop. 
}
}////Dummy.label("com/github/fommil/netlib/Fpcurf",110);
}}              //  Close for() loop. 
}
// c  add contribution of this row to the sum of squares of residual
// c  right hand sides.
////label120:
//// Dummy.label("com/github/fommil/netlib/Fpcurf",120);
fp.val = (fp.val+(yi*yi));
////Dummy.label("com/github/fommil/netlib/Fpcurf",130);
}              //  Close for() loop. 
}
if ((ier.val == (-2))) {
    fp0 = fp.val;
}
    fpint[(n.val-(1))+ _fpint_offset] = fp0;
fpint[((n.val-1)-(1))+ _fpint_offset] = fpold;
nrdata[(n.val-(1))+ _nrdata_offset] = nplus;
// c  backward substitution to obtain the b-spline coefficients.
com.github.fommil.netlib.Fpback.fpback(a,_a_offset,z,_z_offset,nk1,k1,c,_c_offset,nest);
// c  test whether the approximation sinf(x) is an acceptable solution.
if ((iopt < 0)) {
    return; ////Dummy.go_to("com/github/fommil/netlib/Fpcurf",440);
}
    fpms = (fp.val-s);
if ((Math.abs(fpms) < acc)) {
	return;    ////Dummy.go_to("com/github/fommil/netlib/Fpcurf",440);
}
    // c  if f(p=inf) < s accept the choice of knots.
if ((fpms < 0.0e0)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",250);
}
    // c  if n = nmax, sinf(x) is an interpolating spline.
if ((n.val == nmax)) {
	ier.val = -1;
	return;  ////  Dummy.go_to("com/github/fommil/netlib/Fpcurf",430);
}
if(true) throw new RuntimeException("not-expected");
    // c  increase the number of knots.
// c  if n=nest we cannot increase the number of knots because of
// c  the storage capacity limitation.
if ((n.val == nest)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",420);
}
    // c  determine the number of knots nplus we are going to add.
if ((ier.val == 0)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",140);
}
    nplus = 1;
ier.val = 0;
Dummy.go_to("com/github/fommil/netlib/Fpcurf",150);
label140:
   Dummy.label("com/github/fommil/netlib/Fpcurf",140);
npl1 = (nplus*2);
rn = (double)(nplus);
if (((fpold-fp.val) > acc)) {
    npl1 = (int)(((rn*fpms)/((fpold-fp.val))));
}
    nplus = Math.min((nplus*2), Util.max(npl1, (nplus/2), 1));
label150:
   Dummy.label("com/github/fommil/netlib/Fpcurf",150);
fpold = fp.val;
// c  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
// c  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
fpart = 0.0e0;
i = 1;
l = k2;
New = 0;
{
for (it = 1; it <= m; it++) {
if (((x[(it-(1))+ _x_offset] < t[(l-(1))+ _t_offset]) || (l > nk1))) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",160);
}
    New = 1;
l = (l+1);
label160:
   Dummy.label("com/github/fommil/netlib/Fpcurf",160);
term = 0.0e0;
l0 = (l-k2);
{
for (j = 1; j <= k1; j++) {
l0 = (l0+1);
term = (term+(c[(l0-(1))+ _c_offset]*q[(it-(1))+(j-(1)) * (m)+ _q_offset]));
Dummy.label("com/github/fommil/netlib/Fpcurf",170);
}              //  Close for() loop. 
}
term = ( Math.pow(((w[(it-(1))+ _w_offset]*((term-y[(it-(1))+ _y_offset])))), 2));
fpart = (fpart+term);
if ((New == 0)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",180);
}
    store = (term*half);
fpint[(i-(1))+ _fpint_offset] = (fpart-store);
i = (i+1);
fpart = store;
New = 0;
Dummy.label("com/github/fommil/netlib/Fpcurf",180);
}              //  Close for() loop. 
}
fpint[(nrint-(1))+ _fpint_offset] = fpart;
{
for (l = 1; l <= nplus; l++) {
// c  add a new knot.
//// com.github.fommil.netlib.Fpknot.fpknot(x,_x_offset,m,t,_t_offset,n.val,fpint,_fpint_offset,nrdata,_nrdata_offset,nrint,nest,1);
	throw new RuntimeException("not-supported");
// c  if n=nmax we locate the knots as for interpolation.
//// if ((n.val == nmax)) {
////     Dummy.go_to("com/github/fommil/netlib/Fpcurf",10);
//// }
    // c  test whether we cannot further increase the number of knots.
//// if ((n.val == nest)) {
////     Dummy.go_to("com/github/fommil/netlib/Fpcurf",200);
//// }
////     Dummy.label("com/github/fommil/netlib/Fpcurf",190);
}              //  Close for() loop. 
}
// c  restart the computations with the new set of knots.
Dummy.label("com/github/fommil/netlib/Fpcurf",200);
}              //  Close for() loop. 
}
// c  test whether the least-squares kth degree polynomial is a solution
// c  of our approximation problem.
label250:
   Dummy.label("com/github/fommil/netlib/Fpcurf",250);
if ((ier.val == (-2))) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",440);
}
    // cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// c  part 2: determination of the smoothing spline sp(x).                c
// c  ***************************************************                 c
// c  we have determined the number of knots and their position.          c
// c  we now compute the b-spline coefficients of the smoothing spline    c
// c  sp(x). the observation matrix a is extended by the rows of matrix   c
// c  b expressing that the kth derivative discontinuities of sp(x) at    c
// c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
// c  ponding weights of these additional rows are set to 1/p.            c
// c  iteratively we then have to determine the value of p such that      c
// c  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
// c  the least-squares kth degree polynomial corresponds to p=0, and     c
// c  that the least-squares spline corresponds to p=infinity. the        c
// c  iteration process which is proposed here, makes use of rational     c
// c  interpolation. since f(p) is a convex and strictly decreasing       c
// c  function of p, it can be approximated by a rational function        c
// c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
// c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
// c  to calculate the new value of p such that r(p)=s. convergence is    c
// c  guaranteed by taking f1>0 and f3<0.                                 c
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// c  evaluate the discontinuity jump of the kth derivative of the
// c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
//// com.github.fommil.netlib.Fpdisc.fpdisc(t,_t_offset,n.val,k2,b,_b_offset,nest);
	throw new RuntimeException("not-supported");
// c  initial value for p.
/**** p1 = 0.0e0;
f1 = (fp0-s);
p3 = (-(one));
f3 = fpms;
p = (double)(0.f);
{
for (i = 1; i <= nk1; i++) {
p = (p+a[(i-(1))+(1-(1)) * (nest)+ _a_offset]);
Dummy.label("com/github/fommil/netlib/Fpcurf",255);
}              //  Close for() loop.
}
rn = (double)(nk1);
p = (rn/p);
ich1 = 0;
ich3 = 0;
n8 = (n.val-nmin);
// c  iteration process to find the root of f(p) = s.
{
for (iter = 1; iter <= maxit; iter++) {
// c  the rows of matrix b with weight 1/p are rotated into the
// c  triangularised observation matrix a which is stored in g.
pinv = (one/p);
{
for (i = 1; i <= nk1; i++) {
c[(i-(1))+ _c_offset] = z[(i-(1))+ _z_offset];
g[(i-(1))+(k2-(1)) * (nest)+ _g_offset] = 0.0e0;
{
for (j = 1; j <= k1; j++) {
g[(i-(1))+(j-(1)) * (nest)+ _g_offset] = a[(i-(1))+(j-(1)) * (nest)+ _a_offset];
Dummy.label("com/github/fommil/netlib/Fpcurf",260);
}              //  Close for() loop. 
}
Dummy.label("com/github/fommil/netlib/Fpcurf",260);
}              //  Close for() loop. 
}
{
for (it = 1; it <= n8; it++) {
// c  the row of matrix b is rotated into triangle by givens transformation
{
for (i = 1; i <= k2; i++) {
h[(i-(1))] = (b[(it-(1))+(i-(1)) * (nest)+ _b_offset]*pinv);
Dummy.label("com/github/fommil/netlib/Fpcurf",270);
}              //  Close for() loop. 
}
yi = 0.0e0;
{
for (j = it; j <= nk1; j++) {
piv = h[(1-(1))];
// c  calculate the parameters of the givens transformation.
////com.github.fommil.netlib.Fpgivs.fpgivs(piv,g,(j-(1))+(1-(1)) * (nest)+ _g_offset,cos,sin);
	throw new RuntimeException("not-supported");
// c  transformations to right hand side.
//// com.github.fommil.netlib.Fprota.fprota(cos,sin,yi,c,(j-(1))+ _c_offset);
	throw new RuntimeException("not-supported");
if ((j == nk1)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",300);
}
    i2 = k1;
if ((j > n8)) {
    i2 = (nk1-j);
}
    {
for (i = 1; i <= i2; i++) {
// c  transformations to left hand side.
i1 = (i+1);
////com.github.fommil.netlib.Fprota.fprota(cos,sin,h,(i1-(1)),g,(j-(1))+(i1-(1)) * (nest)+ _g_offset);
	throw new RuntimeException("not-supported");
h[(i-(1))] = h[(i1-(1))];
Dummy.label("com/github/fommil/netlib/Fpcurf",280);
}              //  Close for() loop. 
}
h[((i2+1)-(1))] = 0.0e0;
Dummy.label("com/github/fommil/netlib/Fpcurf",290);
}              //  Close for() loop. 
}
Dummy.label("com/github/fommil/netlib/Fpcurf",300);
}              //  Close for() loop. 
}
// c  backward substitution to obtain the b-spline coefficients.
com.github.fommil.netlib.Fpback.fpback(g,_g_offset,c,_c_offset,nk1,k2,c,_c_offset,nest);
// c  computation of f(p).
fp.val = 0.0e0;
l = k2;
{
for (it = 1; it <= m; it++) {
if (((x[(it-(1))+ _x_offset] < t[(l-(1))+ _t_offset]) || (l > nk1))) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",310);
}
    l = (l+1);
label310:
   Dummy.label("com/github/fommil/netlib/Fpcurf",310);
l0 = (l-k2);
term = 0.0e0;
{
for (j = 1; j <= k1; j++) {
l0 = (l0+1);
term = (term+(c[(l0-(1))+ _c_offset]*q[(it-(1))+(j-(1)) * (m)+ _q_offset]));
Dummy.label("com/github/fommil/netlib/Fpcurf",320);
}              //  Close for() loop. 
}
fp.val = (fp.val+( Math.pow(((w[(it-(1))+ _w_offset]*((term-y[(it-(1))+ _y_offset])))), 2)));
Dummy.label("com/github/fommil/netlib/Fpcurf",330);
}              //  Close for() loop. 
}
// c  test whether the approximation sp(x) is an acceptable solution.
fpms = (fp.val-s);
if ((Math.abs(fpms) < acc)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",440);
}
    // c  test whether the maximal number of iterations is reached.
if ((iter == maxit)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",400);
}
    // c  carry out one more step of the iteration process.
p2 = p;
f2 = fpms;
if ((ich3 != 0)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",340);
}
    if ((((f2-f3)) > acc)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",335);
}
    // c  our initial choice of p is too large.
p3 = p2;
f3 = f2;
p = (p*con4);
if ((p <= p1)) {
    p = ((p1*con9)+(p2*con1));
}
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",360);
label335:
   Dummy.label("com/github/fommil/netlib/Fpcurf",335);
if ((f2 < 0.0e0)) {
    ich3 = 1;
}
    label340:
   Dummy.label("com/github/fommil/netlib/Fpcurf",340);
if ((ich1 != 0)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",350);
}
    if ((((f1-f2)) > acc)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",345);
}
    // c  our initial choice of p is too small
p1 = p2;
f1 = f2;
p = (p/con4);
if ((p3 < 0.f)) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",360);
}
    if ((p >= p3)) {
    p = ((p2*con1)+(p3*con9));
}
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",360);
label345:
   Dummy.label("com/github/fommil/netlib/Fpcurf",345);
if ((f2 > 0.0e0)) {
    ich1 = 1;
}
    // c  test whether the iteration process proceeds as theoretically
// c  expected.
label350:
   Dummy.label("com/github/fommil/netlib/Fpcurf",350);
if (((f2 >= f1) || (f2 <= f3))) {
    Dummy.go_to("com/github/fommil/netlib/Fpcurf",410);
}
    // c  find the new value for p.
//// p = com.github.fommil.netlib.Fprati.fprati(p1,f1,p2,f2,p3,f3);
	throw new RuntimeException("not-supported");
Dummy.label("com/github/fommil/netlib/Fpcurf",360);
}              //  Close for() loop. 
}
// c  error codes and messages.
label400:
   Dummy.label("com/github/fommil/netlib/Fpcurf",400);
ier.val = 3;
Dummy.go_to("com/github/fommil/netlib/Fpcurf",440);
label410:
   Dummy.label("com/github/fommil/netlib/Fpcurf",410);
ier.val = 2;
Dummy.go_to("com/github/fommil/netlib/Fpcurf",440);
label420:
   Dummy.label("com/github/fommil/netlib/Fpcurf",420);
ier.val = 1;
Dummy.go_to("com/github/fommil/netlib/Fpcurf",440);
label430:
   Dummy.label("com/github/fommil/netlib/Fpcurf",430);
ier.val = -1;
label440:
   Dummy.label("com/github/fommil/netlib/Fpcurf",440);
Dummy.go_to("com/github/fommil/netlib/Fpcurf",999999);
Dummy.label("com/github/fommil/netlib/Fpcurf",999999);
return; ****/
   }
} // End class.
