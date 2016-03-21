/*
 *  Produced by f2java.  f2java is part of the Fortran-
 *  -to-Java project at the University of Tennessee Netlib
 *  numerical software repository.
 *
 *  FITPACK source is found in SciPy
 *  https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 *
 *  Fortran input file: fprota.f
 *  f2java version: 0.8.1
 *
 */

package com.github.fommil.netlib;
import java.lang.*;
import org.netlib.util.*;



public class Fprota {

// c  subroutine fprota applies a givens rotation to a and b.
// c  ..
// c  ..scalar arguments..
// c ..local scalars..
// c  ..

public static void fprota (double cos,
double sin,
doubleW a,
doubleW b)  {

double stor1= 0.0d;
double stor2= 0.0d;
stor1 = a.val;
stor2 = b.val;
b.val = ((cos*stor2)+(sin*stor1));
a.val = ((cos*stor1)-(sin*stor2));
/****
Dummy.go_to("com/github/fommil/netlib/Fprota",999999);
Dummy.label("com/github/fommil/netlib/Fprota",999999);
return;
****/
   }
} // End class.
