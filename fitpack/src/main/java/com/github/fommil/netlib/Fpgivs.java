/*
 *  Produced by f2java.  f2java is part of the Fortran-
 *  -to-Java project at the University of Tennessee Netlib
 *  numerical software repository.
 *
 *  FITPACK source is found in SciPy
 *  https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 *
 *  Fortran input file: fpgivs.f
 *  f2java version: 0.8.1
 *
 */

package com.github.fommil.netlib;
import java.lang.*;
import org.netlib.util.*;



public class Fpgivs {

// c  subroutine fpgivs calculates the parameters of a givens
// c  transformation .
// c  ..
// c  ..scalar arguments..
// c  ..local scalars..
// c  ..function references..
// c  ..

public static void fpgivs (double piv,
doubleW ww,
doubleW cos,
doubleW sin)  {

double dd= 0.0d;
double one= 0.0d;
double store= 0.0d;
one = (double)(0.1e+01f);
store = Math.abs(piv);
if ((store >= ww.val)) {
    dd = (store*Math.sqrt((one+( Math.pow(((ww.val/piv)), 2)))));
}
    if ((store < ww.val)) {
    dd = (ww.val*Math.sqrt((one+( Math.pow(((piv/ww.val)), 2)))));
}
    cos.val = (ww.val/dd);
sin.val = (piv/dd);
ww.val = dd;
/****
Dummy.go_to("com/github/fommil/netlib/Fpgivs",999999);
Dummy.label("com/github/fommil/netlib/Fpgivs",999999);
return;
****/
   }
} // End class.
