/*
 *  maths_base.h
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 23/06/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifndef __MATHS_BASE_H
#define __MATHS_BASE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if 0
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif
 
#ifdef __PLANCK__
#include "HL2_likely/tools/errorlist.h"
#include "HL2_likely/tools/io.h"
#else
#include "errorlist.h"
#include "io.h"
#endif


#define math_base                -300
#define math_negative            -1 + math_base
#define math_singularValue       -2 + math_base
#define math_tooManySteps        -3 + math_base
#define math_underflow	         -4 + math_base
#define math_infnan	         -5 + math_base
#define math_wrongValue	         -6 + math_base
#define math_alloc	         -7 + math_base
#define math_interpoloutofrange	 -8 + math_base
#define math_interpol2small      -9 + math_base
#define math_interpol2big       -10 + math_base
#define math_stackTooSmall      -11 + math_base
#define math_overflow	        -12 + math_base
#define math_unknown            -13 + math_base


typedef double my_complex[2];

/* Mathematical constants */
#define pi     3.14159265358979323846
#define pi_sqr 9.86960440108935861883
#define twopi  6.28318530717958647693
#define ln2    0.69314718
#define ln2pi  1.837877066409
#define arcmin 2.90888208665721580e-4
#define arcsec 4.84813681e-6

/* Confidence levels, = erf({1,2,3}/sqrt(2)) */
#define conf_68 0.6827
#define conf_90 0.9000
#define conf_95 0.9545
#define conf_99 0.9973

/* Small numbers */
#define EPSILON  1.0e-5
#define EPSILON1 1.0e-8
#define EPSILON2 1.0e-18

/* Square of the absolute value of x (=|x|^2) */
#define ABSSQR(x) ((x)[0]*(x)[0] + (x)[1]*(x)[1])

/* Scalar product (x,y) */
#define SP(x,y) ((x)[0]*(y)[0]+(x)[1]*(y)[1])


#endif

