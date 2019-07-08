/************************************************/
/*                                              */
/*  CMATH.  Copyright (c) 1989 Design Software  */
/*                                              */
/************************************************/


/* This C code written by ...  Peter & Nigel,
   ----------------------      Design Software,
                               42 Gubberley St,
                               Kenmore, 4069,
                               Australia.         */

/*-----------------------------------------------------------------*/



#ifndef INC_2LAB_CDECOMP_H
#define INC_2LAB_CDECOMP_H



#define  PROTOTYPE  1

#define  STDLIBH  1

#define  STRINGH  1


/* Machine epsilon :
   ---------------
   This is the smallest real value that your machine can distinguish.
   You may estimate it using the following piece of code.

   one = 1.0;
   two = 2.0;
   epsilon = one;
   do { epsilon /= two; } while (one+epsilon > one);
   epsilon *= two;

   We have selected a value that should allow the CMATH routines to
   work on a large range of machines but may not make the best use of
   a machine offering higher precision arithmetic.  Note that the
   values we supply will allow CMATH to work on ALL of the machine-
   compiler combinations mentioned above.

   */

#define  EPSILON  2.2e-16

/* Machine Overflow :
   ----------------
   The largest real number that can be represented. */

#define  OVRFLOW  1.0e+75

/* Machine Underflow :
   -----------------
   The smallest positive number that can be represented. */

#define  UNDRFLOW  1.0e-75

/*-----------------------------------------------------------------*/


#define  DECOMP_C    101


int decomp (int n, int ndim,
            double *a, double *cond,
            int pivot[], int *flag);

int solve (int n, int ndim,
           double *a, double b[],
           int pivot[]);




#endif //INC_2LAB_CDECOMP_H