#ifndef _MATRIX_MULTIPLY_MACROS_SSE_H_
#define _MATRIX_MULTIPLY_MACROS_SSE_H_

/* 
Jernej Barbic
barbic@cs.cmu.edu
Carnegie Mellon University
2006
Version 1.0

Performs 4x4 matrix vs 4-vector multiply, using MMX/SSE:
Computes wVector = AMatrix * vVector, 
   for  4x4 matrix AMatrix, and 4-vectors vVector and wVector,
   using Intel's Multimedia Extensions SSE. 
SSE was introduced in 1999 with Intel Pentium III, and should
be supported in all Intel Pentium processors thereafter.
This code doesn't use any SSE2 instructions or higher.
Other vector/matrix sizes (other than 4) are not supported.
All matrix/vector entries must be single-precision floats (32 bits).

This code is a C macro, written for Microsoft Visual Studio .NET 2003, 
for inclusion into C code at places where you need a fast matrix-vector 
product, such as for example if transforming three-dimensional vectors 
using rigid body transformations, camera view transformations, etc.
It won't compile e.g. under Linux in gcc since gcc uses a 
different syntax to insert inline assembly into C code
(but should be easy to port). The machine code itself, however, 
is just the same.

Matrix "AMatrix" must be 16-byte aligned and stored row-major:
Ordering of AMatrix entries:
  0  1  2  3
  4  5  6  7
  8  9 10 11
 12 13 14 15
Vectors "vVector" and "wVector" must also be 16-byte aligned.
For non-heap variables (i.e. global variables, local variables),
you can align using __declspec(align(16)) (see below). 
For variables dynamically 
allocated on the heap, you need to use some other method
(for example, allocate a buffer using malloc, then make 
a pointer pointing to the first aligned address within the buffer).
If you don't align, your application will crash. 

I don't claim that this code is the fastest possible, and
you are welcome to let me know if you manage to optimize it
even further.

Also, if using SSE3, the last part (summing the results) could
be significantly simplified using the horizontal 
parallel add instruction HADDPS. The Wikipedia pages on SSE
are a good starting point if you are interested in MMX/SSE.

Example usage:

#include "matrixMultiplyMacros_SSE.h"
#include <stdio.h>

int main()
{
  __declspec(align(16)) float vVector[4] = {3.14, 2.718, -1.5, 2.3};
  __declspec(align(16)) float wVector[4]; // will hold the result
  __declspec(align(16)) float AMatrix[16] = {1.2f, 0.8f, 0.3f, -0.2f,
                                             -0.4f, 0.2f, 0.7f, -1.5f,
                                             7.2f, -3.5f, 2.76f, 1.79f,
                                             1.3f, -0.6f, 1.1f, -2.31f};
  MATRIX_VECTOR_MULTIPLY4X4_SSE(AMatrix,vVector,wVector);
  printf("Result is: %f %f %f %f\n", 
    wVector[0], wVector[1], wVector[2], wVector[3]);

  return 0;
}
*/

#ifdef _WIN32
#define MATRIX_VECTOR_MULTIPLY4X4_SSE(AMatrix,vVector,wVector) MATRIX_VECTOR_MULTIPLY4X4(AMatrix,vVector,wVector)
#endif

#ifdef __GNUC__
#define MATRIX_VECTOR_MULTIPLY4X4_SSE(AMatrix,vVector,wVector)\
\
asm(\
    "xorps %%xmm0, %%xmm0\n\t"\
    "xorps %%xmm1, %%xmm1\n\t"\
    "xorps %%xmm2, %%xmm2\n\t"\
    "xorps %%xmm3, %%xmm3\n\t"\
    "\n\t"\
    /* load vVector into xmm7 */\
    "movaps (%1), %%xmm7\n\t"\
    "\n\t"\ 
    /* multiply row 0 */\
    "movaps (%0), %%xmm0\n\t"\
    "mulps  %%xmm7, %%xmm0\n\t"\
    "\n\t"\
    /* multiply row 1 */\
    "movaps 16(%0), %%xmm1\n\t"\
    "mulps  %%xmm7, %%xmm1\n\t"\
    "\n\t"\
    /* multiply row 2 */\
    "movaps 32(%0), %%xmm2\n\t"\
    "mulps  %%xmm7, %%xmm2\n\t"\
    "\n\t"\
    /* multiply row 3 */\
    "movaps 48(%0), %%xmm3\n\t"\
    "mulps  %%xmm7, %%xmm3\n\t"\
    "\n\t"\
    /* sum the results  */\
    "movaps %%xmm0, %%xmm4\n\t"\
    "shufps $0x44, %%xmm1, %%xmm4\n\t"\
    "movaps %%xmm0, %%xmm5\n\t"\
    "shufps $0xEE, %%xmm1, %%xmm5\n\t"\
    "addps  %%xmm5, %%xmm4\n\t"\
    "\n\t"\
    "movaps %%xmm2, %%xmm6\n\t"\
    "shufps $0x44, %%xmm3, %%xmm6\n\t"\
    "movaps %%xmm2, %%xmm7\n\t"\
    "shufps $0xEE, %%xmm3, %%xmm7\n\t"\
    "addps  %%xmm7, %%xmm6\n\t"\
    "\n\t"\
    "movaps %%xmm4, %%xmm5\n\t"\
    "movaps %%xmm4, %%xmm7\n\t"\
    "shufps $0x88, %%xmm6, %%xmm5\n\t"\
    "shufps $0xDD, %%xmm6, %%xmm7\n\t"\
    "addps  %%xmm7, %%xmm5\n\t"\
    "\n\t"\
    /* store result */\
    "movaps %%xmm5, (%2)\n\t"\
    : /* no output registers */\
    : "r" (AMatrix), "r" (vVector), "r" (wVector)\
    : "%xmm0", "%xmm1", "%xmm2", "%xmm3",\
      "%xmm4", "%xmm5", "%xmm6", "%xmm7", "memory"\
  );
#endif

#endif
