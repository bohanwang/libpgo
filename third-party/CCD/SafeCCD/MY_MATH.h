//**************************************************************************************
//  Copyright (C) 2002 - 2004, Huamin Wang
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//**************************************************************************************
// MY_MATH library.
// Contents:
//		Min and Max functions, float2integer functions, random functions
//		Vector and matrix functions
//**************************************************************************************
#ifndef __MY_MATH_H__
#define __MY_MATH_H__

#define MY_PI		3.14159265358979323846f
#define PI_180		0.0174532925f
#define INV_PI		0.31830988618379067154f
#define INV_TWOPI	0.15915494309189533577f
#define MY_INFINITE 999999999
#define INV_255		0.00392156862745098039f

#include <string.h>
#include <stdio.h>
#include <math.h>

#if defined(__GNUG__)
#define FORCEINLINE //__attribute__((always_inline))
#else
#define FORCEINLINE __forceinline
#endif


#define ADD(a, b, c)	{c[0]=a[0]+b[0]; c[1]=a[1]+b[1]; c[2]=a[2]+b[2];}
#define SUB(a, b, c)	{c[0]=a[0]-b[0]; c[1]=a[1]-b[1]; c[2]=a[2]-b[2];}
#define DOT(a, b)		(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define CROSS(a, b, r)	{r[0]=a[1]*b[2]-a[2]*b[1]; r[1]=a[2]*b[0]-a[0]*b[2]; r[2]=a[0]*b[1]-a[1]*b[0];}
#define	MIN(a,b)		((a)<(b)?(a):(b))
#define	MAX(a,b)		((a)>(b)?(a):(b))
#define CLAMP(a, l, h)  (((a)>(h))?(h):(((a)<(l))?(l):(a)))
#define SIGN(a)			((a)<0?-1:1)




//**************************************************************************************
// Min and Max Functions.
//**************************************************************************************
template<class T> FORCEINLINE 
T Min(const T a, const T b, const T c)
{
	T r=a;
	if(b<r) r=b;
	if(c<r) r=c; 
	return r;
}

template<class T> FORCEINLINE 
T Min(const T a, const T b, const T c, const T d)
{
	T r=a;
	if(b<r) r=b;
	if(c<r) r=c; 
	if(d<r) r=d;
	return r;
}

template<class T> FORCEINLINE 
T Max(const T a, const T b, const T c)	
{
	T r=a;
	if(b>r) r=b;
	if(c>r) r=c; 
	return r;
}

template<class T> FORCEINLINE 
T Max(const T a, const T b, const T c, const T d)
{
	T r=a;
	if(b>r) r=b;
	if(c>r) r=c; 
	if(d>r) r=d; 
	return r;
}

template<class T> FORCEINLINE 
T Max_By_Abs(const T a, const T b)
{
	return fabsf(a)>fabsf(b)?a:b;
}

template<class T> FORCEINLINE
T Min_By_Abs(const T a, const T b)
{
	return fabsf(a)<fabsf(b)?a:b;
}

//**************************************************************************************
// Integer Functions.
//**************************************************************************************
#if (defined(__linux__) && defined(__i386__)) || defined(WIN32)
#define FAST_INT 1
#endif
#define _doublemagicroundeps		(.5-1.4e-11)	//almost .5f = .5f - 1e^(number of exp bit)

FORCEINLINE
int Round(double val) 
{
#ifdef FAST_INT
#define _doublemagic				double(6755399441055744.0)	//2^52 * 1.5,  uses limited precision to floor
	val		= val + _doublemagic;
	return ((long*)&val)[0];
#undef _doublemagic
#else
	return int (val+_doublemagicroundeps);
#endif
}

FORCEINLINE 
int Float_to_Int(double val) 
{
#ifdef FAST_INT
	return (val<0) ?  Round(val+_doublemagicroundeps) :
	Round(val-_doublemagicroundeps);
#else
	return (int)val;
#endif
}

FORCEINLINE 
int Floor(double val) 
{
#ifdef FAST_INT
	return Round(val - _doublemagicroundeps);
#else
	return (int)floorf(val);
#endif
}

FORCEINLINE 
int Ceiling(double val) 
{
#ifdef FAST_INT
	return Round(val + _doublemagicroundeps);
#else
	return (int)ceilf(val);
#endif
}

//**************************************************************************************
// Math Functions.
//**************************************************************************************

template<class T> FORCEINLINE 
void Swap(T &a, T &b)
{
	T c=a; a=b; b=c;
}

//**************************************************************************************
// Random Functions. RandomFloat and RandomUInt
//**************************************************************************************
//  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//	are met:
//		1.	Redistributions of source code must retain the above copyright
//			notice, this list of conditions and the following disclaimer.
//		2.	Redistributions in binary form must reproduce the above copyright
//			notice, this list of conditions and the following disclaimer in the
//			documentation and/or other materials provided with the distribution.
//		3.	The names of its contributors may not be used to endorse or promote
//			products derived from this software without specific prior written
//			permission.
//	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
static unsigned long mt[N];		/* the array for the state vector  */
static int mti=N+1;				/* mti==N+1 means mt[N] is not initialized */

FORCEINLINE 
void init_genrand(unsigned long seed) 
{
	mt[0]= seed & 0xffffffffUL;
	for (mti=1; mti<N; mti++) {
		mt[mti] =
		(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
}
FORCEINLINE 
unsigned long genrand_int32(void)
{
	unsigned long y;
	static unsigned long mag01[2]={0x0UL, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N) { /* generate N words at one time */
		int kk;

		if (mti == N+1)   /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* default initial seed */

		for (kk=0;kk<N-M;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (;kk<N-1;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}
	y = mt[mti++];
	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);
	return y;
}

FORCEINLINE 
float RandomFloat()			/* generates a random number on [0,1)-real-interval */
{
	return genrand_int32()*((float)1.0/(float)4294967296.0);	/* divided by 2^32 */
}

FORCEINLINE 
float RandomFloat2()		/* generates a random number on [0,1]-real-interval */
{
	return genrand_int32()*((float)1.0/(float)4294967295.0);	/* divided by 2^32-1 */
}

FORCEINLINE 
unsigned long RandomUInt() 
{
	return genrand_int32();
}

#undef N
#undef M
#undef MATRIX_A 
#undef UPPER_MASK
#undef LOWER_MASK


//**************************************************************************************
// 2D matrix functions.
//**************************************************************************************
template <class T> FORCEINLINE
void Matrix_Inverse_2(T *A, T *R)
{
	T inv_det=1/(A[0]*A[3]-A[1]*A[2]);
	R[0]= A[3]*inv_det;
	R[1]=-A[1]*inv_det;
	R[2]=-A[2]*inv_det;
	R[3]= A[0]*inv_det;
}

template <class T> FORCEINLINE
void Matrix_Transpose_2(T *A, T *R)
{
	memcpy(R, A, sizeof(T)*4);
	Swap(R[1], R[2]);
}

template <class T> FORCEINLINE
void Matrix_Product_2(T *A, T *B, T *R)
{
	T temp_R[4];
	temp_R[0]=A[0]*B[0]+A[1]*B[2];
	temp_R[1]=A[0]*B[1]+A[1]*B[3];
	temp_R[2]=A[2]*B[0]+A[3]*B[2];
	temp_R[3]=A[2]*B[1]+A[3]*B[3];
	R[0]=temp_R[0];
	R[1]=temp_R[1];
	R[2]=temp_R[2];
	R[3]=temp_R[3];
}

template <class T> FORCEINLINE
void Matrix_Vector_Product_2(T *A, T *x, T *r)	//r=A*x
{
	r[0]=A[0]*x[0]+A[1]*x[1];
	r[1]=A[2]*x[0]+A[3]*x[1];	
}

template <class T> FORCEINLINE
void Matrix_T_Product_2(T *A, T *B, T *R)
{
	for(int i=0; i<2; i++)
	for(int j=0; j<2; j++)
		R[i*2+j]=A[i]*B[j]+A[i+2]*B[j+2];	
}

template <class T> FORCEINLINE
void ED_2(T *G, T *w, T *V) // G-> V'w V
{
	T a=1;
	T b=-(G[0]+G[3]);
	T c=G[0]*G[3]-G[1]*G[2];
	T delta=(b*b-4*c);
	if(delta<0)	{delta=0;}
	else		delta=sqrtf(delta);

	w[0]=(-b+delta)*0.5f;
	w[1]=(-b-delta)*0.5f;
	T inv_length;
	a=G[0]-w[0];
	b=G[1];

	if(fabsf(a)<1e-6f && fabsf(b)<1e-6f)
	{
		V[0]=1;
		V[1]=0;
		V[2]=0;
		V[3]=1;
	}
	else
	{
		inv_length=1/sqrtf(a*a+b*b);
		V[0]= b*inv_length;
		V[1]=-a*inv_length;
		V[2]=-V[1];
		V[3]= V[0];
	}
	if(V[0]<0)
	{
		V[0]=-V[0];
		V[1]=-V[1];
		V[2]=-V[2];
		V[3]=-V[3];
	}	
	if(w[0]<0)	w[0]=0;
	if(w[1]<0)	w[1]=0;
}

//**************************************************************************************
// 3D vector functions.
//**************************************************************************************
template <class T> FORCEINLINE 
T Magnitude(T *x)
{
	return sqrtf(DOT(x,x));
}

template <class T> FORCEINLINE 
T Normalize(T *x)
{
	T m=Magnitude(x);
	if(m<1e-14f)	return m;//{printf("ERROR: vector cannot be normalized.\n"); return m;}
	T inv_m=1/m;
	x[0]*=inv_m;
	x[1]*=inv_m;
	x[2]*=inv_m;
	return m;
}

template <class T> FORCEINLINE 
T Area_Squared(T* V0, T* V1, T* V2)
{
	T E10[3], E20[3], N[3];
	E10[0]=V1[0]-V0[0];
	E10[1]=V1[1]-V0[1];
	E10[2]=V1[2]-V0[2];
	E20[0]=V2[0]-V0[0];
	E20[1]=V2[1]-V0[1];
	E20[2]=V2[2]-V0[2];
	Cross(E10, E20, N);
	return Dot(N, N);
}

template <class T> FORCEINLINE
T Normal(T *p0, T *p1, T *p2, T *normal)
{
	T e0[3], e1[3];
	for(int i=0; i<3; i++)
	{
		e0[i]=p1[i]-p0[i];
		e1[i]=p2[i]-p0[i];
	}
	CROSS(e0, e1, normal);
	return Normalize(normal);
}

//**************************************************************************************
// 3D matrix functions.
//**************************************************************************************

template <class T> FORCEINLINE
void Matrix_Add_3(T *A, T a, T *B, T b, T *R)	//R=aA+bB
{
	for(int i=0; i<9; i++)
		R[i]=A[i]*a+B[i]*b;
}

template <class T> FORCEINLINE
void Matrix_Add_3(T *A, T *B, T *R)				//R=A+B
{
	for(int i=0; i<9; i++)
		R[i]=A[i]+B[i];
}

template <class T> FORCEINLINE
void Matrix_Substract_3(T *A, T *B, T *R)		//R=A-B
{
	for(int i=0; i<9; i++)
		R[i]=A[i]-B[i];
}

template <class T> FORCEINLINE
void Matrix_Inverse_3(T *A, T *R)				//R=inv(A)
{
	R[0]=A[4]*A[8]-A[7]*A[5];
	R[1]=A[7]*A[2]-A[1]*A[8];
	R[2]=A[1]*A[5]-A[4]*A[2];
	R[3]=A[5]*A[6]-A[3]*A[8];
	R[4]=A[0]*A[8]-A[2]*A[6];
	R[5]=A[2]*A[3]-A[0]*A[5];
	R[6]=A[3]*A[7]-A[4]*A[6];
	R[7]=A[1]*A[6]-A[0]*A[7];
	R[8]=A[0]*A[4]-A[1]*A[3];
	T inv_det=1/(A[0]*R[0]+A[3]*R[1]+A[6]*R[2]);
	for(int i=0; i<9; i++)
		R[i]*=inv_det;
}

template <class T> FORCEINLINE					//R=A'
void Matrix_Transpose_3(T *A, T *R)
{
	memcpy(R, A, sizeof(T)*9);
	Swap(R[1], R[3]);
	Swap(R[2], R[6]);
	Swap(R[5], R[7]);
}

template <class T> FORCEINLINE
void Matrix_Factorization_3(T *A, T *R)			//R=chol(A), Chelosky factorization
{
	R[0]=sqrtf(A[0]);
	R[1]=A[1]/R[0];
	R[2]=A[2]/R[0];
	R[3]=0;
	R[4]=sqrtf(A[4]-R[1]*R[1]);
	R[5]=(A[5]-R[1]*R[2])/R[4];
	R[6]=0;
	R[7]=0;
	R[8]=sqrtf(A[8]-R[2]*R[2]-R[5]*R[5]);
}

template <class T> FORCEINLINE
void Matrix_Product_3(T *A, T *B, T *R)		//R=A*B
{
	R[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
	R[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
	R[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
	R[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
	R[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
	R[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
	R[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
	R[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
	R[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}

template <class T> FORCEINLINE
void Matrix_T_Product_3(T *A, T *B, T *R)	//R=A'*B
{

	R[0]=A[0]*B[0]+A[3]*B[3]+A[6]*B[6];
	R[1]=A[0]*B[1]+A[3]*B[4]+A[6]*B[7];
	R[2]=A[0]*B[2]+A[3]*B[5]+A[6]*B[8];
	R[3]=A[1]*B[0]+A[4]*B[3]+A[7]*B[6];
	R[4]=A[1]*B[1]+A[4]*B[4]+A[7]*B[7];
	R[5]=A[1]*B[2]+A[4]*B[5]+A[7]*B[8];
	R[6]=A[2]*B[0]+A[5]*B[3]+A[8]*B[6];
	R[7]=A[2]*B[1]+A[5]*B[4]+A[8]*B[7];
	R[8]=A[2]*B[2]+A[5]*B[5]+A[8]*B[8];
}

template <class T> FORCEINLINE
void Matrix_Vector_Product_3(T *A, T *x, T *r)	//r=A*x
{
	r[0]=A[0]*x[0]+A[1]*x[1]+A[2]*x[2];
	r[1]=A[3]*x[0]+A[4]*x[1]+A[5]*x[2];
	r[2]=A[6]*x[0]+A[7]*x[1]+A[8]*x[2];
}

template <class T> FORCEINLINE
void Matrix_T_Vector_Product_3(T *A, T *x, T *r)//r=A'*x
{
	r[0]=A[0]*x[0]+A[3]*x[1]+A[6]*x[2];
	r[1]=A[1]*x[0]+A[4]*x[1]+A[7]*x[2];
	r[2]=A[2]*x[0]+A[5]*x[1]+A[8]*x[2];
}

//**************************************************************************************
// Arbitrary dimensional functions.
//**************************************************************************************
template <class T> FORCEINLINE T 
Norm(T *x, int number=3)	//infinite norm
{
	T ret=0;
	for(int i=0; i<number; i++)	
		if(ret<fabsf(x[i]))	ret=fabsf(x[i]);
	return ret;
}

template <class T> FORCEINLINE T 
Dot(T *x, T *y, int number)
{
	T ret=0;
	for(int i=0; i<number; i++)	ret+=x[i]*y[i];
	return ret;
}

template <class T> FORCEINLINE
void Matrix_Transpose(T *A, T *R, int nx, int ny)				//R=A'
{
	for(int i=0; i<nx; i++)
	for(int j=0; j<ny; j++)
		R[j*nx+i]=A[i*ny+j];
}

template <class T> FORCEINLINE
void Matrix_Product(T *A, T *B, T *R, int nx, int ny, int nz)	//R=A*B
{
	memset(R, 0, sizeof(T)*nx*nz);
	for(int i=0; i<nx; i++)
	for(int j=0; j<nz; j++)
		for(int k=0; k<ny; k++)
			R[i*nz+j]+=A[i*ny+k]*B[k*nz+j];
}

template <class T> FORCEINLINE
void Matrix_Self_Product(T *A, T *R, int nx, int ny)			//R=A'*A
{
	memset(R, 0, sizeof(T)*ny*ny);
	for(int i=0; i<ny; i++)
	for(int j=i; j<ny; j++)
	{
		for(int k=0; k<nx; k++)
			R[i*ny+j]+=A[k*ny+i]*A[k*ny+j];
		if(i!=j)	R[j*ny+i]=R[i*ny+j];		
	}
}

template <class T>
void Gaussian_Elimination(T *a, int n, T *b)
{
	int* indxc=new int[n];
	int* indxr=new int[n];
	int* ipiv =new int[n];
	int i,icol,irow,j,k,l,ll;
	T big,dum,pivinv,temp;

	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) 
	{ 
		big=0.0;	
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) 
				{
					if (ipiv[k] ==0) 
					{
						if (fabs(a[j*n+k]) >= big) 
						{
							big=fabs(a[j*n+k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);

		if (irow != icol) 
		{
			for (l=0;l<n;l++) {temp=a[irow*n+l]; a[irow*n+l]=a[icol*n+l]; a[icol*n+l]=temp;}
			temp=b[irow]; b[irow]=b[icol]; b[icol]=temp;
		}
		indxr[i]=irow; 
		indxc[i]=icol; 
		if (a[icol*n+icol] == 0.0) printf("Error: Singular Matrix in Gaussian_Elimination.");
		pivinv=1.0/a[icol*n+icol];
		a[icol*n+icol]=1.0;
		for (l=0;l<n;l++) a[icol*n+l] *= pivinv;
		b[icol] *= pivinv;

		for (ll=0;ll<n;ll++) 
			if (ll != icol) 
			{
				dum=a[ll*n+icol];
				a[ll*n+icol]=0.0;
				for (l=0;l<n;l++) a[ll*n+l] -= a[icol*n+l]*dum;
				b[ll] -= b[icol]*dum;
			}
	}

	for (l=n-1;l>1;l--) 
	{
		if (indxr[l] != indxc[l])
		for (k=0;k<n;k++)
		{
			temp=a[k*n+indxr[l]];
			a[k*n+indxr[l]]=a[k*n+indxc[l]];
			a[k*n+indxc[l]]=temp;
		}
	} 
	delete []ipiv;
	delete []indxr;
	delete []indxc;
}



#endif
