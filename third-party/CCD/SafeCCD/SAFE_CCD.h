//**************************************************************************************
//  Copyright (C) 2014 - 2014. Huamin Wang
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
//	Safe Continuous Collision Detection (SAFE_CCD Version 1.0)
//	Please use double precision as the floating point TYPE.
//	Single floating point precision should be used just for testing purposes.
//
//	To use this library, call Set_Coefficients(B) first.
//	B can be different for different VT and EE pairs, but it must stay the 
//	same over the whole animation process. Read the paper or the demo to see
//	how B should be defined.
//
//	After that, call Vertex_Triangle_CCD and Edge_Edge_CCD for detection.
//	The functions can return barycentric coordinates of the collision point
//	as well. But I do not recommend using it for collision response calculation.
//**************************************************************************************
#ifndef     __SAFE_CCD_H__
#define     __SAFE_CCD_H__

#include    "MY_MATH.h"

#include <cfloat>


template <class TYPE>
class SAFE_CCD
{
public:
    TYPE B;			//We use the same B for both vt and ee collisions.
	TYPE epsilon;

	//VV collision parameters	
	TYPE D2;		//D^2 is slightly enlarged to avoid the use of delta
    //VE collision parameters
    TYPE R;
    TYPE S;    
    TYPE rho;
    TYPE psi;
    TYPE lambda2;
    //VT collision parameters
    TYPE vt_mu;
    TYPE vt_tau;
    //EE collision parameters
    TYPE ee_mu;
    TYPE ee_tau;    
    TYPE eta;
	//Collision culling parameter
	TYPE gap;		//Maximum distance when a false positive happens
	TYPE interval;	//The interval to ensure the existence of real t
	TYPE new_vt_mu;
	TYPE new_ee_mu;
	TYPE w_bound;	//Upper bound on the weight if a collision happens
    
	
	SAFE_CCD()
	{}

//**************************************************************************************
//  Set up the toelerance coefficien values using the given bound B.
//**************************************************************************************
    void Set_Coefficients(TYPE _B)
    {
		B=_B;
		//Set up precision values
		TYPE half_epsilon;
        TYPE quad_epsilon;
		if(sizeof(TYPE)==4)
		{
			epsilon=FLT_EPSILON;
			quad_epsilon=0.015625f;
			half_epsilon=quad_epsilon*quad_epsilon;
		}
		if(sizeof(TYPE)==8)
		{
			epsilon=DBL_EPSILON;
			quad_epsilon=0.0001220703125;
			half_epsilon=quad_epsilon*quad_epsilon;
		}
		//VV collision parameters 
        D2			= 1678*B*B*half_epsilon;
		//VE collision parameters
        R			= 6*B*quad_epsilon;
        S			= 3*R;
        rho			= 3*B*B*half_epsilon;
        psi			= 725*B*B*half_epsilon;
        lambda2		= 2100*B*B*half_epsilon;        
		//VT collision parameters
        vt_mu		= 64*B*B*B*half_epsilon*quad_epsilon;
        vt_tau		= 64*B*B*half_epsilon;
		//EE collision parameters
		ee_mu		= 64*B*B*B*half_epsilon*quad_epsilon;
		ee_tau		= 64*B*B*half_epsilon;
        eta			= 4*B*quad_epsilon;
		//The culling parameter
		gap			= 41*B*quad_epsilon*2;
		interval	= 0.001;
		new_vt_mu	= vt_mu+746*B*B*B*epsilon;
		new_ee_mu	= ee_mu+746*B*B*B*epsilon;
		w_bound		= 10*B*B*B*B*interval;	//We use 10 instead 9, to acount for all of the rounding errors.
    }
    
//**************************************************************************************
//  Vertex-vertex CCD test. 
//**************************************************************************************
    bool Vertex_Vertex_CCD(	const TYPE xi_0[], const TYPE xi_1[], 
							const TYPE xj_0[], const TYPE xj_1[], TYPE &t,
							bool need_additional_test=true)
    {
        TYPE xji[3], vji[3], xji_1[3];
		SUB(xj_0, xi_0, xji);
		SUB(xj_1, xi_1, xji_1);
		SUB(xji_1, xji, vji);
		if(need_additional_test)
        {
			TYPE a=-DOT(xji, vji);
			TYPE b= DOT(vji, vji);
			if(a>0 && a<b)
			{
			    t=a/b;
			    xji[0]+=vji[0]*t;
			    xji[1]+=vji[1]*t;
			    xji[2]+=vji[2]*t;
			    if(DOT(xji, xji)<D2) return true;
			}
		}
		//Now test the case when t=1
		t=1;
		if(DOT(xji_1, xji_1)<D2) return true;
        return false;
    }

//**************************************************************************************
//  Vertex-edge CCD test. 
//	t is the collision time estimate.
//	r is the barycentric weight of the collision point: 1-r, r.
//**************************************************************************************
    bool Vertex_Edge_CCD(	const TYPE x0_0[], const TYPE x0_1[], 
							const TYPE xi_0[], const TYPE xi_1[], 
							const TYPE xj_0[], const TYPE xj_1[], TYPE &t, TYPE *r=0, 
							bool need_additional_test=true)
    {       
		//We perform basic planar collision culling before doing actual test.
		//The gap term is introduced here to ensure the consistency.
		//It is equivalent to twice the maximum distance when a false positive happens.

		//Axis plane culling
		if(MAX(x0_0[0], x0_1[0])+gap<MIN(MIN(xi_0[0], xi_1[0]), MIN(xj_0[0], xj_1[0])))		return false;
		if(MAX(x0_0[1], x0_1[1])+gap<MIN(MIN(xi_0[1], xi_1[1]), MIN(xj_0[1], xj_1[1])))		return false;
		if(MAX(x0_0[2], x0_1[2])+gap<MIN(MIN(xi_0[2], xi_1[2]), MIN(xj_0[2], xj_1[2])))		return false;
		if(MIN(x0_0[0], x0_1[0])-gap>MAX(MAX(xi_0[0], xi_1[0]), MAX(xj_0[0], xj_1[0])))		return false;
		if(MIN(x0_0[1], x0_1[1])-gap>MAX(MAX(xi_0[1], xi_1[1]), MAX(xj_0[1], xj_1[1])))		return false;
		if(MIN(x0_0[2], x0_1[2])-gap>MAX(MAX(xi_0[2], xi_1[2]), MAX(xj_0[2], xj_1[2])))		return false;
		//Diagonal plane culling
		if(MAX(x0_0[0]+x0_0[1], x0_1[0]+x0_1[1])+gap<MIN(MIN(xi_0[0]+xi_0[1], xi_1[0]+xi_1[1]), MIN(xj_0[0]+xj_0[1], xj_1[0]+xj_1[1])))	return false;
		if(MAX(x0_0[1]+x0_0[2], x0_1[1]+x0_1[2])+gap<MIN(MIN(xi_0[1]+xi_0[2], xi_1[1]+xi_1[2]), MIN(xj_0[1]+xj_0[2], xj_1[1]+xj_1[2])))	return false;
		if(MAX(x0_0[2]+x0_0[0], x0_1[2]+x0_1[0])+gap<MIN(MIN(xi_0[2]+xi_0[0], xi_1[2]+xi_1[0]), MIN(xj_0[2]+xj_0[0], xj_1[2]+xj_1[0])))	return false;
		if(MAX(x0_0[0]-x0_0[1], x0_1[0]-x0_1[1])+gap<MIN(MIN(xi_0[0]-xi_0[1], xi_1[0]-xi_1[1]), MIN(xj_0[0]-xj_0[1], xj_1[0]-xj_1[1])))	return false;
		if(MAX(x0_0[1]-x0_0[2], x0_1[1]-x0_1[2])+gap<MIN(MIN(xi_0[1]-xi_0[2], xi_1[1]-xi_1[2]), MIN(xj_0[1]-xj_0[2], xj_1[1]-xj_1[2])))	return false;
		if(MAX(x0_0[2]-x0_0[0], x0_1[2]-x0_1[0])+gap<MIN(MIN(xi_0[2]-xi_0[0], xi_1[2]-xi_1[0]), MIN(xj_0[2]-xj_0[0], xj_1[2]-xj_1[0])))	return false;
		if(MIN(x0_0[0]+x0_0[1], x0_1[0]+x0_1[1])-gap>MAX(MAX(xi_0[0]+xi_0[1], xi_1[0]+xi_1[1]), MAX(xj_0[0]+xj_0[1], xj_1[0]+xj_1[1])))	return false;
		if(MIN(x0_0[1]+x0_0[2], x0_1[1]+x0_1[2])-gap>MAX(MAX(xi_0[1]+xi_0[2], xi_1[1]+xi_1[2]), MAX(xj_0[1]+xj_0[2], xj_1[1]+xj_1[2])))	return false;
		if(MIN(x0_0[2]+x0_0[0], x0_1[2]+x0_1[0])-gap>MAX(MAX(xi_0[2]+xi_0[0], xi_1[2]+xi_1[0]), MAX(xj_0[2]+xj_0[0], xj_1[2]+xj_1[0])))	return false;
		if(MIN(x0_0[0]-x0_0[1], x0_1[0]-x0_1[1])-gap>MAX(MAX(xi_0[0]-xi_0[1], xi_1[0]-xi_1[1]), MAX(xj_0[0]-xj_0[1], xj_1[0]-xj_1[1])))	return false;
		if(MIN(x0_0[1]-x0_0[2], x0_1[1]-x0_1[2])-gap>MAX(MAX(xi_0[1]-xi_0[2], xi_1[1]-xi_1[2]), MAX(xj_0[1]-xj_0[2], xj_1[1]-xj_1[2])))	return false;
		if(MIN(x0_0[2]-x0_0[0], x0_1[2]-x0_1[0])-gap>MAX(MAX(xi_0[2]-xi_0[0], xi_1[2]-xi_1[0]), MAX(xj_0[2]-xj_0[0], xj_1[2]-xj_1[0])))	return false;

		//Calculate the vectors    
        TYPE x0i[3], x0j[3], xji[3];
        TYPE v0i[3], vji[3], v0j[3];
		for(int n=0; n<3; n++)
        {
            x0i[n]=x0_0[n]-xi_0[n];
            x0j[n]=x0_0[n]-xj_0[n];
            xji[n]=xj_0[n]-xi_0[n];
            v0i[n]=(x0_1[n]-x0_0[n])-(xi_1[n]-xi_0[n]);
            v0j[n]=(x0_1[n]-x0_0[n])-(xj_1[n]-xj_0[n]);
            vji[n]=(xj_1[n]-xj_0[n])-(xi_1[n]-xi_0[n]);
        }
        TYPE xa[3], xb[3], xc[3], temp[3];
        CROSS(v0i, vji, xa);
        CROSS(x0i, vji, temp);
        CROSS(v0i, xji, xb);
        ADD(temp, xb, xb);
        CROSS(x0i, xji, xc);


		TYPE times[128], time;
        int  time_number=0;
		//CASE 1
		times[time_number++]=1;
						
		//Perform one more culling, by predicting the minimum value of the left term in Equation 1
		if(need_additional_test)
        {
			TYPE left[3]={0, 0, 0};
			if((xa[0]+xb[0]+xc[0])*xc[0]>0 && (-xb[0]/(2*xa[0])<0 || -xb[0]/(2*xa[0])>1)) left[0]=MIN(fabs(xa[0]+xb[0]+xc[0]), fabs(xc[0]));
			if((xa[1]+xb[1]+xc[1])*xc[1]>0 && (-xb[1]/(2*xa[1])<0 || -xb[1]/(2*xa[1])>1)) left[1]=MIN(fabs(xa[1]+xb[1]+xc[1]), fabs(xc[1]));
			if((xa[2]+xb[2]+xc[2])*xc[2]>0 && (-xb[2]/(2*xa[2])<0 || -xb[2]/(2*xa[2])>1)) left[2]=MIN(fabs(xa[2]+xb[2]+xc[2]), fabs(xc[2]));
			if(left[0]+left[1]+left[2]<=S*B) 
			{
				TYPE new_xa[3]={xa[0], xa[1], xa[2]};
				TYPE new_xb[3]={xb[0], xb[1], xb[2]};
				TYPE new_xc[3]={xc[0], xc[1], xc[2]};
				for(int n=0; n<3; n++)
				{
					if(fabs(new_xa[n])<rho)	new_xa[n]=0;
					if(fabs(new_xb[n])<rho)	new_xb[n]=0;
				}
			
				//CASE 3: six possibilities
				TYPE threshold=10*B*epsilon;	//This threshold is for culling purposes.
				TYPE a, b, c;
				//X=Y
				a=-xji[0]+xji[1];
				b= vji[0]-vji[1];
				if(a>0 && b>a || a<0 && b<a)    
				{
					time=a/b;
					temp[0]=xji[0]+time*vji[0];
					temp[1]=xji[1]+time*vji[1];
					temp[2]=xji[2]+time*vji[2];
					if((temp[0]>-threshold && temp[1]>-threshold) || (temp[0]<threshold && temp[1]<threshold))
						if(fabs(temp[0])+threshold>fabs(temp[2]))
							times[time_number++]=time;
				}
				//Y=Z
				a=-xji[1]+xji[2];
				b= vji[1]-vji[2];
				if(a>0 && b>a || a<0 && b<a)
				{
					time=a/b;
					temp[0]=xji[0]+time*vji[0];
					temp[1]=xji[1]+time*vji[1];
					temp[2]=xji[2]+time*vji[2];
					if((temp[1]>-threshold && temp[2]>-threshold) || (temp[1]<threshold && temp[2]<threshold))
						if(fabs(temp[1])+threshold>fabs(temp[0]))
							times[time_number++]=time;  
				}
				//Z=X
				a=-xji[2]+xji[0];
				b= vji[2]-vji[0];
				if(a>0 && b>a || a<0 && b<a)
				{
					time=a/b;
					temp[0]=xji[0]+time*vji[0];
					temp[1]=xji[1]+time*vji[1];
					temp[2]=xji[2]+time*vji[2];
					if((temp[2]>-threshold && temp[0]>-threshold) || (temp[2]<threshold && temp[0]<threshold))
						if(fabs(temp[2])+threshold>fabs(temp[1]))
							times[time_number++]=time;
				}

				//X=-Y
				a=-xji[0]-xji[1];
				b= vji[0]+vji[1];
				if(a>0 && b>a || a<0 && b<a)
				{
					time=a/b;
					temp[0]=xji[0]+time*vji[0];
					temp[1]=xji[1]+time*vji[1];
					temp[2]=xji[2]+time*vji[2];
					if((temp[0]>-threshold && temp[1]<threshold) || (temp[0]<threshold && temp[1]>-threshold))
						if(fabs(temp[0])+threshold>fabs(temp[2]))
							times[time_number++]=time;
				}
				//Y=-Z
				a=-xji[1]-xji[2];
				b= vji[1]+vji[2];
				if(a>0 && b>a || a<0 && b<a)
				{
					time=a/b;
					temp[0]=xji[0]+time*vji[0];
					temp[1]=xji[1]+time*vji[1];
					temp[2]=xji[2]+time*vji[2];
					if((temp[1]>-threshold && temp[2]<threshold) || (temp[1]<threshold && temp[2]>-threshold))
						if(fabs(temp[1])+threshold>fabs(temp[0]))
							times[time_number++]=time;
				}
				//Z=-X
				a=-xji[2]-xji[0];
				b= vji[2]+vji[0];
				if(a>0 && b>a || a<0 && b<a)
				{
					time=a/b;
					temp[0]=xji[0]+time*vji[0];
					temp[1]=xji[1]+time*vji[1];
					temp[2]=xji[2]+time*vji[2];
					if((temp[2]>-threshold && temp[0]<threshold) || (temp[2]<threshold && temp[0]>-threshold))
						if(fabs(temp[2])+threshold>fabs(temp[1]))
							times[time_number++]=time;
				}

				//CASE 4
				TYPE root[3];
				int root_number;
				Quadratic_Solver(new_xa[0], new_xb[0], new_xc[0], root, root_number);
				for(int r=0; r<root_number; r++)    times[time_number++]=root[r];
				Quadratic_Solver(new_xa[1], new_xb[1], new_xc[1], root, root_number);
				for(int r=0; r<root_number; r++)    times[time_number++]=root[r];
				Quadratic_Solver(new_xa[2], new_xb[2], new_xc[2], root, root_number);
				for(int r=0; r<root_number; r++)    times[time_number++]=root[r];
					
				//CASE 2
				for(int sign_x=-1; sign_x<=1; sign_x+=2)
				for(int sign_y=-1; sign_y<=1; sign_y+=2)
				for(int sign_z=-1; sign_z<=1; sign_z+=2)
				{
					a=sign_x*new_xa[0]+sign_y*new_xa[1]+sign_z*new_xa[2];
					b=sign_x*new_xb[0]+sign_y*new_xb[1]+sign_z*new_xb[2];
					c=sign_x*new_xc[0]+sign_y*new_xc[1]+sign_z*new_xc[2];
					if(a<0) continue;	//Avoid symmetric possibilities
				
					TYPE new_b;
					new_b=b-S*vji[0];
					if(new_b<0 && -new_b<2*a)	times[time_number++]=-new_b/(2*a);				
					new_b=b+S*vji[0];
					if(new_b<0 && -new_b<2*a)  	times[time_number++]=-new_b/(2*a);        
					new_b=b-S*vji[1];
					if(new_b<0 && -new_b<2*a)   times[time_number++]=-new_b/(2*a);				
					new_b=b+S*vji[1];
					if(new_b<0 && -new_b<2*a)	times[time_number++]=-new_b/(2*a);            
					new_b=b-S*vji[2];
					if(new_b<0 && -new_b<2*a)	times[time_number++]=-new_b/(2*a);
					new_b=b+S*vji[2];
					if(new_b<0 && -new_b<2*a)	times[time_number++]=-new_b/(2*a);
				}
			}
		}
		Quick_Sort(times, 0, time_number-1);

		//Proximity test here
        TYPE x0it[3], xjit[3];        
        for(int l=0; l<time_number; l++)
        {
            t=times[l];
            xjit[0]=xji[0]+vji[0]*t;
            xjit[1]=xji[1]+vji[1]*t;
            xjit[2]=xji[2]+vji[2]*t;        
            
			//Evaluate F_0ij
            TYPE f=0;
            f+=fabs((xa[0]*t+xb[0])*t+xc[0]);
            f+=fabs((xa[1]*t+xb[1])*t+xc[1]);
            f+=fabs((xa[2]*t+xb[2])*t+xc[2]);
            f-=S*Max(fabs(xjit[0]), fabs(xjit[1]), fabs(xjit[2])); 
            if(f>psi)								continue;
			//Evaluate other proximity conditions
			TYPE xjit_xjit=DOT(xjit, xjit);
			if(xjit_xjit<=lambda2)					continue;
            x0it[0]=x0i[0]+v0i[0]*t;
            x0it[1]=x0i[1]+v0i[1]*t;
            x0it[2]=x0i[2]+v0i[2]*t;   
            TYPE x0it_xjit=DOT(x0it, xjit);
            if(x0it_xjit<0 || x0it_xjit>xjit_xjit)	continue;

			//Calculate the barycentric weight
			if(r)	*r=x0it_xjit/xjit_xjit;
            return true;
		}
        return false;
    }   
    
//**************************************************************************************
//  Vertex-triangle CCD test. 
//	b2 and b3 are the barycentric weights of the collision point: 1-b2-b3, b2, b3.
//**************************************************************************************  
    bool Vertex_Triangle_CCD(TYPE x0_0[], TYPE x0_1[],
      TYPE x1_0[], TYPE x1_1[],
      TYPE x2_0[], TYPE x2_1[],
      TYPE x3_0[], TYPE x3_1[], TYPE &t, TYPE *b2=0, TYPE *b3=0)
    {
		//Axis plane culling
		if(MAX(x0_0[0], x0_1[0])+gap<Min(MIN(x1_0[0], x1_1[0]), MIN(x2_0[0], x2_1[0]), MIN(x3_0[0], x3_1[0])))		return false;
		if(MAX(x0_0[1], x0_1[1])+gap<Min(MIN(x1_0[1], x1_1[1]), MIN(x2_0[1], x2_1[1]), MIN(x3_0[1], x3_1[1])))		return false;
		if(MAX(x0_0[2], x0_1[2])+gap<Min(MIN(x1_0[2], x1_1[2]), MIN(x2_0[2], x2_1[2]), MIN(x3_0[2], x3_1[2])))		return false;
		if(MIN(x0_0[0], x0_1[0])-gap>Max(MAX(x1_0[0], x1_1[0]), MAX(x2_0[0], x2_1[0]), MAX(x3_0[0], x3_1[0])))		return false;
		if(MIN(x0_0[1], x0_1[1])-gap>Max(MAX(x1_0[1], x1_1[1]), MAX(x2_0[1], x2_1[1]), MAX(x3_0[1], x3_1[1])))		return false;
		if(MIN(x0_0[2], x0_1[2])-gap>Max(MAX(x1_0[2], x1_1[2]), MAX(x2_0[2], x2_1[2]), MAX(x3_0[2], x3_1[2])))		return false;
		//Diagonal plane culling
		if(MAX(x0_0[0]+x0_0[1], x0_1[0]+x0_1[1])+gap<Min(MIN(x1_0[0]+x1_0[1], x1_1[0]+x1_1[1]), MIN(x2_0[0]+x2_0[1], x2_1[0]+x2_1[1]), MIN(x3_0[0]+x3_0[1], x3_1[0]+x3_1[1])))	return false;
		if(MAX(x0_0[1]+x0_0[2], x0_1[1]+x0_1[2])+gap<Min(MIN(x1_0[1]+x1_0[2], x1_1[1]+x1_1[2]), MIN(x2_0[1]+x2_0[2], x2_1[1]+x2_1[2]), MIN(x3_0[1]+x3_0[2], x3_1[1]+x3_1[2])))	return false;
		if(MAX(x0_0[2]+x0_0[0], x0_1[2]+x0_1[0])+gap<Min(MIN(x1_0[2]+x1_0[0], x1_1[2]+x1_1[0]), MIN(x2_0[2]+x2_0[0], x2_1[2]+x2_1[0]), MIN(x3_0[2]+x3_0[0], x3_1[2]+x3_1[0])))	return false;
		if(MAX(x0_0[0]-x0_0[1], x0_1[0]-x0_1[1])+gap<Min(MIN(x1_0[0]-x1_0[1], x1_1[0]-x1_1[1]), MIN(x2_0[0]-x2_0[1], x2_1[0]-x2_1[1]), MIN(x3_0[0]-x3_0[1], x3_1[0]-x3_1[1])))	return false;
		if(MAX(x0_0[1]-x0_0[2], x0_1[1]-x0_1[2])+gap<Min(MIN(x1_0[1]-x1_0[2], x1_1[1]-x1_1[2]), MIN(x2_0[1]-x2_0[2], x2_1[1]-x2_1[2]), MIN(x3_0[1]-x3_0[2], x3_1[1]-x3_1[2])))	return false;
		if(MAX(x0_0[2]-x0_0[0], x0_1[2]-x0_1[0])+gap<Min(MIN(x1_0[2]-x1_0[0], x1_1[2]-x1_1[0]), MIN(x2_0[2]-x2_0[0], x2_1[2]-x2_1[0]), MIN(x3_0[2]-x3_0[0], x3_1[2]-x3_1[0])))	return false;
		if(MIN(x0_0[0]+x0_0[1], x0_1[0]+x0_1[1])-gap>Max(MAX(x1_0[0]+x1_0[1], x1_1[0]+x1_1[1]), MAX(x2_0[0]+x2_0[1], x2_1[0]+x2_1[1]), MAX(x3_0[0]+x3_0[1], x3_1[0]+x3_1[1])))	return false;
		if(MIN(x0_0[1]+x0_0[2], x0_1[1]+x0_1[2])-gap>Max(MAX(x1_0[1]+x1_0[2], x1_1[1]+x1_1[2]), MAX(x2_0[1]+x2_0[2], x2_1[1]+x2_1[2]), MAX(x3_0[1]+x3_0[2], x3_1[1]+x3_1[2])))	return false;
		if(MIN(x0_0[2]+x0_0[0], x0_1[2]+x0_1[0])-gap>Max(MAX(x1_0[2]+x1_0[0], x1_1[2]+x1_1[0]), MAX(x2_0[2]+x2_0[0], x2_1[2]+x2_1[0]), MAX(x3_0[2]+x3_0[0], x3_1[2]+x3_1[0])))	return false;
		if(MIN(x0_0[0]-x0_0[1], x0_1[0]-x0_1[1])-gap>Max(MAX(x1_0[0]-x1_0[1], x1_1[0]-x1_1[1]), MAX(x2_0[0]-x2_0[1], x2_1[0]-x2_1[1]), MAX(x3_0[0]-x3_0[1], x3_1[0]-x3_1[1])))	return false;
		if(MIN(x0_0[1]-x0_0[2], x0_1[1]-x0_1[2])-gap>Max(MAX(x1_0[1]-x1_0[2], x1_1[1]-x1_1[2]), MAX(x2_0[1]-x2_0[2], x2_1[1]-x2_1[2]), MAX(x3_0[1]-x3_0[2], x3_1[1]-x3_1[2])))	return false;
		if(MIN(x0_0[2]-x0_0[0], x0_1[2]-x0_1[0])-gap>Max(MAX(x1_0[2]-x1_0[0], x1_1[2]-x1_1[0]), MAX(x2_0[2]-x2_0[0], x2_1[2]-x2_1[0]), MAX(x3_0[2]-x3_0[0], x3_1[2]-x3_1[0])))	return false;


		TYPE times[1024];
        int  time_number=0;        
        
        TYPE x01[3], x02[3], x21[3], x31[3];
        TYPE v01[3], v02[3], v21[3], v31[3];
        for(int n=0; n<3; n++)
        {
            x01[n]=x0_0[n]-x1_0[n];
            x02[n]=x0_0[n]-x2_0[n];
            x21[n]=x2_0[n]-x1_0[n];
            x31[n]=x3_0[n]-x1_0[n];      
            v01[n]=(x0_1[n]-x0_0[n])-(x1_1[n]-x1_0[n]);
            v02[n]=(x0_1[n]-x0_0[n])-(x2_1[n]-x2_0[n]);
            v21[n]=(x2_1[n]-x2_0[n])-(x1_1[n]-x1_0[n]);
            v31[n]=(x3_1[n]-x3_0[n])-(x1_1[n]-x1_0[n]);
        }        

        TYPE a, b, c, d, temp[3];
        CROSS(v21, v31, temp);
        a=DOT(v01, temp);
        b=DOT(x01, temp);
        CROSS(v21, x31, temp);
        b+=DOT(v01, temp);
        c =DOT(x01, temp);
        CROSS(x21, v31, temp);
        b+=DOT(v01, temp);
        c+=DOT(x01, temp);
        CROSS(x21, x31, temp);
        c+=DOT(v01, temp);
        d =DOT(x01, temp);
		        
        TYPE root[3];
        int  root_number;
        Cubic_Solver(a, b, c, d, vt_mu, root, root_number);
        for(int i=0; i<root_number; i++)    times[time_number++]=root[i];

		bool need_additional_test=false;

        for(int i=0; i<time_number; i++)
        {
            t=times[i];

			//First, determine whether t is sufficiently accurate
			TYPE l_t=MAX(t-interval, 0);
			TYPE r_t=MIN(t+interval, 1);
			TYPE l_value=((a*l_t+b)*l_t+c)*l_t+d;
			TYPE r_value=((a*r_t+b)*r_t+c)*r_t+d;
			if(fabs(l_value)<=new_vt_mu || fabs(r_value)<=new_vt_mu)	need_additional_test=true;
						
			TYPE x01t[3], x02t[3], x21t[3], x31t[3];
            for(int n=0; n<3; n++)
            {
                x01t[n]=x01[n]+v01[n]*t;
                x02t[n]=x02[n]+v02[n]*t;
                x21t[n]=x21[n]+v21[n]*t;
                x31t[n]=x31[n]+v31[n]*t;            
            }
            
            TYPE N[3];
            CROSS(x21t, x31t, N);
            TYPE n_norm     = Max(fabs(N[0]), fabs(N[1]), fabs(N[2]));
            TYPE x21_norm   = Max(fabs(x21t[0]), fabs(x21t[1]), fabs(x21t[2]));
            TYPE x31_norm   = Max(fabs(x31t[0]), fabs(x31t[1]), fabs(x31t[2]));
            TYPE x32_norm   = Max(fabs(x31t[0]-x21t[0]), fabs(x31t[1]-x21t[1]), fabs(x31t[2]-x21t[2]));
            
			if(n_norm<vt_tau)									{need_additional_test=true;	continue;}
			if(eta*Max(x21_norm, x31_norm, x32_norm)>n_norm)	{need_additional_test=true;	continue;}
            TYPE nn=DOT(N, N);
            CROSS(x21t, x01t, temp);
            TYPE weight2101=DOT(N, temp);
            if(weight2101<0)
			{
				if(weight2101>=-w_bound)				need_additional_test=true;
				continue;
			}
			if(weight2101>nn)
			{
				if(weight2101<=nn+w_bound)				need_additional_test=true;
				continue;
			}
            CROSS(x01t, x31t, temp);
			TYPE weight0131=DOT(N, temp);
            if(weight0131<0)
			{
				if(weight0131>=-w_bound)				need_additional_test=true;
				continue;
			}
			if(weight0131>nn)
			{
				if(weight0131<=nn+w_bound)				need_additional_test=true;
				continue;
			}
			
			if(weight2101+weight0131<0)
			{
				if(weight2101+weight0131>=-w_bound)		need_additional_test=true;
				continue;
			}
			if(weight2101+weight0131>nn)
			{
				if(weight2101+weight0131<=nn+w_bound)	need_additional_test=true;
				continue;
			}
			
			//Calculate the barycentric weights
			if(b2)	*b2=weight0131/nn;
			if(b3)	*b3=weight2101/nn;
            return true;
        }
        

		if(Vertex_Edge_CCD(x0_0, x0_1, x1_0, x1_1, x2_0, x2_1, t, b2, need_additional_test))
		{
			if(b3)	*b3=0;
			return true;
		}
		if(Vertex_Edge_CCD(x0_0, x0_1, x2_0, x2_1, x3_0, x3_1, t, b3, need_additional_test))
		{
			if(b2)	*b2=1-*b3;
			return true;
		}
		if(Vertex_Edge_CCD(x0_0, x0_1, x1_0, x1_1, x3_0, x3_1, t, b3, need_additional_test))
		{
			if(b2)	*b2=0;
			return true;
		}
		if(Vertex_Vertex_CCD(x0_0, x0_1, x1_0, x1_1, t, need_additional_test))
		{
			if(b2)	*b2=0;
			if(b3)	*b3=0;
			return true;
		}
		if(Vertex_Vertex_CCD(x0_0, x0_1, x2_0, x2_1, t, need_additional_test))
		{
			if(b2)	*b2=1;
			if(b3)	*b3=0;
			return true;
		}
		if(Vertex_Vertex_CCD(x0_0, x0_1, x3_0, x3_1, t, need_additional_test))
		{
			if(b2)	*b2=0;
			if(b3)	*b3=1;
			return true;
		}
		
        
        return false;    
    }    
//**************************************************************************************
//  Vertex-triangle CCD test.
//	r and s are the barycentric weights of the collision point: 1-r, r; 1-s, s.
//**************************************************************************************  
    bool Edge_Edge_CCD(	TYPE x0_0[], TYPE x0_1[], 
						TYPE x1_0[], TYPE x1_1[], 
						TYPE x2_0[], TYPE x2_1[], 
						TYPE x3_0[], TYPE x3_1[], TYPE &t, TYPE *r=0, TYPE *s=0)
    {
		//Axis plane culling
		if(MAX(MAX(x0_0[0], x0_1[0]), MAX(x1_0[0], x1_1[0]))+gap<MIN(MIN(x2_0[0], x2_1[0]), MIN(x3_0[0], x3_1[0])))		return false;
		if(MAX(MAX(x0_0[1], x0_1[1]), MAX(x1_0[1], x1_1[1]))+gap<MIN(MIN(x2_0[1], x2_1[1]), MIN(x3_0[1], x3_1[1])))		return false;
		if(MAX(MAX(x0_0[2], x0_1[2]), MAX(x1_0[2], x1_1[2]))+gap<MIN(MIN(x2_0[2], x2_1[2]), MIN(x3_0[2], x3_1[2])))		return false;
		if(MIN(MIN(x0_0[0], x0_1[0]), MIN(x1_0[0], x1_1[0]))-gap>MAX(MAX(x2_0[0], x2_1[0]), MAX(x3_0[0], x3_1[0])))		return false;
		if(MIN(MIN(x0_0[1], x0_1[1]), MIN(x1_0[1], x1_1[1]))-gap>MAX(MAX(x2_0[1], x2_1[1]), MAX(x3_0[1], x3_1[1])))		return false;
		if(MIN(MIN(x0_0[2], x0_1[2]), MIN(x1_0[2], x1_1[2]))-gap>MAX(MAX(x2_0[2], x2_1[2]), MAX(x3_0[2], x3_1[2])))		return false;
		//Diagonal plane culling
		if(MAX(MAX(x0_0[0]+x0_0[1], x0_1[0]+x0_1[1]), MAX(x1_0[0]+x1_0[1], x1_1[0]+x1_1[1]))+gap<MIN(MIN(x2_0[0]+x2_0[1], x2_1[0]+x2_1[1]), MIN(x3_0[0]+x3_0[1], x3_1[0]+x3_1[1])))	return false;
		if(MAX(MAX(x0_0[1]+x0_0[2], x0_1[1]+x0_1[2]), MAX(x1_0[1]+x1_0[2], x1_1[1]+x1_1[2]))+gap<MIN(MIN(x2_0[1]+x2_0[2], x2_1[1]+x2_1[2]), MIN(x3_0[1]+x3_0[2], x3_1[1]+x3_1[2])))	return false;
		if(MAX(MAX(x0_0[2]+x0_0[0], x0_1[2]+x0_1[0]), MAX(x1_0[2]+x1_0[0], x1_1[2]+x1_1[0]))+gap<MIN(MIN(x2_0[2]+x2_0[0], x2_1[2]+x2_1[0]), MIN(x3_0[2]+x3_0[0], x3_1[2]+x3_1[0])))	return false;
		if(MAX(MAX(x0_0[0]-x0_0[1], x0_1[0]-x0_1[1]), MAX(x1_0[0]-x1_0[1], x1_1[0]-x1_1[1]))+gap<MIN(MIN(x2_0[0]-x2_0[1], x2_1[0]-x2_1[1]), MIN(x3_0[0]-x3_0[1], x3_1[0]-x3_1[1])))	return false;
		if(MAX(MAX(x0_0[1]-x0_0[2], x0_1[1]-x0_1[2]), MAX(x1_0[1]-x1_0[2], x1_1[1]-x1_1[2]))+gap<MIN(MIN(x2_0[1]-x2_0[2], x2_1[1]-x2_1[2]), MIN(x3_0[1]-x3_0[2], x3_1[1]-x3_1[2])))	return false;
		if(MAX(MAX(x0_0[2]-x0_0[0], x0_1[2]-x0_1[0]), MAX(x1_0[2]-x1_0[0], x1_1[2]-x1_1[0]))+gap<MIN(MIN(x2_0[2]-x2_0[0], x2_1[2]-x2_1[0]), MIN(x3_0[2]-x3_0[0], x3_1[2]-x3_1[0])))	return false;
		if(MIN(MIN(x0_0[0]+x0_0[1], x0_1[0]+x0_1[1]), MIN(x1_0[0]+x1_0[1], x1_1[0]+x1_1[1]))-gap>MAX(MAX(x2_0[0]+x2_0[1], x2_1[0]+x2_1[1]), MAX(x3_0[0]+x3_0[1], x3_1[0]+x3_1[1])))	return false;
		if(MIN(MIN(x0_0[1]+x0_0[2], x0_1[1]+x0_1[2]), MIN(x1_0[1]+x1_0[2], x1_1[1]+x1_1[2]))-gap>MAX(MAX(x2_0[1]+x2_0[2], x2_1[1]+x2_1[2]), MAX(x3_0[1]+x3_0[2], x3_1[1]+x3_1[2])))	return false;
		if(MIN(MIN(x0_0[2]+x0_0[0], x0_1[2]+x0_1[0]), MIN(x1_0[2]+x1_0[0], x1_1[2]+x1_1[0]))-gap>MAX(MAX(x2_0[2]+x2_0[0], x2_1[2]+x2_1[0]), MAX(x3_0[2]+x3_0[0], x3_1[2]+x3_1[0])))	return false;
		if(MIN(MIN(x0_0[0]-x0_0[1], x0_1[0]-x0_1[1]), MIN(x1_0[0]-x1_0[1], x1_1[0]-x1_1[1]))-gap>MAX(MAX(x2_0[0]-x2_0[1], x2_1[0]-x2_1[1]), MAX(x3_0[0]-x3_0[1], x3_1[0]-x3_1[1])))	return false;
		if(MIN(MIN(x0_0[1]-x0_0[2], x0_1[1]-x0_1[2]), MIN(x1_0[1]-x1_0[2], x1_1[1]-x1_1[2]))-gap>MAX(MAX(x2_0[1]-x2_0[2], x2_1[1]-x2_1[2]), MAX(x3_0[1]-x3_0[2], x3_1[1]-x3_1[2])))	return false;
		if(MIN(MIN(x0_0[2]-x0_0[0], x0_1[2]-x0_1[0]), MIN(x1_0[2]-x1_0[0], x1_1[2]-x1_1[0]))-gap>MAX(MAX(x2_0[2]-x2_0[0], x2_1[2]-x2_1[0]), MAX(x3_0[2]-x3_0[0], x3_1[2]-x3_1[0])))	return false;


        TYPE times[1024];
        int  time_number=0;        
        
        TYPE x10[3], x20[3], x32[3];
        TYPE v10[3], v20[3], v32[3];
        for(int n=0; n<3; n++)
        {
            x10[n]=x1_0[n]-x0_0[n];
            x20[n]=x2_0[n]-x0_0[n];
            x32[n]=x3_0[n]-x2_0[n];

            v10[n]=(x1_1[n]-x1_0[n])-(x0_1[n]-x0_0[n]);
            v20[n]=(x2_1[n]-x2_0[n])-(x0_1[n]-x0_0[n]);
            v32[n]=(x3_1[n]-x3_0[n])-(x2_1[n]-x2_0[n]);
        }        

        TYPE a, b, c, d, temp[3];
        CROSS(v10, v32, temp);
        a=DOT(v20, temp);
        b=DOT(x20, temp);
        CROSS(v10, x32, temp);
        b+=DOT(v20, temp);
        c =DOT(x20, temp);
        CROSS(x10, v32, temp);
        b+=DOT(v20, temp);
        c+=DOT(x20, temp);
        CROSS(x10, x32, temp);
        c+=DOT(v20, temp);
        d =DOT(x20, temp);
		
        TYPE root[3];
        int  root_number;
        Cubic_Solver(a, b, c, d, ee_mu, root, root_number);
        for(int i=0; i<root_number; i++)    times[time_number++]=root[i];
	
		bool need_additional_test=false;
		bool tag=false;
		       
        for(int i=0; i<time_number; i++)
        {
            t=times[i];

			//First, determine whether t is sufficiently accurate
			TYPE l_t=MAX(t-interval, 0);
			TYPE r_t=MIN(t+interval, 1);
			TYPE l_value=((a*l_t+b)*l_t+c)*l_t+d;
			TYPE r_value=((a*r_t+b)*r_t+c)*r_t+d;
			if(fabs(l_value)<=new_vt_mu || fabs(r_value)<=new_vt_mu)	{tag=true; need_additional_test=true;}
			
			TYPE x10t[3], x20t[3], x32t[3];
            for(int n=0; n<3; n++)
            {
                x10t[n]=x10[n]+v10[n]*t;
                x20t[n]=x20[n]+v20[n]*t;
                x32t[n]=x32[n]+v32[n]*t;      
            }
            
            TYPE N[3];
            CROSS(x10t, x32t, N);
            TYPE n_norm     = Max(fabs(N[0]), fabs(N[1]), fabs(N[2]));
            TYPE x10_norm   = Max(fabs(x10t[0]), fabs(x10t[1]), fabs(x10t[2]));
            TYPE x32_norm   = Max(fabs(x32t[0]), fabs(x32t[1]), fabs(x32t[2]));
            
			if(n_norm<ee_tau)								{need_additional_test=true;		continue;}
			if(eta*MAX(x10_norm, x32_norm)>n_norm)			{need_additional_test=true;		continue;}
            TYPE nn=DOT(N, N);
            CROSS(x20t, x10t, temp);
            TYPE weight2010=DOT(N, temp);
            if(weight2010<0)
			{
				if(weight2010>=-w_bound)	need_additional_test=true;
				continue;
			}
			if(weight2010>nn)
			{
				if(weight2010<=nn+w_bound)	need_additional_test=true;
				continue;
			}
            CROSS(x20t, x32t, temp);
			TYPE weight2032=DOT(N, temp);
            if(weight2032<0)
			{
				if(weight2032>=-w_bound)	need_additional_test=true;
				continue;
			}
			else if(weight2032>nn)			
			{
				if(weight2032<=nn+w_bound)	need_additional_test=true;
				continue;
			}
			
			if(r)	*r=weight2032/nn;
			if(s)	*s=weight2010/nn;
            return true;
        }

 		if(Vertex_Edge_CCD(x2_0, x2_1, x0_0, x0_1, x1_0, x1_1, t, r, need_additional_test))
		{
			if(s)	*s=0;				
			return true;
		}
		if(Vertex_Edge_CCD(x3_0, x3_1, x0_0, x0_1, x1_0, x1_1, t, r, need_additional_test))
		{
			if(s)	*s=1;
			return true;
		}
		if(Vertex_Edge_CCD(x0_0, x0_1, x2_0, x2_1, x3_0, x3_1, t, s, need_additional_test))
		{
			if(r)	*r=0;
			return true;
		}
		if(Vertex_Edge_CCD(x1_0, x1_1, x2_0, x2_1, x3_0, x3_1, t, s, need_additional_test))
		{
			if(r)	*r=1;
			return true;
		}
		if(Vertex_Vertex_CCD(x0_0, x0_1, x2_0, x2_1, t, need_additional_test))
		{
			if(r)	*r=0;
			if(s)	*s=0;
			return true;
		}
		if(Vertex_Vertex_CCD(x0_0, x0_1, x3_0, x3_1, t, need_additional_test))
		{
			if(r)	*r=0;
			if(s)	*s=1;
			return true;
		}
		if(Vertex_Vertex_CCD(x1_0, x1_1, x2_0, x2_1, t, need_additional_test))
		{
			if(r)	*r=1;
			if(s)	*s=0;
			return true;
		}
		if(Vertex_Vertex_CCD(x1_0, x1_1, x3_0, x3_1, t, need_additional_test))
		{
			if(r)	*r=1;
			if(s)	*s=1;
			return true;
		}
		

        return false;    
    }    
//**************************************************************************************
//  This analytic quadratic solver finds roots within (0, 1) only.
//**************************************************************************************
	void Quadratic_Solver(TYPE a, TYPE b, TYPE c, TYPE root[], int &root_number)
	{
		if(a<0)     {a=-a; b=-b; c=-c;}
		root_number=0;
		TYPE delta=b*b-4*a*c;
	
		//if(delta<0) return;		//Clamp delta to 0 if negative
		if(delta<=0)
		{
			if(-b>0 && -b<2*a)  root[root_number++]=-b/(2*a);
			return;
		}
		
		if(b<=0)
		{
			TYPE temp=-b+sqrt(delta);
			TYPE twice_c=2*c;
			TYPE twice_a=2*a;
			if(twice_c>0 && twice_c<temp)	root[root_number++]=twice_c/temp;
			if(temp<twice_a)				root[root_number++]=temp/twice_a;
		}		
		else
		{
			TYPE temp=-b-sqrt(delta);
			TYPE twice_c=2*c;
			TYPE twice_a=2*a;
			if(twice_a<temp)				root[root_number++]=temp/twice_a;
			if(twice_c<0 && temp<twice_c)	root[root_number++]=twice_c/temp;
		}
	}
//**************************************************************************************
//  This cubic solver uses the Newton-Bisection method to find roots within [0, 1].
//	It is based on the implementation in the book "Numerical Recipe in C", but there 
//	are several differences. Please see the paper for more details.
//**************************************************************************************
    void Cubic_Solver(TYPE a, TYPE b, TYPE c, TYPE d, TYPE mu, TYPE root[], int &root_number)
    {
        root_number=0;
        
		//Build the intervals. There are 4 nodes at most.
		TYPE    nodes[4]={0, 0, 0, 0};
        int     node_number=1;
        TYPE    min_max[2];
        int     min_max_number;
        Quadratic_Solver(3*a, 2*b, c, min_max, min_max_number);
        if(min_max_number==2 && min_max[0]>min_max[1])  
		{
			nodes[node_number++]=min_max[1];
			nodes[node_number++]=min_max[0];		
		}
		else
		{
			for(int i=0; i<min_max_number; i++)
				nodes[node_number++]=min_max[i];
		}
        nodes[node_number++]=1;
		        
        //Detect a root in every interval.
        for(int i=0; i<node_number-1; i++)
        {
            TYPE x1=nodes[i];
            TYPE x2=nodes[i+1];
			
			//Obtain lower and upper nodes and their function values.
            TYPE fl=((a*x1+b)*x1+c)*x1+d;
            if(fabs(fl)<mu)     {root[root_number++]=x1; continue;}
            TYPE fh=((a*x2+b)*x2+c)*x2+d;
            if(fabs(fh)<mu)     {root[root_number++]=x2; continue;}
            if(fl>0 && fh>0 || fl<0 && fh<0)    continue;			
			TYPE xl, xh;
            if(fl<0)    {xl=x1; xh=x2;}
            else        {xh=x1; xl=x2;}
			
			//Start with bisection
            TYPE dxold	= fabs(x2-x1);
            TYPE dx		= dxold;
			TYPE rts	= (x1+x2)*0.5;
			TYPE f		= ((a*rts+b)*rts+c)*rts+d;
            TYPE df		= (3*a*rts+2*b)*rts+c;
			
            int j=0;
            for(j=0; j<1024; j++)
            {
                if(((rts-xh)*df-f)*((rts-xl)*df-f)<0 && fabs(2.0*f)<fabs(dxold*df))	//Try Newton first
                {					
                    dxold=dx;
                    dx=f/df;
                    rts=rts-dx;
                    if(rts>=xh && rts>=xl || rts<=xh && rts<=xl) //Switch back to bisection if out of range 	
                    {
                        dxold=dx;
                        dx=0.5*(xh-xl);
                        rts=(xl+xh)*0.5;
                    }
                }
				else //Now do bisection
                {
                    dxold=dx;
                    dx=0.5*(xh-xl);
                    rts=(xl+xh)*0.5;                    
                }
                //Prepare for the next iteration
                f=((a*rts+b)*rts+c)*rts+d;
                if(fabs(f)<mu)     {root[root_number++]=rts; break;}
				df=(3*a*rts+2*b)*rts+c;
                if(f<0) xl=rts;
                else    xh=rts;

            }
            // for performance, I comment it out
            if(j==1024)     printf("ERROR: Fails to converge within 1000 iterations...\n");       
        }
    }
//**************************************************************************************
//  Quick sort. (Use it to sort the potential collision time list before verification)
//**************************************************************************************  
	int Quick_Sort_Partition(TYPE* t, int p, int r)
	{
		TYPE x=t[r];
		int i=p-1;
		for(int j=p; j<=r-1; j++)
			if( t[j]<=x)
			{
				i++;
				Swap(t[i], t[j]);

			}
		Swap(t[i+1], t[r]);
		return i+1;
	}
	 
	void Quick_Sort(TYPE* t, int p, int r)
	{
		if(p<r)
		{
			int q=Quick_Sort_Partition(t, p, r);
			Quick_Sort(t, p, q-1);
			Quick_Sort(t, q+1, r);
		}
	}
};

                      
#endif




