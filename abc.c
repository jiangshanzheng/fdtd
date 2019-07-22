#include<stdio.h>
#include<math.h>
#include"fdtd.h"
#include"abc.h"
//for b to far side
//six elements for N (num of Y) row in y (L R boundary)
#define L(M,N,Q)  *(B0_012_01 + N * 6 + Q * 3 + M  )
#define R(g,M,N,Q)  *(B0_012_01 + N * 6 + Q * 3 + M  )

//six elements for M (num of X) col in x (B T boundary)
#define B(g,M,N,Q)  *(B0_012_01 + M * 6 + Q * 3 + N  )
#define T(g,M,N,Q)  *(B0_012_01 + M * 6 + Q * 3 + N  )


//2ed order abc
static double coef[3];
//inner two points and one point on boundary for past q and q-1 I0p I1p I2p I0pp I1pp I2pp 
//L R B T
static double * B0_012_01;
static double * B1_012_01;
static double * B2_012_01;
static double * B3_012_01;


void BoundaryInit(Grid *g){

/*
 L0 = g->sizeY;
 L1 = g->sizeY;
 L2 = g->sizeX;
 L3 = g->sizeX;
*/	
	ALLOC_2D(B0_012_01,g->L0,6,double);
	ALLOC_2D(B1_012_01,g->L1,6,double);
	ALLOC_2D(B2_012_01,g->L2,6,double);
	ALLOC_2D(B3_012_01,g->L3,6,double);

	//Sc
	double tmp1 = g->cdtd;
	double tmp2 = - 1./( 1./tmp1 + 2. + tmp1);
	coef[0] = (1./tmp1 - 2. + tmp1)*tmp2;
	coef[1] = 2.*(tmp1 - 1./tmp1)*tmp2;
	coef[2] = -4.*(1./tmp1 + tmp1)*tmp2;
}

int abc(Grid *g){

	int L0 = g->sizeY;
	int L1 = g->sizeY;
	int L2 = g->sizeX;
	int L3 = g->sizeX;

	int mm,nn;
	switch (g->type){
		case OneD:
			exit(127);
			break;
		case TMz:
			//L
			for(nn=0;nn<g->sizeY;nn++){
				Ez(g,0,nn) = coef[0]*(Ez(g,2,nn) + L(0,nn,1) )+ coef[1]* (L(0,nn,0) + L(2,nn,0) - Ez(g,1,nn) - L(1,nn,1)) + coef[2] * L(1,nn,0) - L(2,nn,1);
				for(mm=0;mm<3;mm++){
					L(mm,nn,1) = L(mm,nn,0);
					L(mm,nn,0) = Ez(g,mm,nn);

				}
			}
/*
			//R
			for(nn=0;nn<g->sizeY;nn++){
				Ez(g,L1-1,nn) = coef[0]*(Ez(g,L1-3,nn) + R(0,nn,1)) + coef[1]* (R(0,nn,0) + R(2,nn,0) - Ez(g,L1-2,nn) - R(1,nn,1)) + coef[2] * R(1,nn,0) - R(2,nn,1);
				for(mm=0;mm<3;mm++){
					R(mm,nn,1) = R(mm,nn,0);
					R(mm,nn,0) = Ez(g,(L1-1)-mm,nn);

				}
			}
			//B
			for(mm=0;mm<g->sizeX;mm++){
				Ez(g,mm,0) = coef[0]*(Ez(g,mm,2) + B(mm,0,1) )+ coef[1]* (B(mm,0,0) + B(mm,2,0) - Ez(g,mm,1) - B(mm,1,1)) + coef[2] * B(mm,1,0) - B(mm,2,1);
				for(nn=0;nn<3;nn++){
					B(mm,nn,1) = B(mm,nn,0);
					B(mm,nn,0) = Ez(g,mm,nn);

				}

			}

			//T
			for(mm=0;mm<g->sizeX;mm++){
				Ez(mm,L3-1) = coef[0]*(Ez(mm,L3-3) + T(mm,0,1) )+ coef[1]* (T(mm,0,0) + T(mm,2,0) - Ez(g,mm,L3-2) - T(mm,1,1)) + coef[2] * T(mm,1,0) - T(mm,2,1);
				for(nn=0;nn<3;nn++){
					T(mm,nn,1) = T(mm,nn,0);
					T(mm,nn,0) = Ez(g,mm,(L3-1) - nn);

				}
			}
*/
			break;
		case TEz:
			exit(127); 
			break;

	}


	return 0;

}

