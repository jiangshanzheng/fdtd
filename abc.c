#include<stdio.h>
#include"fdtd.h"
#include"abc.h"
//for b to far side
//six elements for N (num of Y) row in y (L R boundary)
#define L(g,M,N,Q)  *(B0_012_01 + N * 6 + Q * 3 + M  )
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


	int L0 = g->sizeY;
	int L1 = g->sizeY;
	int L2 = g->sizeX;
	int L3 = g->sizeX;

	ALLOC_2D(B0_012_01,L0,6,double);
	ALLOC_2D(B1_012_01,L1,6,double);
	ALLOC_2D(B2_012_01,L2,6,double);
	ALLOC_2D(B3_012_01,L3,6,double);

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
		case TMz:
			//L
			for(nn=0;nn<g->sizeY;nn++){
				Ez(g,0,nn) = coef[0]*(Ez(g,2,nn) + L(g,0,nn,1) )+ coef[1]* (L(g,0,nn,0) + L(g,2,nn,0) - Ez(g,1,nn) - L(g,1,nn,1)) + coef[2] * L(g,1,nn,0) - L(g,2,nn,1);
				for(mm=0;mm<3;mm++){
					L(g,mm,nn,1) = L(g,mm,nn,0);
					L(g,mm,nn,0) = Ez(g,mm,nn);

				}

			}

			//R
			for(nn=0;nn<g->sizeY;nn++){
				Ez(g,L1-1,nn) = coef[0]*(Ez(g,L1-3,nn) + R(g,0,nn,1)) + coef[1]* (R(g,0,nn,0) + R(g,2,nn,0) - Ez(g,L1-2,nn) - R(g,1,nn,1)) + coef[2] * R(g,1,nn,0) - R(g,2,nn,1);
				for(mm=0;mm<3;mm++){
					L(g,mm,nn,1) = L(g,mm,nn,0);
					L(g,mm,nn,0) = Ez(g,(L1-1)-mm,nn);

				}
			}
			//B
			for(mm=0;mm<g->sizeX;mm++){
				Ez(g,mm,0) = coef[0]*(Ez(g,mm,2) + B(g,mm,0,1) )+ coef[1]* (B(g,mm,0,0) + B(g,mm,2,0) - Ez(g,mm,1) - B(g,mm,1,1)) + coef[2] * B(g,mm,1,0) - B(g,mm,2,1);
				for(nn=0;nn<3;nn++){
					L(g,mm,nn,1) = L(g,mm,nn,0);
					L(g,mm,nn,0) = Ez(g,mm,nn);

				}

			}

			//T
			for(mm=0;mm<g->sizeX;mm++){
				Ez(g,mm,L3-1) = coef[0]*(Ez(g,mm,L3-3) + L(g,mm,0,1) )+ coef[1]* (L(g,mm,0,0) + L(g,mm,2,0) - Ez(g,mm,L3-2) - L(g,mm,1,1)) + coef[2] * L(g,mm,1,0) - L(g,mm,2,1);
				for(nn=0;nn<3;nn++){
					T(g,mm,nn,1) = T(g,mm,nn,0);
					T(g,mm,nn,0) = Ez(g,mm,(L3-1) - nn);

				}
			}

	}


	return 0;

}

