#include<stdio.h>
#include<math.h>
#include"fdtd.h"
#include"abc.h"
//for b to far side
//six elements for N (num of Y) row in y (L R boundary)
#define EzL(M,N,Q)  *(TMzBL_012_01 + N * 6 + Q * 3 + M  )
#define EzR(M,N,Q)  *(TMzBR_012_01 + N * 6 + Q * 3 + M  )

//six elements for M (num of X) col in x (B T boundary)
#define EzB(M,N,Q)  *(TMzBB_012_01 + M * 6 + Q * 3 + N  )
#define EzT(M,N,Q)  *(TMzBT_012_01 + M * 6 + Q * 3 + N  )

#define EyL(M,N,Q)  *(TEzBL_012_01 + N * 6 + Q * 3 + M  )
#define EyR(M,N,Q)  *(TEzBR_012_01 + N * 6 + Q * 3 + M  )

//six elements for M (num of X) col in x (B T boundary)
#define ExB(M,N,Q)  *(TEzBB_012_01 + M * 6 + Q * 3 + N  )
#define ExT(M,N,Q)  *(TEzBT_012_01 + M * 6 + Q * 3 + N  )



//2ed order abc
static double coef[3];
//inner two points and one point on boundary for past q and q-1 I0p I1p I2p I0pp I1pp I2pp 
//L R B T
static double * TMzBL_012_01;
static double * TMzBR_012_01;
static double * TMzBB_012_01;
static double * TMzBT_012_01;

static double * TEzBL_012_01;
static double * TEzBR_012_01;
static double * TEzBB_012_01;
static double * TEzBT_012_01;



void ABCInit(Grid *g){
//alloc big is not a problem
	switch (g->type){
		case OneD:
			exit(127);
			break;
		case TMz:
			ALLOC_2D(TMzBL_012_01,g->LL,6,double);
			ALLOC_2D(TMzBR_012_01,g->LR,6,double);
			ALLOC_2D(TMzBB_012_01,g->LB,6,double);
			ALLOC_2D(TMzBT_012_01,g->LT,6,double);
			break;
		case TEz:
			ALLOC_2D(TEzBL_012_01,g->LL,6,double);
			ALLOC_2D(TEzBR_012_01,g->LR,6,double);
			ALLOC_2D(TEzBB_012_01,g->LB,6,double);
			ALLOC_2D(TEzBT_012_01,g->LT,6,double);
			break;
		}
	//Sc
	double tmp1 = g->cdtds;
	double tmp2 = - 1./( 1./tmp1 + 2. + tmp1);
	coef[0] = (1./tmp1 - 2. + tmp1)*tmp2;
	coef[1] = 2.*(tmp1 - 1./tmp1)*tmp2;
	coef[2] = -4.*(1./tmp1 + tmp1)*tmp2;
}

int abc(Grid *g){

	int mm,nn;

	int bmm,bnn;
	switch (g->type){
		case OneD:
			exit(127);
			break;
		case TMz:
			//L
			for(nn=g->BLb;nn<g->BLt;nn++){
				bnn = nn - g->BLb;
				Ez(g,0,nn) = coef[0]*(Ez(g,2,nn) + EzL(0,bnn,1) )+ coef[1]* (EzL(0,bnn,0) + EzL(2,bnn,0) - Ez(g,1,nn) - EzL(1,bnn,1)) + coef[2] * EzL(1,bnn,0) - EzL(2,bnn,1);
				for(mm=0;mm<3;mm++){
					EzL(mm,bnn,1) = EzL(mm,bnn,0);
					EzL(mm,bnn,0) = Ez(g,mm,nn);

				}
			}

			//R
			for(nn=g->BRb;nn<g->BRt;nn++){
				bnn = nn - g->BRb;
				Ez(g,g->sizeX-1,nn) = coef[0]*(Ez(g,g->sizeX-3,nn) + EzR(0,bnn,1)) + coef[1]* (EzR(0,bnn,0) + EzR(2,bnn,0) - Ez(g,g->sizeX-2,nn) - EzR(1,bnn,1)) + coef[2] * EzR(1,bnn,0) - EzR(2,bnn,1);
				for(mm=0;mm<3;mm++){
					EzR(mm,bnn,1) = EzR(mm,bnn,0);
					EzR(mm,bnn,0) = Ez(g,(g->sizeX-1)-mm,nn);

				}
			}
			//B
			for(mm=g->BBl;mm<g->BBr;mm++){
				bmm = mm - g->BBl;
				Ez(g,mm,0) = coef[0]*(Ez(g,mm,2) + EzB(bmm,0,1) )+ coef[1]* (EzB(bmm,0,0) + EzB(bmm,2,0) - Ez(g,mm,1) - EzB(bmm,1,1)) + coef[2] * EzB(bmm,1,0) - EzB(bmm,2,1);
				for(nn=0;nn<3;nn++){
					EzB(bmm,nn,1) = EzB(bmm,nn,0);
					EzB(bmm,nn,0) = Ez(g,mm,nn);

				}

			}

			//T
			for(mm=g->BTl;mm<g->BTr;mm++){
				bmm = mm - g->BTl;
				Ez(g,mm,g->sizeY-1) = coef[0]*(Ez(g,mm,g->sizeY-3) + EzT(bmm,0,1) )+ coef[1]* (EzT(bmm,0,0) + EzT(bmm,2,0) - Ez(g,mm,g->sizeY-2) - EzT(bmm,1,1)) + coef[2] * EzT(bmm,1,0) - EzT(bmm,2,1);
				for(nn=0;nn<3;nn++){
					EzT(bmm,nn,1) = EzT(bmm,nn,0);
					EzT(bmm,nn,0) = Ez(g,mm,(g->sizeY-1) - nn);

				}
			}
			break;


		case TEz:
			//L
			for(nn=g->BLb;nn<g->BLt;nn++){
				bnn = nn - g->BLb;
				Ey(g,0,nn) = coef[0]*(Ey(g,2,nn) + EyL(0,bnn,1) )+ coef[1]* (EyL(0,bnn,0) + EyL(2,bnn,0) - Ey(g,1,nn) - EyL(1,bnn,1)) + coef[2] * EyL(1,bnn,0) - EyL(2,bnn,1);
				for(mm=0;mm<3;mm++){
					EyL(mm,bnn,1) = EyL(mm,bnn,0);
					EyL(mm,bnn,0) = Ey(g,mm,nn);

				}
			}

			//R
			for(nn=g->BRb;nn<g->BRt;nn++){
				bnn = nn - g->BRb;
				Ey(g,g->sizeX-1,nn) = coef[0]*(Ey(g,g->sizeX-3,nn) + EyR(0,bnn,1)) + coef[1]* (EyR(0,bnn,0) + EyR(2,bnn,0) - Ey(g,g->sizeX-2,nn) - EyR(1,bnn,1)) + coef[2] * EyR(1,bnn,0) - EyR(2,bnn,1);
				for(mm=0;mm<3;mm++){
					EyR(mm,bnn,1) = EyR(mm,bnn,0);
					EyR(mm,bnn,0) = Ey(g,(g->sizeX-1)-mm,nn);

				}
			}
			//B
			for(mm=g->BBl;mm<g->BBr;mm++){
				bmm = mm - g->BBl;
				Ex(g,mm,0) = coef[0]*(Ex(g,mm,2) + ExB(bmm,0,1) )+ coef[1]* (ExB(bmm,0,0) + ExB(bmm,2,0) - Ex(g,mm,1) - ExB(bmm,1,1)) + coef[2] * ExB(bmm,1,0) - ExB(bmm,2,1);
				for(nn=0;nn<3;nn++){
					ExB(bmm,nn,1) = ExB(bmm,nn,0);
					ExB(bmm,nn,0) = Ex(g,mm,nn);

				}

			}

			//T
			for(mm=g->BTl;mm<g->BTr;mm++){
				bmm = mm - g->BTl;
				Ex(g,mm,g->sizeY-1) = coef[0]*(Ex(g,mm,g->sizeY-3) + ExT(bmm,0,1) )+ coef[1]* (ExT(bmm,0,0) + ExT(bmm,2,0) - Ex(g,mm,g->sizeY-2) - ExT(bmm,1,1)) + coef[2] * ExT(bmm,1,0) - ExT(bmm,2,1);
				for(nn=0;nn<3;nn++){
					ExT(bmm,nn,1) = ExT(bmm,nn,0);
					ExT(bmm,nn,0) = Ex(g,mm,(g->sizeY-1) - nn);

				}
			}
			break;

	}


	return 0;

}

