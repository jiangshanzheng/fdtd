/*************************************************************************
	> File Name: fdtd.h
	> Author: 
	> Mail: 
	> Created Time: Sat Jul 20 15:41:28 2019
 ************************************************************************/

#ifndef _FDTD_H
#define _FDTD_H
#include<stdio.h>
#include<stdlib.h>
enum GRID_TYPE {OneD, TMz, TEz,ThreeD};

struct Grid{
	//TMz
	double * hx, *hy, * ez;
	double *ezx, *ezy;
	double * chxe, *chye, * cezxe, * cezye;
	double * chxh, *chyh, * cezxh, * cezyh;
	//double * sig_my,  * sig_mx, *sig_x, * sig_y;
	//TEz
	double * ex, *ey, * hz;
	double *hzx, *hzy;
	double * cexe, *ceye, * chzxe, * chzye;
	double * cexh, *ceyh, * chzxh, * chzyh;
//boundary
	int PML_Layers;
	int LL, BLb, BLt;
	int LR, BRb, BRt;
	int LB, BBl, BBr;
	int LT, BTl, BTr;


	int sizeX,sizeY,sizeZ;
	int Ti,sizeT;
	int output,DUMP;
	int type;
	double imp0, c, cdtds, dt, dx;
};

typedef struct Grid Grid;

//define alloc macro
#define ALLOC_1D(PNTR, NUM, TYPE)                     \
	PNTR =  (TYPE * )calloc((NUM), sizeof(TYPE));   \
	if(!PNTR){				      \
		perror("ALLOC_1D");		      \
		fprintf(stderr,"Allocation error for " #PNTR ".\n");	\
		exit(-1);				\
	}						\

#define ALLOC_2D(PNTR, NUMX, NUMY, TYPE)                     \
	PNTR =  (TYPE * )calloc( (NUMX) * (NUMY) , sizeof(TYPE));   \
	if(!PNTR){				      \
		perror("ALLOC_2D");		      \
		fprintf(stderr,"Allocation error for " #PNTR ".\n");	\
		exit(-1);				\
	}						\

//index macro,M row,N col
/* TMz */
//Hy M=x-1 N=y
//fuck
#define Hx(G,M,N) *(G->hx + (M) * (G->sizeY -1) + (N))
#define Hy(G,M,N) *(G->hy + (M) * (G->sizeY) + (N))
#define Ez(G,M,N) *(G->ez + (M) * (G->sizeY) + (N))
#define Ezx(G,M,N) *(G->ezx + (M) * (G->sizeY) + (N))
#define Ezy(G,M,N) *(G->ezy + (M) * (G->sizeY) + (N))

#define cHxE(G,M,N) *(G->chxe + (M) * (G->sizeY -1) + (N))
#define cHxH(G,M,N) *(G->chxh + (M) * (G->sizeY -1) + (N))

#define cHyE(G,M,N) *(G->chye + (M) * (G->sizeY) + (N))
#define cHyH(G,M,N) *(G->chyh + (M) * (G->sizeY) + (N))

#define cEzxE(G,M,N) *(G->cezxe + (M) * (G->sizeY) + (N))
#define cEzxH(G,M,N) *(G->cezxh + (M) * (G->sizeY) + (N))

#define cEzyE(G,M,N) *(G->cezye + (M) * (G->sizeY) + (N))
#define cEzyH(G,M,N) *(G->cezyh + (M) * (G->sizeY) + (N))




/* TEz */
//Hz M=x-1 N=y-1
//N changes fastest
#define Ex(G,M,N) *(G->ex + (M) * (G->sizeY) + (N))
#define Ey(G,M,N) *(G->ey + (M) * (G->sizeY-1) + (N))
#define Hz(G,M,N) *(G->hz + (M) * (G->sizeY-1) + (N))
#define Hzx(G,M,N) *(G->hzx + (M) * (G->sizeY-1) + (N))
#define Hzy(G,M,N) *(G->hzy + (M) * (G->sizeY-1) + (N))

#define cExE(G,M,N) *(G->cexe + (M) * (G->sizeY) + (N))
#define cExH(G,M,N) *(G->cexh + (M) * (G->sizeY) + (N))

#define cEyE(G,M,N) *(G->ceye + (M) * (G->sizeY-1) + (N))
#define cEyH(G,M,N) *(G->ceyh + (M) * (G->sizeY-1) + (N))

#define cHzxE(G,M,N) *(G->chzxe + (M) * (G->sizeY-1) + (N))
#define cHzxH(G,M,N) *(G->chzxh + (M) * (G->sizeY-1) + (N))

#define cHzyE(G,M,N) *(G->chzye + (M) * (G->sizeY-1) + (N))
#define cHzyH(G,M,N) *(G->chzyh + (M) * (G->sizeY-1) + (N))

#endif
