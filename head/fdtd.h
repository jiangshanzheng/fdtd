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
	//TEz
	double * ex, *ey, * hz;

	int sizeX,sizeY,sizeZ;
	int Ti,sizeT;
	int type;
	double cdtd;
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

/* TEz */
//Hz M=x-1 N=y-1
#define Ex(G,M,N) *(G->hx + (M) * (G->sizeY) + (N))
#define Ey(G,M,N) *(G->hy + (M) * (G->sizeY-1) + (N))
#define Hz(G,M,N) *(G->ez + (M) * (G->sizeY-1) + (N))
#endif
