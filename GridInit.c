/*************************************************************************
	> File Name: gridinit.c
	> Author: 
	> Mail: 
	> Created Time: Sat Jul 20 16:39:37 2019
 ************************************************************************/

#include<stdio.h>
#include<math.h>
#include"fdtd.h"
#include"GridInit.h"
int GridInit(Grid * g){

	g->type = TMz;
	g->sizeX = 101;
	g->sizeY = 101;
	g->sizeT = 300;
	g->cdtd = 1/sqrt(2.0);
	g->BLs = 0;
	g->BLe = 101;
       	g->BRs = 0;
	g->BRe = 0;
	g->BBs = 0;
	g->BBe = 0;
	g->BTs = 0;
	g->BTe = 0;
	g->L0 = g->BLe - g->BLs;
	g->L1 = g->BRe - g->BRs;
	g->L2 = g->BBe - g->BBs;
	g->L3 = g->BTe - g->BTs;

	switch (g->type){
		case OneD:
			break;
		case TMz:
			ALLOC_2D( g->hx, g->sizeX , g->sizeY -1, double);
			ALLOC_2D( g->hy, g->sizeX -1 , g->sizeY, double);
			ALLOC_2D( g->ez, g->sizeX , g->sizeY, double);
			break;
		case TEz:
			ALLOC_2D( g->ex, g->sizeX-1 , g->sizeY, double);
			ALLOC_2D( g->ey, g->sizeX , g->sizeY-1, double);
			ALLOC_2D( g->hz, g->sizeX-1 , g->sizeY-1, double);
			break;
	}

	return 0;
}

