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

	ALLOC_2D( g->hx, g->sizeX , g->sizeY -1, double);
	ALLOC_2D( g->hy, g->sizeX -1 , g->sizeY, double);
	ALLOC_2D( g->ez, g->sizeX , g->sizeY, double);

	return 0;
}

