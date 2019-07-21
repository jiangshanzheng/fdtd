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

	g->type = TEz;
	g->sizeX = 101;
	g->sizeY = 101;
	g->sizeT = 300;
	g->cdtd = 1/sqrt(2.0);
	
	switch (g->type){
		case TMz:
			ALLOC_2D( g->hx, g->sizeX , g->sizeY -1, double);
			ALLOC_2D( g->hy, g->sizeX -1 , g->sizeY, double);
			ALLOC_2D( g->ez, g->sizeX , g->sizeY, double);
		case TEz:
			ALLOC_2D( g->ex, g->sizeX-1 , g->sizeY, double);
			ALLOC_2D( g->ey, g->sizeX , g->sizeY-1, double);
			ALLOC_2D( g->hz, g->sizeX-1 , g->sizeY-1, double);
	}

	return 0;
}

