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
	g->sizeX = 1280;
	g->sizeY = 10240;
	g->sizeT = 3000;
	g->output = 1;
	g->DUMP	= 100;
	g->cdtd = 1/sqrt(2.0);
	
	//TEz finial points must less then num it E components
	//macro will not check over boundary
	//allocation for boundy fields is large then its real size
	g->BLb = 0;
	g->BLt = g->sizeY -1 ;
	g->BLt = 0;
       	
	g->BRb = 0;
	g->BRt = g->sizeY -1;
	g->BRt = 0;
	
	g->BBl = 0;
	g->BBr = g->sizeX-1;

	g->BTl = 0;
	g->BTr = g->sizeX-1;
	g->LL = g->BLt - g->BLb;
	g->LR = g->BRt - g->BRb;
	g->LB = g->BBr - g->BBl;
	g->LT = g->BTr - g->BTl;

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

