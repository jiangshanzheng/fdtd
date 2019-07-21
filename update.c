/*************************************************************************
	> File Name: update.c
	> Author: 
	> Mail: 
	> Created Time: Sat Jul 20 16:54:59 2019
 ************************************************************************/
#include"fdtd.h"
#include"update.h"
#include<stdio.h>
//q in E or q+1/2 in H
int updateH(Grid *g){
	int mm,nn;
	switch (g->type){
	case OneD:
		exit(127);
	case TMz:
		for(mm=0; mm < g->sizeX;mm++){
			for(nn=0;nn < g->sizeY -1;nn++)
				//printf("mm=%d,nn=%d\n",mm,nn); fflush(stderr);
				Hx(g,mm,nn) = Hx(g,mm,nn) - (Ez(g,mm,nn+1) - Ez(g,mm,nn)) * g->cdtd;
		}
		for(mm=0; mm < g->sizeX-1;mm++){
			for(nn=0;nn < g->sizeY;nn++)
			Hy(g,mm,nn) = Hy(g,mm,nn) + (Ez(g,mm+1,nn) - Ez(g,mm,nn)) * g->cdtd;
		}
	default: 
		break;
	}
	return 0;

}


int updateE(Grid *g){
	int mm,nn;
	switch (g->type){
	case OneD:
		exit(127);
	case TMz:
		for(mm=1; mm < g->sizeX-1;mm++){
			for(nn=1;nn < g->sizeY-1;nn++)
			Ez(g,mm,nn) = Ez(g,mm,nn) + ( ( Hy(g,mm,nn) - Hy(g,mm-1,nn) ) - (Hx(g,mm,nn) - Hx(g,mm,nn-1) ) ) * g->cdtd;
		}
	default:
		 break;
	}
	return 0;

}
