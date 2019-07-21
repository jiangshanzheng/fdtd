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
//all H can be calc
int updateH(Grid *g){
	int mm,nn;
	switch (g->type){
	case OneD:
		exit(127);
	case TMz:
		for(mm=0; mm < g->sizeX;mm++){
			for(nn=0;nn < g->sizeY -1;nn++)
				//printf("mm=%d,nn=%d\n",mm,nn); fflush(stderr);
				//wrong index in fdtd.h
				Hx(g,mm,nn) = Hx(g,mm,nn) - (Ez(g,mm,nn+1) - Ez(g,mm,nn)) * g->cdtd;
		}
		for(mm=0; mm < g->sizeX-1;mm++){
			for(nn=0;nn < g->sizeY;nn++)
				Hy(g,mm,nn) = Hy(g,mm,nn) + (Ez(g,mm+1,nn) - Ez(g,mm,nn)) * g->cdtd;
		}
	case TEz:
		for(mm=0; mm < g->sizeX-1;mm++){
			for(nn=0;nn < g->sizeY-1;nn++)
			//mm nn in H for half, in Ey mm for int, nn for half; in Ex mm for half, nn for int
			//the number of hz is M-1 * N-1 same to definition, so no cut-off on Hz, what is cut off!!!left start with mm nn, num depend on left
			Hz(g,mm,nn) = Hz(g,mm,nn) + ( ( Ex(g,mm,nn+1) - Ex(g,mm,nn) ) - (Ey(g,mm+1,nn) - Ey(g,mm,nn) ) ) * g->cdtd;
		}
	default: 
		break;
	}
	return 0;

}

//some components of E can not be calculated, so that's cut off?
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
	case TEz:
		//Ex at nn=0 nn=sizeX-1, all H start n-1
		for(mm=0; mm < g->sizeX-1;mm++){
			for(nn=1;nn < g->sizeY -1;nn++)
				Ex(g,mm,nn) = Ex(g,mm,nn) + (Hz(g,mm,nn) - Hz(g,mm,nn-1)) * g->cdtd;
		}
		//Ey at nn=0 nn=sizeX-1,H start m-1
		for(mm=1; mm < g->sizeX-1;mm++){
			for(nn=0;nn < g->sizeY-1;nn++)
				Ey(g,mm,nn) = Ey(g,mm,nn) - (Hz(g,mm,nn) - Hz(g,mm-1,nn)) * g->cdtd;
		}
	default:
		 break;
	}
	return 0;

}
