/*************************************************************************
  > File Name: update.c
  > Author: 
  > Mail: 
  > Created Time: Sat Jul 20 16:54:59 2019
 ************************************************************************/
#include"fdtd.h"
#include"update_pml.h"
#include<stdio.h>
//q in E or q+1/2 in H
//all H can be calc
//we access them for mm then nn, so the first x=0,y=* data will be storged
int updateH_pml(Grid *g){
	int mm,nn;
	double PMLs = g-> PML_Layers;
	switch (g->type){
		case OneD:
			exit(127);
			break;
		case TMz:
			for(mm=0; mm < g->sizeX;mm++){
				for(nn=0;nn < g->sizeY -1;nn++){
					//printf("mm=%d,nn=%d\n",mm,nn); fflush(stderr);
					//wrong index in fdtd.h
					Hx(g,mm,nn) = cHxH(g,mm,nn) * Hx(g,mm,nn) - cHxE(g,mm,nn) * (Ez(g,mm,nn+1) - Ez(g,mm,nn));
				}
			}
			for(mm=0; mm < g->sizeX-1;mm++){
				for(nn=0;nn < g->sizeY;nn++){
					Hy(g,mm,nn) = cHyH(g,mm,nn) * Hy(g,mm,nn) + cHyE(g,mm,nn)*(Ez(g,mm+1,nn) - Ez(g,mm,nn));
				}
			}
			break;
		case TEz:
			for(mm=0; mm < g->sizeX-1;mm++){
				for(nn=0;nn < g->sizeY-1;nn++){
					//mm nn in H for half, in Ey mm for int, nn for half; in Ex mm for half, nn for int
					//the number of hz is M-1 * N-1 same to definition, so no cut-off on Hz, what is cut off!!!left start with mm nn, num depend on left
					//Hz(g,mm,nn) = Hz(g,mm,nn) + ( ( Ex(g,mm,nn+1) - Ex(g,mm,nn) ) - (Ey(g,mm+1,nn) - Ey(g,mm,nn) ) ) * g->cdtds;
					Hzx(g,mm,nn) = cHzxH(g,mm,nn)*Hzx(g,mm,nn) + cHzxE(g,mm,nn)*( - ( Ey(g,mm,nn) - Ey(g,mm-1,nn) ));
					Hzy(g,mm,nn) = cHzyH(g,mm,nn)*Hzy(g,mm,nn) + cHzyE(g,mm,nn)*( ( Ex(g,mm,nn) - Ex(g,mm,nn-1) ));
					Hz(g,mm,nn) = Hzx(g,mm,nn) + Hzy(g,mm,nn);
				}
			}
			break;
		default: 
			break;
	}
	return 0;

}

//some components of E can not be calculated, so that's cut off?
int updateE_pml(Grid *g){
	int mm,nn;
	switch (g->type){
		case OneD:
			exit(127);
			break;
		case TMz:
			for(mm=1; mm < g->sizeX-1;mm++){
				for(nn=1;nn < g->sizeY-1;nn++){
					//printf("%.1f-%.1f-%.1f-%.1f\n",cEzxE(g,mm,nn),cEzxH(g,mm,nn),cEzyE(g,mm,nn),cEzyH(g,mm,nn));
					Ezx(g,mm,nn) = cEzxE(g,mm,nn)*Ezx(g,mm,nn) + cEzxH(g,mm,nn)*( ( Hy(g,mm,nn) - Hy(g,mm-1,nn) ));
					Ezy(g,mm,nn) = cEzyE(g,mm,nn)*Ezy(g,mm,nn) + cEzyH(g,mm,nn)*(- ( Hx(g,mm,nn) - Hx(g,mm,nn-1) ));
					Ez(g,mm,nn) = Ezx(g,mm,nn) + Ezy(g,mm,nn);
					//Ez(g,mm,nn) = Ez(g,mm,nn) + ( ( Hy(g,mm,nn) - Hy(g,mm-1,nn) ) - (Hx(g,mm,nn) - Hx(g,mm,nn-1) ) ) * g->cdtds;
				}
			}
			break;
		case TEz:
			//Ex at nn=0 nn=sizeX-1, all H start n-1
			for(mm=0; mm < g->sizeX-1;mm++){
				for(nn=1;nn < g->sizeY -1;nn++)
					Ex(g,mm,nn) = cExE(g,mm,nn)*Ex(g,mm,nn) + cExH(g,mm,nn) * (Hz(g,mm,nn) - Hz(g,mm,nn-1));
			}
			//Ey at nn=0 nn=sizeX-1,H start m-1
			for(mm=1; mm < g->sizeX-1;mm++){
				for(nn=0;nn < g->sizeY-1;nn++)
					Ey(g,mm,nn) = cEyE(g,mm,nn) * Ey(g,mm,nn) - cEyH(g,mm,nn) * (Hz(g,mm,nn) - Hz(g,mm-1,nn));
			}
			break;
		default:
			break;
	}
	return 0;

}
