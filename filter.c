/*************************************************************************
  > File Name: filter.c
  > Author: 
  > Mail: 
  > Created Time: Sun 18 Aug 2019 12:02:37 AM UTC
 ************************************************************************/

#include<stdio.h>
#include"fdtd.h"
double filter(Grid *g, int mm, int nn);
int filterH(Grid *g){
	return 0;

}

//some components of E can not be calculated, so that's cut off?
int filterE(Grid *g){
	int mm,nn;
	double filter_coef;
	switch (g->type){
		case OneD:
			exit(127);
			break;
		case TMz:
			for(mm=1; mm < g->sizeX-1;mm++){
				for(nn=1;nn < g->sizeY-1;nn++){
					filter_coef = filter(g, mm, nn);
					Ez(g,mm,nn) = filter_coef * Ez(g,mm,nn) ;		}
			}
			break;
		case TEz:
			//Ex at nn=0 nn=sizeX-1, all H start n-1
			for(mm=0; mm < g->sizeX-1;mm++){
				for(nn=1;nn < g->sizeY -1;nn++){
					//filter_coef = filter(g, mm, nn);
					//Ex(g,mm,nn) = filter_coef * Ex(g,mm,nn);
					//				Ex(g,mm,nn) =  Ex(g,mm,nn);
				}
			}
			//Ey at nn=0 nn=sizeX-1,H start m-1
			for(mm=1; mm < g->sizeX-1;mm++){
				for(nn=0;nn < g->sizeY-1;nn++){
					filter_coef = filter(g, mm, nn);
					Ey(g,mm,nn) = filter_coef* Ey(g,mm,nn) ;
				}
			}
			break;
		default:
			break;
	}
	return 0;

}
inline double filter(Grid *g, int mm, int nn){
	double ret;
	if ( (nn <  g->WG_LEN) &&  (nn > 1) ){
		if( mm%g->WG_WIDTH <=2 ){
			ret = 0.;
			//			printf("xm=%d\n",mm);
		}
		else 
			ret = 1.;
	}
	else
		ret = 1.;
	return ret;
}
