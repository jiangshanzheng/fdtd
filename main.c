//hard source excitation FDTD
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fdtd.h"
#include"inc.h"
#include"update.h"
#include"GridInit.h"
#include"abc.h"

#define mu0 1
#define epsi0 1
#define DUMP 1
#define output 1
int main(){
	Grid *g;
	ALLOC_1D(g,1,Grid);
	GridInit(g);
	BoundaryInit(g);
	FILE * tmpEz = NULL; 
	FILE * tmpBx = NULL; 
	FILE * tmpBy = NULL;

	FILE * tmpBz = NULL; 
	FILE * tmpEx = NULL; 
	FILE * tmpEy = NULL; 
	switch (g->type){
		case OneD:
			exit(127);
			break;
		case TMz:
			tmpEz = fopen("tmpEz.dat","wb");
			tmpBx = fopen("tmpBx.dat","wb");
			tmpBy = fopen("tmpBy.dat","wb");
			break;
		case TEz:
			tmpBz = fopen("tmpBz.dat","wb");
			tmpEx = fopen("tmpEx.dat","wb");
			tmpEy = fopen("tmpEy.dat","wb");
			break;
	}
	for(g->Ti=0; g->Ti < g->sizeT; g->Ti ++){
		if (output &&  0 == g->Ti % DUMP){
			switch (g->type){
				case OneD:
					exit(127);
					break;
				case TMz:
					fwrite(g->hx, sizeof(double), (g->sizeX)*(g->sizeY-1), tmpBx);
					fwrite(g->hy, sizeof(double), (g->sizeX-1)*(g->sizeY), tmpBy); 
					fwrite(g->ez, sizeof(double), (g->sizeX)*(g->sizeY), tmpEz); 
					break;
				case TEz:
					fwrite(g->ex, sizeof(double), (g->sizeX-1)*(g->sizeY), tmpEx);
					fwrite(g->ey, sizeof(double), (g->sizeX)*(g->sizeY-1), tmpEy); 
					fwrite(g->hz, sizeof(double), (g->sizeX-1)*(g->sizeY-1), tmpBz);
					break;
			}
		}
		updateH(g);
		updateE(g);
		int mm = g->sizeX/2;
		int nn;
		for(nn=0;nn< g->sizeY-1;nn++){
			Ey(g,mm,nn) = 0;
		}
		abc(g);
		printf("steps=%d\n",g->Ti);
		int x0=2.;		
		int x1=13.;
		int num_wg = 2;
		int width_wg = 16;
		int SLOC=100;
		int loc;
		double phs;
		int i_wg;
		for(i_wg=0;i_wg<num_wg;i_wg++){
			int offset = i_wg * width_wg;
			for(loc=x0+offset;loc<x1+offset;loc++){
			//	printf("seting");
				phs= num_wg / (loc)/(width_wg) *2*M_PI;
				Ey(g,loc,SLOC) += sin(1.*g->Ti - phs);
			}
		}
//		Ey(g,50,50) = sin(0.1*g->Ti);
//		Ey(g,g->sizeX/2,g->sizeY/2) = inc(g,0.0);
	}
		switch (g->type){	
			case OneD:
				exit(127);
				break;
			case TMz:
				fclose(tmpEz);
				fclose(tmpBy);
				fclose(tmpBx);
				break;
			case TEz:
				fclose(tmpBz);
				fclose(tmpEy);
				fclose(tmpEx);
				break;
		}
		return 0;
	}

