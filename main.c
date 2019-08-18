//hard source excitation FDTD
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"fdtd.h"
#include"inc.h"
#include"inc.h"
#include"filter.h"
#include"update.h"
#include"update_pml.h"
#include"GridInit.h"
#include"abc.h"

#define mu0 1
#define epsi0 1
int main(){
	Grid *g;
	ALLOC_1D(g,1,Grid);
	GridInit(g);
	//BoundaryInit(g);
	PMLInit(g);
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
	clock_t begin, end;
	for(g->Ti=0; g->Ti < g->sizeT; g->Ti ++){

		//	storge
		if (g->output &&  0 == g->Ti % g->DUMP){
			begin = clock();
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
			end=clock();
			double Tstorage = (double)(end - begin) / CLOCKS_PER_SEC;
			printf("outputing...,storge time = %f\n",Tstorage);
		}

		begin = clock();

		//printf("%f\n",Hz(g,g->sizeX/2,g->sizeY/2));
		//updateH_pml(g);
		updateH(g);
		//Hz(g,g->sizeX/2,g->sizeY/2) += 10.0*sin(0.1*g->Ti);

		//updateE_pml(g);
		updateE(g);
		abc(g);
		filterE(g);
		int width = g->WG_WIDTH;
		int mm,nn;		
		double omega = 0.11;
		double amp = 0.55;
		int LOC = 256;
		double  coef =  (  1 -  exp (   (  -1 *         (  g->Ti * omega * 1. ) ) ) ) ;
		for(mm=0;mm<g->sizeX;mm++){
			double  kperp =   (       (  1 *         floor (         (  mm / width ) ) ) / (g->sizeX/width) ) ;
			double phs=(  g->Ti * omega -     (  2 *  (  M_PI * kperp ) ) ) ;
			//printf("%f\n",phs);			
					if( ((int)mm) % ((int)width) ==(int)(width/2) )
					Ex(g,mm,LOC) += amp *  sin ( phs ) *  coef;
			//printf("%f\n",Ez(g,mm,80));			
		}
		//Ez(g,g->sizeX/2,g->sizeY/2) += sin(0.1*g->Ti);

		end=clock();
		double Tfield = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("steps=%d,time=%f\n",g->Ti,Tfield);

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

