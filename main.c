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
//	BoundaryInit(g);
	FILE * tmpEz = NULL; 
	FILE * tmpBx = NULL; 
	FILE * tmpBy = NULL;

	FILE * tmpBz = NULL; 
	FILE * tmpEx = NULL; 
	FILE * tmpEy = NULL; 
	printf("%d",g->type);
	switch (g->type){
		case TMz:
			tmpEz = fopen("tmpEz.dat","wb");
			tmpBx = fopen("tmpBx.dat","wb");
			tmpBy = fopen("tmpBy.dat","wb");
		case TEz:
			tmpBz = fopen("tmpBz.dat","wb");
			tmpEx = fopen("tmpEx.dat","wb");
			tmpEy = fopen("tmpEy.dat","wb");
	}
	for(g->Ti=0; g->Ti < g->sizeT; g->Ti ++){
		if (output &&  0 == g->Ti % DUMP){
			switch (g->type){
				case TMz:
					fwrite(g->hx, sizeof(double), (g->sizeX)*(g->sizeY-1), tmpBx);
					fwrite(g->hy, sizeof(double), (g->sizeX-1)*(g->sizeY), tmpBy); 
					fwrite(g->ez, sizeof(double), (g->sizeX)*(g->sizeY), tmpEz); 
				case TEz:
					fwrite(g->ex, sizeof(double), (g->sizeX-1)*(g->sizeY), tmpEx);
					fwrite(g->ey, sizeof(double), (g->sizeX)*(g->sizeY-1), tmpEy); 
					fwrite(g->hz, sizeof(double), (g->sizeX-1)*(g->sizeY-1), tmpBz);
			}
		}
		updateH(g);
		updateE(g);
	//	abc(g);
	//	printf("steps=%d\n",g->Ti);
		int loc;
		double phs;
		for(loc=20;loc<80;loc++){
			phs=loc/60.*2*M_PI;
			Ez(g,loc,g->sizeY/2) = sin(0.1*g->Ti)*cos(phs);
		}
	}
	switch (g->type){
		case TMz:
			fclose(tmpEz);
			fclose(tmpBy);
			fclose(tmpBx);
		case TEz:
			fclose(tmpBz);
			fclose(tmpEy);
			fclose(tmpEx);
	}
	return 0;
}

