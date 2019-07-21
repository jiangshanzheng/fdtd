//hard source excitation FDTD
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fdtd.h"
#include"inc.h"
#include"update.h"
#include"GridInit.h"

#define mu0 1
#define epsi0 1
#define DUMP 1
#define output 1
int main(){
	Grid *g;
	ALLOC_1D(g,1,Grid);
	GridInit(g);
	FILE * tmpEz = NULL; 
	FILE * tmpBx = NULL; 
	FILE * tmpBy = NULL;

	FILE * tmpBz = NULL; 
	FILE * tmpEx = NULL; 
	FILE * tmpEy = NULL; 
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
					fwrite(g->hx, sizeof(double), (g->sizeX-1)*(g->sizeY), tmpEx);
					fwrite(g->hy, sizeof(double), (g->sizeX)*(g->sizeY-1), tmpEy); 
					fwrite(g->ez, sizeof(double), (g->sizeX-1)*(g->sizeY-1), tmpBz);
			}
		}
		updateH(g);
		updateE(g);
		printf("steps=%d\n",g->Ti);
		Hz(g,g->sizeX/2,g->sizeY/2) = inc(g,0.0);
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

