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
	FILE * tmpEz = fopen("tmpEz.dat","wb");
	FILE * tmpBx = fopen("tmpBx.dat","wb");
	FILE * tmpBy = fopen("tmpBy.dat","wb");
	for(g->Ti=0; g->Ti < g->sizeT; g->Ti ++){
		if (output &&  0 == g->Ti % DUMP){
		fwrite(g->hx, sizeof(double), (g->sizeX)*(g->sizeY-1), tmpBx);
		fwrite(g->hy, sizeof(double), (g->sizeX-1)*(g->sizeY), tmpBy); 
		fwrite(g->ez, sizeof(double), (g->sizeX)*(g->sizeY), tmpEz); 
		}
//	printf("fuck");
	updateH(g);
	updateE(g);
	printf("steps=%d\n",g->Ti);
	Ez(g,g->sizeX/2,g->sizeY/2) = inc(g,0.0);
	}
	fclose(tmpEz);
	fclose(tmpBy);
	fclose(tmpBx);
	return 0;
}

