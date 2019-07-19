//hard source excitation FDTD
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define SIZE 200
#define STEP 2500
#define dt 1
#define dx 1.0
#define c 1
#define mu0 1
#define epsi0 1
#define im0 1
int main(){
	double * ez, * by;
	//boundary condition EB = 0. PEC PMC
	ez=(double *)calloc(sizeof(double), SIZE); 
	by=(double *)calloc(sizeof(double), SIZE);
	FILE * tmpE = fopen("tmpE","wb");
	FILE * tmpB = fopen("tmpB","wb");
	for(int ti=0;ti<STEP;ti++){
		fwrite(ez, sizeof(double), SIZE, tmpE); fwrite(by, sizeof(double), SIZE, tmpB);
		//push B, this is done on half grid and half time
		//*(by+SIZE-1) = *(by + SIZE-2);
		for(int gi=0;gi<SIZE-1;gi++){
			*(by + gi ) = *(by + gi )  +  ( *(ez + gi + 1 ) - *(ez + gi ) ) * dt/dx / im0;
		}

		*(by+49) -= exp(- (ti - 30.) *  (ti - 30.) /100.) /im0 ;
		//push E, this is done on integer grid and time 
		//*ez = *(ez+1);
		for(int gi=1;gi<SIZE;gi++){
			*(ez + gi ) = *(ez + gi )  + 1/mu0/epsi0*  ( *(by + gi ) - *(by + gi -1 ) ) * dt/dx * im0;
			*(ez+gi) = exp(- (ti - gi - 30.) *  (ti -gi - 30.) /100.) ;
		}
		//set source every time steps(ti - 30) 
		//*(ez+50) += exp(- (ti +0.5 -(-0.5) - 30.) *  (ti +0.5 - (- 0.5) - 30.) /100.) ;
		//*(ez+50) +=  sin(0.1*ti) ;
	}
	fclose(tmpE);
	fclose(tmpB);
	return 0;
}

