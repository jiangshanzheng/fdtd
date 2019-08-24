/*************************************************************************
  > File Name: gridinit.c
  > Author: 
  > Mail: 
  > Created Time: Sat Jul 20 16:39:37 2019
 ************************************************************************/

#include<stdio.h>
#include<math.h>
#include"fdtd.h"
#include"GridInit.h"
void coef_Ex(Grid *g,int smm, int emm, int snn, int enn, int bond);
void coef_Ey(Grid *g,int smm, int emm, int snn, int enn, int bond);
void coef_Ezx(Grid *g,int smm, int emm, int snn, int enn, int bond);
void coef_Ezy(Grid *g,int smm, int emm, int snn, int enn, int bond);
void coef_Hx(Grid *g,int smm, int emm, int snn, int enn, int bond);
void coef_Hy(Grid *g,int smm, int emm, int snn, int enn, int bond);
void coef_Hzx(Grid *g,int smm, int emm, int snn, int enn, int bond);
void coef_Hzy(Grid *g,int smm, int emm, int snn, int enn, int bond);

double sig_x(Grid * g, int d);
double sig_y(Grid * g, int d);
double sig_mx(Grid * g, int d);
double sig_my(Grid * g, int d);
static double mur=1.0;
static double epsir=1.0;

int GridInit(Grid * g){

	g->type = TMz;
	g->sizeX = 256;
	g->sizeY = 256;
	g->sizeT = 500;
	g->output = 1;
	//per DUMP
	g->DUMP	= 1;
	g->dt  = 1./sqrt(2);
	g->dx = 1.;
	g->c = 1.;
	g->imp0 = 1.;
	g->cdtds = g->c *g->dt/ g->dx;

	//TEz finial points must less then num it E components
	//macro will not check over boundary
	//allocation for boundy fields is large then its real size
	g->PML_Layers = 20;
	g->sig_max = 8;
	g->BLb = 0;
	g->BLt = g->sizeY -1 ;
	g->BLt = 0;

	g->BRb = 0;
	g->BRt = g->sizeY -1;
	g->BRt = 0;

	g->BBl = 0;
	g->BBr = g->sizeX-1;
	g->BBr = 0;

	g->BTl = 0;
	g->BTr = g->sizeX-1;
	g->BTr = 0;
	g->LL = g->BLt - g->BLb;
	g->LR = g->BRt - g->BRb;
	g->LB = g->BBr - g->BBl;
	g->LT = g->BTr - g->BTl;

	int mm,nn;
	switch (g->type){
		case OneD:
			break;
		case TMz:
			ALLOC_2D( g->hx, g->sizeX , g->sizeY -1, double);
			ALLOC_2D( g->hy, g->sizeX -1 , g->sizeY, double);
			ALLOC_2D( g->ez, g->sizeX , g->sizeY, double);
			//coef
			ALLOC_2D( g->ezx, g->sizeX , g->sizeY, double);
			ALLOC_2D( g->ezy, g->sizeX , g->sizeY, double);
			//should check allocation
			
			ALLOC_2D( g->chxe, g->sizeX , g->sizeY -1, double);
			ALLOC_2D( g->chxh, g->sizeX , g->sizeY -1, double);

			ALLOC_2D( g->chye, g->sizeX -1 , g->sizeY, double);
			ALLOC_2D( g->chyh, g->sizeX -1 , g->sizeY, double);

			ALLOC_2D( g->cezxe, g->sizeX , g->sizeY, double);
			ALLOC_2D( g->cezxh, g->sizeX , g->sizeY, double);

			ALLOC_2D( g->cezye, g->sizeX , g->sizeY, double);
			ALLOC_2D( g->cezyh, g->sizeX , g->sizeY, double);
			//initializing 
			for(mm=0; mm < g->sizeX;mm++){
				for(nn=0;nn < g->sizeY -1;nn++){
					cHxE(g,mm,nn) = g->imp0 * g->cdtds;
					cHxH(g,mm,nn) = 1.0;
				}
			}
			for(mm=0; mm < g->sizeX-1;mm++){
				for(nn=0;nn < g->sizeY;nn++){
					cHyE(g,mm,nn) = g->imp0 * g->cdtds;
					cHyH(g,mm,nn) = 1.0;
				}
			}
			for(mm=0; mm < g->sizeX;mm++){
				for(nn=0;nn < g->sizeY;nn++){
					cEzxH(g,mm,nn) = g->imp0 * g->cdtds;
					cEzxE(g,mm,nn) = 1.0;
					cEzyH(g,mm,nn) = g->imp0 * g->cdtds;
					cEzyE(g,mm,nn) = 1.0;
				}
			}

			break;
		case TEz:
			ALLOC_2D( g->ex, g->sizeX-1 , g->sizeY, double);
			ALLOC_2D( g->ey, g->sizeX , g->sizeY-1, double);
			ALLOC_2D( g->hz, g->sizeX-1 , g->sizeY-1, double);
			//coef		
			ALLOC_2D( g->hzx, g->sizeX-1 , g->sizeY-1, double);
			ALLOC_2D( g->hzy, g->sizeX-1 , g->sizeY-1, double);
			//should check allocation
			ALLOC_2D( g->ceye, g->sizeX , g->sizeY -1, double);
			ALLOC_2D( g->ceyh, g->sizeX , g->sizeY -1, double);

			ALLOC_2D( g->cexe, g->sizeX -1 , g->sizeY, double);
			ALLOC_2D( g->cexh, g->sizeX -1 , g->sizeY, double);

			ALLOC_2D( g->chzxe, g->sizeX-1 , g->sizeY-1, double);
			ALLOC_2D( g->chzxh, g->sizeX-1 , g->sizeY-1, double);

			ALLOC_2D( g->chzye, g->sizeX-1 , g->sizeY-1, double);
			ALLOC_2D( g->chzyh, g->sizeX-1 , g->sizeY-1, double);
			for(mm=0; mm < g->sizeX;mm++){
				for(nn=0;nn < g->sizeY -1;nn++){
					cEyH(g,mm,nn) = g->imp0 * g->cdtds;
					cEyE(g,mm,nn) = 1.0;
				}
			}
			for(mm=0; mm < g->sizeX-1;mm++){
				for(nn=0;nn < g->sizeY;nn++){
					cExH(g,mm,nn) = g->imp0 * g->cdtds;
					cExE(g,mm,nn) = 1.0;
				}
			}
			for(mm=0; mm < g->sizeX-1;mm++){
				for(nn=0;nn < g->sizeY-1;nn++){
					cHzxE(g,mm,nn) = g->imp0 * g->cdtds;
					cHzxH(g,mm,nn) = 1.0;
					cHzyE(g,mm,nn) = g->imp0 * g->cdtds;
					cHzyH(g,mm,nn) = 1.0;
				}
			}


			break;
	}

	return 0;
}

int PMLInit(Grid * g){
	//	extern double sig_my,sig_mx,sig_x,sig_y;
	//	extern double mur=1.0;
	//	extern double epsir=1.0;
	int smm,emm,snn,enn;
	int mm,nn;
	switch (g->type){
		case OneD:
			break;
		case TMz:
			//ridges vertical to Y
			//sig_mx=0.; sig_x=0.;
			//sig_my=VAL; sig_y=VAL;

			//Bottom 
			//-------
			//.......
			//.......
			smm=0; emm=g->sizeX;
			snn=0; enn=g->PML_Layers-1;
			coef_Hx(g,smm,emm,snn,enn,enn);

			snn=0;enn=g->PML_Layers;
			coef_Ezy(g,smm,emm,snn,enn,enn);

			//Top
			//......
			//......
			//------
			smm=0; emm=g->sizeX;
			snn=g->sizeY - g-> PML_Layers-1; enn=g->sizeY-1; 
			coef_Hx(g,smm,emm,snn,enn,snn);

			snn = g-> sizeY - g-> PML_Layers; enn=g->sizeY;
			coef_Ezy(g,smm,emm,snn,enn,snn);

			//ridges vertical to X
			//sig_my=0.0; sig_y=0.0;
			//sig_x=VAL; sig_mx=VAL; 

			//Left 
			//..|
			//..|
			//..|
			//..|
			//..|
			//nn set entire
			snn=0; enn=g->sizeY;


			smm=0; emm=g->PML_Layers-1;
			coef_Hy(g,smm,emm,snn,enn,emm);

			smm=0;emm=g->PML_Layers;
			coef_Ezx(g,smm,emm,snn,enn,emm);

			//ridge R	
			//  |..//
			//  |..//
			//  |..//
			//  |..//
			//  |..//
			smm=g->sizeX - g-> PML_Layers -1; emm=g->sizeX-1;
			coef_Hy(g,smm,emm,snn,enn,smm);

			smm=g->sizeX - g-> PML_Layers; emm=g->sizeX;
			coef_Ezx(g,smm,emm,snn,enn,smm);
	
			break;
		case TEz:
			//ridges vertical to Y
			//sig_mx=0.; sig_x=0.;
			//sig_my=VAL; sig_y=VAL;

			//Bottom 
			//-------
			//.......
			//.......
			//mm set entire
			smm=0; emm=g->sizeX;

			snn=0; enn=g->PML_Layers-1;
			coef_Ex(g,smm,emm,snn,enn,enn);

			snn=0;enn=g->PML_Layers-1;
			coef_Hzy(g,smm,emm,snn,enn,enn);
			
			//Top
			//......
			//......
			//------
			snn=g->sizeY - g-> PML_Layers-1; enn=g->sizeY-1;
			coef_Ex(g,smm,emm,snn,enn,snn);

			snn=g->sizeY - g-> PML_Layers-1; enn=g->sizeY-1;
			coef_Hzy(g,smm,emm,snn,enn,snn);

			//ridges vertical to X
			//sig_my=0.0; sig_y=0.0;
			//sig_x=VAL; sig_mx=VAL; 

			//Left 
			//..|
			//..|
			//..|
			//..|
			//..|
			//nn set entire
			snn=0; enn=g->sizeY;
			smm=0; emm=g->PML_Layers-1;
			coef_Ey(g,smm,emm,snn,enn,emm);

			snn=0; enn=g->sizeY -1;
			smm=0; emm=g->PML_Layers-1;
			coef_Hzx(g,smm,emm,snn,enn,emm);
			

			//ridge R	
			//  |..//
			//  |..//
			//  |..//
			//  |..//
			//  |..//

			snn=0; enn=g->sizeY;
			smm=g->sizeX - g-> PML_Layers-1; emm=g->sizeX-1;
			coef_Ey(g,smm,emm,snn,enn,smm);

			snn=0; enn=g->sizeY - 1;
			smm=g->sizeX - g-> PML_Layers-1; emm=g->sizeX-1;
			coef_Hzx(g,smm,emm,snn,enn,smm);
			break;
	}
	return 0;
}
//TMz sig_m x Hy
inline void coef_Hy(Grid *g,int smm, int emm, int snn, int enn, int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cHyE(g,mm,nn) = g->imp0 * g->cdtds /mur / (1.0 + sig_mx(g,mm-bond) * g->dt / 2/ mur);  
			cHyH(g,mm,nn) = (1. - sig_mx(g,mm-bond) * g->dt/2/mur)/(1. + sig_mx(g,mm-bond) * g->dt/2/mur);      
		}
	}
}

//TMz sig_m y Hx
inline void coef_Hx(Grid *g,int smm, int emm, int snn, int enn, int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cHxE(g,mm,nn) = g->imp0 * g->cdtds /mur / (1.0 + sig_my(g,nn-bond) * g->dt / 2/ mur)  ;  
			cHxH(g,mm,nn) = (1. - sig_my(g,nn-bond) * g->dt/2/mur)/(1. + sig_my(g,nn-bond) * g->dt/2/mur);
		}
	}
}
//TMz sig x Ezx
inline void coef_Ezx(Grid *g,int smm, int emm, int snn, int enn, int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cEzxH(g,mm,nn) = g->imp0 * g->cdtds /epsir / (1.0 + sig_x(g,mm-bond) * g->dt / 2/ epsir); 
			cEzxE(g,mm,nn) = (1. - sig_x(g,mm-bond) * g->dt/2/epsir)/(1. + sig_x(g,mm-bond) * g->dt/2/epsir);   
		}
	}
}
//TMz sig y Ezy
inline void coef_Ezy(Grid *g,int smm, int emm, int snn, int enn, int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cEzyH(g,mm,nn) = g->imp0 * g->cdtds /epsir / (1.0 + sig_y(g,nn-bond) * g->dt / 2/ epsir); 
			cEzyE(g,mm,nn) = (1. - sig_y(g,nn-bond) * g->dt/2/epsir)/(1. + sig_y(g,nn-bond) * g->dt/2/epsir);   

		}
	}
}

//TEz
inline void coef_Ey(Grid *g,int smm, int emm, int snn, int enn, int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cEyH(g,mm,nn) = g->imp0 * g->cdtds /epsir / (1.0 + sig_x(g,mm-bond) * g->dt / 2/ epsir);  
			cEyE(g,mm,nn) = (1. - sig_x(g,mm-bond) * g->dt/2/epsir)/(1. + sig_x(g,mm-bond) * g->dt/2/epsir);      
		}
	}
}
inline void coef_Ex(Grid *g,int smm, int emm, int snn, int enn, int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cExH(g,mm,nn) = g->imp0 * g->cdtds /epsir / (1.0 + sig_y(g,nn-bond) * g->dt / 2/ epsir);  
			cExE(g,mm,nn) = (1. - sig_y(g,nn-bond) * g->dt/2/epsir)/(1. + sig_y(g,nn-bond) * g->dt/2/epsir);      
		}
	}
}
inline void coef_Hzx(Grid *g,int smm, int emm, int snn, int enn, int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cHzxE(g,mm,nn) = g->imp0 * g->cdtds /mur / (1.0 + sig_mx(g,mm-bond) * g->dt / 2/ mur); 
			cHzxH(g,mm,nn) = (1. - sig_mx(g,mm-bond) * g->dt/2/mur)/(1. + sig_mx(g,mm-bond) * g->dt/2/mur);   
		}
	}
}

inline void coef_Hzy(Grid *g,int smm, int emm, int snn, int enn,int bond){
	int mm,nn;
	for(mm=smm; mm < emm; mm++){
		for(nn= snn ;nn<enn ; nn++){
			cHzyE(g,mm,nn) = g->imp0 * g->cdtds /mur / (1.0 + sig_my(g,nn-bond) * g->dt / 2/ mur); 
			cHzyH(g,mm,nn) = (1. - sig_my(g,nn-bond) * g->dt/2/mur)/(1. + sig_my(g,nn-bond) * g->dt/2/mur);   
		}
	}
}
inline double sig_x(Grid * g, int d){
	return g->sig_max * fabs(d)/(double)g->PML_Layers;
	//return 0.;
}
inline double sig_mx(Grid * g, int d){
	return g->sig_max * fabs(d)/(double)g->PML_Layers;
	//return 0.;
}
inline double sig_y(Grid * g, int d){
	return g->sig_max * fabs(d)/(double)g->PML_Layers;
	//return 0.;
}
inline double sig_my(Grid * g, int d){
	return g->sig_max * fabs(d)/(double)g->PML_Layers;
	//return 0.;
}
//mur epsi visible
