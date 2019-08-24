/*************************************************************************
  > File Name: inc.c
  > Author: 
  > Mail: 
  > Created Time: Sat Jul 20 18:55:25 2019
 ************************************************************************/

#include"fdtd.h"
#include<math.h>
#include"inc.h"
#define ppw 20
double inc(Grid *g, double location){
	//double phs;
	double coef;
	coef = M_PI *  ( (1.1 * g->Ti - location)/ppw - 1.);
	coef=coef*coef;
	return ( 1. - 2.* coef )* exp(-coef);
}

double point(Grid *g,double omega, double amp, double phs, int LOCX,int LOCY){
	int mm,nn;		
	int LOC = 128;
	double  coef =  (  1 -  exp (   (  -1 *         ( g->Ti*omega * 1. ) ) ) ) ;
	return amp *  sin ( g->Ti*omega + phs ) *  coef;
}

