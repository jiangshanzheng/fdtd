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
	coef = M_PI *  ( (g->cdtd * g->Ti - location)/ppw - 1.);
	coef=coef*coef;
	return ( 1. - 2.* coef )* exp(-coef);
}
