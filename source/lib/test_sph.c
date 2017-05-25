#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "math2.h"

double Work[9],Work2[9];

void Sph_Cart2Sc(double *Cart, double *Sc){
	Work[0]=sqrt(Cart[0]*Cart[0]+Cart[1]*Cart[1]+Cart[2]*Cart[2]);
	Work[1]=Work[0]>0?acos(Cart[2]/Work[0]):0;
	Work[2]=atan2(Cart[1],Cart[0]);
	if(Work[2]<0) Work[2]+=2.0*PI;
	Sc[0]=Work[0];
	Sc[1]=Work[1];
	Sc[2]=Work[2];
	return;
}
void main(){
	double pt0[3]={0,0,0};
	double pt1[3]={250,250,5000};
	double pt2[3]={-250,250,5000};
	double pt3[3]={-250,-250,5000};
	double pt4[4]={200,200,4000};
	double Sc1[3];
	double Sc2[3];
	double Sc3[3];
	double Sc4[3];
		
	Sph_Cart2Sc(pt1,Sc1);
	Sph_Cart2Sc(pt2,Sc2);
	Sph_Cart2Sc(pt3,Sc3);
	Sph_Cart2Sc(pt4,Sc4);

	printf("Sc1[0]=%f, Sc1[1]=%f, Sc1[2]=%f\n", Sc1[0], Sc1[1], Sc1[2]);
	printf("Sc2[0]=%f, Sc2[1]=%f, Sc2[2]=%f\n", Sc2[0], Sc2[1], Sc2[2]);
	printf("Sc3[0]=%f, Sc3[1]=%f, Sc3[2]=%f\n", Sc3[0], Sc3[1], Sc3[2]);
	printf("Sc4[0]=%f, Sc4[1]=%f, Sc4[2]=%f\n", Sc4[0], Sc4[1], Sc4[2]);

	return;
}

