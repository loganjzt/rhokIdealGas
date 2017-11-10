
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <string>
#include <math.h>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]){

	char f[40];
    FILE *data;
    int i,j;

	double dr = 0.1;
	int natom = 100;
	double r,x;

	double p[int(natom/dr)];
	
	for(i=0;i<natom/dr;i++) p[i] = 0.0;

    for(i=0;i<natom/dr;i++){
		p[i] = 0.0;
		r = dr*i;
		for(j = 0 ; j < natom*100 ; j++ ){
			x = j*dr;
			p[i] += ( cyl_bessel_j(1,r*x) * pow( cyl_bessel_j(0,x),natom) + cyl_bessel_j(1,r*(x+dr)) * pow( cyl_bessel_j(0,x+dr) ,natom ) )*dr/2.0;
		}
		p[i] = p[i] * r;	
	}

	sprintf(f,"p_A_n%d.dat",natom);
	data = fopen(f,"w");

	for(i=0;i<natom/dr;i++){
		fprintf(data,"%e\t%e\t%e\n",i*dr,p[i],-log(p[i]));
	}

    return 0;
}
