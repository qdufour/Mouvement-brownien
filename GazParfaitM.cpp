#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

using namespace std;


int coord(int L, int r) {

    double x = rand() % (L - 2*r) + r;
    return x;

}


int place(double a[], double b[], double x, double y, int r, int n) {

    int o = 0;
    int j;
    for (j=0; j<n; j++){
        if (sqrt((x-a[j])*(x-a[j]) + (y-b[j])*(y-b[j])) < r*2){
            o = 1;
            break;
        }
    }
    return o;

}


int adsorption (double a[], double b[], int L, int r, int N, int n) {

    int o = 0;
    double c = coord(L,r);
    double d = coord(L,r);
   
    int p = place(a, b, c, d, r, n);
    if (p == 1){
        o = 1;
    }
    else {
    a[n] = c;
    b[n] = d;
   
    }
    return o;

}  


int main() {

    int r = 10;
    int L = 800;
    double m = 5;
    srand(time(0));
   
    const int N = 100;
    double x[N];
    double y[N];
    double vx[N];
    double vy[N];
    double dt =1;
    double dt0 = 0;
    int tmax = 100;

    int n = 0;
    int j;
    int k;
    double t;
    int collision=0;
    double vquad = 0;
    ofstream f("fichier.res");

   
    while (n <= N) {
        int q = adsorption(x, y, L, r, N, n);
            if (q == 1) {
                continue;
            }
            else {
                vx[n] = (rand() % 10) - 5;
                vy[n] = (rand() % 10) - 5;              
                n++;
            }
    }
    for (k=0; k<100; k++) {
    L = k*L + L;
    for (t=0; t<=tmax; t+=dt) { 
	vquad = 0;   
        for (j=0; j<=N; j++) {
           
            if (L-r-(x[j]+vx[j]*dt) < 0) {
                collision++;
                dt0 = (abs((x[j]-L))/abs(vx[j]));
                vx[j] = -vx[j];
                x[j] = L + (dt-dt0)*vx[j];    
         }          
            else if ((x[j]+vx[j]*dt)-r < 0) {
                collision++;
                dt0 = x[j]/abs(vx[j]);
                vx[j] = -vx[j];
                x[j] = vx[j]*(dt-dt0);
             
            }
            else if (L-r-(y[j]+vy[j]*dt) < 0) {
                collision++;
                dt0 = (abs((y[j]-L))/abs(vy[j]));
                vy[j] = -vy[j];
                y[j] = L + (dt-dt0)*vy[j];               
            }
            else if ((y[j]+vy[j]*dt)-r < 0) {
                collision++;
                dt0 = y[j]/abs(vy[j]);
                vy[j] = -vy[j];
                y[j] = vy[j]*(dt-dt0);
            }         
            else {      
                x[j] += vx[j]*dt;
                y[j] += vy[j]*dt;
            }      
        vquad +=  (vx[j]*vx[j]+vy[j]*vy[j])*(vx[j]*vx[j]+vy[j]*vy[j])/N;
        }      
    }

    double V = L*L*10e-18;
    double kb = 1.38e-23;
    double collisiont = collision/(tmax);
    double vitesse_quadratique_moyenne = sqrt(vquad)*10e-10;

    double tempe = (m*vitesse_quadratique_moyenne*vitesse_quadratique_moyenne)/(3*kb);
    double pressiont = (N*m*vitesse_quadratique_moyenne*vitesse_quadratique_moyenne)/(3*V);
    double pressione = (collisiont*2*m*vitesse_quadratique_moyenne*10e9)/(4*L);
    double u = (N*kb)/(V);
    double energie = 0.5*m*vitesse_quadratique_moyenne*vitesse_quadratique_moyenne;


 
    f << V << " " << pressiont/tempe << " " << pressione/tempe << " " << kb*tempe << " " << energie << endl;

    }
    f.close();
    return 0;
}
