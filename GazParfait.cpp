#include <cstdlib>
#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <ctime>
#include <fstream>

using namespace std;
using namespace sf;


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
    double t;

    CircleShape cercle(r); 
    cercle.setFillColor(Color::Cyan); 
    cercle.setOrigin(r,r); 
    RenderWindow window(VideoMode(L, L), "Gaz Parfait");
    window.setFramerateLimit(50); 
   
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

    for (t=0; t<=tmax; t+=dt) { 
        window.clear(Color::White);
	
        for (j=0; j<N; j++) {
            cercle.setPosition(x[j],y[j]);
            window.draw(cercle); 
        }
        window.display();
	  
        for (j=0; j<=N; j++) {
           
            if (L-r-(x[j]+vx[j]*dt) < 0) {
                dt0 = (abs((x[j]-L))/abs(vx[j]));
                vx[j] = -vx[j];
                x[j] = L + (dt-dt0)*vx[j];    
         }          
            else if ((x[j]+vx[j]*dt)-r < 0) {
                dt0 = x[j]/abs(vx[j]);
                vx[j] = -vx[j];
                x[j] = vx[j]*(dt-dt0);
             
            }
            else if (L-r-(y[j]+vy[j]*dt) < 0) {
                dt0 = (abs((y[j]-L))/abs(vy[j]));
                vy[j] = -vy[j];
                y[j] = L + (dt-dt0)*vy[j];               
            }
            else if ((y[j]+vy[j]*dt)-r < 0) {
                dt0 = y[j]/abs(vy[j]);
                vy[j] = -vy[j];
                y[j] = vy[j]*(dt-dt0);
            }         
            else {      
                x[j] += vx[j]*dt;
                y[j] += vy[j]*dt;
            }      
        } 
    cout << "N*kb/V = " << u << endl;
    cout << "pthéorique/T = " << pressiont/tempe << endl;
    cout << "pexperimental/T = " << pressione/tempe << endl;
    cout << "la température du système est de " << tempe << endl;
    cout << "Le nombre de particules est de "<< N << endl;
    cout << "L'aire du système est de "<< V << endl;
    cout << "La pression théorique est de "<< pressiont << endl;
    cout << "La pression expérimentale est de "<< pressione << endl;
    cout << "energie théorique d'équirépartition = " << kb*tempe << endl;
    cout << "energie expérimentale = " << energie << endl;     
    }

    return 0;
}
