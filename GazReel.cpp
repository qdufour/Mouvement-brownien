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
 
    int L;
    cout << " Longueur du côté de la fenêtre (inférieur à 1000) =  ";
    cin >> L;

    int Nmax = (40*L*L)/(100*M_PI*400);
    int N;
    cout << " Nombre de particules (inférieur à " << Nmax << ") = " ;
    cin >> N;
    while ( N > Nmax) {
    cout << " Nombre de particules (inférieur à " << Nmax << ") = " ;
    cin >> N;
    }
    
    int nbespece;
    cout << " Nombre d'espèces de particules (inférieur à " << N << ") = " ;
    cin >> nbespece;
    while ( nbespece > N) {
    cout << " Nombre d'espèces de particules (inférieur à " << N << ") = " ;
    cin >> nbespece;
    }

    int k = 0;
    int particule;
    int total = 0;
    int nbparticules[nbespece];
    while ( k < nbespece) {
    cout << " Nombre de particules de l'espèce " << k+1 << " = ";
    cin >> particule;
    nbparticules[k] = particule;
    k++;
    total += particule;
    cout << " " << total << " particules au total." << endl;
    if ( k == nbespece && total != N) {
    cout << " Nombre de particules erroné." << endl;
    k = 0;
    total = 0;
    }
    }

    double masse[nbespece];
    double mmax = 1;
    double masseparticule;
    k = 0;
    while ( k < nbespece) {
    cout << " Masse des particules de l'espèce " << k+1 << " (inférieure à " << mmax << ") = ";
    cin >> masseparticule;
    if ( masseparticule > mmax) {
    continue;
    }
    else {
    	masse[k] = masseparticule;
    	k++;
    }
    }
    
    double rayon[nbespece];
    double rmax = 20;
    double rayonparticule;
    k = 0;
    while ( k < nbespece) {
    cout << " Rayon des particules de l'espèce " << k+1 << " (inférieure à " << rmax << ") = ";
    cin >> rayonparticule;
    if ( rayonparticule > rmax) {
    continue;
    }
    else {
	rayon[k] = rayonparticule;
        k++;
    }
    }

    double m[N];
    double r[N];
    
    srand(time(0));
       
    double x[N];
    double y[N];
    double vx[N];
    double vy[N];
    double dt = 1.0;
    double t;
    double dt0;
    double tmax = 200;

    int n=0;
    int j;
    int i;

    double theta;
    double distance;
    double tcoll[N];
    double t0;
    double x1;
    double y1;
    double x2;
    double y2;

    CircleShape cercle(r); 
    cercle.setFillColor(Color::Cyan); 
    cercle.setOrigin(r,r); 
    RenderWindow window(VideoMode(L, L), "Gaz Réel");
    window.setFramerateLimit(50); 

    while (n < N) { 
        int q = adsorption(x, y, L, r[n], N, n);
            if (q == 1) {
                continue;
            }  
	    else {
		vx[n] = (rand() % 10) - 5;
                vy[n] = (rand() % 10) - 5;
                n++; 
	    }  
    }
    
    int npartiel = 0;
    for (j=0; j<N; j++) {
	for (k=0; k<nbespece; k++) {
		npartiel += nbparticules[k];
		if (j<npartiel) {
			m[j] = masse[k];
			r[j] = rayon[k];
			npartiel = 0;
			break;
		}
		else {
			continue;
		}
	}
    }

    for (t=0; t<=tmax; t+=dt) {
        window.clear(Color::White);
	
        for (j=0; j<N; j++) {
            cercle.setPosition(x[j],y[j]);
            window.draw(cercle); 
	    
            tcoll[j] = 0.0;
        }
        window.display();
        
	for (j=0; j<N; j++) {

		for (i=0; i<N; i++) {
			if (tcoll[j] != 0) {  
                		break;
                	}
			else if (i == j or tcoll[i] != 0) {
				continue;
			}

                	distance = sqrt(((x[j]+vx[j]*dt)-(x[i]+vx[i]*dt))*((x[j]+vx[j]*dt)-(x[i]+vx[i]*dt)) + ((y[j]+vy[j]*dt)-(y[i]+vy[i]*dt))*((y[j]+vy[j]*dt)-(y[i]+vy[i]*dt)));

                	if (distance < r[j]+r[i]) {
		        	for (t0=0;t0<=dt;t0+=dt/100) {
		    			if (sqrt(((x[j]+vx[j]*t0)-(x[i]+vx[i]*t0))*((x[j]+vx[j]*t0)-(x[i]+vx[i]*t0)) + ((y[j]+vy[j]*t0)-(y[i]+vy[i]*t0))*((y[j]+vy[j]*t0)-(y[i]+vy[i]*t0))) < (r[j]+r[i])+0.1) {
						tcoll[j] = t0; // Détermination des positions et du t lors de la collision 
						tcoll[i] = t0;
						x1 = x[j]+(vx[j]*t0);
						y1 = y[j]+(vy[j]*t0);
						x2 = x[i]+(vx[i]*t0);
						y2 = y[i]+(vy[i]*t0);
						break;
					
					}
		 	         }
				 
				 double vx1 = vx[j];
				 double vy1 = vy[j];
				 double vx2 = vx[i];
				 double vy2 = vy[i];

				 double un1 = x2-x1; // Définition du vecteur unitaire normal à la collision
				 double un2 = y2-y1;
				 double un = sqrt(un1*un1+un2*un2);
				 un1 = un1/n;
				 un2 = un2/n;

				 double ut1 = -un2; // Définition du vecteur unitaire tangentielle à la collision
				 double ut2 = un1;
			         
			         if( (x2-x1) > 0 && (y2-y1) > 0) { // Définition de l'angle entre les vecteurs vitesses
				 	theta = acos((x2-x1)/(r[j]+r[i]));  
			         }
			         if( (x2-x1) < 0 && (y2-y1) < 0) {
				 	theta = acos((x1-x2)/(r[j]+r[i]));
			         }
			         if( (x2-x1) > 0 && (y2-y1) < 0) {
				 	theta = acos((x1-x2)/(r[j]+r[i]));
			         }
			         if( (x2-x1) < 0 && (y2-y1) > 0) {
				 	theta = acos((x2-x1)/(r[j]+r[i]));
			         }

			         double vn1 = cos(theta)*(vx2)+sin(theta)*(vy2); // Projection des vecteurs vitesses sur la base (un,ut)
			         double vt1 = cos(theta)*(vy1)-sin(theta)*(vx1);
			         double vn2 = cos(theta)*(vx1)+sin(theta)*(vy1);
			         double vt2 = cos(theta)*(vy2)-sin(theta)*(vx2);
			         
			         double vt10 = vt1; // Vitesses tangentielles après impact
			         double vt20 = vt2;

				 double vn10 = (vn1*(m[j]-m[i])+2*m[i]*vn2)/(m[i]+m[j]); // Vitesses normales après impact
				 double vn20 = (vn2*(m[i]-m[j])+2*m[j]*vn1)/(m[i]+m[j]);

				 vx[j] = cos(theta)*vn1 - sin(theta)*vt1; // Composantes des vitesses après impact dans la base (x,y)
				 vy[j] = sin(theta)*vn1 + cos(theta)*vt1;
				 vx[i] = cos(theta)*vn2 - sin(theta)*vt2;
				 vy[i] = sin(theta)*vn2 + cos(theta)*vt2;  

				 x[j] = x1 + vx[j]*(dt-t0); // Attribution des nouvelles positions à la fin de l'intervalle dt
				 y[j] = y1 + vy[j]*(dt-t0);
				 x[i] = x2 + vx[i]*(dt-t0);
				 y[i] = y2 + vy[i]*(dt-t0);
		                 break;
		        }
		    }

		    if (tcoll[j] != 0) {
			continue;
		    } 
		    else if (L-r[j]-(x[j]+vx[j]*dt) < 0) { 
		        
		        dt0 = ((x[j]-L)/vx[j]);
		        vx[j] = -vx[j]; 
		        x[j] = L + (dt-dt0)*vx[j];
		    }
		    else if ((x[j]+vx[j]*dt)-r[j] < 0) {
		        
		        dt0 = -x[j]/vx[j];
		        vx[j] = -vx[j]; 
		        x[j] = vx[j]*(dt-dt0);
		    }
		    else if (L-r[j]-(y[j]+vy[j]*dt) < 0) {
		        
		        dt0 = ((y[j]-L)/vy[j]);
		        vy[j] = -vy[j]; 
		        y[j] = L + (dt-dt0)*vy[j];
		    }
		    else if ((y[j]+vy[j]*dt)-r[j] < 0) {
		        
		        dt0 = -y[j]/vy[j];
		        vy[j] = -vy[j]; 
		        y[j] = vy[j]*(dt-dt0);
		    }
		    else {      
		        x[j] += vx[j]*dt;
		        y[j] += vy[j]*dt;
		    }
		}
	    }

    return 0;
}
