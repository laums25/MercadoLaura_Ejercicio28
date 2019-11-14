#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


const double g = 9.8;
const double k = 0.9;
const double delta = 0.01;
const int N = 25;
const double ang = 34.0;
const double v0 = 22.0;
const double v0x = v0* cos(ang*M_PI/180.0);
const double v0y = v0* sin(ang*M_PI/180.0);

void rk(double t, double delta, double p[]);
double pro(double t, double p[], int i);


int main(void)
{
  //ofstream fout("datos-array.dat");
  double time;
  double p[N] = {1, 0};
  for(time = 0; time <= 10; time += delta) {
    cout << time << "\t" << p[0] << "\t" << p[1] << endl;
    rk(time, delta, p);
  }
  //fout.close();

  return 0;
}

double pro(double t, double p[], double i)
{
  if (0 == i) {
    return p[1];
  }
  else if (1 == i) {
    double vt= sqrt(v0x*v0x + v0y*v0y);  
    double dx=- k*v0x*v0x/vt; 
    double dy=- g -k *v0y*v0y / vt;
    return dy;
  }
  else {
    std::cerr << "ERROR. Asking unexistent function f, with i = " <<  i << std::endl;
  }
}


void rk(double t, double delta, double p[]){
    double k1[N], k2[N], k3[N], k4[N];
    double yt[N];
    int i;
    for (i=0; i<N; ++i)
    {
        k1[i] = delta * pro(t,p,i);
    }
    for (i=0; i<N; ++i)
    {
        yt[i] = p[i] + k1[i]/2;
    }
    for (i=0; i<N; ++i)
    {
        k2[i] = delta * pro(t + delta/2, p,i);
    }
    for (i=0; i<N; ++i)
    {
        yt[i] = p[i] + k2[i]/2;
    }
    for (i=0; i<N; ++i)
    {
        k3[i] = delta * pro(t + delta/2, p,i);
    }
    for (i=0; i<N; ++i)
    {
        yt[i] = p[i] + k3[i];
    }  
    for (i=0; i<N; ++i)
    {
        k4[i] = delta * pro(t + delta, p,i);
    }
    for (i=0; i<N; ++i)
    {
        p[i] = p[i] + (1.0/6.0)*(k1[i]+ 2*k2[i] + 2*k3[i] + k4[i]);
    } 
}
