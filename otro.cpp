#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double g = 9.8;
double k = 0.9;
double delta = 0.01;
int N = 25;
double ang = 34.0;
double v0 = 22.0;
double v0x = v0* cos(ang*M_PI/180.0);
double v0y = v0* sin(ang*M_PI/180.0);

double f0(double t, double y0, double y1); 
double f1(double t, double y0, double y1); 
void rk4(double t, double h, double & y0, double & y1); 

int main(void)
{
  ofstream fout("datos-array.txt");
  double x = 1, v = 0;
  double time;
  double y[N] = {1, 0};
  for(time = 0; time <= 10; time += DeltaT) {
    cout << time << "\t" << x << "\t" << v << endl;
    rk4(time, DeltaT, x, v);
    fout << time << "\t" << y[0] << "\t" << y[1] << endl;
    rk4_vector(time, DeltaT, y);
  }
  fout.close();

  return 0;
}


double dx(double t, double y0, double y1)
{
  return y1;
}

double dx1(double t, double y0, double y1)
{
  return - k*v0x*v0x/vt; 

}

double dy(double t, double T0, double T1)
{
  return y1;
}

double dy1(double t, double T0, double T1)
{
  return - g -k *T0*T0 / vt;
}





void rk4(double t, double h, double & y0, double & y1, double & T0, double & T1) // metodo de runge kutta 4 orden
{
  double k10, k11, k20, k21, k30, k31, k40, k41, K1, k1y, k2, k2y, k3, k3y, k4, k4y;
  k10 = h*dx(t, y0, y1);
  k11 = h*dx1(t, y0, y1);
  k20 = h*dx(t+h/2, y0 + k10/2, y1 + k11/2);
  k21 = h*dx1(t+h/2, y0 + k10/2, y1 + k11/2);
  k30 = h*dx(t+h/2, y0 + k20/2, y1 + k21/2);
  k31 = h*dx1(t+h/2, y0 + k20/2, y1 + k21/2);
  k40 = h*dx(t + h, y0 + k30, y1 + k31);
  k41 = h*dx1(t + h, y0 + k30, y1 + k31);
  
  k1 = h*dy(t, T0, T1);
  k1y = h*dy1(t, T0, T1);
  k2 = h*dy(t+h/2, T0 + k1/2, T1 + k1y/2);
  k2y = h*dy1(t+h/2, T0 + k1/2, T1 + k1y/2);
  k3 = h*dy(t+h/2, T0 + k2/2, T1 + k2y/2);
  k3y = h*dy1(t+h/2, T0 + k2/2, T1 + k2y/2);
  k4 = h*dy(t + h, T0 + k3, T1 + k3y);
  k4y = h*dy1(t + h, T0 + k3, T1 + k3y);

  

  y0 = y0 + (1.0/6.0)*(k10 + 2*k20 + 2*k30 + k40);
  y1 = y1 + (1.0/6.0)*(k11 + 2*k21 + 2*k31 + k41);
}