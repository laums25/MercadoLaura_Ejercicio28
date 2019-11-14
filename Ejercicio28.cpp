#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double g = 9.8;
double k = 0.9;
double delta = 0.01;
double ang = 34.0;
double v0 = 22.0;
double v0x = v0* cos(ang*M_PI/180.0);
double v0y = v0* sin(ang*M_PI/180.0);
double vt=sqrt(v0x*v0x+v0y*v0y);

double dx(double t, double y0, double y1); 
double dx1(double t, double y0, double y1);
double dy(double t, double T0, double T1); 
double dy1(double t, double T0, double T1); 

void rk4(double t, double h, double & y0, double & y1, double & T0, double & T1); 

int main(void)
{
  ofstream fout("Ejercicio28.dat");
  double y1 = 0, vex = v0x, vey = v0y, T1=0;
  double time;
  for(time = 0; time <= 1.9; time += delta) {
    fout << time << "\t" << y1 << "\t" << vex << "\t" << T1 << "\t" << vey << "\t" <<endl;
    rk4(time, delta, y1, vex, T1, vey);
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
  return T1;
}

double dy1(double t, double T0, double T1)
{
  return - g -k *v0y*v0y / vt;
}


void rk4(double t, double h, double & y0, double & y1, double & T0, double & T1) // metodo de runge kutta 4 orden
{
    
  double k10, k11, k20, k21, k30, k31, k40, k41, k10y, k1y, k20y, k2y, k30y, k3y, k40y, k4y;
    
  k10 = h*dx(t, y0, y1);
  k11 = h*dx1(t, y0, y1);
  k20 = h*dx(t+h/2, y0 + k10/2, y1 + k11/2);
  k21 = h*dx1(t+h/2, y0 + k10/2, y1 + k11/2);
  k30 = h*dx(t+h/2, y0 + k20/2, y1 + k21/2);
  k31 = h*dx1(t+h/2, y0 + k20/2, y1 + k21/2);
  k40 = h*dx(t + h, y0 + k30, y1 + k31);
  k41 = h*dx1(t + h, y0 + k30, y1 + k31);
  
  k10y = h*dy(t, T0, T1);
  k1y = h*dy1(t, T0, T1);
  k20y = h*dy(t+h/2, T0 + k10y/2, T1 + k1y/2);
  k2y = h*dy1(t+h/2, T0 + k10y/2, T1 + k1y/2);
  k30y = h*dy(t+h/2, T0 + k20y/2, T1 + k2y/2);
  k3y = h*dy1(t+h/2, T0 + k20y/2, T1 + k2y/2);
  k40y = h*dy(t + h, T0 + k30y, T1 + k3y);
  k4y = h*dy1(t + h, T0 + k30y, T1 + k3y);

  y0 = y0 + (1.0/6.0)*(k10 + 2*k20 + 2*k30 + k40);
  y1 = y1 + (1.0/6.0)*(k11 + 2*k21 + 2*k31 + k41);
    
  T0 = T0 + (1.0/6.0)*(k10y + 2*k20y + 2*k30y + k40y);
  T1 = T1 + (1.0/6.0)*(k1y + 2*k2y + 2*k3y + k4y);

}