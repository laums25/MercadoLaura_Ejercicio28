#include <iostream>
#include <fstream>
#include <cmath>

float g = 9.8;
float k = 0.9;
float delta = 0.01;
int N = 25;
float ang = 34.0;
float v0 = 22.0;
float v0x = v0* cos(ang*M_PI/180.0);
float v0y = v0* sin(ang*M_PI/180.0);

int main(void)
{
  ofstream fout("datos-array.txt");
  float x = 0.0, y = 0.0;
  float time;
  double p[N] = {1, 0};
  for(time = 0; time <= 10; time += DeltaT) {
    fout << time << "\t" << y[0] << "\t" << y[1] << endl;
    rk4_vector(time, DeltaT, y);
  }
  fout.close();

  return 0;
}

void pro(double t, double y0, double y1){
    float vt= sqrt(v0x*v0x + v0y*v0y);
    float dx=- k*v0x*v0x/vt;
    float dy=- g -k *v0y*v0y / vt;
}


float rk(float t, float delta, float y[]){
    float k1[N], k2[N], k3[N], k4[N];
    float yt[N];
    int i;
    for (i=0; i<N; ++i)
    {
        k1[i] = delta * pro(t,y,i);
    }
    for (i=0; i<N; ++i)
    {
        yt[i] = y[i] + k1[i]/2;
    }
    for (i=0; i<N; ++i)
    {
        k2[i] = delta * pro(t + h/2, y,i);
    }
    for (i=0; i<N; ++i)
    {
        yt[i] = y[i] + k2[i]/2;
    }
    for (i=0; i<N; ++i)
    {
        k3[i] = delta * pro(t + h/2, y,i);
    }
    for (i=0; i<N; ++i)
    {
        yt[i] = y[i] + k3[i];
    }  
    for (i=0; i<N; ++i)
    {
        k4[i] = delta * pro(t + h, y,i);
    }
    for (i=0; i<N; ++i)
    {
        y[i] = y[i] + (1.0/6.0)*(k1[i]+ 2*k2[i] + 2*k3[i] + k4[i]);
    } 
}
