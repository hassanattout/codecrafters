#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define nmax 100

int n,t1;
float k1x,k2x,k3x,k4x;
float k1v,k2v,k3v,k4v;
float tmax = 10.0;
float t[nmax+1];
float xan[nmax+1],van[nmax+1];
float u1[nmax+1];
float u2[nmax+1];
float xnum[nmax+1],vnum[nmax+1];

//Question 1

float xana(float t){
  return exp(-3*t/8)*(-2*cos(sqrt(311)*t/8)-6/sqrt(311)*sin(sqrt(311)*t/8))+2;
}

float vana(float t){
  return 80*exp(-3*t/8)*sqrt(311)*sin(sqrt(311)*t/8)/311;
}

float f1(float u1, float u2){
  return u2;
}

float f2(float u1, float u2){
  return -5*u1-3.0/4.0*u2+10;
}

int main(){
  float dt = tmax/nmax;
  FILE*output = fopen("diagramme","w");
  if(output == NULL){
    puts("Erreur");
    return 0;
  }
  for(t1=0;t1<tmax;t1++){
    fprintf(output,"%f\t%f\n",xana(t1),vana(t1));
    printf("xana: %f\tvana: %f\n",xana(t1),vana(t1));
  }
  fclose(output);
//Question 2
  FILE*out = fopen("xana","w");
  if(out == NULL){
    puts("Erreur");
    return 0;
  }
  for(n=0;n<nmax;n++){
    t[n] = n*dt;
    xan[n] = xana(t[n]);
    van[n] = vana(t[n]);
    k1x = f1(u1[n],u2[n]);
    u1[n+1] = u1[n]+dt*k1x;
    k1v = f2(u1[n],u2[n]);
    u2[n+1] = u2[n]+dt*k1v;
    xnum[n+1] = u1[n+1];                      //solution numerique approchee
    vnum[n+1] = u2[n+1];                     //solution numerique approchee
  }
  float Eeuler = fabs(xnum[nmax]-xana(tmax));
  fprintf(out,"%f",xnum[n]);
  printf("Erreur nmax: %f\n",Eeuler);
  fclose(out);
  
//Question 4: En augmentant nmax, l'erreur devient de plus en plus precise. Elle tend vers 0. De plus pour nmax = 50 000, le x range est invalide.
  
//Question 5 Méthode d'Euler Modifiée
  FILE*out1 = fopen("erreur.eul","w");
  if(out1 == NULL){
    puts("Erreur");
    return 0;
  }
  for(n=0;n<nmax;n++){
    t[n] = n*dt;
    xan[n] = xana(t[n]);
    van[n] = vana(t[n]);
    k1x = f1(u1[n],u2[n]);
    k1v = f2(u1[n],u2[n]);
    k2x = f1(u1[n]+dt/2*k1x,u2[n]+dt/2*k1v);
    k2v = f2(u1[n]+dt/2*k1x,u2[n]+dt/2*k1v);
    u1[n+1] = u1[n]+dt*k2x;
    u2[n+1] = u2[n]+dt*k2v;
    xnum[n+1] = u1[n+1];
    vnum[n+1] = u2[n+1];
  }
  Eeuler = fabs(xnum[nmax]-xana(tmax));
  fprintf(out,"%f\t%f\n",Eeuler,dt);
  printf("Erreur euler modifie RK2: %f\n",Eeuler);
  fclose(out);
//Question 6 Méthode RK4
  FILE*rk4 = fopen("erreur.eul","w");
  if(rk4 == NULL){
    puts("Erreur");
    return 0;
  }
  for(n=0;n<nmax;n++){
    t[n] = n*dt;
    xan[n] = xana(t[n]);
    van[n] = vana(t[n]);
    k1x = f1(u1[n],u2[n]);
    k1v = f2(u1[n],u2[n]);
    k2x = f1(u1[n]+dt/2*k1x,u2[n]+dt/2*k1v);
    k2v = f2(u1[n]+dt/2*k1x,u2[n]+dt/2*k1v);
    k3x = f1(u1[n]+dt/2*k2x,u2[n]+dt/2*k2v);
    k3v = f2(u1[n]+dt/2*k2x,u2[n]+dt/2*k2v);
    k4x = f1(u1[n]+dt*k3x,u2[n]+dt*k3v);
    k4v = f2(u1[n]+dt*k3x,u2[n]+dt*k3v);
    u1[n+1] = u1[n]+dt/6*(k1x+2*k2x+2*k3x+k4x);
    u2[n+1] = u2[n]+dt/6*(k1v+2*k2v+2*k3v+k4v);
    xnum[n+1] = u1[n+1];
    vnum[n+1] = u2[n+1];
  }
  float Erreur = fabs(xnum[nmax]-xana(tmax));
  fprintf(out,"%f\t%f\n",Eeuler,dt);
  printf("Erreur RK4: %f",Erreur);
  fclose(rk4);
}