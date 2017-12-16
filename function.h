#ifndef FUNCTION_H
#define FUNCTION_H
#include"Global_params.h"
 
void ReadInputs();
void Allocate();
void Calc_Mu_Seg(char ,float c[3], float *, float);
void EquilibriumMatrixTau();
void SystemSetup();
void MatrixSetup(int, double);
int derivative_Location(int , int ,int ,int ,int );	
int CheckGrad_Mu(int ,int );
//double Composition(int NComponents,int ,int ,double Mu_stencil[5*NComponents]);
float Composition(char ,int ,float *,float);
void PopulateSeeds(int ,char ,int );
void ProbabilityNucleation(char, int);
float Free_Energy(char ,float *,float);
float GP(char phase,float *,float temp,int gidx,int gidy);
int Location(int ,int ,int);
double H(int ,int , struct node *List);
//double Der_H(int ,int);
float Derivative_H(int ,int ,int, int, int);
void Derivative_mu(char ,float * ,float);
float tau(int,int,int);
double int_energy(int );
double Laplacian(int ,int, int);
void StencilGeneration(int i,int j);
void Calc_Flux(int,int);
int Periodic(int, int);
void PFSolver(int i,int j,int t);
void Updation(int i,int j);
void Visualization(int t);
float AtpertoWtper(float c[2],int );
//void Superimpose(int );
void Volume_Fraction_Austenite(int t);
void Check_Phases(int );
void FreeMemory();
#endif
