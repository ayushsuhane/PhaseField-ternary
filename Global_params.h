#ifndef GLOBAL_PARAMS_H
#define GLOBAL_PARAMS_H

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_math.h"
gsl_rng *URN;
gsl_rng *PDF;

//Readinput
double *c_init,*K,*InterD,*m_slope,*c_ref,*c_eq;
double tau_constant;
float *mu_seg_low,*mu_seg_high;
//double T_ref,T_Init,T_Hold,Timesec,Lenx,Leny,aspectratio,eps,Init_strainenergy,se,Diffusivity,Vm;
float *Comp_Initial,*Comp_Hold;
float T_Init,T_Hold,Time_Sec,RUNTIME;
int Grid_x,Grid_y,Num_Grains;
int move_matrix[5][2];
float Heating_Rate,Vol_Fraction,Aspect_Ratio,Molar_Volume,Epsilon,Strain_Energy,Surface_Energy,D_Max;
int NComponents,NPhases,Num_Grains,Seeds,Seeds_Ferrite,N_Cut,Random_Seed;
float E_scaling,L_scaling,t_scaling;
float dx_nd,dx_real,dt_nd,dt_real,Len_x,Len_y,dy_nd;
int N_Timesteps;
int track_x,track_y;
float *Mu_stencil;
float *mu_initial,*mu_equilibrium,*mu_bulk;
char phase_init,**Elements;
float *Molecular_Weight,*Bulk_Composition,**Fit_Constant,**Diffusion_Matrix,*Bulk_Comp_Mole;
double C_pearlite;
double volfrac_pearlite;
//Initialize
double *mole_T,*mole_ref,*mole_init;
float R;
//Allocate
float *Mu_Matrix,*New_Mu_Matrix,*Comp;
double *derivative_c_matrix,*derivative_mu_matrix;
double *tau_computed;
char phase_array[3];
struct node
{
	int grainnum;
	float phi;
	char phase;
	int flag;
};
struct node *Matrix,*New_Matrix,*stencil,*PFstencil,*list,*update_list,blank_matrix;

//Solver
int Lenlist;
float *Div_Flux,*Source;
float *updated_mu;
#endif


