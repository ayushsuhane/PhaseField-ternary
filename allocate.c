#include"function.h"
#include"Global_params.h"
#include<stdlib.h>
#include"setup_fn.h"
#include<math.h>
#include<gsl/gsl_linalg.h>

void EquilibriumMatrixTau(){
 extern int NComponents,NPhases;
 extern double *tau_computed;
 tau_computed = (double *)malloc((NPhases*(NPhases))*sizeof(double)); 
/*
NPHases = 3
alpha-alpha alpha-gamma alpha-theta
gamma-alpha gamma-gamma gamma-theta
theta-alpha theta-gamma theta-theta

NPhases = 2
alpha-alpha alpha-gamma
gamma-alpha gamma-gamma
*/
 extern float Epsilon,*mu_equilibrium,T_Hold;
 double derivative_mu_eq[(NComponents-1)*(NComponents-1)],derivative_c_eq[(NComponents-1)*(NComponents-1)];
 double Diff_Mat[(NComponents-1)*(NComponents-1)];
 extern float **Diffusion_Matrix;
 double Inv_Diff_Mat[(NComponents-1)*(NComponents-1)];
 int i=0,j=0,k=0,l=0;
 double temp = T_Hold;
 int s;
 for(i=0;i<(NComponents-1)*(NComponents-1);i++){
  derivative_mu_eq[i] = 0.0;
  derivative_c_eq[i] = 0.0;
  Diff_Mat[i] = 0.0;
  Inv_Diff_Mat[i] = 0.0;
 }
 for(i=0;i<NComponents-1;i++){
  for(j=0;j<NComponents-1;j++){
    if(i==j) derivative_mu_eq[j*(NComponents-1) + i] = 1.0/(Composition('g',NComponents-1,mu_equilibrium,temp)) + 1.0/Composition('g',i,mu_equilibrium,temp);
    else derivative_mu_eq[j*(NComponents-1) + i] = 1.0/Composition('g',NComponents-1,mu_equilibrium,temp);   
   }
  }
/*
 gsl_matrix_view m = gsl_matrix_view_array(derivative_mu_eq,(NComponents-1),(NComponents-1));
 gsl_matrix_view inv = gsl_matrix_view_array(derivative_c_eq,(NComponents-1),(NComponents-1));
 gsl_permutation * p = gsl_permutation_alloc((NComponents-1));
 gsl_linalg_LU_decomp(&m.matrix,p,&s);
 gsl_linalg_LU_invert(&m.matrix,p,&inv.matrix);
 gsl_permutation_free(p);
*/
  derivative_c_eq[0] = Composition('g',0,mu_equilibrium,temp)*(1.0 - Composition('g',0,mu_equilibrium,temp));
  derivative_c_eq[1] = -Composition('g',0,mu_equilibrium,temp)*Composition('g',1,mu_equilibrium,temp);
  derivative_c_eq[2] = -Composition('g',0,mu_equilibrium,temp)*Composition('g',1,mu_equilibrium,temp);
  derivative_c_eq[3] = Composition('g',1,mu_equilibrium,temp)*(1.0 - Composition('g',1,mu_equilibrium,temp));

 for(i=0;i<(NComponents-1);i++){
  for(k=0;k<(NComponents-1);k++){ 
   for(l=0;l<NComponents-1;l++){
    Diff_Mat[i*(NComponents-1) + k] += Diffusion_Matrix[1][i*(NComponents-1) + l]*derivative_c_eq[l*(NComponents-1) + k];
   }
  }
 }

 gsl_matrix_view m2 = gsl_matrix_view_array(Diff_Mat,(NComponents-1),(NComponents-1));
 gsl_matrix_view inv2 = gsl_matrix_view_array(Inv_Diff_Mat,(NComponents-1),(NComponents-1));
 gsl_permutation * q = gsl_permutation_alloc((NComponents-1));
 gsl_linalg_LU_decomp(&m2.matrix,q,&s);
 gsl_linalg_LU_invert(&m2.matrix,q,&inv2.matrix);
 gsl_permutation_free(q);

 extern double tau_constant; 
// tau_constant = 3.5*((Composition('g',0,mu_equilibrium,T_Hold) - Composition('a',0,mu_equilibrium,T_Hold))*(Composition('g',0,mu_equilibrium,T_Hold) - Composition('a',0,mu_equilibrium,T_Hold)))*Inv_Diff_Mat[0*(NComponents-1) + 0];
 double Comp_Diff[NComponents-1];
 extern char phase_array[3];
 
 double Mul_AB[NComponents-1];
 double val;
 for(i=0;i<NPhases;i++){
  for(j=i;j<NPhases;j++){
   if(i!=j){
    val = 0;
    for(k=0;k<(NComponents-1);k++){
     Comp_Diff[k] = Composition(phase_array[i],k,mu_equilibrium,T_Hold) - Composition(phase_array[j],k,mu_equilibrium,T_Hold);
     Mul_AB[k] = 0;
    }
    for(l=0;l<NComponents-1;l++){
     for(k=0;k<NComponents-1;k++) Mul_AB[l] += Comp_Diff[k]*Inv_Diff_Mat[k*(NComponents-1) + l];
    }
    for(k=0;k<(NComponents-1);k++) val += Mul_AB[k]*Comp_Diff[k];
    tau_computed[j*(NPhases) + i] = val;
   }
  tau_computed[j*(NPhases)+i] *= Epsilon*(0.063828+0.158741);
  tau_computed[i*(NPhases)+j] = tau_computed[j*(NPhases)+i]; 
  }
 }
 tau_constant = 1.5*tau_computed[0*NPhases + 1];
 printf("TAU MATRIX\n");
 for(i=0;i<NPhases;i++) tau_computed[i*(NPhases) + i] = tau_constant;
 for(i=0;i<NPhases;i++){  
  for(j=0;j<NPhases;j++){
   printf("%f ",tau_computed[j*(NComponents-1) + i]);
  }
  printf("\n");
 }
 /*******************************/
/****Comment out...only for testing******/
 /*for(i=0;i<NPhases*NPhases;i++) tau_computed[i] *= 0.1;
 /******************************/

 printf("tau_constant = %f\n", tau_constant);
}


void Allocate(){

 extern int Grid_x,Grid_y,
	    NComponents,
	    N_Cut;
 extern struct node *Matrix,*New_Matrix;
 extern float *Mu_Matrix,
	      *New_Mu_Matrix,
	      *Comp; 
 extern float *mu_initial,
	      *mu_equilibrium;
 extern float *mu_bulk;
 int i=0;

 extern struct node blank_matrix;
 /*definition of blank matrix*/
 blank_matrix.grainnum = 0;
 blank_matrix.phi = 0.0;
 blank_matrix.phase = 'n';
 blank_matrix.flag = 0;
	      
 Matrix = (struct node*)malloc((Grid_x*Grid_y*N_Cut)*sizeof(struct node));  //free in memory.c
 New_Matrix = (struct node*)malloc((Grid_x*Grid_y*N_Cut)*sizeof(struct node)); //free in memory.c
 for(i=0;i<Grid_x*Grid_y*N_Cut;i++){
  Matrix[i] = blank_matrix;
  New_Matrix[i] = blank_matrix;
 }


 Mu_Matrix = (float *)malloc((Grid_x*Grid_y*(NComponents-1))*sizeof(float));  //free in memory.c
 New_Mu_Matrix = (float *)malloc((Grid_x*Grid_y*(NComponents-1))*sizeof(float));  //free in memory.c
 Comp = (float *)malloc((Grid_x*Grid_y*(NComponents))*sizeof(float));  //free in memory.c
 for(i=0;i<(Grid_x*Grid_y*(NComponents-1));i++){
  Mu_Matrix[i] = 0.0;
  New_Mu_Matrix[i] = 0.0;
 }
 for(i=0;i<(Grid_x*Grid_y*(NComponents));i++) Comp[i] = 0.0;

 printf("Memory for matrix allocated\n");

 mu_initial = (float*)malloc((NComponents-1)*sizeof(float));  //free in memory.c
 mu_equilibrium = (float*)malloc((NComponents-1)*sizeof(float));  //free in memory.c
 mu_bulk = (float*)malloc((NComponents-1)*sizeof(float));   //free in memory.c


 extern float *Bulk_Composition,*Molecular_Weight;
 extern float *Bulk_Comp_Mole;

 Bulk_Comp_Mole = (float *)malloc(NComponents*sizeof(float));  //free in memory.c
 Wt2Mole(Bulk_Composition,Molecular_Weight,Bulk_Comp_Mole);
 Scaling();
 Mu_Calc(mu_initial,T_Init);
 Mu_Calc(mu_equilibrium,T_Hold);
 extern float R;
 R = 1.0;
 extern int track_x,track_y;

 track_x = 0.375*Grid_x;
 track_y = 0.5*Grid_y;
 float *mu_alpha = (float*)malloc((NComponents-1)*sizeof(float));
// mu_alpha = (float*)malloc((NComponents-1)*sizeof(float));

// track_x = 6;
// track_y = 76;
 extern double *derivative_c_matrix;
 derivative_c_matrix = (double *)malloc(Grid_x*Grid_y*NPhases*(NComponents-1)*(NComponents-1)*sizeof(double));
 extern double *derivative_mu_matrix;
 derivative_mu_matrix = (double *)malloc(Grid_x*Grid_y*NPhases*(NComponents-1)*(NComponents-1)*sizeof(double));
 for(i=0;i<Grid_x*Grid_y*NPhases*(NComponents-1)*(NComponents-1);i++){
  derivative_c_matrix[i] = 0.0;
  derivative_mu_matrix[i] = 0.0;
 } 
 EquilibriumMatrixTau();

 for(i=0;i<NComponents-1;i++) mu_bulk[i] = R*T_Init*log(Bulk_Comp_Mole[i]/Bulk_Comp_Mole[NComponents-1]);
 for(i=0;i<NComponents-1;i++) {
  printf("mu bulk %d = %lf\n",i,mu_bulk[i]);
  printf("mu_equilibrium %d = %lf\n",i,mu_equilibrium[i]);
  printf("mu_initial %d = %lf\n",i,mu_initial[i]);
 }

/* GP Check */
 extern float *Comp_Initial;
 Mu_Calc_Alpha(Comp_Initial,mu_alpha,T_Init);
 printf("at T_initial : GP_alpha = %lf, GP_beta = %lf\n",GP('a',mu_initial,T_Init,0,0),GP('g',mu_initial,T_Init,0,0));
 printf("at T_initial : G_alpha = %lf, G_beta = %lf\n",Free_Energy('a',mu_initial,T_Init),Free_Energy('g',mu_initial,T_Init));
 for(i=0;i<NComponents-1;i++) printf("%d at T_initial : c_alpha = %lf, c_beta = %lf\n",i,Composition('a',i,mu_initial,T_Init),Composition ('g',i,mu_initial,T_Init));
/*          */
 free(mu_alpha);
}

