#include"function.h"
#include"Global_params.h"
void FreeMemory(){

extern float *Molecular_Weight;
extern char **Elements;
extern float *Bulk_Composition; 
extern float **Fit_Constant;
extern float **Diffusion_Matrix;
extern float *Comp_Initial,*Comp_Hold;
extern float *Mu_Matrix,*New_Mu_Matrix,*Comp;
extern float *mu_initial,*mu_equilibrium,*Bulk_Comp_Mole,*mu_bulk;
extern double *derivative_c_matrix,*derivative_mu_matrix,*tau_computed;
extern int NComponents,NPhases;
extern struct node *Matrix,*New_Matrix;
int i=0;


 free(Molecular_Weight);
 free(Bulk_Composition);
 for(i=0;i<NComponents;i++){
  free(Fit_Constant[i]);
  free(Elements[i]);
 }
 free(Fit_Constant);
 for(i=0;i<NPhases;i++) free(Diffusion_Matrix[i]);
 free(Diffusion_Matrix);
 free(Comp_Initial);
 free(Comp_Hold);
 free(Mu_Matrix);
 free(New_Mu_Matrix);
 free(Matrix);
 free(New_Matrix);
 free(Comp);
 free(mu_initial);
 free(mu_equilibrium);
 free(mu_bulk);
 free(Bulk_Comp_Mole);
free(derivative_c_matrix);
free(derivative_mu_matrix);
free(tau_computed);
}
