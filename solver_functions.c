#include"function.h"
#include"Global_params.h"
#include<math.h>
#include<gsl/gsl_linalg.h>

float Composition(char phase,int n,float *mu,float temp){
 extern float **Fit_Constant,R;
 extern int NComponents;
 int i=0;
 float sum = 0.0;
 float B[NComponents-1];
 if(phase=='a'){
  for(i=0;i<NComponents-1;i++) B[i] = Fit_Constant[i][0]*(temp - Fit_Constant[i][1])/(Fit_Constant[i][1]) - Fit_Constant[NComponents-1][0]*(temp - Fit_Constant[NComponents-1][1])/(Fit_Constant[NComponents-1][1]);
 }
 if(phase=='g'){
  for(i=0;i<NComponents-1;i++) B[i] = 0.0;
 }
 for(i=0;i<NComponents-1;i++) sum += exp((mu[i] - B[i])/(R*temp));
 if(n == NComponents-1) return(1.0/(1.0+sum));
 else return(exp((mu[n] - B[n])/(R*temp))/(1.0 + sum));
/*
if(phase == 'p') 
	{
	extern double C_pearlite;
	return(C_pearlite);
	}
*/
}



float Free_Energy(char phase,float *mu,float temp){
 extern float R;
 extern int NComponents;
 extern float **Fit_Constant;
 int i=0;
 float F=0;
 float sum=0.0,sum_R = 0.0;
 for(i=0;i<NComponents;i++) sum_R += R*temp*(Composition(phase,i,mu,temp)*log(Composition(phase,i,mu,temp)));
 if(phase == 'a'){
  for(i=0;i<NComponents;i++){
   sum += (Fit_Constant[i][0]*(temp - Fit_Constant[i][1])/Fit_Constant[i][1])*Composition(phase,i,mu,temp);
  }
  F = sum_R + sum;
 }
 else if(phase == 'g') F = sum_R;
 else if(phase == 'p'){
  extern double G0,L;
 }
 return(F);
}

float GP(char phase,float *mu,float temp,int gidx,int gidy){
 float grand_potential = 0.0;
 float sum = 0.0;
 extern int NComponents;
 extern int track_x,track_y;
 int i =0,j=0;
 for(i=0;i<NComponents-1;i++) sum = sum + mu[i]*Composition(phase,i,mu,temp);
 grand_potential = Free_Energy(phase,mu,temp)-sum;
/*
if(gidx == track_x && gidy == track_y) 
   {
   printf("phase = %c Free energy =%f ",phase,Free_Energy(phase,mu,temp));
   for(j=0;j<NComponents-1;j++) printf("mu_%d = %f Composition = %f ",j,mu[j],Composition(phase,j,mu,temp));
   printf("Composition_Fe = %lf ",Composition(phase,NComponents-1,mu,temp));
   printf("sum = %f Grand potential = %f\n",sum,grand_potential);
   }
*/
 return(grand_potential);
}




float tau(int Lenlist,int gidx,int gidy){
 extern struct node *Matrix,*list;
 int i=0,j=0,k=0,l=0;
 extern double tau_constant;
 double val=0,val_denom=0.0,sum = 0,sum_phi = 0,val_num=0.0;
 char phase_i,phase_j;
 extern char phase_array[3];
 extern double *tau_computed;
 float TAU=0.0;
 extern int track_x,track_y;
 extern int Grid_x,Grid_y;
 int key_i,key_j;
 double TAU_num,TAU_denom;
 for(i=0;i<Lenlist;i++){
  phase_i = list[i].phase;
  for(k=0;k<3;k++){
   if(phase_i == phase_array[k]) key_i = k;
  }
  for(j=i;j<Lenlist;j++){
   phase_j = list[j].phase;
   for(k=0;k<3;k++){
    if(phase_j == phase_array[k]) key_j = k;
   }
   if(phase_i != phase_j){ 
    TAU_num += tau_computed[key_j*(NPhases) + key_i]*list[i].phi*list[j].phi;
    TAU_denom += list[i].phi*list[j].phi;
   }
  }
 }
 if(TAU_denom != 0.0 && TAU_num != 0.0) TAU = TAU_num/TAU_denom;
 else TAU = tau_constant;
 
 if(track_x == gidx && track_y == gidy){
  printf("TAU = %f denominator = %f numerator = %f\n",TAU,TAU_denom,TAU_num);
  printf("\n");
  if(TAU == 0) printf("Tau is Zero\n");
 }
 if(TAU<0) TAU *= -1.0;
 return(TAU);
}

int CheckGrad_Mu(int gidx,int gidy){
 extern float *Mu_Matrix;
 extern int move_matrix[5][2];
 extern int NComponents;
 int i=0,j=0,k=0,n=0,x,y;
 float val;
 int flag = 0;
 for(k=0;k<NComponents-1;k++){
  val = Mu_Matrix[Location(gidx,gidy,k)];
  for(n=0;n<5;n++){
   x = gidx+move_matrix[n][0];
   y = gidy+move_matrix[n][1]; 
   if(val != Mu_Matrix[Location(x,y,k)]){
    flag = 1;
    return(flag);
   }
  }
 }
return(flag);
}

int Location(int x,int y,int k){
 extern int Grid_x,Grid_y;
 return(k*Grid_x*Grid_y + y*Grid_x + x);
}

double H(int k,int Lenlist,struct node *sol_list){
 int i =0,j=0;
 double val = 0;
 if(Lenlist>2){
  for(i=0;i<Lenlist;i++){
   for(j=i+1;j<Lenlist;j++) {
    if((j!=k) && i!=k && sol_list[i].flag!=0) val = val + (sol_list[i].phi)*(sol_list[j].phi);
   }
  }
 }
 else val = 0.0;
 val = (sol_list[k].phi)*(sol_list[k].phi)*(3.0 - 2.0*sol_list[k].phi) + 2.0*val*sol_list[k].phi; 
 return(val);
}

/*
float Der_H(int k,int Lenlist){
extern struct node *PFstencil;
float val= 0;
int i=0,j=0;
if(Lenlist>2)
   {
   for(i=0;i<Lenlist;i++) 
       {
       for(j=i+1;j<Lenlist;j++) 
           {
           if((j!=k) && i!=k) val = val + (PFstencil[2*Lenlist + i].phi)*(PFstencil[2*Lenlist + j].phi);
           }
       }
   }
val = 6.0*PFstencil[2*Lenlist + k].phi*(1.0 - PFstencil[2*Lenlist + k].phi) + 2.0*val;
return(val);
}
*/

float Derivative_H(int k,int l,int Lenlist,int gidx, int gidy){
 extern struct node *Matrix,*list;
 double val= 0.0;
 int i=0,j=0;
 if(Lenlist>2){
  if(k==l){
   for(i=0;i<Lenlist;i++){ 
    for(j=i+1;j<Lenlist;j++){ 
     if((j!=k) && i!=k) val = val + (list[i].phi*list[j].phi);
     }
    }
   }
  else{
   for(i=0;i<Lenlist;i++){ 
    if((i!=k) && (i!=l)) val = val + list[i].phi;
   }
  }
 }
 else val = 0;
 if(k==l) val = 6.0*list[k].phi*(1.0 - list[k].phi) + 2.0*val;
 else val = 2.0*(val)*list[k].phi;
 return(val);
}


double Laplacian(int k, int gidx, int gidy){
 extern struct node *Matrix,*list;
 extern int move_matrix[5][2];
 extern int Lenlist,N_Cut;
 extern float dx_nd,dy_nd;
 extern int track_x,track_y;
 double value=0.0;
 double phi_value[5] = {0.0};
 int n=0,i=0;
 int x=0,y=0;
 for(n=0;n<5;n++){
  x = gidx+move_matrix[n][0];
  y = gidy+move_matrix[n][1];
  for(i=0;i<N_Cut;i++){ 
   if(list[k].grainnum == Matrix[Location(x,y,i)].grainnum){
    phi_value[n] = Matrix[Location(x,y,i)].phi;
    break;
   }
  }
 }
 value = (phi_value[3] + phi_value[1] - 2.0*phi_value[2])/(dx_nd*dx_nd) + (phi_value[0] + phi_value[4] - 2.0*phi_value[2])/(dy_nd*dy_nd); 
 return(value);
}

void Calc_Flux(int gidx,int gidy){
 extern int NComponents;
 extern struct node *Matrix;
 int i=0,j=0,n=0,k=0,l=0;
 float Mobility[5*(NComponents-1)*(NComponents-1)];
 float M[(NComponents-1)*(NComponents-1)];
 extern float T_Hold,**Diffusion_Matrix,*Mu_Matrix, Epsilon,dx_nd,dy_nd;
 extern int Lenlist,track_x,track_y,move_matrix[5][2];
 extern double *derivative_c_matrix;
 float mu[NComponents-1];	
 int val=0;  //Identifies the phase at the point
 int x,y;
 int key_phase;
 extern char phase_array[3];
 for(k=0;k<(5*(NComponents-1)*(NComponents-1));k++) Mobility[k] = 0.0;
 for(n=0;n<5;n++){
  x = gidx+move_matrix[n][0];
  y = gidy+move_matrix[n][1];
  for(l=0;l<NComponents-1;l++) mu[l] = Mu_Matrix[Location(x,y,l)];
  for(k=0;k<Lenlist;k++){
   for(l=0;l<NPhases;l++){
    if(Matrix[Location(x,y,k)].phase == phase_array[l]) key_phase = l;
   }
   for(i=0;i<((NComponents-1)*(NComponents-1));i++) M[i] = 0.0;
   for(i=0;i<NComponents-1;i++){
    for(j=0;j<NComponents-1;j++){
     for(l=0;l<NComponents-1;l++) M[i*(NComponents-1) + j] = M[i*(NComponents-1) + j] + Diffusion_Matrix[key_phase][i*(NComponents-1) + l]*derivative_c_matrix[derivative_Location(key_phase,j,l,x,y)]*Matrix[Location(x,y,k)].phi;
     Mobility[n*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j] = Mobility[n*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j] + M[i*(NComponents-1) + j];
    }
   }
  }
 }

 if(gidx==track_x && gidy == track_y){
  for(l=0;l<NComponents-1;l++) mu[l] = Mu_Matrix[Location(gidx,gidy,l)];
  printf("\n");
  printf("Phase Mobility\n");
  for(i=0;i<(NComponents-1);i++){ 
   for(j=0;j<(NComponents-1);j++)printf("%f ",M[i*(NComponents-1)+j]);
  }
  printf("\n");
  printf("Mobility values\n");
  for(n=0;n<5;n++){
   for(i=0;i<NComponents-1;i++){
    for(j=0;j<NComponents-1;j++){
     printf("%f ",Mobility[n*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j]);
    }
   }
  }
 }
 float *Flux_x,*Flux_y;
 Flux_x = (float *)malloc((2*(NComponents-1))*sizeof(float));
 Flux_y = (float *)malloc((2*(NComponents-1))*sizeof(float));
 for(i=0;i<(2*(NComponents-1));i++) Flux_x[i] = 0.0;
 for(i=0;i<(2*(NComponents-1));i++) Flux_y[i] = 0.0;
 double div_flux_x[NComponents-1],div_flux_y[NComponents-1];
 for(i=0;i<NComponents-1;i++){
  div_flux_x[i] = 0.0;
  div_flux_y[i] = 0.0;
 }
 extern float *Div_Flux;
 Div_Flux = (float *)malloc((NComponents-1)*sizeof(float));
 for(i=0;i<NComponents-1;i++) Div_Flux[i] = 0.0;
 for(i=0;i<NComponents-1;i++){
  for(j=0;j<NComponents-1;j++){
   Flux_x[0*(NComponents-1) + i] = Flux_x[0*(NComponents-1) + i] + (0.5)*(Mobility[1*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j] + Mobility[2*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j])*(Mu_Matrix[Location(gidx,gidy,j)] - Mu_Matrix[Location(gidx+move_matrix[1][0],gidy+move_matrix[1][1],j)])/(dx_nd);
   Flux_x[1*(NComponents-1) + i] = Flux_x[1*(NComponents-1) + i] + (0.5)*(Mobility[2*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j] + Mobility[3*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j])*(Mu_Matrix[Location(gidx+move_matrix[3][0],gidy+move_matrix[3][1],j)] - Mu_Matrix[Location(gidx,gidy,j)])/(dx_nd);
   Flux_y[0*(NComponents-1) + i] = Flux_y[0*(NComponents-1) + i] + (0.5)*(Mobility[4*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j] + Mobility[2*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j])*(Mu_Matrix[Location(gidx,gidy,j)] - Mu_Matrix[Location(gidx+move_matrix[4][0],gidy+move_matrix[4][1],j)])/(dy_nd);
   Flux_y[1*(NComponents-1) + i] = Flux_y[1*(NComponents-1) + i] + (0.5)*(Mobility[2*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j] + Mobility[0*(NComponents-1)*(NComponents-1) + i*(NComponents-1) + j])*(Mu_Matrix[Location(gidx+move_matrix[0][0],gidy+move_matrix[0][1],j)] - Mu_Matrix[Location(gidx,gidy,j)])/(dy_nd);
  }
 }
 for(i=0;i<NComponents-1;i++){
  div_flux_x[i] = (Flux_x[1*(NComponents-1) + i] - Flux_x[0*(NComponents-1) + i])/(dx_nd);
  div_flux_y[i] = (Flux_y[1*(NComponents-1) + i] - Flux_y[0*(NComponents-1) + i])/(dy_nd);
     //if(div_flux_x[i]!=0) printf("%f %f \n",div_flux_x[i],div_flux_y[i]);
 }
 for(i=0;i<NComponents-1;i++) {
  Div_Flux[i] = div_flux_x[i] + div_flux_y[i];
     //printf("%f ",Div_Flux[i]);
 }

 free(Flux_x);
 free(Flux_y);
}


