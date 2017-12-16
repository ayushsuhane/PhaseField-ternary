#include"function.h"
#include"Global_params.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_linalg.h>

double tor(int k,int x,int y){
 extern int Lenlist;
 extern struct node *Matrix;
 int i=0,j=0;
 double sum = 0;
 for(i=0;i<Lenlist;i++){	
  for(j=i+1;j<Lenlist;j++){
   if(i!=k && (j!=k)) sum = sum + Matrix[Location(x,y,j)].phi*Matrix[Location(x,y,i)].phi;
  }
 }
 return(sum);
}
	
double Sumphi(int k,int x,int y){
 extern int Lenlist;
 extern struct node *Matrix;
 int i=0;
 double value =0.0;
 for(i=0;i<Lenlist;i++){
  if(i!=k) value = value + Matrix[Location(x,y,i)].phi;
 }
 return(value);
}

void PFSolver(int gidx,int gidy,int t){
 extern struct node *list,*update_list,*Matrix,blank_matrix;
 extern int N_Cut, Lenlist,Grid_x,Grid_y, move_matrix[5][2];
 extern int track_x,track_y;
 extern float Epsilon,Surface_Energy, dt_nd, T_Hold,T_Init;
 extern float *Mu_Matrix, *updated_mu, *mu_equilibrium, *Comp;
 extern double tau_constant,*derivative_c_matrix;
 extern char phase_array[3];	
 int i=0,j=0,k=0,l=0,n=0, val = 0,value_phi=1,x=0,y=0,key_phase,s,flag=1;
 double Chemical[Lenlist],Interface[Lenlist],Net[Lenlist],lambda = 0.0,p =0;
 float TAU=0,B[NComponents-1];
 double Der_Ht = 0.0;
 double Kai[(NComponents-1)*(NComponents-1)], Inv_Kai[(NComponents-1)*(NComponents-1)]; 
 
 update_list = (struct node*)malloc((N_Cut)*sizeof(struct node));
 updated_mu = (float *)malloc((NComponents-1)*sizeof(float));

 for(n=0;n<NComponents;n++){ 
  updated_mu[n] = 0.0;
  Comp[n*Grid_x*Grid_y + gidy*Grid_x + gidx] = 0.0;
 }
 for(k=0;k<N_Cut;k++) update_list[k] = blank_matrix;
 for(k=0;k<Lenlist;k++){
 Chemical[k] = 0.0;
 Interface[k] = 0.0;
 Net[k] = 0.0;
 }
 for(n=0;n<NComponents-1;n++){
   for(k=0;k<NComponents-1;k++) Kai[n*(NComponents-1) + k] = 0.0;
  }
 if(t<1000){
  value_phi = 0;
  TAU = 0.05*tau_constant;
 }
 else value_phi = 1;

 float dummy_mu[NComponents-1];
/************************************************************************************************************/
 if(Lenlist > 1){
  for(k=0;k<NComponents-1;k++) dummy_mu[k] = Mu_Matrix[Location(gidx,gidy,k)];
  if(value_phi == 1) TAU = tau(Lenlist,gidx,gidy);
  for(k=0;k<Lenlist;k++){
   Interface[k] = 2.0*Epsilon*Epsilon*Laplacian(k,gidx,gidy) - ((16.0/((3.14)*(3.14)))*Sumphi(k,gidx,gidy))-40.0*tor(k,gidx,gidy);
   for(l=0;l<Lenlist;l++) Chemical[k] = Chemical[k] + Epsilon*value_phi*(-1.0)*(GP(list[l].phase,dummy_mu,T_Hold,gidx,gidy)*Derivative_H(l,k,Lenlist,gidx,gidy));
   lambda += Interface[k] + Chemical[k];
  }
  for(k=0;k<Lenlist;k++){
   Net[k] = (Interface[k] + Chemical[k]) - (lambda/(Lenlist));
   update_list[k] = list[k];
   if(TAU != 0) update_list[k].phi  = list[k].phi + (Surface_Energy*dt_nd/(TAU*Epsilon*Epsilon))*(Net[k]);
  }
  if(gidx == track_x && gidy == track_y){
   for(k=0;k<Lenlist;k++) printf("Net Force %c %d = %lf Chemical = %lf Interface = %lf\n",list[k].phase,k,Net[k], Chemical[k],Interface[k]);
  }
  lambda = 0.0;
  for(k=0;k<Lenlist;k++){
   if(update_list[k].phi > 1.0) update_list[k].phi = 1.0;
   if(update_list[k].phi < 0.0) update_list[k].phi = 0.0;
   lambda = lambda + update_list[k].phi;	 
  }
  if(lambda != 1.0){
   for(k=0;k<Lenlist;k++) update_list[k].phi = update_list[k].phi/lambda;
  }
 }
 else{
  for(k=0;k<N_Cut;k++) update_list[k] = list[k];
 }


 
 if(t > 100 ) flag = CheckGrad_Mu(gidx,gidy);
//
 if(flag == 1 && t>50){
  Calc_Flux(gidx,gidy);
  extern float *Source,*Div_Flux;
  int key_grain;
  int phi_val = 0.0;
  Source = (float *)malloc((NComponents-1)*sizeof(float));
  for(k=0;k<NComponents-1;k++){
   Source[k] = 0.0;
   dummy_mu[k] = Mu_Matrix[Location(gidx,gidy,k)];
  }
  for(n=0;n<NComponents-1;n++){
   for(k=0;k<Lenlist;k++){
    if(update_list[k].grainnum == list[k].grainnum){
     Der_Ht = 0.0;
     for(i=0;i<Lenlist;i++) Der_Ht = Der_Ht + ((update_list[i].phi-list[i].phi)/(dt_nd))*Derivative_H(k,i,Lenlist,gidx,gidy); 
     Source[n] = Source[n] + (1.0)*Composition(list[k].phase,n,dummy_mu,T_Hold)*Der_Ht;
    }
    else Source[n] = Source[n] + 0.0;
   }
  }
//  for(i=0;i<(NComponents-1)*(NComponents-1);i++) Kai[i] = 0.0;
  for(i=0;i<NComponents-1;i++){
   for(j=0;j<NComponents-1;j++){
    for(k=0;k<Lenlist;k++){	
     phi_val = list[k].phi;
     for(n=0;n<3;n++){
      if(list[k].phase == phase_array[n]) key_phase = n; 
     }
     Kai[i*(NComponents-1) +j] = Kai[i*(NComponents-1)+j] + derivative_c_matrix[derivative_Location(key_phase,j,i,gidx,gidy)]*H(k,Lenlist,list);
    }
   }
  }
  
  /********Direct value of inverse matrix for 3 components*******/
  float product;
  product = Kai[0]*Kai[3] - Kai[1]*Kai[2];
  Inv_Kai[0] = Kai[3]/product;
  Inv_Kai[1] = -Kai[1]/product;
  Inv_Kai[2] = -Kai[2]/product;
  Inv_Kai[3] = Kai[0]/product;
 /****************************************************************/
/*
  gsl_matrix_view m = gsl_matrix_view_array(Kai,NComponents-1,NComponents-1);
  gsl_matrix_view inv = gsl_matrix_view_array(Inv_Kai,NComponents-1,NComponents-1);
  gsl_permutation * pr = gsl_permutation_alloc(NComponents-1);
  gsl_linalg_LU_decomp(&m.matrix,pr,&s);
  gsl_linalg_LU_invert(&m.matrix,pr,&inv.matrix);
  gsl_permutation_free (pr);
*/
  for(k=0;k<Lenlist;k++){
    if((list[k].phase == 'p')) val = val+1;
   }
  if(val==0){
   for(i=0;i<NComponents-1;i++) B[i] = (Div_Flux[i]) - Source[i];
   if(gidx == track_x && gidy == track_y){
    for(k=0;k<NComponents-1;k++) printf("Div_Flux = %lf, Source = %lf\n",Div_Flux[k],Source[k]);
   }
   for(i=0;i<NComponents-1;i++){
    for(j=0;j<NComponents-1;j++) updated_mu[i] += (Inv_Kai[i*(NComponents-1)+j]*B[j]);
   }
//   for(i=0;i<NComponents-1;i++) updated_mu[i] = value_phi*dt_nd*updated_mu[i] + dummy_mu[i];
     for(i=0;i<NComponents-1;i++) updated_mu[i] = dt_nd*updated_mu[i] + dummy_mu[i];
  }
  else{
   for(k=0;k<NComponents-1;k++) updated_mu[k] = dummy_mu[k];
  }
 
 free(Div_Flux);
 free(Source);
 }
 else{
 for(k=0;k<NComponents-1;k++) updated_mu[k] = Mu_Matrix[Location(gidx,gidy,k)];
 }
for(n=0;n<NComponents-1;n++){ 
 for(k=0;k<Lenlist;k++) Comp[n*Grid_x*Grid_y + gidy*Grid_x + gidx] = Comp[n*Grid_x*Grid_y + gidy*Grid_x +gidx] + Composition(list[k].phase,n,updated_mu,T_Hold)*H(k,Lenlist,update_list);
  }
 if(gidx == track_x && gidy == track_y){
  for(k=0;k<NComponents-1;k++) printf("Mu %d: updated - %lf, previous = %lf\n",k,updated_mu[k],Mu_Matrix[Location(gidx,gidy,k)]);
 }
free(list);

}


void Updation(int gidx,int gidy){
 extern float *updated_mu,*New_Mu_Matrix;
 extern struct node *update_list,*New_Matrix,blank_matrix;
 extern int Lenlist,N_Cut,NComponents,Grid_x,Grid_y;
 double min = 0.000001;
 int k=0,n=0,count = 0;
/****Mu_Update*****/
 for(n=0;n<NComponents-1;n++) New_Mu_Matrix[n*Grid_x*Grid_y + gidy*Grid_x + gidx] = updated_mu[n];
/***********************/
/****Phi Update *********/
 for(k=0;k<Lenlist;k++){
  if((update_list[k].phi)>(min) && update_list[k].flag!=0){
   New_Matrix[Location(gidx,gidy,count)] = update_list[k];
   count = count + 1;
  }
 }
 for(k=count;k<N_Cut;k++){
  New_Matrix[Location(gidx,gidy,k)] = blank_matrix;
 }
/**************************/ 
 free(update_list);
 free(updated_mu);
}
