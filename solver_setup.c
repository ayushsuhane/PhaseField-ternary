#include"Global_params.h"
#include"function.h"
#include<math.h>
#include<gsl/gsl_linalg.h>


int CheckMuDiff(int index_i,int index_j){
 extern float *Mu_Matrix, *New_Mu_Matrix;
 extern int NComponents,Grid_x,Grid_y;
 int i=0;

 for(i=0;i<NComponents-1;i++){
  if(Mu_Matrix[i*Grid_x*Grid_y + index_j*Grid_x + index_i] != New_Mu_Matrix[i*Grid_x*Grid_y + index_j*Grid_x + index_i]) return(1);
 }
 return(0);
}

int derivative_Location(int phase, int matrix_i,int matrix_j,int index_i,int index_j){
 extern int Grid_x,Grid_y,NComponents;
 return(matrix_i + matrix_j*(NComponents-1) + index_i*(NComponents-1)*(NComponents-1) + index_j*(NComponents-1)*(NComponents-1)*Grid_x + phase*(NComponents-1)*(NComponents-1)*Grid_x*Grid_y); 
}

void Calculate_derivatives(int index_i,int index_j, double time){
 extern double *derivative_mu_matrix;
 extern double *derivative_c_matrix;
 extern float *New_Mu_Matrix,*Mu_Matrix;
 extern int track_x,track_y;
 char phase;
 float mu[NComponents-1];
 extern float T_Hold;
 float temp;
 temp = T_Hold;
 extern char phase_array[3];
 int i=0,j=0,k=0;
 
 if(time!= 0){
  for(k=0;k<NComponents-1;k++) mu[k] = New_Mu_Matrix[Location(index_i,index_j,k)];
 }
 else{
  for(k=0;k<NComponents-1;k++) mu[k] = Mu_Matrix[Location(index_i,index_j,k)];
 }
 if(index_i == track_x && index_j == track_y){
  for(k=0;k<NComponents-1;k++) printf("mu[%d] = %lf ",k,mu[k]);
//  printf("T_Hold = %lf Comp = %lf \t",T_Hold,Composition('a',0,mu,temp));
 }
 for(k=0;k<NPhases;k++){
  phase = phase_array[k];
  for(i=0;i<NComponents-1;i++){
   for(j=0;j<NComponents-1;j++){
    if(i==j) derivative_mu_matrix[derivative_Location(k,i,j,index_i,index_j)] = 1.0/(Composition(phase,NComponents-1,mu,temp)) + 1.0/Composition(phase,i,mu,temp);
    else derivative_mu_matrix[derivative_Location(k,i,j,index_i,index_j)] = 1.0/Composition(phase,NComponents-1,mu,temp);
//    if(index_i == track_x && index_j == track_y) printf("Composition = %lf dmu/dc = %lf Loc = %d phase = %c\n",Composition(phase,NComponents-1,mu,temp),derivative_mu_matrix[derivative_Location(k,i,j,index_i,index_j)],derivative_Location(k,i,j,index_i,index_j),phase);   
   }
  }
 }
 double dummy_matrix[(NComponents-1)*(NComponents-1)],dummy_inv_matrix[(NComponents-1)*(NComponents-1)];
 int s;
 for(k=0;k<NPhases;k++){
  phase = phase_array[k];
  for(i=0;i<(NComponents-1)*(NComponents-1);i++){
   dummy_matrix[i] = 0.0;
   dummy_inv_matrix[i] = 0.0;
  }
  for(i=0;i<NComponents-1;i++){
   for(j=0;j<NComponents-1;j++){
    dummy_matrix[j*(NComponents-1) + i] =  derivative_mu_matrix[derivative_Location(k,i,j,index_i,index_j)];
   }
  }

/*Directly inverse value for 3 component system******/
 float product;
 product = dummy_matrix[0]*dummy_matrix[3] - dummy_matrix[1]*dummy_matrix[2];
 if(index_i == track_x && index_j == track_y) printf("Product = %lf\n",product);
 /* 
  dummy_inv_matrix[0] = dummy_matrix[3]/product;
  dummy_inv_matrix[1] = -dummy_matrix[1]/product;
  dummy_inv_matrix[2] = -dummy_matrix[2]/product;
  dummy_inv_matrix[3] = dummy_matrix[1]/product;
  */
  dummy_inv_matrix[0] = Composition(phase,0,mu,temp)*(1.0 - Composition(phase,0,mu,temp));
  dummy_inv_matrix[1] = -Composition(phase,0,mu,temp)*Composition(phase,1,mu,temp);
  dummy_inv_matrix[2] = -Composition(phase,0,mu,temp)*Composition(phase,1,mu,temp);
  dummy_inv_matrix[3] = Composition(phase,1,mu,temp)*(1.0 - Composition(phase,1,mu,temp));
/*
  dummy_inv_matrix[0] = 1/dummy_matrix[0];
  dummy_inv_matrix[1] = 1/dummy_matrix[2];
  dummy_inv_matrix[2] = 1/dummy_matrix[1];
  dummy_inv_matrix[3] = 1/dummy_matrix[3];	
*/
/****************************************************/


/*
  gsl_matrix_view m = gsl_matrix_view_array(dummy_matrix,(NComponents-1),(NComponents-1));
  gsl_matrix_view inv = gsl_matrix_view_array(dummy_inv_matrix,(NComponents-1),(NComponents-1));
  gsl_permutation * p = gsl_permutation_alloc((NComponents-1));
  gsl_linalg_LU_decomp(&m.matrix,p,&s);
  gsl_linalg_LU_invert(&m.matrix,p,&inv.matrix);
  gsl_permutation_free(p);
*/
  for(i=0;i<NComponents-1;i++){
   for(j=0;j<(NComponents-1);j++)derivative_c_matrix[derivative_Location(k,i,j,index_i,index_j)] = dummy_inv_matrix[j*(NComponents-1) + i];
  }
 }
 if(index_i == track_x && index_j == track_y){
  for(i=0;i<NComponents-1;i++){
   for(j=0;j<NComponents-1;j++){
    printf("Derivative _mu %d %d= %lf Derivative_C = %lf \t",i,j,derivative_mu_matrix[derivative_Location(0,i,j,index_i,index_j)],derivative_c_matrix[derivative_Location(0,i,j,index_i,index_j)]);
   }
   printf("\n");
  }
// printf("Composition  = %lf\n",Composition('g',NComponents-1,mu,temp));
 }
}




void Superimpose(int BC, double time){
/************Superimpose Potential*********/

 extern float *Mu_Matrix,*New_Mu_Matrix;
 extern struct node *Matrix;
 
 /****************************/
 int i=0,j=0,k=0;
 extern int Grid_x,Grid_y,N_Cut,NComponents;
 if(time!= 0.0){
  for(j=1;j<Grid_y-1;j++){
   for(i=1;i<Grid_x-1;i++){
    for(k=0;k<NComponents-1;k++){
     Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i];
     New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = 0.0;
    }
    for(k=0;k<N_Cut;k++){
     Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = blank_matrix;
     if(New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag != 0) Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i];
			/*********************Clear New_Matrix ******************/
     New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i]= blank_matrix;
    }
   }
  }
 }
 
 if(BC == 0){ //PBC
  for(i=1;i<Grid_x-1;i++){
   for(k=0;k<NComponents-1; k++){
    Mu_Matrix[k*Grid_x*Grid_y + (0)*Grid_x + i] =  Mu_Matrix[k*Grid_x*Grid_y + (Grid_y-2)*Grid_x + i];
    Mu_Matrix[k*Grid_x*Grid_y + (Grid_y-1)*Grid_x + i] =  Mu_Matrix[k*Grid_x*Grid_y + (1)*Grid_x + i];
   }
   for(k=0;k<N_Cut;k++){
    Matrix[k*Grid_x*Grid_y + (0)*Grid_x + i] =  Matrix[k*Grid_x*Grid_y + (Grid_y-2)*Grid_x + i];
    Matrix[k*Grid_x*Grid_y + (Grid_y-1)*Grid_x + i] =  Matrix[k*Grid_x*Grid_y + (1)*Grid_x + i];
   }
  }
  for(j=1;j<Grid_y-1;j++){
   for(k=0;k<NComponents-1; k++){
    Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + 0] =  Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + (Grid_x-2)];
    Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + (Grid_x-1)] =  Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + 1];
   }
   for(k=0;k<N_Cut;k++){
    Matrix[k*Grid_x*Grid_y + (j)*Grid_x + 0] =  Matrix[k*Grid_x*Grid_y + (j)*Grid_x + (Grid_x-2)];
    Matrix[k*Grid_x*Grid_y + (j)*Grid_x + (Grid_x-1)] = Matrix[k*Grid_x*Grid_y + (j)*Grid_x + 1];
   }
  }
 }
 if(BC == 1){ //No Flux
  for(i=1;i<Grid_x-1;i++){
   for(k=0;k<NComponents-1; k++){
    Mu_Matrix[k*Grid_x*Grid_y + (0)*Grid_x + i] =  Mu_Matrix[k*Grid_x*Grid_y + (1)*Grid_x + i];
    Mu_Matrix[k*Grid_x*Grid_y + (Grid_y-1)*Grid_x + i] =  Mu_Matrix[k*Grid_x*Grid_y + (Grid_y - 2)*Grid_x + i];
   }
   for(k=0;k<N_Cut;k++){
    Matrix[k*Grid_x*Grid_y + (0)*Grid_x + i] =  Matrix[k*Grid_x*Grid_y + (1)*Grid_x + i];
    Matrix[k*Grid_x*Grid_y + (Grid_y-1)*Grid_x + i] =  Matrix[k*Grid_x*Grid_y + (Grid_y - 2)*Grid_x + i];
   }
  }
  for(j=1;j<Grid_y-1;j++){
   for(k=0;k<NComponents-1; k++){
    Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + 0] =  Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + (1)];
    Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + (Grid_x-1)] =  Mu_Matrix[k*Grid_x*Grid_y + (j)*Grid_x + Grid_x-2];
   }
   for(k=0;k<N_Cut;k++){
    Matrix[k*Grid_x*Grid_y + (j)*Grid_x + 0] =  Matrix[k*Grid_x*Grid_y + (j)*Grid_x + 1];
    Matrix[k*Grid_x*Grid_y + (j)*Grid_x + (Grid_x-1)] = Matrix[k*Grid_x*Grid_y + (j)*Grid_x + Grid_x-2];
   }
  }
 }
}


void MatrixSetup(int BC, double time){

 int i=0,j=0;
 int flag = 1;
 extern int Grid_x,Grid_y;
 
 for(i=1;i<Grid_x-1;i++){
  for(j=1;j<Grid_y-1;j++){
   if(time != 0.0) flag = CheckMuDiff(i,j);
   if(flag == 1) Calculate_derivatives(i,j,time);
  }
 }
 Superimpose(BC,time);
}
