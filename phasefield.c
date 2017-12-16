#include"stdio.h"
#include"function.h"
#include"Global_params.h"

int main(){
 
 ReadInputs(); 

 Allocate();

 SystemSetup();

 extern int Grid_x,Grid_y;
 extern float RUNTIME,Time_Sec,dt_real;
 int t=0,i=0,j=0;
 printf("Gridx = %d, Gridy = %d \n",Grid_x,Grid_y);
 RUNTIME = 0.0;


 while(RUNTIME<Time_Sec){
  if(t%1 == 0)printf("Timestep = %d\n",t);
  MatrixSetup(0,RUNTIME);/*Populate Values of Phi, Mu, dmu/d(c), dc/d(mu)*/
  for(i=1;i<(Grid_x-1);i++){
   for(j=1;j<(Grid_y-1);j++){
    StencilGeneration(i,j);
    PFSolver(i,j,t);
    Updation(i,j);
   }
  }
//      Superimpose(0); /*Copying calculated values to updated matrix*/
  if(t%(10000)==0) Visualization(t);
      //if((t+1)%100 == 0) Volume_Fraction_Austenite(t);
  t = t + 1;
  RUNTIME =RUNTIME + dt_real;
  printf("Runtime = %lf,total time = %lf,dt = %lf\n",RUNTIME,Time_Sec,dt_real);
 }
 FreeMemory();
 return(0);
}

