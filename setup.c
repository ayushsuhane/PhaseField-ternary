#include"function.h"
#include"Global_params.h"
//#include"setup_fn.h"
#include<math.h>
#include<stdlib.h>

void SystemSetup(){

extern float Aspect_Ratio,
	      *mu_initial,
	      *mu_equilibrium,
              *mu_bulk;
extern int Grid_x,Grid_y,  
	   Num_Grains,
	   Seeds;
extern float *mu_seg_low,*mu_seg_high;
extern struct node *Matrix;
UniformRNGen();
int option;
int i=0,j=0;
option = 1;
int cell,l,k;


 if(option == 1) { //FlatInterface
  printf("Creating a flat profile with three phases\n");
  for(i=0;i<Grid_x;i++){
   for(j=0;j<Grid_y;j++){
    if(i<(0.375*Grid_x-1)){
        //if(j<(0.5*Ly-1))
        //if(i<=0.5*(Lx) && j<=((0.5*Ly)))
     cell = 0*Grid_x*Grid_y + j*Grid_x + i;
     if(Matrix[cell].flag == 0){
      Matrix[cell].grainnum  = 1;
      Matrix[cell].phi = 1.0;
      Matrix[cell].phase = 'a';
      Matrix[cell].flag = 1;
      for(l=0;l<(NComponents-1);l++) Mu_Matrix[l*Grid_x*Grid_y + j*Grid_x + i] = mu_initial[l];
     }
    }
     //else
    else if(i>(0.625*Grid_x-1)){
      //else if((i>0.5*(Lx)) && (j<(0.5*Ly)))
     cell = 0*Grid_x*Grid_y + j*Grid_x + i;
     if(Matrix[cell].flag == 0){
      Matrix[cell].grainnum  = 1;
      Matrix[cell].phi = 1.0;
      Matrix[cell].phase = 'a';
      Matrix[cell].flag = 1;
      for(l=0;l<(NComponents-1);l++) Mu_Matrix[l*Grid_x*Grid_y + j*Grid_x + i] = mu_initial[l];
     }
    }
    else{
     cell = 0*Grid_x*Grid_y + j*Grid_x + i;
     if(Matrix[cell].flag == 0){
      Matrix[cell].grainnum  = 2;
      Matrix[cell].phi = 1.0;
      Matrix[cell].phase = 'g';
      Matrix[cell].flag = 1;
      for(l=0;l<(NComponents-1);l++) Mu_Matrix[l*Grid_x*Grid_y + j*Grid_x + i] = mu_initial[l];
     }
    }
//    if(j>0.5*Grid_x-1) Mu_Matrix[1*Grid_x*Grid_y + j*Grid_x + i] = mu_seg_high[1];	
   }
  }
 }
if(option == 2) //Spherical
	{
	for(i=0;i<Grid_x;i++)
                {
                for(j=0;j<Grid_y;j++)
                        {
                        if(((i-(Grid_x/2))*(i-Grid_x/2) + (j-(Grid_y/2))*(j-(Grid_y/2))) <= 100.0 )
                                {
                                cell = 0*Grid_x*Grid_y + j*Grid_x + i;
                                if(Matrix[cell].flag == 0)
                                        {
                                        Matrix[cell].grainnum  = 1;
                                        Matrix[cell].phi = 1.0;
                                        Matrix[cell].phase = 'a';
                                        Matrix[cell].flag = 1;
                                        for(l=0;l<(NComponents-1);l++) Mu_Matrix[l*Grid_x*Grid_y + j*Grid_x + i] = mu_initial[l];
                                        }
                                }
                        else
                                {
                                cell = 0*Grid_x*Grid_y + j*Grid_x + i;
                                if(Matrix[cell].flag == 0)
                                        {
                                        Matrix[cell].grainnum  = 2;
                                        Matrix[cell].phi = 1.0;
                                       // Matrix[cell].phase = 'g';
                                        Matrix[cell].phase = 'g';
                                        Matrix[cell].flag = 1;
                                        for(l=0;l<(NComponents-1);l++) Mu_Matrix[l*Grid_x*Grid_y + j*Grid_x + i] = mu_initial[l];
                                        //for(l=0;l<NComponents;l++) Mu_Matrix[l*Lx*Ly + j*Lx + i] = mu_initial[l];
                                        }

                                }

                        }
		}

	}

 double radius = 11.0;
 int rx=0,ry=0;
 int nx[Num_Grains],ny[Num_Grains];
 double dmin=0.0,d;
 extern double volfrac_pearlite;
 double vf=0.0;
 int grain=0;
 double xdis=0,ydis=0;
 extern float Aspect_Ratio; 
 if(option == 6){ //Voronoi
  printf("Polycrystal with three phases\n"); 
/*****************************************************************************/
/*******************Voronoi**************************************************/
  for(i=0;i<Num_Grains;i++) {
   nx[i] = (gsl_rng_uniform(URN))*Grid_x;
   ny[i] = (gsl_rng_uniform(URN))*Grid_y;
  }
  for(i=0;i<Grid_y;i++){
   for(j=0;j<Grid_x;j++){
    dmin = pow(((Aspect_Ratio)*(Aspect_Ratio)*(Grid_x-1)*(Grid_x-1) + (Grid_y-1)*(Grid_y-1)),0.5);
    l = -1;
    for(k=0;k<Num_Grains;k++){
     xdis = nx[k] - j;
     if(fabs(nx[k]-j)>Grid_x/2) xdis = Grid_x - fabs(nx[k] - j); 
     ydis = ny[k] - i;
     if(fabs(ny[k]-i)>Grid_y/2) ydis = Grid_y - fabs(ny[k] - i); 
     d = pow(((Aspect_Ratio)*(Aspect_Ratio)*(xdis)*(xdis) + (ydis)*(ydis)),0.5);
     if (d<dmin) {
      dmin = d;
      l = k;
     }
    }
    cell = 0*Grid_x*Grid_y + i*Grid_x + j;
    Matrix[cell].grainnum = l+1;
    Matrix[cell].phi = 1.0;
    Matrix[cell].phase = 'g';//to change in case of ferrite 
    Matrix[cell].flag = 1;
    for(l=0;l<NComponents;l++) Mu_Matrix[l*Grid_x*Grid_y + i*Grid_x + j] = mu_seg_low[l];
   }
  }
  printf("Constructed initial homogeneous microstructure\n");
	
  int count = 0;
  while(vf<volfrac_pearlite){
   grain = gsl_rng_uniform(URN)*Num_Grains;
   printf("Grain number selected = %d\n",grain);
   for(i=0;i<Grid_y;i++){
    count = 0;
    for(j=0;j<Grid_x;j++){
     cell = 0*Grid_x*Grid_y + i*Grid_x + j;
     if(Matrix[cell].grainnum==grain){
      if(Matrix[cell].phase =='p'){
       count = count + 1;
       break;
      }
      else{
       Matrix[cell].phase = 'p';
       vf = vf + 1.0/(Grid_x*Grid_y);
       printf("vf = %lf\n",vf);
      }
     }
    }
    if(count>0) break;
   }
  }
  printf("Full initial microstructuree Created\n");

/****************************************************************************/	
/****************************************************************************/
/***************Site Saturated Nucleation***********************************/
  ProbabilityNucleation('a',Num_Grains);
	//PopulateSeeds(Seeds,'a',Num_Grains);
	//PopulateSeeds(Seeds_Ferrite,'a',Num_Grains+Seeds);
  for(i=0;i<Grid_x;i++){
   for(j=0;j<Grid_y;j++){
    if(j>0.2*Grid_x-1 && j<0.4*Grid_x-1) Mu_Matrix[1*Grid_x*Grid_y + j*Grid_x + i] = mu_seg_high[1];	
    if(j>0.6*Grid_x-1 && j<0.8*Grid_x-1) Mu_Matrix[1*Grid_x*Grid_y + j*Grid_x + i] = mu_seg_high[1];	
   }
  }
 }
}
