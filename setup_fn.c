#include"Global_params.h"
#include"setup_fn.h"
#include<math.h>

void Calc_Mu_Seg(char phase, float c[3], float *m, float temp){
 extern float **Fit_Constant;
 extern float T_Hold;
 extern float R; 
 int i=0;
 float sum =0;
 extern int NComponents;
 R = 8.314;
 for(i=0;i<3;i++) printf("C : %lf\n", c[i]);
 if(phase == 'a'){
 for(i=0;i<NComponents-1;i++) m[i] = (-1)*Fit_Constant[NComponents-1][0]*((temp/T_Hold) - Fit_Constant[NComponents-1][1])/Fit_Constant[NComponents-1][1] + Fit_Constant[i][0]*((T_Hold/T_Hold) - Fit_Constant[i][1])/Fit_Constant[i][1] + (temp/T_Hold)*log(c[i]/c[2]);
 }
 else if(phase == 'g'){
  for(i=0;i<NComponents-1;i++){
   m[i] = (temp/T_Hold)*log(c[i]/(c[2]));
   printf("%f %d %f",R,i,m[i]);
  }
 }

}

void Wt2Mole(float *Bulk_C,float *Molecular_Weight,float *Comp_Mole){
 extern int NComponents;
 int i=0,j=0;
 float sum =0.0;
 for(i=0;i<NComponents;i++) sum += Bulk_C[i]/Molecular_Weight[i];
 for(i=0;i<NComponents;i++){
  Comp_Mole[i] = (Bulk_C[i]/Molecular_Weight[i])/sum;
  printf("%f ",Comp_Mole[i]);
  printf("SUM = %lf\n",sum);
 }
 printf("\n");
}

void Scaling(){

 extern float E_scaling,L_scaling,t_scaling;
 extern float T_Hold,Molar_Volume,Surface_Energy,D_Max,Time_Sec,T_Init;
 extern float Epsilon;
 extern float dx_nd,dy_nd,dt_nd;
 extern float dx_real,dt_real;
 extern int N_Timesteps,Grid_x,Grid_y;
 extern float Len_x,Len_y;

 float K = 200.0;

 E_scaling = 8.314*T_Hold/Molar_Volume;
 L_scaling = Surface_Energy/E_scaling;
 t_scaling = L_scaling*L_scaling/D_Max;


 printf("Energy Scaling = %e\n",E_scaling);
 printf("Length Scaling = %e\n",L_scaling);
 printf("time Scaling = %e\n",t_scaling);

 dx_nd = (22.0/7.0)*(22.0/7.0)*1.414*Epsilon/(4.0*10.0);
 dy_nd = (22.0/7.0)*(22.0/7.0)*1.414*Epsilon/(4.0*10.0);
 dx_real = dx_nd*L_scaling;

 dt_nd = dx_nd*dx_nd/K;
 dt_real = dt_nd*t_scaling;

 N_Timesteps =(int)(Time_Sec/dt_real);

 Len_x = Grid_x*dx_real;
 Len_y = Grid_y*dx_real;

 T_Init /= T_Hold;
 T_Hold /= T_Hold;
 printf("Domain Length = %e X %e\n",Len_x,Len_y);
 printf("Number of timesteps for given time of %f : %d\n",Time_Sec,N_Timesteps);
 printf("dt = %e, dx = %e\n",dt_real,dx_real);
 printf("Non dimensional : dt = %f, dx = %f\n",dt_nd,dx_nd);
 printf("Holding Temperature : %f, Initial Temperature : %f\n",T_Hold,T_Init);
}


void Mu_Calc(float *mu,float T){
 extern int NComponents;
 extern float **Fit_Constant;
 extern float *Comp_Initial,*Comp_Hold;
 double sum = 0.0;

 int i=0,j=0;
 if(T==T_Init){
  for(i=0;i<(NComponents-1);i++) sum += Comp_Initial[2*i+1];
  for(i=0;i<(NComponents-1);i++){
   mu[i] = (T)*log(Comp_Initial[2*i+1]/(1.0 - sum));
  }
  for(j=0;j<NComponents-1;j++) printf("C : %lf\n",Comp_Initial[2*j+1]);
  printf("C : %lf\n",1.0-sum); 
  for(j=0;j<NComponents-1;j++)printf("Initial Equilibrium Potential: Gamma:  mu_%d = %f\n",j,mu[j]);

  sum = 0.0;
  for(i=0;i<(NComponents-1);i++) sum += Comp_Initial[2*i];
  for(i=0;i<(NComponents-1);i++){
   mu[i] = (1.0)*(Fit_Constant[i][0]*(T - Fit_Constant[i][1])/Fit_Constant[i][1] - Fit_Constant[NComponents-1][0]*(T - Fit_Constant[NComponents-1][1])/Fit_Constant[NComponents-1][1] + T*log(Comp_Initial[2*i]/(1.0 - sum)));
  }
  for(j=0;j<NComponents-1;j++)printf("C : %lf\n",Comp_Initial[2*j]);
  printf("C : %lf\n",1.0-sum);
  for(j=0;j<NComponents-1;j++) printf("Initial Equilibrium Potential: Alpha :mu_%d = %f\n",j,mu[j]);
   
 }
 else if(T==T_Hold){
  sum = 0.0;
  for(i=0;i<(NComponents-1);i++) sum += Comp_Hold[2*i+1];
  for(i=0;i<(NComponents-1);i++){
   mu[i] = (T)*log(Comp_Hold[2*i+1]/(1.0 - sum));
   printf("Holding Equilibrium Potential: Gamma:  mu_%d = %f\n",i,mu[i]);
  }

  sum = 0.0;
  for(i=0;i<(NComponents-1);i++) sum += Comp_Hold[2*i];
  for(i=0;i<(NComponents-1);i++){
   mu[i] = (1.0)*(Fit_Constant[i][0]*(T - Fit_Constant[i][1])/Fit_Constant[i][1] - Fit_Constant[NComponents-1][0]*(T - Fit_Constant[NComponents-1][1])/Fit_Constant[NComponents-1][1] + T*log(Comp_Hold[2*i]/(1.0 - sum)));
   printf("Holding Equilibrium Potential: Alpha :mu_%d = %f\n",i,mu[i]);
  }
 }
}

void Mu_Calc_Alpha(float *Comp_Mole,float *mu,float T){
int i = 0;
double sum=0.0;
double add=0.0;
extern int NComponents;
extern float **Fit_Constant;
sum = Fit_Constant[i][0]*(T - Fit_Constant[i][1])/Fit_Constant[i][1] - Fit_Constant[NComponents-1][0]*(T - Fit_Constant[NComponents-1][1])/Fit_Constant[NComponents-1][1];

for(i=0;i<(NComponents-1);i++) add += Comp_Hold[2*i];
for(i=0;i<NComponents-1;i++) mu[i] = T*log(Comp_Mole[2*i]/(1-add)) + sum;
}

void Mu_Calc_Beta(float *Comp_Mole,float *mu,float T){
int i = 0;
double add=0.0;
extern int NComponents;

for(i=0;i<(NComponents-1);i++) add += Comp_Hold[2*i+1];
for(i=0;i<NComponents-1;i++) mu[i] = T*log(Comp_Mole[2*i+1]/(1-add));
}


void PopulateSeeds(int seeds,char phase,int totalgrains){
int k=0,i=0,j=0,l=0;
extern struct node *Matrix;
extern float *mu_initial,*mu_bulk;
double radius = 10.0;
int rx,ry;
int cell;
/* For unequal probability of nucleation in Mn rich and Mn lean regions */
///*
double prob = 0.3;
int flag = 0,nuc_flag = 0;
double random;

//*/
extern int Grid_x,Grid_y;
for(k=0;k<seeds;k++)
    {
    do
     {
     nuc_flag = 0;
     rx = (gsl_rng_uniform(URN))*Grid_x;
     ry = (gsl_rng_uniform(URN))*Grid_y;
     random = gsl_rng_uniform(URN);
     if(((j>0.2*Grid_x-1 && j<0.4*Grid_x-1) || (j>0.6*Grid_x-1 && j<0.8*Grid_x-1))) flag = 1;
     else flag = 0;
     if(flag == 1 && random<=prob) nuc_flag = 1;
     if(flag == 0 && random>=prob) nuc_flag = 1;	
     }while(Matrix[0*Grid_x*Grid_y + ry*Grid_x + rx].phase == 'a' || nuc_flag == 0);
     for(i=0;i<Grid_x;i++)
         {
         for(j=0;j<Grid_y;j++)
             {
             if((i-rx)*(i-rx) + (j-ry)*(j-ry) <= radius*radius)
                {
                cell = 0*Grid_x*Grid_y + j*Grid_x + i;
                Matrix[cell].grainnum  = totalgrains + k + 1;
                Matrix[cell].phi = 1.0;
                Matrix[cell].phase = phase;
                Matrix[cell].flag = 1;
                for(l=0;l<(NComponents-1);l++) Mu_Matrix[l*Grid_x*Grid_y + j*Grid_x + i] = mu_initial[l];//to change in case of austenite
                }
	     if(j>0.2*Grid_x-1 && j<0.4*Grid_x-1) Mu_Matrix[1*Grid_x*Grid_y + j*Grid_x + i] = mu_bulk[1];	
             if(j>0.6*Grid_x-1 && j<0.8*Grid_x-1) Mu_Matrix[1*Grid_x*Grid_y + j*Grid_x + i] = mu_bulk[1];	
             }
          }
     }
}

void ProbabilityNucleation(char phase,int totalgrains){
 extern struct node *Matrix;
 extern float *mu_initial;
 double radius = 7.0;
 int i=0,j=0,k=0,l=0,m=0;
 extern int Grid_x,Grid_y;

 int interfacepoints[Grid_x*Grid_y];
 for(i=0;i<Grid_x*Grid_y-1;i++) interfacepoints[i] = 0;
 int cell_n[4],cell;
 for(i=1;i<Grid_x-1;i++){
  for(j=1;j<Grid_y-1;j++){
   if(interfacepoints[j*Grid_x + i] != 1){
    cell = 0*Grid_x*Grid_y + j*Grid_x + i;
    cell_n[0] = 0*Grid_x*Grid_y + j*Grid_x + i+1;
    cell_n[2] = 0*Grid_x*Grid_y + j*Grid_x + i-1;
    cell_n[3] = 0*Grid_x*Grid_y + (j+1)*Grid_x + i;
    cell_n[1] = 0*Grid_x*Grid_y + (j-1)*Grid_x + i;
    for(k=0;k<4;k++){
     if(Matrix[cell].grainnum != Matrix[cell_n[k]].grainnum){
      interfacepoints[j*Grid_x + i] = 1;
      interfacepoints[cell_n[k]] = 1;
      break;
     }
    }
   }
  }
 }
/*
for(i=0;i<Grid_x;i++)
 for(j=0;j<Grid_y;j++)
  if(interfacepoints[j*Grid_x+i]==1) printf("interface %d %d %d\n",i,j,Matrix[j*Grid_x+i].grainnum);
*/

/******Look for some justification for probability of nucleation*******************/

 double prob_ferrite = 0.002,prob_pearlite=0.005,prob_austenite = 0.005,prob;
 double rand;
 extern float *mu_seg_low;
 int count=1;
 int x,y;
 int pass = 0;
 for(i=1;i<Grid_x-1;i++){
  for(j=1;j<Grid_y-1;j++){
   pass = 0;
   if(interfacepoints[j*Grid_x+i] == 1 && (Matrix[j*Grid_x+i].phase == 'a') && phase == 'g'){
    prob = prob_ferrite;
    pass = 1;
   }
   else if(interfacepoints[j*Grid_x + i] ==1 && Matrix[j*Grid_x+i].phase =='p' && phase == 'g'){
    prob = prob_pearlite;
    pass = 1;
   }
   else if(interfacepoints[j*Grid_x + i] == 1 && Matrix[j*Grid_x+i].phase == 'g' && phase == 'a'){
    prob = prob_austenite;
    pass = 1;
   }
   else pass = 0;
   if(pass == 1){
    rand = (gsl_rng_uniform(URN));
    if(rand<=prob){
     for(l=0;l<Grid_x;l++){
      for(m=0;m<Grid_y;m++){
                /*
 *  *                 if((l-i) > Lx/2) x = l-Lx/2;
 *   *                                 if((i-l) > Lx/2) x = l+Lx/2;
 *    *                                                 if((m-j) > Ly/2) y = m-Ly/2;
 *     *                                                                 if((j-m) > Ly/2) y = m+Ly/2;
 *      *                                                                                 */
       if((l-i)*(l-i) + (m-j)*(m-j) <= radius*radius){
        cell = 0*Grid_y*Grid_x + m*Grid_x + l;
        Matrix[cell].grainnum = totalgrains + count;
        Matrix[cell].phi = 1.0;
        Matrix[cell].phase = phase;
        Matrix[cell].flag = 1.0;
        for(k=0;k<NComponents-1;k++) Mu_Matrix[k*Grid_x*Grid_y + m*Grid_x + l] = mu_seg_low[k];
       }
      }
     }
     count = count + 1;
     printf("Total Grains = %d\n",count+totalgrains);
    }
   }
  }
 }
}

