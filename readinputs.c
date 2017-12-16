#include"function.h"
#include<stdio.h>
#include"Global_params.h"
#include<stdlib.h>

void ReadInputs(){


 extern int NPhases,NComponents,Grid_x,Grid_y,Num_Grains,Seeds,Seeds_Ferrite,N_Cut,Random_Seed;
 extern float *Molecular_Weight;
 extern char **Elements;
 extern float *Bulk_Composition;
 extern float Heating_Rate,Time_Sec,Aspect_Ratio,Molar_Volume,Strain_Energy,Surface_Energy,D_Max;
 extern float T_Init,T_Hold,Epsilon;
 extern float *mu_seg_low,*mu_seg_high;
 extern char phase_array[3]; 
 phase_array[0] = 'a';
 phase_array[1] = 'g';
 phase_array[2] = 'p';  
 

 int i=0,j=0,k=0;
 int flag_seg = 0;
 char phase_seg;
 float c_seg[3],c_seg_mole[3];

 FILE *input = fopen("input.txt","r");
 if(input == NULL){
  printf("Couldnot find input file");
 }


 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&NComponents);
 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&NPhases);
 fscanf(input,"%*[^\n]\n",NULL);

 Elements = (char**)malloc(NComponents*sizeof(char*)); //free in memory.c
 for(i=0;i<NComponents;i++){
  Elements[i] = (char*)malloc(2*sizeof(char));
 }
 Molecular_Weight = (float*)malloc(NComponents*sizeof(float));  //free in memory.c
 Bulk_Composition = (float*)malloc(NComponents*sizeof(float));  //free in memory.c

 for(i=0;i<NComponents;i++){
  fscanf(input,"%s %f\n",Elements[i],&Molecular_Weight[i]);
  printf("%s %f\n",Elements[i],Molecular_Weight[i]);
 }

 fscanf(input,"%*[^\n]\n",NULL);

 for(i=0;i<NComponents;i++){
  fscanf(input,"%s %f\n",Elements[i],&Bulk_Composition[i]);
  printf("%s %f\n",Elements[i],Bulk_Composition[i]);
 }

fscanf(input,"%*[^\n]\n",NULL);
fscanf(input,"%f\n",&Heating_Rate);
printf("Heating Rate:%f\n",Heating_Rate);


 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&T_Init);
 printf("T Initial:%f\n",T_Init);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&T_Hold);
 printf("Holding Temperature:%f\n",T_Hold);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&Time_Sec);
 printf("Time(sec) :%f\n",Time_Sec);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&Grid_x);
 fscanf(input,"%d\n",&Grid_y);
 printf("Grid Dimensions:%dx%d\n",Grid_x,Grid_y);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&Num_Grains);
 printf("Total Number of Grains:%d\n",Num_Grains);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&Vol_Fraction);
 printf("Volume fraction of pearlite :%f\n",Vol_Fraction);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&Seeds);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&Seeds_Ferrite);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&Aspect_Ratio);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&N_Cut);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&Molar_Volume);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&Epsilon);
 printf("Epsilon:%f\n",Epsilon);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&Strain_Energy);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&Surface_Energy);
 printf("Surface Energy:%f\n",Surface_Energy);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%f\n",&D_Max);
 printf("Maximum Diffusivity:%e\n",D_Max);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&Random_Seed);

 fscanf(input,"%*[^\n]\n",NULL);
 fscanf(input,"%d\n",&flag_seg);

 


/***********************************************************************************************/
/****Reading constants evaluated from free energy fitting************************************/
 FILE *constant = fopen("constants.inp","r");

 extern float **Fit_Constant; //Free in memory.c
 Fit_Constant = (float**)malloc(NComponents*sizeof(float*));
 for(i=0;i<NComponents;i++) Fit_Constant[i] = (float*)malloc(2*sizeof(float));

 for(i=0;i<NComponents;i++){
  for(j=0;j<2;j++)fscanf(constant,"%f ",&Fit_Constant[i][j]);
  fscanf(constant,"\n");
 }

/********Fit constants divide by T for scaling***********/
 for(i=0;i<NComponents;i++)
  for(j=0;j<NComponents;j++)
   Fit_Constant[i][j] /= T_Hold;
/*******************************************************/


 for(i=0;i<NComponents;i++){
  printf("%s ",Elements[i]);
  for(j=0;j<2;j++) printf("%f ",Fit_Constant[i][j]);
  printf("\n");
 }

 fclose(constant);

 if(flag_seg == 1){
  mu_seg_low = (float*)malloc((NComponents-1)*sizeof(float));
  mu_seg_high = (float*)malloc((NComponents-1)*sizeof(float));
  for(i=0;i<2;i++){
   fscanf(input,"%*[^\n]\n",NULL);
   fscanf(input,"%f %f %f %c\n",&c_seg[0],&c_seg[1],&c_seg[2],&phase_seg);
   Wt2Mole(c_seg,Molecular_Weight,c_seg_mole);
   if(i==0) Calc_Mu_Seg('g',c_seg_mole,mu_seg_low, T_Init);
   else Calc_Mu_Seg('g',c_seg_mole,mu_seg_high,T_Init);
  }
 }
 printf("\n");
 printf("Segregation mu :\n LOW - %f %f\n HIGH - %f %f\n",mu_seg_low[0],mu_seg_low[1],mu_seg_high[0],mu_seg_high[1]);
 fclose(input);
/*******************************************************************************************/
/*************************Diffusivity File************************************************/
 FILE *diffuse = fopen("diffusivity.inp","r");

 extern float **Diffusion_Matrix;  //free in memory.c
 Diffusion_Matrix = (float**)malloc(NPhases*sizeof(float*));
 for(i=0;i<NPhases;i++) Diffusion_Matrix[i] = (float*)malloc(((NComponents-1)*(NComponents-1))*sizeof(float));

 for(i=0;i<5;i++) fscanf(diffuse,"%*[^\n]\n",NULL);
 for(i=0;i<NPhases;i++){
  for(k=0;k<NComponents-1;k++){
   for(j=0;j<NComponents-1;j++) fscanf(diffuse,"%f ",&Diffusion_Matrix[i][k*(NComponents-1) + j]); 
    fscanf(diffuse,"\n");
  }
  fscanf(diffuse,"%*[^\n]\n",NULL);
 }

 for(i=0;i<NPhases;i++){
  printf("Phase : %d\n",i);
  for(k=0;k<(NComponents-1);k++){
   for(j=0;j<(NComponents-1);j++){
    printf("%f ",Diffusion_Matrix[i][k*(NComponents-1) + j]);
    printf("\n");
   }
  }
 }

/***********************************************************************/
 extern float *Comp_Initial; //Free in memory.c
 Comp_Initial = (float*)malloc(2*(NComponents-1)*sizeof(float));
 extern float *Comp_Hold; //Free in memory.c
 Comp_Hold = (float*)malloc(2*(NComponents-1)*sizeof(float));
 FILE *fc = fopen("init_comp.inp","r");
 for(i=0;i<(2*(NComponents-1));i++) fscanf(fc,"%e ",&Comp_Initial[i]);
 fscanf(fc,"\n");
 for(i=0;i<(2*(NComponents-1));i++) fscanf(fc,"%e ",&Comp_Hold[i]);
 fclose(fc);
 
 
 printf("Completed reading Files\n");
/**************************************************************************************************/

}
