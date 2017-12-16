#include"function.h"
#include<stdio.h>
#include<stdlib.h>
#include"Global_params.h"

float AtpertoWtper(float c[2], int component){
 extern float *Molecular_Weight;
 int i=0,j=0;
 float sum = 0,comp_sum = 0;
 extern int NComponents;
 for(i=0;i<NComponents-1;i++){
  sum += Molecular_Weight[i]*c[i];
  comp_sum += c[i];
 }
 sum = sum + (1.0 - comp_sum)*Molecular_Weight[NComponents-1];
 return (100*Molecular_Weight[component]*c[component]/sum);
}

void Visualization(int t){
 printf("Output in file for t = %d\n",t);
 extern struct node *New_Matrix;
 extern float *New_Mu_Matrix,*Comp;
 
 extern int Grid_x,Grid_y,N_Cut,NComponents,NPhases;
 int i,j,k,l,cell,count; 
 char phase_filename[50],phase[50],ind_phase_alpha[50],ind_phase_gamma[50],ind_phase_cem[50],line_filename[50],comp[50],grain_file[50];
 sprintf(ind_phase_alpha,"outphase_alpha_%d.txt",t);
 sprintf(ind_phase_gamma,"outphase_gamma_%d.txt",t);
 sprintf(ind_phase_cem,"outphase_cem_%d.txt",t);
 sprintf(line_filename,"line_%d.txt",t);
 sprintf(comp,"comp_%d.txt",t);
 sprintf(phase,"phase_%d.txt",t);
 sprintf(grain_file,"grain_%d.txt",t);
 FILE *fo_alpha = fopen(ind_phase_alpha,"w");
 FILE *fo_gamma = fopen(ind_phase_gamma,"w");
 FILE *fo_cem = fopen(ind_phase_cem,"w");
 FILE *fout = fopen(phase,"w");
 FILE *fp = fopen(line_filename,"w");
 FILE *fo_comp = fopen(comp,"w");
 FILE *fo_grain = fopen(grain_file,"w");
/********************************************/
/*****************Differentiate between grain boundaries and different phase in single 2D plot ***********/
 double grain_val = 0.0;
 for(i=1;i<Grid_x-1;i++){
  for(j=1;j<Grid_y-1;j++){
   count = 0;
   grain_val = 0.0;
   for(k=0;k<N_Cut;k++){
    cell = k*Grid_x*Grid_y + j*Grid_x + i;
    if(New_Matrix[cell].grainnum != 0) count = count + 1;
    grain_val = grain_val + New_Matrix[cell].phi*10*New_Matrix[cell].grainnum;
   }
   if(count > 1) fprintf(fout,"%d %d 3\n",i,j);
   else{
    cell = 0*Grid_x*Grid_y +  j*Grid_x + i;
    if(New_Matrix[cell].phase == 'a') fprintf(fout,"%d %d 1\n",i,j);
    if(New_Matrix[cell].phase == 'g') fprintf(fout,"%d %d 2\n",i,j);
    if(New_Matrix[cell].phase == 'p') fprintf(fout,"%d %d 0\n",i,j);
   }
   fprintf(fo_grain,"%d %d %lf\n",i,j,grain_val);
  }	
  fprintf(fout,"\n");
  fprintf(fo_grain,"\n");	
 }
 fclose(fout);
 fclose(fo_grain);
/**********************************************************************************************************/
/********Writing individual phi contribution of phases in a file********************/
 double dummy_phi_alpha=0.0,dummy_phi_gamma=0.0,dummy_phi_cem=0.0;
 extern int track_y;
 int n=0;
 for(i=1;i<Grid_x-1;i++){
  for(j=1;j<Grid_y-1;j++){
   //for(n=0;n<NComponents-1;n++) printf("%lf ",Comp[n*Grid_x*Grid_y + j*Grid_x + i]);
   dummy_phi_alpha = 0.0; 
   dummy_phi_gamma = 0.0;
   dummy_phi_cem = 0.0;
   for(k=0;k<N_Cut;k++){
    cell = k*Grid_x*Grid_y + j*Grid_x + i;
    if(New_Matrix[cell].phase == 'a') dummy_phi_alpha = dummy_phi_alpha + New_Matrix[cell].phi;
    if(New_Matrix[cell].phase == 'g') dummy_phi_gamma = dummy_phi_gamma + New_Matrix[cell].phi;
    if(New_Matrix[cell].phase == 'p') dummy_phi_cem = dummy_phi_cem + New_Matrix[cell].phi;
   }
   fprintf(fo_alpha,"%d %d %lf\n",i,j,dummy_phi_alpha);
   fprintf(fo_gamma,"%d %d %lf\n",i,j,dummy_phi_gamma);
   fprintf(fo_cem,"%d %d %lf\n",i,j,dummy_phi_cem);
   fprintf(fo_comp,"%d %d %lf ",i,j,Comp[0*Grid_x*Grid_y + j*Grid_x + i]);
   for(k=1;k<NComponents-1;k++) fprintf(fo_comp,"%lf ",Comp[k*Grid_x*Grid_y + j*Grid_x + i]);
   fprintf(fo_comp,"\n");
		
		/******Writing along a slice******/
   if(j==track_y){
    fprintf(fp,"%d %lf %lf %lf ",i,dummy_phi_alpha,dummy_phi_gamma,dummy_phi_cem); 
    for(k=0;k<(NComponents-1);k++){
     cell = k*Grid_x*Grid_y + j*Grid_x + i;
     fprintf(fp,"%lf %lf ",New_Mu_Matrix[cell],Comp[cell]);
    }
    fprintf(fp,"\n");	
   }
		/*********************************/

  }
  fprintf(fo_alpha,"\n");
  fprintf(fo_gamma,"\n");
  fprintf(fo_cem,"\n");
  fprintf(fo_comp,"\n");
 }
 fclose(fo_alpha);
 fclose(fo_gamma);
 fclose(fo_cem);
 fclose(fp);
 fclose(fo_comp);
/********************************************************************************/
/****Average Properties*****************************************/
 extern char phase_array[3];
 double dummy_comp[NComponents-1],dummy_comp_gamma[NComponents-1],dummy_comp_alpha[NComponents-1];
 double total_comp[(NComponents-1)],total_comp_alpha[NComponents-1],total_comp_gamma[NComponents-1];
 int key_phase,flag;
 int grain_num_ferrite[1000],grain_num_austenite[1000];
 double vol_fraction[NPhases];
 struct node dummy_list[N_Cut];
 float dummy_mu_list[NComponents-1];
 char phase_id;
 extern float T_Hold;
 float temp,total_vf = 0.0,c[NComponents-1];
 extern float dt_real;
 float Effective_vf[NPhases];
 FILE *avg_out = fopen("average_output.txt","a");
 if(t==0) fprintf(avg_out,"Time VF_ferrite VF_austenite Composition_C_atper Composition_Mn_atper, Composition_C_wtper, Composition_Mn_wtper, Composition_C_alpha_avg_atper, Composition_Mn_alpha_avg_atper, Composition_C_gamma_avg_atper,  Composition_Mn_gamma_avg_atper, Composition_C_alpha_avg_wtper, Composition_Mn_alpha_avg_wtper, Composition_C_gamma_avg_wtper,  Composition_Mn_gamma_avg_wtper\n");
fprintf(avg_out,"%lf,",t*dt_real);
 temp = T_Hold;
 for(i=0;i<NPhases;i++) vol_fraction[i] = 0.0;
 /***Volume Fraction*/
 for(i=1;i<Grid_x-1;i++){
  for(j=1;j<Grid_y-1;j++){
   for(l=0;l<N_Cut;l++){
    cell = l*Grid_x*Grid_y  + j*(Grid_x) + i;
    for(k=0;k<NPhases;k++){
     if(New_Matrix[cell].phase == phase_array[k]) vol_fraction[k] += New_Matrix[cell].phi/((Grid_x-2)*(Grid_y-2));
    }
   }
  }
 }
 for(i=0;i<NPhases;i++){
  fprintf(avg_out,"%lf,",vol_fraction[i]);
  total_vf += vol_fraction[i];
 }
 fprintf(avg_out,"%lf,",total_vf);
 
 //fprintf(avg_out,"\n");
 /*******************/
 for(l=0;l<NComponents-1;l++){
  total_comp[l] = 0.0;
  dummy_comp[l] = 0.0;
  dummy_comp_alpha[l] = 0.0;
  dummy_comp_gamma[l] = 0.0;
 }

 for(i=1;i<Grid_x-1;i++){
  for(j=1;j<Grid_y-1;j++){
   for(l=0;l<N_Cut;l++) dummy_list[l] = New_Matrix[Location(i,j,l)];
   for(l=0;l<NComponents-1;l++) dummy_mu_list[l] = New_Mu_Matrix[Location(i,j,l)]; 
   for(l=0;l<N_Cut;l++){ 
    cell =  l*Grid_x*Grid_y + j*Grid_x + i;
    if(Matrix[cell].flag !=0){
     phase_id = New_Matrix[cell].phase;
     for(k=0;k<NPhases;k++){
      if(phase_id == phase_array[k]) key_phase = k;
     }
     for(k=0;k<NComponents-1;k++){
       dummy_comp[k] += Composition(phase_id,k,dummy_mu_list,temp)*H(l,N_Cut,dummy_list);
       if(phase_id == 'a') dummy_comp_alpha[k] += Composition(phase_id,k,dummy_mu_list,temp)*H(l,N_Cut,dummy_list);
       if(phase_id == 'g') dummy_comp_gamma[k] += Composition(phase_id,k,dummy_mu_list,temp)*H(l,N_Cut,dummy_list);        
     }
    }
   }
   for(k=0;k<NComponents-1;k++){
    total_comp[k] += dummy_comp[k]; 
    dummy_comp[k] = 0.0;
   }
  }
 }
 for(i=0;i<NPhases;i++) Effective_vf[i] = (Grid_x-2)*(Grid_y-2)*vol_fraction[i];
 //for(i=0;i<NComponents-1;i++) fprintf(avg_out,"%lf ",dummy_comp_alpha[i]/Effective_vf[0]);
 //for(i=0;i<NComponents-1;i++) fprintf(avg_out,"%lf ",dummy_comp_gamma[i]/Effective_vf[1]);
 c[0] = dummy_comp_gamma[0]/Effective_vf[1]*vol_fraction[1] + (dummy_comp_alpha[0]/Effective_vf[0])*vol_fraction[0];
 c[1] = dummy_comp_gamma[1]/Effective_vf[1]*vol_fraction[1] + (dummy_comp_alpha[1]/Effective_vf[0])*vol_fraction[0];
 for(i=0;i<NComponents-1;i++) fprintf(avg_out,"%lf,",(dummy_comp_gamma[i]/Effective_vf[1])*vol_fraction[1] + (dummy_comp_alpha[i]/Effective_vf[0])*vol_fraction[0]);
 for(i=0;i<NComponents-1;i++) fprintf(avg_out,"%lf,",AtpertoWtper(c,i));
 c[0] =  dummy_comp_alpha[0]/Effective_vf[0];
 c[1] =  dummy_comp_alpha[1]/Effective_vf[0];

 fprintf(avg_out,"%lf,%lf,",dummy_comp_alpha[0]/Effective_vf[0],dummy_comp_alpha[1]/Effective_vf[0]);
 for(i=0;i<NComponents-1;i++) fprintf(avg_out,"%lf,",AtpertoWtper(c,i));
 c[0] =  dummy_comp_gamma[0]/Effective_vf[1];
 c[1] =  dummy_comp_gamma[1]/Effective_vf[1];
 fprintf(avg_out,"%lf,%lf,",dummy_comp_gamma[0]/Effective_vf[1],dummy_comp_gamma[1]/Effective_vf[1]);
 for(i=0;i<NComponents-1;i++) fprintf(avg_out,"%lf,",AtpertoWtper(c,i));
/*******************************************************************/
/*********Grain Sizes************************************************/



 fprintf(avg_out,"\n");
 fclose(avg_out);

 
/****************************************************************/


}
