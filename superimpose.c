#include"Global_params.h"
#include"function.h"


void Superimpose(int flag){
/************Superimpose Potential*********/

extern float *Mu_Matrix,*New_Mu_Matrix;
int i=0,j=0,k=0;
extern int Grid_x,Grid_y;
if(flag == 0)
   {
   for(k=0;k<(NComponents-1);k++)
       {
       for(j=0;j<Grid_y;j++)
	   {
	   for(i=0;i<Grid_x;i++) 
	       {
	       Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i];
	       New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = 0.0;
	       }
	   }
	}
   }
else if(flag == 1)
	{
	for(k=0;k<(NComponents-1);k++)
		{
		for(j=1;j<(Grid_y-1);j++)
		    {
		    for(i=0;i<Grid_x;i++) 
			{
			Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i];
			New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = 0.0;
			}
		    }
		j = 0;
		for(i=0;i<Grid_x;i++)
		    {
		    Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = Mu_Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i];
		    New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = 0.0;
		    }
		j = Grid_y-1;
		for(i=0;i<Grid_x;i++)
		    {
		    Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = Mu_Matrix[k*Grid_x*Grid_y + (j-1)*Grid_x + i];
		    New_Mu_Matrix[k*Grid_x*Grid_y + j*Grid_x + i] = 0.0;
		    }
		}
	}
/****************************************/
/**************Superimpose Matrix*********/

if(flag==0)
   {
   for(k=0;k<N_Cut;k++)
       {
       for(j=0;j<Grid_y;j++)
	   {
	   for(i=0;i<Grid_x;i++) 
	       {
	       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = 0;
	       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = 0.0;
	       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = 'n';
	       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = 0;
	       if(New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag != 0)
		  {
		  Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum;
		  Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi;
		  Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase;
		  Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag;
		  }
			/*********************Clear New_Matrix ******************/
	       New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = 0;
	       New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = 0.0;
	       New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = 'n';
	       New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = 0;
				
			/******************************************/
	       }
	   }
	}
   }
else if(flag == 1)
	{
	for(k=0;k<N_Cut;k++)
	    {
	    for(j=1;j<(Grid_y-1);j++)
		{
		for(i=0;i<Grid_x;i++) 
		    {
		    Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = 0;
		    Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = 0.0;
		    Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = 'n';
		    Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = 0;
		    if(New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag != 0)
		       {
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum;
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi;
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase;
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag;
		       }
			/*********************Clear New_Matrix ******************/
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = 0;
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = 0.0;
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = 'n';
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = 0;
			/******************************************/
		    }
		}
	    }
	for(k=0;k<N_Cut;k++)
	    {
	    j = 0;
	    for(i=0;i<Grid_x;i++) 
		{
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = 0;
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = 0.0;
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = 'n';
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = 0;
			//if(Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i].flag != 0)
			//	{
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i].grainnum;
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i].phi;
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i].phase;
		Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i].flag;
				//printf("Matrix : %c %lf, matrix +1 : %c %lf\n",Matrix[k*Grid_x*Grid_y +j*Grid_x + i].phase,Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi,Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i].phase,Matrix[k*Grid_x*Grid_y + (j+1)*Grid_x + i].phi);	
			//	}
			/*********************Clear New_Matrix ******************/
		New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = 0;
		New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = 0.0;
		New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = 'n';
		New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = 0;
				
			/******************************************/
		}
		j = Grid_y-1;
		for(i=0;i<Grid_x;i++) 
		    {
		    if(Matrix[k*Grid_x*Grid_y + (j-1)*Grid_x + i].flag != 0)
		       {
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = Matrix[k*Grid_x*Grid_y + (j-1)*Grid_x + i].grainnum;
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = Matrix[k*Grid_x*Grid_y + (j-1)*Grid_x + i].phi;
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = Matrix[k*Grid_x*Grid_y + (j-1)*Grid_x + i].phase;
		       Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = Matrix[k*Grid_x*Grid_y + (j-1)*Grid_x + i].flag;
		       }
			/*********************Clear New_Matrix ******************/
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].grainnum = 0;
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phi = 0.0;
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].phase = 'n';
		    New_Matrix[k*Grid_x*Grid_y + j*Grid_x + i].flag = 0;
				
			/******************************************/
		    }
	     }
	}
}
