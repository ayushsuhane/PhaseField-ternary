#include"function.h"
#include"Global_params.h"


void StencilGeneration(int gidx,int gidy){
 extern struct node *list,
                   *Matrix;
 extern int Grid_x,Grid_y,N_Cut;
 extern int Lenlist;
 extern int track_x,track_y;
 extern int move_matrix[5][2];

 move_matrix[0][0] = 0;
 move_matrix[0][1] = 1;
 move_matrix[1][0] = -1;
 move_matrix[1][1] = 0;
 move_matrix[2][0] = 0;
 move_matrix[2][1] = 0;
 move_matrix[3][0] = 1;
 move_matrix[3][1] = 0;
 move_matrix[4][0] = 0;
 move_matrix[4][1] = -1;

//move_matrix[][2] = {0, 1, -1, 0, 0, 0, 1, 0, 0, -1};
 int k=0;
 int x=gidx,y=gidy;
 list = (struct node*)malloc((N_Cut)*sizeof(struct node));
 int point = 0,l=0,loc = 0,n=0,count=0;
 for(k=0;k<N_Cut;k++){
  if(Matrix[Location(x,y,k)].flag != 0) {
   list[count] = Matrix[Location(x,y,k)];
   count = count + 1;
  }
 }
 for(n=0;n<5;n++){   
  x = gidx+move_matrix[n][0];
  y = gidy+move_matrix[n][1];
  for(k=0;k<N_Cut;k++){
   point = 0;
   if((Matrix[Location(x,y,k)].flag != 0)){
    for(l=0;l<count;l++){
     if(list[l].grainnum == Matrix[Location(x,y,k)].grainnum) point = point + 1;
    }
    if(point == 0){
     list[count] = Matrix[Location(x,y,k)];
     list[count].phi = 0.0;
     count = count + 1;
    }
   }
  }
  if(count == N_Cut){
   printf("ERROR in cut size\n");
   printf("gidx = %d, gidy = %d\n",gidx,gidy);
  }
 }
 int i=0,j=0;
 Lenlist = count;
 if(gidx == (track_x) && gidy == track_y){ 
//if((gidx == (track_x) && gidy ==track_y) || (gidx == (Grid_x - track_x-1) && gidy ==(Grid_y-track_y-1)))
   /*
   for(j=0;j<Grid_y;j++) 
    for(i=0;i<Grid_x;i++)
      printf("%d ",Matrix[Location(i,j,0)].grainnum);
   printf("\n");
   */
  printf("gidx = %d gidy = %d Lenlist = %d\n",gidx,gidy,Lenlist);
  printf("List:\n");
  for(k=0;k<Lenlist;k++) printf("%d %lf %c\n",list[k].grainnum,list[k].phi,list[k].phase);
  printf("Neighbours\n");
  for(n=0;n<5;n++) {
   for(l=0;l<Lenlist;l++){
    x = gidx+move_matrix[n][0];
    y = gidy+move_matrix[n][1];
    if(Matrix[Location(x,y,l)].flag!=0)  printf("%d %lf %c\t",Matrix[Location(x,y,l)].grainnum,Matrix[Location(x,y,l)].phi,Matrix[Location(x,y,l)].phase);
   }
   printf("\n");
  }
 }
}
	

