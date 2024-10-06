#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double **gauss_seidel(double **T,int m,int n,double B,double w,double e,int *iteration)
{FILE *fpr=fopen("norm_residue_vs_iteration_q1a.txt","w");
 double residue_sqr_sum;
 do
 {residue_sqr_sum=0;
  for(int i=1;i<n-1;i++)
  { for(int j=1;j<m-1;j++)
    {double TI=T[i][j];
     T[i][j]=(T[i][j+1]+T[i][j-1]+pow(B,2)*(T[i+1][j]+T[i-1][j]))/2/(1+pow(B,2));
     T[i][j]=TI+w*(T[i][j]-TI);
     double residue=TI*2*(1+pow(B,2))-pow(B,2)*T[i+1][j]-pow(B,2)*T[i-1][j]- T[i][j+1]-T[i][j-1]; 
    
     residue_sqr_sum+=pow(residue,2);
     //printf("%lf ",residue);
    }
  }
  *iteration=*iteration+1;
 fprintf(fpr,"%d  %lf \n",*iteration,sqrt(residue_sqr_sum));  
 }
 while(sqrt(residue_sqr_sum)>e);
 fclose(fpr);
 
 return T;
}


double **prepare_Temp(int m,int n,double Tu,double Td,double Tl,double Tr,double Ti)
{
  
  double **temp=malloc(sizeof(double *)*n);
  for(int j=0;j<n;j++)
  {
   temp[j]=malloc(sizeof(double)*m);
  }
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<m;j++) 
     {
      temp[i][j]=Ti;
     } 
  
  }
  for(int i=1;i<n-1;i++)
  {temp[i][0]=Tl;
   temp[i][m-1]=Tr;
  }
  for(int i=1;i<m-1;i++)
  {temp[0][i]=Tu;
   temp[n-1][i]=Td;
  }
  temp[0][0]=(Tu+Tl)/2;
  temp[0][m-1]=(Tu+Tr)/2;
  temp[n-1][0]=(Td+Tl)/2;
  temp[n-1][m-1]=(Td+Tr)/2;
  
  return temp;
}



int main()
{FILE *fp;
 fp=fopen("input1.txt","r");
 
 double X,Y,dx,dy,K,diffusivity,Td,Tu,Tl,Tr,Ti,accuracy,relaxation;
 fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&X,&Y,&dx,&dy,&K,&diffusivity,&Td,&Tu,&Tl,&Tr,&Ti,&accuracy,&relaxation);
 printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",X,Y,dx,dy,K,diffusivity,Td,Tu,Tl,Tr,Ti,accuracy,relaxation);
 int m=X/dx+1;
 int n=Y/dy+1;
 double B=dx/dy;
 double **temp=prepare_Temp(m,n,Tu,Td,Tl,Tr,Ti);
 fclose(fp);
 
 
 
 int iter;
 fp=fopen("relaxation_vs_iteration_1a.txt","w");
 relaxation=0.6;
 while(relaxation<1.5)
 {iter=0;
 temp=gauss_seidel(temp,m,n,B,relaxation,accuracy,&iter);
 temp=prepare_Temp(m,n,Tu,Td,Tl,Tr,Ti);
 fprintf(fp,"%lf %d \n",relaxation,iter);
 relaxation+=0.2;
 }
 fclose(fp);
 
 iter=0;
 temp=prepare_Temp(m,n,Tu,Td,Tl,Tr,Ti);
 temp=gauss_seidel(temp,m,n,B,relaxation=1,accuracy,&iter);
 fp=fopen("output1a.txt","w");
 printf("%d",iter);
 for(int i=0;i<n;i++)
 { for(int j=0;j<m;j++)
   {
    fprintf(fp,"%lf ",temp[i][j]);
   }
  fprintf(fp,"\n");
 } 
 
 
 
 }
