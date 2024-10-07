#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double **crank_nicholson_2d(double **T,int m,int n,double B,double C,double A,double w,double e,int *iteration)
{FILE *fpr=fopen("timehistory_midpt.txt","w");
 double T_n[n][m];
 for(int i=0;i<n;i++)
  { for(int j=0;j<m;j++)
    {
     T_n[i][j]=T[i][j];
     
    
    }
  }
 double residue_sqr_sum;
 do
 {residue_sqr_sum=0;
  
  for(int i=1;i<n-1;i++)
  { for(int j=1;j<m-1;j++)
    {
     T_n[i][j]=((2-B)*T[i][j]+C*(T_n[i][j+1]+T_n[i][j-1]+T[i][j+1]+T[i][j-1])+A*(T_n[i+1][j]+T_n[i-1][j]+T[i+1][j]+T[i-1][j]))/B;
     
     double residue= T_n[i][j]-T[i][j];
    
     residue_sqr_sum+=pow(residue,2);
     //printf("%lf ",residue);
    }
  }
  
  for(int i=1;i<n-1;i++)
  { for(int j=1;j<m-1;j++)
    {
     T[i][j]=T_n[i][j];
     
    
    }
  }
  
  
  
  
  
  *iteration=*iteration+1;
 fprintf(fpr,"%f  %lf \n",*iteration*0.2,T[n/2][m/2]);  //dt=0.2
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
 fp=fopen("input2.txt","r");
 
 double X,Y,dx,dy,dt,K,diffusivity,Td,Tu,Tl,Tr,Ti,accuracy,relaxation;
 fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&X,&Y,&dx,&dy,&dt,&K,&diffusivity,&Td,&Tu,&Tl,&Tr,&Ti,&accuracy,&relaxation);
 printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",X,Y,dx,dy,dt,K,diffusivity,Td,Tu,Tl,Tr,Ti,accuracy,relaxation);
 int m=X/dx+1;
 int n=Y/dy+1;
 double B=1+2*diffusivity*dt/pow(dx,2)/2+2*diffusivity*dt/pow(dy,2)/2,C=diffusivity*dt/pow(dx,2)/2,A=diffusivity*dt/pow(dy,2)/2;
 double **temp=prepare_Temp(m,n,Tu,Td,Tl,Tr,Ti);
 fclose(fp);
 
 
 
 
 int iter=0;
 temp=crank_nicholson_2d(temp,m,n,B,C,A,relaxation=1,accuracy,&iter);
 fp=fopen("output2c.txt","w");
 printf("%d",iter);
 for(int i=0;i<n;i++)
 { for(int j=0;j<m;j++)
   {
    fprintf(fp,"%lf ",temp[i][j]);
   }
  fprintf(fp,"\n");
 } 
 
 
 
 }
