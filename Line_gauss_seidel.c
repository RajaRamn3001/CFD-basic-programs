#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//-------------------------------------------------------------------
void LU_factorise(double b1,double b2,double b3,double l[],double v[],double a[],int n)
{
 float d[n],b[n];
 for(int i=0;i<n;i++)
 {d[i]=b2;
  a[i]=b3;
  b[i]=b1;
 }
 v[0]=d[0];
 for(int i=1;i<n;i++)
 {l[i]=b[i]/v[i-1];
  v[i]=d[i]-l[i]*a[i-1];
 }


}
//------------------------------------------------------------------------
void solve_LU(double C[],double g[],double h[],double k[],int n)
{double y[n],x[n];
 y[0]=C[0];
 for(int i=1;i<n;i++)
 {y[i]=C[i]-g[i]*y[i-1];
 }
 x[n-1]=y[n-1]/h[n-1];
 for(int i=n-2;i>=0;i--)
 {x[i]=(y[i]-k[i]*x[i+1])/h[i];
 }
 
 for(int i=0;i<n;i++)
 {C[i]=x[i];
 }

}
//-----------------------------------------------------------------------
void prepare_C(double C[],double **T,double s,int n,int i)
 { C[0]=pow(s,2)*(T[i+1][1]+T[i-1][1])+T[i][0];
   C[n-3]=pow(s,2)*(T[i+1][n-2]+T[i-1][n-2])+T[i][n-1];

   for(int j=1;j<n-3;j++)
      {C[j]=pow(s,2)*(T[i+1][j+1]+T[i-1][j+1]);}

 }
//-----------------------------------------------------------------------------
double **line_gauss_seidel(double **T,int m,int n,double b1,double b2,double b3,double s,double w,double e,int *iteration)
{FILE *fpr=fopen("norm_residue_vs_iteration_q1b.txt","w");
 //-----sweep in x_direction
 double C[m-2],l[m-2],v[m-2],a[m-2];
 LU_factorise(b1,b2,b3,l,v,a,m-2);
 double residue_sqr_sum;
 
 do{
     residue_sqr_sum=0;
     for(int i=1;i<n-1;i++)
       { prepare_C(C,T,s,m,i);
         solve_LU(C,l,v,a,sizeof(C)/sizeof(C[0]));
        for(int j=1;j<m-1;j++)
         {double T_old=T[i][j]; 
           T[i][j]=C[j-1];
           T[i][j]=T_old-w*(T_old-T[i][j]);
          double residue=T[i][j]*2*(1+pow(s,2))-pow(s,2)*T[i][j+1]-pow(s,2)*T[i][j-1]- T[i+1][j]-T[i-1][j]; 
     residue_sqr_sum+=pow(residue,2);
    }
   
  }
  *iteration=*iteration+1;
  fprintf(fpr,"%d  %lf \n",*iteration,sqrt(residue_sqr_sum)); 
  }while(sqrt(residue_sqr_sum)>e); 
  
   
 
 
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
 //---------------------------------------
 //sweeping in x_direction
 double bi,di,ai;
 bi=-1;
 di=2*(1+pow(B,2));
 ai=-1;
 
 
 
 fclose(fp);
 fp=fopen("output1b.txt","w");
 
 int iter;
 fp=fopen("relaxation_vs_iteration_1b.txt","w");
 relaxation=0.6;
 while(relaxation<1.5)
 {iter=0;
  temp=line_gauss_seidel(temp,m,n,bi,di,ai,B,relaxation,accuracy,&iter);
  temp=prepare_Temp(m,n,Tu,Td,Tl,Tr,Ti);
  fprintf(fp,"%lf %d \n",relaxation,iter);
  relaxation+=0.2;
 }
 
fclose(fp);
fp=fopen("output1b.txt","w"); 
fprintf(fp,"\n");
 iter=0;
 temp=prepare_Temp(m,n,Tu,Td,Tl,Tr,Ti);
 temp=line_gauss_seidel(temp,m,n,bi,di,ai,B,relaxation=1,accuracy,&iter);
 for(int i=0;i<n;i++)
 { for(int j=0;j<m;j++)
   {
    fprintf(fp,"%lf ",temp[i][j]);
   }
  fprintf(fp,"\n");
 } 
 
 
 
 }
