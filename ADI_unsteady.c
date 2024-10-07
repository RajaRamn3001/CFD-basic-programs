#include<stdio.h>
#include<stdlib.h>
#include<math.h>


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
void prepare_C_x_sweep(double C[],double **T,double A,double B,int n,int i)
 { C[0]=B*(T[i+1][1]+T[i-1][1])+(1-2*B)*T[i][1]+A*T[i][0];
   C[n-3]=B*(T[i+1][n-2]+T[i-1][n-2])+(1-2*B)*T[i][n-2]+A*T[i][n-1];

   for(int j=1;j<n-3;j++)
      {C[j]=B*(T[i+1][j+1]+T[i-1][j+1])+(1-2*B)*T[i][j+1];}

 }
 //-----------------------------------------------------------------------
 void prepare_C_y_sweep(double C[],double **T,double A,double B,int n,int j)
 { C[0]=A*(T[1][j+1]+T[1][j-1])+(1-2*A)*T[1][j]+B*T[0][j];
   C[n-3]=A*(T[n-2][j+1]+T[n-2][j-1])+(1-2*A)*T[n-2][j]+B*T[n-1][j];

   for(int i=1;i<n-3;i++)
      {C[i]=A*(T[i+1][j-1]+T[i+1][j+1])+(1-2*A)*T[i+1][j];}

 }
double **ADI_2d(double **T,int m,int n,double b_x1,double b_x2,double b_x3,double b_y1,double b_y2,double b_y3,double A,double B,double w,double e,int *iteration)
{FILE *fpr=fopen("timehistory_midpt_ADI.txt","w");
 double T_n[n][m];
 for(int i=0;i<n;i++)
  { for(int j=0;j<m;j++)
    {
     T_n[i][j]=T[i][j];
     
    
    }
  }
 double C_x[m-2],C_y[n-2],l_x[m-2],v_x[m-2],a_x[m-2],l_y[n-2],v_y[n-2],a_y[n-2];
 LU_factorise(b_x1,b_x2,b_x3,l_x,v_x,a_x,m-2);
 LU_factorise(b_y1,b_y2,b_y3,l_y,v_y,a_y,n-2);
 double residue_sqr_sum;
 do
 {residue_sqr_sum=0;
  
  for(int i=1;i<n-1;i++)
  {prepare_C_x_sweep(C_x,T,A,B,m,i);
  solve_LU(C_x,l_x,v_x,a_x,sizeof(C_x)/sizeof(C_x[0]));
    for(int k=1;k<m-1;k++)
         {double T_old=T[i][k]; 
           T_n[i][k]=C_x[k-1];
           T_n[i][k]=T_old-w*(T_old-T_n[i][k]);
         }
  
  }
  //--------------------------------
  for(int i=1;i<n-1;i++)
  { for(int j=1;j<m-1;j++)
    {
     T[i][j]=T_n[i][j];
    }
  }
  //----------------------------------------- 
  
  
  for(int j=1;j<m-1;j++)
    {prepare_C_y_sweep(C_y,T,A,B,n,j);
     solve_LU(C_y,l_y,v_y,a_y,sizeof(C_y)/sizeof(C_y[0]));
     
     for(int k=1;k<n-1;k++)
         {double T_old=T[k][j]; 
           T_n[k][j]=C_y[k-1];
           T_n[k][j]=T_old-w*(T_old-T_n[k][j]);
     
     
     double residue= T_n[k][j]-T[k][j];
    
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
 double A=diffusivity*dt/pow(dx,2)/2,B=diffusivity*dt/pow(dy,2)/2;
  double b_xi,d_xi,a_xi,b_yi,d_yi,a_yi;
 b_xi=-A;
 d_xi=(1+2*A);
 a_xi=-A;
 b_yi=-B;
 d_yi=(1+2*B);
 a_yi=-B;
 double **temp=prepare_Temp(m,n,Tu,Td,Tl,Tr,Ti);
 fclose(fp);
 
 
 
 
 int iter=0;
 temp=ADI_2d(temp,m,n,b_xi,d_xi,a_xi,b_yi,d_yi,a_yi,A,B,relaxation=1,accuracy,&iter);
 fp=fopen("output2d.txt","w");
 printf("%d",iter);
 for(int i=0;i<n;i++)
 { for(int j=0;j<m;j++)
   {
    fprintf(fp,"%lf ",temp[i][j]);
   }
  fprintf(fp,"\n");
 } 
 
 
 
 }
