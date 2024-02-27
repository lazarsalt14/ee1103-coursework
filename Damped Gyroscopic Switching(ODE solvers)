/******************************************************************************************************************************************************************

Roll no.:EE23B049
Name: Nishanth Senthil Kumar
Date: 3 November 2023
Version: 1
Description: To solve given ordinary differential equations using Euler's, Huens and RK45 and compare their accuracies with the help of R^2 values
Input: theta_start, theta_stop, alpha, delta_t
Output: alpha, delta_t, R^2_Euler, R^2_Heun

********************************************************************************************************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define gamma 17.6*pow(10,10) //defining the value of gamma(it is 4*pi*10^-3*17.6*10^6)
#define H pow(10,-10) //defining the value of the magnetising field applied to cause the change
#define R 1 //defining the radius of the planet
#define Hk 0 //definfing the value of Hk as one tenth of the H value


int euler(float theta_start,float h,float alpha,float theta_end);

int huens(float theta_start,float h,float alpha,float theta_end);

int rk45(float theta_start,float h,float alpha,float theta_end);

int goldstandard(float theta_start,float alpha,float theta_end);

int main(int argc, char const *argv[])
{   
    //INPUTS
    float theta_start=atof(argv[1]); //starting of the interval
    float theta_stop=atof(argv[2]); //ending of the interval
    float alpha=atof(argv[3]); //constant in the equation
    float delt=atof(argv[4]); //stepsize

    if(theta_stop>=3.141592 || theta_start<=0){
        printf("Theta_start is too low(<0) or too high (>3.141592), please enter valid input");
        return 1;
    }

    if(alpha>0.3){
        printf("Please input a smaller alpha value(i.e an alpha value less than 0.2)");
        return 2;
    }

    if(delt<0.009){
        printf("Please input a higher step size (i.e >0.01)");
        return 3;
    }

    //the return values of all the functions are the number of elements in each test case for each method of solving
    
    int elements2=huens(theta_start,delt,alpha,theta_stop);
    int elements1=euler(theta_start,delt,alpha,theta_stop);
    int elements4=goldstandard(theta_start,alpha,theta_stop);
    float *gold=malloc(elements4*sizeof(float));
    FILE *fptrmain=fopen("rk45thetagold.txt","r");
    float mean_rk_theta=0;
  
    //finding the mean of the "gold standard" points and to also store the individual points which will later be used to find R^2

    for(int i=0;i<elements4;i++){
        fscanf(fptrmain,"%*f %f",&gold[i]);
        mean_rk_theta=mean_rk_theta+gold[i];
    }
    
    //finding mean
    mean_rk_theta=mean_rk_theta/elements4;

    int times=(delt+0.00001)/0.005; //used to iterate through the arrays at different speeds 
   
    //dynamically alloting data to the arrays using "least" value
    float *eulertheta=malloc(elements1*sizeof(float));
    float *heuntheta=malloc(elements2*sizeof(float));

    //files which will be used to temporarily store the the required values 
    FILE *fptr=fopen("eulertheta.txt","r"); 
    FILE *fptr3=fopen("huenstheta.txt","r");
    
    //will be used to calculate the respective mean values
    float mean_eu_theta=0;
    float mean_hu_theta=0;

    for(int i=0;i<(elements1)-1;i++){
        //extracting the data from the files and putting them into respective arrays
        fscanf(fptr,"%*f %f",&eulertheta[i]);
        mean_eu_theta=mean_eu_theta+eulertheta[i];

    }

    for(int i=0;i<(elements2)-1;i++){

        fscanf(fptr3,"%*f %f",&heuntheta[i]);
        //calculating the summation which will later be used to find the mean
        mean_hu_theta=mean_hu_theta+heuntheta[i];

    }

    //finding the means
    mean_eu_theta=mean_eu_theta/(elements1-1);
    mean_hu_theta=mean_hu_theta/(elements2-1);

    //will be used to find the values of R^2
    float ssr_eu_theta=0;
    float sst_eu_theta=0;
    float ssr_hu_theta=0;
    float sst_hu_theta=0;

    for(int i=0;i<(elements1)-1;i++){

        //finding ssr and sst for each case which will be used to find the value of R^2
        ssr_eu_theta=ssr_eu_theta+pow(eulertheta[i]-gold[(times*i)+1],2);
        sst_eu_theta=sst_eu_theta+pow(mean_rk_theta-eulertheta[i],2);

    }

    for(int i=0;i<(elements2)-1;i++){

        ssr_hu_theta=ssr_hu_theta+pow(heuntheta[i]-gold[(times*i)+1],2);
        sst_hu_theta=ssr_hu_theta+pow(heuntheta[i]-mean_rk_theta,2);
       
    }

    //finding R^2 values
    float rsq_eu_theta=1-(ssr_eu_theta/sst_eu_theta);
    float rsq_hu_theta=1-(ssr_hu_theta/sst_hu_theta);
 
    //printing the required output,alpha,step size,r^2 euler,r^2 huens
    printf("%f %f %f %f",alpha,delt,rsq_eu_theta,rsq_hu_theta);
   
    //closing the files
    fclose(fptr);
    fclose(fptr3);
    fclose(fptrmain);

    //freeing dynamically alloted data
    free(eulertheta);
    free(heuntheta);
    free(gold);

    return 0;
}

int euler(float theta_start,float h,float alpha,float theta_end){

    //opening files to store relevant data and defining the inital values
    float theta=theta_start;
    float phi=-M_PI;
    FILE *fptr=fopen("eulertheta.txt","w");
    FILE *fptr2=fopen("eulerphi.txt","w");
    FILE *fptr3=fopen("eulerplot.txt","w");    
    double t;

    for(t=0;theta<theta_end;t=t+h){

        //converting into cartesian coordinates which is used to plot the trajectory
        float x=R*sin(theta)*cos(phi);
        float y=R*sin(theta)*sin(phi);
        float z=R*cos(theta);

        //printing the relevant values into the files
        fprintf(fptr3,"%f %f %f\n",x,y,z);
        fprintf(fptr,"%f %f\n",t,theta);
        fprintf(fptr2,"%f %f\n",t,phi);
        theta=theta+((gamma*alpha)/((alpha*alpha)+1)*(H*sin(theta)-Hk*sin(theta)*cos(theta)))*h;
      
        //resetting the value of phi after every cycle
        if(phi>M_PI){
            phi=-M_PI;
        }

        phi=phi+((gamma*alpha)/((alpha*alpha)+1)*H*sin(theta))/(alpha*sin(theta))*h;
       
    }

    //closing files
    fclose(fptr);
    fclose(fptr2);
    fclose(fptr3);

    //returning the number of values in the plot
    return (t-h)/h;
    
}

int huens(float theta_start,float h,float alpha,float theta_end){

    //opening files to store relevant data and defining the inital values
    float theta=theta_start;
    float phi=-M_PI;
    FILE *fptr=fopen("huenstheta.txt","w");
    FILE *fptr2=fopen("huensphi.txt","w");
    FILE *fptr3=fopen("huensplot.txt","w");
    float temp,temp2;

    //this will be used find the number of data points in each file and this is also the return value of the function
    int elements=0;
    double t;

    for(t=0;theta<theta_end;t=t+h){

         //converting into cartesian coordinates which is used to plot the trajectory
         float x=R*sin(theta)*cos(phi);
         float y=R*sin(theta)*sin(phi);
         float z=R*cos(theta);

         //printing the relevant values into the files
         fprintf(fptr3,"%f %f %f\n",x,y,z);
         fprintf(fptr,"%f %f\n",t,theta);
         fprintf(fptr2,"%f %f\n",t,phi);

         //finding relevant quantities
         temp=theta;
         temp2=theta+((gamma*alpha)/((alpha*alpha)+1)*(H*sin(theta)-Hk*sin(theta)*cos(theta)))*h;
         if(temp2>=theta_end){
            temp2=theta_end;
         }
         theta=theta+((((gamma*alpha)/((alpha*alpha)+1)*(H*sin(temp2)-Hk*sin(temp2)*cos(temp2))) + ((gamma*alpha)/((alpha*alpha)+1)*(H*sin(temp)-Hk*sin(theta)*cos(theta)))))*h/2;
     
        //resetting the value of phi after each cycle
         if(phi>M_PI) {
           phi=-M_PI;
            }

         phi=phi+((gamma)*H/((alpha*alpha)+1))*h;
         elements++;

     }

     //closing files
     fclose(fptr);
     fclose(fptr2);
     fclose(fptr3);

     //returning the number of elements in the plot
     return (t-h)/h;
}

int goldstandard(float theta_start,float alpha,float theta_end){
    
        //opening files to store relevant data and defining the inital values
        float theta=theta_start;
        float phi=-M_PI;
        float temp;
        float h=0.005;
        FILE *fptr=fopen("rk45thetagold.txt","w");
        

        //will be used to find k1,k2,k3,k4
        float temp1,temp2,temp3,temp4;
        int elements=0;
        double t=0;

        for(t=0;theta<theta_end;t=t+0.005){
             
                //converting into cartesian coordinates which is used to plot the trajectory
                float x=R*sin(theta)*cos(phi);
                float y=R*sin(theta)*sin(phi);
                float z=R*cos(theta);

                //printing the relevant values into the files
                fprintf(fptr,"%f %f\n",t,theta);

                //finding k1,k2,k3,k4 and applying the find the values
                temp1=((gamma*alpha)/((alpha*alpha)+1)*(H*sin(theta)-Hk*sin(theta)*cos(theta)));
                temp2=((gamma*alpha)/((alpha*alpha)+1)*(H*sin(theta+ 0.5*temp1*h)-Hk*sin(theta +0.5*temp1*h)*cos(theta+0.5*temp1*h)));
                temp3=((gamma*alpha)/((alpha*alpha)+1)*(H*sin(theta+ 0.5*temp2*h)-Hk*sin(theta +0.5*temp2*h)*cos(theta+0.5*temp2*h)));
                temp4=((gamma*alpha)/((alpha*alpha)+1)*(H*sin(theta+temp3*h)-Hk*sin(theta + temp3*h)*cos(theta+temp3*h))); 
                theta=theta+((temp1+2*temp2+2*temp3+temp4)/6)*h;

                 //resetting the values of phi
                 if(phi>M_PI) {
                    phi=-M_PI;
                }

                //finding for phi using rk45
                temp1=((gamma)/((alpha*alpha)+1)*H);
                temp2=((gamma)/((alpha*alpha)+1)*H);
                temp3=((gamma)/((alpha*alpha)+1)*H);
                temp4=((gamma)/((alpha*alpha)+1)*H);
                phi=phi+((temp1+2*temp2+2*temp3+temp4)/6)*h;
                elements++;
               
    }

    //closing files
    fclose(fptr);

    //returning the number of elements in the plot
    return (t-h)/h;
    
}
