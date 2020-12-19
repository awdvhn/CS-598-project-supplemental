#include <stdio.h>
#include <math.h>
int main(){
float x[7];
float y[7];
float vx[7];
float ax[7];
float vy[7];
float ay[7];
FILE *fp;
fp=fopen("/home/alan/Desktop/2en30.txt","w");
for(int i=0;i<7;i++){
	x[i]=0.0;
	vx[i]=0.0;
	ax[i]=0.0;	
	y[i]=0.0;
	vy[i]=0.0;
	ay[i]=0.0;
}
float vx0=0.05;
int j;
int i;
float m=0.0245;
float cg1=0.0219;
float cg2=0.0359;
float kg1=470.25;
float kg2=813.2;
float ke=121.17;
float ce;
float deltax;
float deltax2;
float ux;
float uy;
float cnl=0.0161;
float knl=230000000.0;
float EE[1024];
float EA[1024];
//kg2=1300.0;
for(int k=0;k<1;k++){
printf("\b\b\b\b%d",k);
vx[0]=.3;//vx0;
//vx[0]=0.3;
//fprintf(fp,"%f",vx0);
vx0+=0.000025;
for(j=0; j<1024; j++){
	for(i=0;i<7;i++){
	ux+=x[i];
	uy+=y[i];
	}
	EE[j]=0.0;
	EA[j]=0.0;
//	fprintf(fp,", %f",ux,uy);	
		ax[0]=(-1.0/m)*(cg2*vx[0]+kg1*x[0]+ke*(x[0]-y[0])+cnl*(2.0*vx[0]-vx[6]-vx[1])+
			knl*((x[0]-x[6])*(x[0]-x[6])*(x[0]-x[6])+(x[0]-x[1])*(x[0]-x[1])*(x[0]-x[1])));
		ay[0]=(-1.0/m)*(cg1*vy[0]+kg1*y[0]+ke*(y[0]-x[0])+cnl*(2.0*vy[0]-vy[6]-vy[1])+
			knl*((y[0]-y[6])*(y[0]-y[6])*(y[0]-y[6])+(y[0]-y[1])*(y[0]-y[1])*(y[0]-y[1])));
	for(i=1;i<6;i++){
	//	deltax1=x[j]-x[j-1]-x[j+1];
		ax[i]=(-1.0/m)*(cg2*vx[i]+kg2*x[i]+ke*(x[i]-y[i])+cnl*(2.0*vx[i]-vx[i-1]-vx[i+1])+
			knl*((x[i]-x[i-1])*(x[i]-x[i-1])*(x[i]-x[i-1])+(x[i]-x[i+1])*(x[i]-x[i+1])*(x[i]-x[i+1])));
		ay[i]=(-1.0/m)*(cg1*vy[i]+kg1*y[i]+ke*(y[i]-x[i])+cnl*(2.0*vy[i]-vy[i-1]-vy[i+1])+
			knl*((y[i]-y[i-1])*(y[i]-y[i-1])*(y[i]-y[i-1])+(y[i]-y[i+1])*(y[i]-y[i+1])*(y[i]-y[i+1])));
	}
		ax[6]=(-1.0/m)*(cg2*vx[6]+kg2*x[6]+ke*(x[6]-y[6])+cnl*(2.0*vx[6]-vx[5]-vx[0])+
			knl*((x[6]-x[0])*(x[6]-x[0])*(x[6]-x[0])+(x[6]-x[5])*(x[6]-x[5])*(x[6]-x[5])));
		ay[6]=(-1.0/m)*(cg1*vy[0]+kg1*y[6]+ke*(y[0]-x[0])+cnl*(2.0*vy[0]-vy[6]-vy[1])+
			knl*((y[6]-y[0])*(y[6]-y[0])*(y[6]-y[0])+(y[6]-y[5])*(y[6]-y[5])*(y[6]-y[5])));
	EE[j]=0.5*(kg1-kg2)*x[i]*x[i];
	for(i=0;i<7;i++){
		vx[i]+=0.001*ax[i];
		x[i]+=0.001*vx[i];
		vy[i]+=0.001*ay[i];
		y[i]+=0.001*vy[i];
		EE[j]+=0.5*m*vx[i]*vx[i]+0.5*kg2*x[i]*x[i]+0.25*knl*(x[i]-x[(i+1)%7])*(x[i]-x[(i+1)%7])*(x[i]-x[(i+1)%7])*(x[i]-x[(i+1)%7]);
		EA[j]+=0.5*m*vy[i]*vy[i]+0.5*kg1*y[i]*y[i]+0.25*knl*(y[i]-y[(i+1)%7])*(y[i]-y[(i+1)%7])*(y[i]-y[(i+1)%7])*(y[i]-y[(i+1)%7]);


	}
//	fprintf(fp,", %f",ax[0]);		
	ux=0.0;
	uy=0.0;
}
for(i=0;i<1024;i++) fprintf(fp," %f",EE[i]);
fprintf(fp,"\n");
for(i=0;i<1024;i++) fprintf(fp," %f",EA[i]);
}
fclose(fp);

return 420;
}
