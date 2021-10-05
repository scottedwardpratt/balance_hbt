#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <boost/math/special_functions.hpp>
#include "randy.h"
#include "randy.cc"

typedef double FourVector[4];

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;
using namespace boost::math;

CRandy *randy;

double triangleQ(double M,double m1,double m2){
	double Q=pow(M,4)+pow(m1,4)+pow(m2,4)-2.0*M*M*m1*m1-2.0*M*M*m2*m2-2.0*m1*m1*m2*m2;
	if(Q<0.0){
		printf("XXXXXXX  OOPS!  XXXXXXXX  Q^2=%g <0, M=%g, m1+m2=%g, m1=%g, m2=%g\n",Q,M,m1+m2,m1,m2);
		Q=0.0;
	}
	else{
		Q=sqrt(Q/(4.0*M*M));
	}
	double check=sqrt(m1*m1+Q*Q)+sqrt(m2*m2+Q*Q);
	if(fabs(check-M)>1.0E-6){
		printf("XXXXXXX  FAILED CHECK!  XXXXXXXX  check==%g=?%g, m1+m2=%g, m1=%g, m2=%g, Q=%g\n",check,M,m1+m2,m1,m2,Q);
		exit(1);
	}
	return Q;
}

void Boost(FourVector &u,FourVector &p){
	FourVector n={1.0,0.0,0.0,0.0};
	double ndotu,ndotp,udotp=u[0]*p[0]-u[1]*p[1]-u[2]*p[2]-u[3]*p[3];
	double check=p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
	ndotu=u[0];
	ndotp=p[0];
	for(int alpha=0;alpha<4;alpha++){
		p[alpha]=p[alpha]+2.0*ndotp*u[alpha]-(u[alpha]+n[alpha])*(ndotp+udotp)/(1.0+ndotu);
	}
	check=check-p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
	if(fabs(check)>1.0E-6){
		printf("Boost Check Failed = %g\n",check);
		printf("u=(%g,%g,%g,%g)\n",u[0],u[1],u[2],u[3]);
		printf("u^2=%g\n",u[0]*u[0]-u[1]*u[1]-u[2]*u[2]-u[3]*u[3]);
		exit(1);
	}
}

void GenerateP(double M,double m1,double m2,double m3,double m4,double M12,double q12,double M34,double q34,double Qmag,FourVector &p1,FourVector &p2,FourVector &p3,FourVector &p4){
	FourVector u,Q;
	double ctheta,stheta,phi,check;
	int alpha;
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	p1[3]=q12*ctheta;
	p1[1]=q12*stheta*cos(phi);
	p1[2]=q12*stheta*sin(phi);
	p1[0]=sqrt(p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3]+m1*m1);
	p2[3]=-p1[3];
	p2[1]=-p1[1];
	p2[2]=-p1[2];
	p2[0]=sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3]+m2*m2);
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	p3[3]=q34*ctheta;
	p3[1]=q34*stheta*cos(phi);
	p3[2]=q34*stheta*sin(phi);
	p3[0]=sqrt(p3[1]*p3[1]+p3[2]*p3[2]+p3[3]*p3[3]+m3*m3);
	p4[3]=-p3[3];
	p4[1]=-p3[1];
	p4[2]=-p3[2];
	p4[0]=sqrt(p4[1]*p4[1]+p4[2]*p4[2]+p4[3]*p4[3]+m4*m4);
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	Q[3]=Qmag*ctheta;
	Q[1]=Qmag*stheta*cos(phi);
	Q[2]=Qmag*stheta*sin(phi);
	Q[0]=sqrt(Qmag*Qmag+M12*M12);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=Q[alpha]/M12;
	Boost(u,p1);
	Boost(u,p2);
	
	Q[3]=-Q[3];
	Q[1]=-Q[1];
	Q[2]=-Q[2];
	Q[0]=sqrt(Qmag*Qmag+M34*M34);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=Q[alpha]/M34;
	Boost(u,p3);
	Boost(u,p4);
	
	check=p1[0]+p2[0]+p3[0]+p4[0];
	if(fabs(check-M)>1.0E-6){
		printf("Generate P Check Failed for alpha=0, check=%g, Qmag=%g\n",check,Qmag);
		exit(1);
	}
	for(alpha=1;alpha<4;alpha++){
		check=p1[alpha]+p2[alpha]+p3[alpha]+p4[alpha];
		if(fabs(check)>1.0E-6){
			printf("Generate P Check Failed for alpha=%d, check=%g\n",alpha,check);
			exit(1);
		}
	}
	
}

double GetW(double m1,double m2,double m3,double m4,double M,double wmax,double M12,double M34,double &q12,double &q34,double &Q){
	double w;
	q12=triangleQ(M12,m1,m2);
	q34=triangleQ(M34,m3,m4);
	Q=triangleQ(M,M12,M34);
	w=q12*q34*Q/wmax;
	return w;
}

int main(int argc,char *argv[]){
	FourVector p1,p2,p3,p4;
  double m1=10.0,m2=10.0,m3=0.0,m4=0.0,M=1000.0;
	double M12,M34,Q,q12,q34,M12guess,M34guess,fracguess;
	double KE12,KE34,KEtotmax;
	double w,wmax,wmaxmax=0.0,wbar=0.0;
	unsigned int iq,NQ=100,im,NM=100000;
	long long unsigned int Nwbar=0;
	bool inputmasses=false;
	randy=new CRandy(-time(NULL));
	//printf("Enter fracguess: ");
	//scanf("%lf",&fracguess);
	fracguess=0.35;
	if(inputmasses)
		NM=1;
	for (im=0;im<NM;im++){
		wmax=1.0;
		if(inputmasses){
			printf("Enter m1,m2,m3,m4: ");
			scanf("%lf %lf %lf %lf",&m1,&m2,&m3,&m4);
		}
		else{
			m1=0.9999*M*randy->ran();
			m2=(0.9999*M-m1)*randy->ran();
			m3=(0.9999*M-m1-m2)*randy->ran();
			m4=(0.9999*M-m1-m2-m3)*randy->ran();
		}
		wmax=1.0;
		M12guess=m1+m2+fracguess*(M-m1-m2-m3-m4);
		M34guess=m3+m4+fracguess*(M-m1-m2-m3-m4);
		w=GetW(m1,m2,m3,m4,M,wmax,M12guess,M34guess,q12,q34,Q);
		wmax=w*1.27;
		for(iq=0;iq<NQ;iq++){
			KEtotmax=M-m1-m2-m3-m4;
			do{
				KE12=KEtotmax*randy->ran();
				KE34=KEtotmax*randy->ran();
			}while(KE12+KE34>KEtotmax);
			M12=m1+m2+KE12;
			M34=m3+m4+KE34;
			w=GetW(m1,m2,m3,m4,M,wmax,M12,M34,q12,q34,Q);
			wbar+=w;
			Nwbar+=1;
			if(w>wmaxmax){
				wmaxmax=w;
				printf("wmaxmax=%g: m_i=(%g,%g,%g,%g), M12=%g, M34=%g\n",wmaxmax,m1,m2,m3,m4,M12,M34);
			}
			GenerateP(M,m1,m2,m3,m4,M12,q12,M34,q34,Q,p1,p2,p3,p4);			
		}
	}
	printf("wmaxmax=%g, <w>=%g\n",wmaxmax,wbar/double(Nwbar));
	return 0;
}


