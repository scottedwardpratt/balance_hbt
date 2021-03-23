#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	int iphi,iy;
	FILE *fptr;
	double phi,y,dely=0.1,delqinv=5.0,qinv,delphi=10.0;
	vector<double> bfread(6);
	int irun,ib,iqinv,NPHI=18,NY=20,NQINV=100,NRUNS;
	char filename[120];
	vector<vector<double>> bf_phi;
	vector<vector<double>> bf_qinv;
	vector<vector<double>> bf_y;
	bf_phi.resize(6);
	bf_y.resize(6);
	bf_qinv.resize(6);
	printf("Enter NRUNS: ");
	scanf("%d",&NRUNS);
	for(ib=0;ib<6;ib++){
		bf_phi[ib].resize(NPHI);
		bf_qinv[ib].resize(NQINV);
		bf_y[ib].resize(NY);
		bf_phi[ib].assign(NPHI,0.0);
		bf_qinv[ib].assign(NQINV,0.0);
		bf_y[ib].assign(NY,0.0);
	}
	
	for(irun=0;irun<NRUNS;irun++){
		sprintf(filename,"results/bf%d_phi.dat",irun);
		fptr=fopen(filename,"r");
		for(iphi=0;iphi<NPHI;iphi++){
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&phi,
			&bfread[0],&bfread[1],&bfread[2],&bfread[3],&bfread[4],&bfread[5]);
			for(ib=0;ib<6;ib++){
				bf_phi[ib][iphi]+=bfread[ib]/double(NRUNS);
			}
		}
		fclose(fptr);
	}
		
	sprintf(filename,"results/bf_phi.dat");
	fptr=fopen(filename,"w");
	for(iphi=0;iphi<NPHI;iphi++){
		phi=(0.5+iphi)*delphi;
		fprintf(fptr,"%6.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",phi,
		bf_phi[0][iphi],bf_phi[1][iphi],bf_phi[2][iphi],bf_phi[3][iphi],bf_phi[4][iphi],bf_phi[5][iphi]);
	}
	fclose(fptr);
	
	for(irun=0;irun<NRUNS;irun++){
		sprintf(filename,"results/bf%d_y.dat",irun);
		fptr=fopen(filename,"r");
		for(iy=0;iy<NY;iy++){
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&y,
			&bfread[0],&bfread[1],&bfread[2],&bfread[3],&bfread[4],&bfread[5]);
			for(ib=0;ib<6;ib++){
				bf_y[ib][iy]+=bfread[ib]/double(NRUNS);
			}
		}
		fclose(fptr);
	}
	sprintf(filename,"results/bf_y.dat");
	fptr=fopen(filename,"w");
	for(iy=0;iy<NY;iy++){
		y=(0.5+iy)*dely;
		fprintf(fptr,"%6.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",y,
		bf_y[0][iy],bf_y[1][iy],bf_y[2][iy],bf_y[3][iy],bf_y[4][iy],bf_y[5][iy]);
	}
	fclose(fptr);
	
	for(irun=0;irun<NRUNS;irun++){
		sprintf(filename,"results/bf%d_qinv.dat",irun);
		fptr=fopen(filename,"r");
		for(iqinv=0;iqinv<NQINV;iqinv++){
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&qinv,
			&bfread[0],&bfread[1],&bfread[2],&bfread[3],&bfread[4],&bfread[5]);
			for(ib=0;ib<6;ib++){
				if(bfread[ib]!=bfread[ib])
					bfread[ib]=0.0;
				bf_qinv[ib][iqinv]+=bfread[ib]/double(NRUNS);
			}
		}
		fclose(fptr);
	}
	sprintf(filename,"results/bf_qinv.dat");
	fptr=fopen(filename,"w");
	for(iqinv=0;iqinv<NQINV;iqinv++){
		qinv=(0.5+iqinv)*delqinv;
		fprintf(fptr,"%6.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",qinv,
		bf_qinv[0][iqinv],bf_qinv[1][iqinv],bf_qinv[2][iqinv],bf_qinv[3][iqinv],bf_qinv[4][iqinv],bf_qinv[5][iqinv]);
	}
	fclose(fptr);
	return 0;
}


