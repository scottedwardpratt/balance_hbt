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
	const int NRUNS=24,NQ=100,NY=20,NPHI=18;
	int iq,irun,ipair,iy,iphi;
	FILE *fptr;
	char filename[120],dummy[120];
	string pairname[12]={"pipluspiplus","pipluspiminus","piplusKplus","piplusKminus","piplusp","pipluspbar",
	"KplusKplus","KplusKminus","Kplusp","Kpluspbar","pp","ppbar"};
	double cy,cphi,cout,cside,clong,cinv,denomy,denomphi,denomout,denomside,denomlong,denominv;
	double errory,errorphi;
	vector<double> q,dely,delphi;
	vector<double> cf_y,cf_phi,error_y,error_phi,cf_out,cf_side,cf_long,cf_inv;
	vector<double> denom_y,denom_phi,denom_out,denom_side,denom_long,denom_inv;
	q.resize(NQ,0.0);
	dely.resize(NY,0.0);
	delphi.resize(NPHI,0.0);
	cf_y.resize(NY,0.0);
	cf_phi.resize(NPHI,0.0);
	cf_out.resize(NQ,0.0);
	cf_side.resize(NQ,0.0);
	cf_long.resize(NQ,0.0);
	cf_inv.resize(NQ,0.0);
	denom_y.resize(NY,0.0);
	denom_phi.resize(NPHI,0.0);
	denom_out.resize(NQ,0.0);
	denom_side.resize(NQ,0.0);
	denom_long.resize(NQ,0.0);
	denom_inv.resize(NQ,0.0);
	error_y.resize(NY,0.0);
	error_phi.resize(NPHI,0.0);
	
	for(ipair=0;ipair<12;ipair++){
		
		for(iq=0;iq<NQ;iq++){
			cf_out[iq]=cf_side[iq]=cf_long[iq]=cf_inv[iq]=denom_out[iq]=denom_side[iq]=denom_long[iq]=denom_inv[iq]=0.0;
		}
		for(iy=0;iy<NY;iy++){
			cf_y[iy]=denom_y[iy]=error_y[iy]=0.0;
		}
		for(iphi=0;iphi<NPHI;iphi++){
			cf_phi[iphi]=denom_phi[iphi]=error_phi[iphi]=0.0;
		}
		for(irun=0;irun<NRUNS;irun++){
			
			sprintf(filename,"results/%s/cf%d_y.dat",pairname[ipair].c_str(),irun);
			//printf("will read from %s\n",filename);
			fptr=fopen(filename,"r");
			fgets(dummy,120,fptr);
			for(iy=0;iy<NY;iy++){
				fscanf(fptr,"%lf",&dely[iy]);
				fscanf(fptr,"%lf %lf %lf",&cy,&denomy,&errory);
				if(cy!=cy)
					cy=0.0;
				cf_y[iy]+=cy*denomy;
				denom_y[iy]+=denomy;
				error_y[iy]+=errory;
			}
			fclose(fptr);
			
			sprintf(filename,"results/%s/cf%d_phi.dat",pairname[ipair].c_str(),irun);
			//printf("will read from %s\n",filename);
			fptr=fopen(filename,"r");
			fgets(dummy,120,fptr);
			for(iphi=0;iphi<NPHI;iphi++){
				fscanf(fptr,"%lf",&delphi[iphi]);
				fscanf(fptr,"%lf %lf %lf",&cphi,&denomphi,&errorphi);
				if(cphi!=cphi)
					cphi=0.0;
				cf_phi[iphi]+=cphi*denomphi;
				denom_phi[iphi]+=denomphi;
				error_phi[iphi]+=errorphi;
			}
			fclose(fptr);
			
			sprintf(filename,"results/%s/cf%d_outsidelong.dat",pairname[ipair].c_str(),irun);
			//printf("will read from %s\n",filename);
			fptr=fopen(filename,"r");
			fgets(dummy,120,fptr);
			for(iq=0;iq<NQ;iq++){
				fscanf(fptr,"%lf",&q[iq]);
				fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf",
				&cout,&denomout,&cside,&denomside,&clong,&denomlong,&cinv,&denominv);
				if(cout!=cout)
					cout=0.0;
				if(cside!=cside)
					cside=0.0;
				if(clong!=clong)
					clong=0.0;
				if(cinv!=cinv)
					cinv=0.0;
				cf_out[iq]+=cout*denomout;
				denom_out[iq]+=denomout;
				cf_side[iq]+=cside*denomside;
				denom_side[iq]+=denomside;
				cf_long[iq]+=clong*denomlong;
				denom_long[iq]+=denomlong;
				cf_inv[iq]+=cinv*denominv;
				denom_inv[iq]+=denominv;
			}
			fclose(fptr);
		}
		
		sprintf(filename,"results/%s/cf_y.dat",pairname[ipair].c_str());
		//printf("will write to %s\n",filename);
		fptr=fopen(filename,"w");
		for(iy=0;iy<NY;iy++){
			cf_y[iy]=cf_y[iy]/denom_y[iy];
			error_y[iy]=error_y[iy]/double(NRUNS);
			if(cf_y[iy]!=cf_y[iy]){
				cf_y[iy]=0.0;
			}
			fprintf(fptr,"%7.2f %12.8f %12.8f %12.8f\n",dely[iy],cf_y[iy],denom_y[iy],error_y[iy]);
		}
		fclose(fptr);
		
		sprintf(filename,"results/%s/cf_phi.dat",pairname[ipair].c_str());
		//printf("will write to %s\n",filename);
		fptr=fopen(filename,"w");
		for(iphi=0;iphi<NPHI;iphi++){
			cf_phi[iphi]=cf_phi[iphi]/denom_phi[iphi];
			error_phi[iphi]=error_phi[iphi]/double(NRUNS);
			if(cf_phi[iphi]!=cf_phi[iphi]){
				cf_phi[iphi]=0.0;
			}
			fprintf(fptr,"%7.2f %12.8f  %12.8f %12.8f\n",delphi[iphi],cf_phi[iphi],denom_phi[iphi],error_phi[iphi]);
		}
		fclose(fptr);

		sprintf(filename,"results/%s/cf_outsidelong.dat",pairname[ipair].c_str());
		//printf("will write to %s\n",filename);
		fptr=fopen(filename,"w");
		for(iq=0;iq<NQ;iq++){
			cf_out[iq]=cf_out[iq]/denom_out[iq];
			if(cf_out[iq]!=cf_out[iq]){
				cf_out[iq]=0.0;
			}
			cf_side[iq]=cf_side[iq]/denom_side[iq];
			if(cf_side[iq]!=cf_side[iq]){
				cf_side[iq]=0.0;
			}
			cf_long[iq]=cf_long[iq]/denom_long[iq];
			if(cf_long[iq]!=cf_long[iq]){
				cf_long[iq]=0.0;
			}
			cf_inv[iq]=cf_inv[iq]/denom_inv[iq];
			if(cf_inv[iq]!=cf_inv[iq]){
				cf_inv[iq]=0.0;
			}
			fprintf(fptr,"%7.2f %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f %12.8f %12.8f\n",q[iq],
			cf_out[iq],denom_out[iq],cf_side[iq],denom_side[iq],cf_long[iq],denom_long[iq],cf_inv[iq],denom_inv[iq]);
		}
		fclose(fptr);
	}
	
	return 0;
}


