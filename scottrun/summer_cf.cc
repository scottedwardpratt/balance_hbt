#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

//const double PI=4.0*atan(1.0);
//const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	const int NRUNS=24,NQ=100,NY=20,NPHI=18;
	int iq,irun,ipair,iy,iphi;
	FILE *fptr;
	char filename[120],dummy[120];
	string dirname="results_direct";
	string pairname[12]={"pipluspiplus","pipluspiminus","piplusKplus","piplusKminus","piplusp","pipluspbar",
	"KplusKplus","KplusKminus","Kplusp","Kpluspbar","pp","ppbar"};
	double cy,cphi,cout,cside,clong,cinv,denomy,denomphi,denomout,denomside,denomlong,denominv;
	vector<double> q,dely,delphi;
	vector<double> cf_y,cf_phi,error_y,error_phi,cf_out,cf_side,cf_long,cf_inv;
	vector<double> error_out,error_side,error_long,error_inv;
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
	error_out.resize(NQ,0.0);
	error_side.resize(NQ,0.0);
	error_long.resize(NQ,0.0);
	error_inv.resize(NQ,0.0);
	
	for(ipair=0;ipair<12;ipair++){
		
		for(iq=0;iq<NQ;iq++){
			cf_out[iq]=cf_side[iq]=cf_long[iq]=cf_inv[iq]=denom_out[iq]=denom_side[iq]=denom_long[iq]=denom_inv[iq]=0.0;
			error_out[iq]=error_side[iq]=error_long[iq]=error_inv[iq]=0.0;
		}
		for(iy=0;iy<NY;iy++){
			cf_y[iy]=denom_y[iy]=error_y[iy]=0.0;
		}
		for(iphi=0;iphi<NPHI;iphi++){
			cf_phi[iphi]=denom_phi[iphi]=error_phi[iphi]=0.0;
		}
		for(irun=0;irun<NRUNS;irun++){
			//READ IN RESULTS
			sprintf(filename,"%s/%s/cf%d_y.dat",dirname.c_str(),pairname[ipair].c_str(),irun);
			//printf("will read from %s\n",filename);
			fptr=fopen(filename,"r");
			fgets(dummy,120,fptr);
			for(iy=0;iy<NY;iy++){
				fscanf(fptr,"%lf",&dely[iy]);
				fscanf(fptr,"%lf %lf",&cy,&denomy);
				if(cy!=cy)
					cy=0.0;
				cf_y[iy]+=cy*denomy;
				denom_y[iy]+=denomy;
				error_y[iy]+=cy*cy;
			}
			fclose(fptr);
			
			sprintf(filename,"%s/%s/cf%d_phi.dat",dirname.c_str(),pairname[ipair].c_str(),irun);
			//printf("will read from %s\n",filename);
			fptr=fopen(filename,"r");
			fgets(dummy,120,fptr);
			for(iphi=0;iphi<NPHI;iphi++){
				fscanf(fptr,"%lf",&delphi[iphi]);
				fscanf(fptr,"%lf %lf",&cphi,&denomphi);
				if(cphi!=cphi)
					cphi=0.0;
				cf_phi[iphi]+=cphi*denomphi;
				denom_phi[iphi]+=denomphi;
				error_phi[iphi]+=cphi*cphi;
			}
			fclose(fptr);
			
			sprintf(filename,"%s/%s/cf%d_outsidelong.dat",dirname.c_str(),pairname[ipair].c_str(),irun);
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
				error_out[iq]+=cout*cout;
				cf_side[iq]+=cside*denomside;
				denom_side[iq]+=denomside;
				error_side[iq]+=cside*cside;
				cf_long[iq]+=clong*denomlong;
				denom_long[iq]+=denomlong;
				error_long[iq]+=clong*clong;
				cf_inv[iq]+=cinv*denominv;
				denom_inv[iq]+=denominv;
				error_inv[iq]+=cinv*cinv;
			}
			fclose(fptr);
		}
		
		// AVERAGE RESULTS
		for(iy=0;iy<NY;iy++){
			cf_y[iy]=cf_y[iy]/denom_y[iy];
			if(cf_y[iy]!=cf_y[iy]){
				cf_y[iy]=0.0;
			}
		}
		for(iphi=0;iphi<NPHI;iphi++){
			cf_phi[iphi]=cf_phi[iphi]/denom_phi[iphi];
			if(cf_phi[iphi]!=cf_phi[iphi]){
				cf_phi[iphi]=0.0;
			}
		}
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
		}
		
		// WRITE RESULTS
		
		sprintf(filename,"%s/%s/cf_y.dat",dirname.c_str(),pairname[ipair].c_str());
		//printf("will write to %s\n",filename);
		fptr=fopen(filename,"w");
		for(iy=0;iy<NY;iy++){
			error_y[iy]=(error_y[iy]/double(NRUNS))-cf_y[iy]*cf_y[iy];
			error_y[iy]=sqrt(error_y[iy]/double(NRUNS));
			fprintf(fptr,"%7.2f %12.8f %12.8f %12.8f\n",dely[iy],cf_y[iy],denom_y[iy],error_y[iy]);
		}
		fclose(fptr);
		
		
		sprintf(filename,"%s/%s/cf_phi.dat",dirname.c_str(),pairname[ipair].c_str());
		//printf("will write to %s\n",filename);
		fptr=fopen(filename,"w");
		for(iphi=0;iphi<NPHI;iphi++){
			error_phi[iphi]=(error_phi[iphi]/double(NRUNS))-cf_phi[iphi]*cf_phi[iphi];
			error_phi[iphi]=sqrt(error_phi[iphi]/double(NRUNS));
			fprintf(fptr,"%7.2f %12.8f %12.8f %12.8f\n",delphi[iphi],cf_phi[iphi],denom_phi[iphi],error_phi[iphi]);
		}
		fclose(fptr);

		sprintf(filename,"%s/%s/cf_outsidelong.dat",dirname.c_str(),pairname[ipair].c_str());
		//printf("will write to %s\n",filename);
		fptr=fopen(filename,"w");
		for(iq=0;iq<NQ;iq++){
			error_out[iq]=(error_out[iq]/double(NRUNS))-cf_out[iq]*cf_out[iq];
			error_side[iq]=(error_side[iq]/double(NRUNS))-cf_side[iq]*cf_side[iq];
			error_long[iq]=(error_long[iq]/double(NRUNS))-cf_long[iq]*cf_long[iq];
			error_inv[iq]=(error_inv[iq]/double(NRUNS))-cf_inv[iq]*cf_inv[iq];
			error_out[iq]=sqrt(error_out[iq]/double(NRUNS));
			error_side[iq]=sqrt(error_side[iq]/double(NRUNS));
			error_long[iq]=sqrt(error_long[iq]/double(NRUNS));
			error_inv[iq]=sqrt(error_inv[iq]/double(NRUNS));
			fprintf(fptr,"%7.2f %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f %12.8f %12.8f\n",
			q[iq],
			cf_out[iq],denom_out[iq],error_out[iq],cf_side[iq],denom_side[iq],error_side[iq],
			cf_long[iq],denom_long[iq],error_long[iq],cf_inv[iq],denom_inv[iq],error_inv[iq]);
		}
		fclose(fptr);
	}
	return 0;
}


