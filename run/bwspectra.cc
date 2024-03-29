#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
#include "balhbt.h"

using namespace std;

void WriteSpectra(vector<double> &spectra,string filename){
	unsigned int ipt,Nspectra=spectra.size();
	double pt;
	FILE *fptr;
	fptr=fopen(filename.c_str(),"w");
	for(ipt=0;ipt<Nspectra;ipt++){
		if(ipt<15){
			pt=500+(ipt+0.5)*100.0;
		}
		else{
			pt=2000+(ipt-15+0.5)*200.0;
		}
		fprintf(fptr,"%4.0f %g\n",pt,spectra[ipt]);
	}
	fclose(fptr);
}

int GetIPt_PHENIX(double pt){
	int ipt=-1;
	if(pt>500 && pt<20000){
		ipt=floorl(pt/100)-6;
	}
	if(pt>=2.0 && pt<3.0){
		ipt=floorl(pt/200)-10+5;
	}
	return ipt;
}

double GetPTDPTDY_PHENIX(int ipt){
	double ptdptdy,ptplus,ptminus;
	if(ipt<15){
		ptminus=500.0+100.0*ipt;
		ptplus=ptminus+100.0;
		ptdptdy=PI*(ptplus*ptplus-ptminus*ptminus);
	}
	else if(ipt>=15 && ipt<20){
		ptminus=2000.0+(ipt-15)*200.0;
		ptplus=ptminus+200.0;
		ptdptdy=PI*(ptplus*ptplus-ptminus*ptminus);
	}
	else{
		printf("ipt is out of range, =%d\n",ipt);
		exit(1);
	}
	return ptdptdy;
}

double  CalcChiSquared_p(vector<double> &spectra_p, vector<double> &spectra_pi, vector<double> &spectra_K, vector<double> &error_p, vector<double> &error_pi, vector<double> &error_K){
	FILE *fptr=fopen("phenix_data/phenix_proton.txt","r");
	//char dummy[120];
	double ptbar,ptlow,pthigh,spectra_phenix,stathigh,statlow,syshigh,syslow,sigma2,chi2=0.0,chipi=0.0,chiK=0.0;
	int ipt;
	for(ipt=1;ipt<15;ipt++){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf",&ptbar,&ptlow,&pthigh,&spectra_phenix,&stathigh,&statlow,&syshigh,&syslow);
		sigma2=0.25*(stathigh-statlow)*(stathigh-statlow)+0.25*(syshigh-syslow)*(syshigh-syslow);
		sigma2+=error_p[ipt]*error_p[ipt];
		//printf("pt=%g, phenix_spectra=%g, bw_spectra=%g, sigam2=%g\n",ptbar,spectra_phenix,spectra_p[ipt],sigma2);
		chi2+=pow(spectra_p[ipt]-spectra_phenix,2)/sigma2;
	}
	fclose(fptr);
	FILE *fptr1=fopen("phenix_data/phenix_pion.txt","r");
	for(ipt=0;ipt<8;ipt++){
		fscanf(fptr1,"%lf %lf %lf %lf %lf %lf %lf %lf",&ptbar,&ptlow,&pthigh,&spectra_phenix,&stathigh,&statlow,&syshigh,&syslow);
		sigma2=0.25*(stathigh-statlow)*(stathigh-statlow)+0.25*(syshigh-syslow)*(syshigh-syslow);
		sigma2+=error_pi[ipt]*error_pi[ipt];
		//printf("pt=%g, phenix_spectra=%g, bw_spectra=%g, sigam2=%g\n",ptbar,spectra_phenix,spectra_p[ipt],sigma2);
		chipi+=pow(spectra_pi[ipt]-spectra_phenix,2)/sigma2;
	}
	fclose(fptr1);
	FILE *fptr2=fopen("phenix_data/phenix_kaon.txt","r");	
	for(ipt=0;ipt<12;ipt++){
		fscanf(fptr2,"%lf %lf %lf %lf %lf %lf %lf %lf",&ptbar,&ptlow,&pthigh,&spectra_phenix,&stathigh,&statlow,&syshigh,&syslow);
		sigma2=0.25*(stathigh-statlow)*(stathigh-statlow)+0.25*(syshigh-syslow)*(syshigh-syslow);
		sigma2+=error_K[ipt]*error_K[ipt];
		//printf("pt=%g, phenix_spectra=%g, bw_spectra=%g, sigam2=%g\n",ptbar,spectra_phenix,spectra_p[ipt],sigma2);
		chiK+=pow(spectra_K[ipt]-spectra_phenix,2)/sigma2;
	}	
	return chi2+chipi+chiK;
}

int main(int argc,char *argv[]){
	int run_number=0;
	double A1, A2, A3, TF, UPERP;
	if (argc != 6) {
		printf("Usage: balance_hbt A TF UPERP\n");
		exit(-1);
	}
	else{
		A1=atof(argv[1]);
		A2=atof(argv[2]);
		A3=atof(argv[3]);
		TF=atof(argv[4]);
		UPERP=atof(argv[5]);
	}
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/respars.txt");
	parmap.ReadParsFromFile("parameters/bfpars.txt");
	parmap.ReadParsFromFile("parameters/bwpars.txt");
	parmap.set("BW_T",TF);
	parmap.set("BW_UPERP",UPERP);
	CBalHBT *balhbt=new CBalHBT(&parmap,run_number);
	
	double Tchem=150.0,taumax=100.0,strangecontent,udcontent,pt;
	double ptbar_pi=0.0,ptbar_K=0.0,ptbar_p=0.0;
	vector<vector<double>> bfnorm;
	vector<CStableInfo *> stablevec;
	unsigned int id,id1,i,iprod,ipt,Nspectrap=15, Nspectrapi=8, NspectraK=12;
	long long unsigned int NK=0,Np=0,Npi=0;
	vector<double> spectra_pi,spectra_K,spectra_p;
	vector<double> error_pi,error_K,error_p;
	spectra_pi.resize(Nspectrapi);
	spectra_K.resize(NspectraK);
	spectra_p.resize(Nspectrap);
	error_pi.resize(Nspectrapi);
	error_K.resize(NspectraK);
	error_p.resize(Nspectrap);
	for(ipt=0;ipt<Nspectrap;ipt++){
		spectra_p[ipt]=0.0;
		error_p[ipt]=0.0;
	}
	for(ipt=0;ipt<Nspectrapi;ipt++){
		spectra_pi[ipt]=0.0;
		error_pi[ipt]=0.0;
	}
	for(ipt=0;ipt<NspectraK;ipt++){
		spectra_K[ipt]=0.0;
	    error_K[ipt]=0.0;
	}
	CHBTPart *part;
	long long int imc,pid;
	int NMC=parmap.getI("BW_NMC",10000);

	balhbt->reslist->Tf=Tchem;
	balhbt->reslist->CalcEoSandChiandQdens(balhbt->reslist->Tf,balhbt->reslist->Pf,balhbt->reslist->epsilonf,balhbt->reslist->nf,balhbt->reslist->densityf,
	balhbt->reslist->maxweightf,balhbt->reslist->chif,strangecontent,udcontent);
	balhbt->reslist->FindFinalProducts(taumax);
	balhbt->bw=new CblastWave(&parmap,balhbt->randy,balhbt->reslist);
	balhbt->GetStableInfo(balhbt->reslist,taumax,stablevec,bfnorm);
	//balhbt->InitHBT(stablevec,"parameters/hbtpars.txt");
	
	vector<CHBTPart *> partvec(1);
	for(id=0;id<1;id++){
		partvec[id]=new CHBTPart();
	}
	vector<vector<CHBTPart *>> productvec(1);
	
	for(imc=0;imc<NMC;imc++){
		balhbt->GetPart(stablevec,id1);
		partvec[0]->resinfo=stablevec[id1]->resinfo;
		if(balhbt->randy->ran()<0.5){
			partvec[0]->PartAntipart();
		}
		
		balhbt->bw->GetXP(partvec);
		balhbt->GetDecayProducts(partvec[0],productvec[0]);
		
		for(int anti=0;anti<2;anti++){
			if(anti==1){
				partvec[0]->PartAntipart();
				for(unsigned long int iprod=0;iprod<productvec[0].size();iprod++){
					productvec[0][iprod]->PartAntipart();
				}
			}
		}
		
		for(i=0;i<1;i++){
			for(iprod=0;iprod<productvec[i].size();iprod++){
				part=productvec[i][iprod];
				pid=abs(part->resinfo->code);
				if(pid==211 || pid==2212 || pid==321){
					pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
					ipt=GetIPt_PHENIX(pt);
					if(ipt>=0){
						if(pid==211){
							if(ipt<Nspectrapi)
								spectra_pi[ipt]+=1;
							ptbar_pi+=pt;
							Npi+=1;
						}
						if(pid==321){
							if(ipt<NspectraK)
								spectra_K[ipt]+=1;	
							ptbar_K+=pt;
							NK+=1;
						}	
						if(pid==2212){
							if(ipt<Nspectrap)
								spectra_p[ipt]+=1;
							ptbar_p+=pt;
							Np+=1;
						}				
					}
				}
				delete part;
			}
			productvec[i].clear();
		}
		//if((imc+1)%(NMC/10)==0)
		//	printf("finished %ld percent\n",lrint(100.0*imc/double(NMC)));
	}
	for(ipt=0;ipt<Nspectrapi;ipt++){
		error_pi[ipt]=sqrt(spectra_pi[ipt]);		
		error_pi[ipt]=1.0E6*(A1/double(NMC))*error_pi[ipt]/GetPTDPTDY_PHENIX(ipt);		
		spectra_pi[ipt]=1.0E6*(A1/double(NMC))*spectra_pi[ipt]/GetPTDPTDY_PHENIX(ipt);
	}
	for(ipt=0;ipt<Nspectrap;ipt++){
		error_p[ipt]=sqrt(spectra_p[ipt]);
		error_p[ipt]=1.0E6*(A3/double(NMC))*error_p[ipt]/GetPTDPTDY_PHENIX(ipt);
		spectra_p[ipt]=1.0E6*(A3/double(NMC))*spectra_p[ipt]/GetPTDPTDY_PHENIX(ipt);
		
	}
	for(ipt=0;ipt<NspectraK;ipt++){
		error_K[ipt]=sqrt(spectra_K[ipt]);
		error_K[ipt]=1.0E6*(A2/double(NMC))*error_K[ipt]/GetPTDPTDY_PHENIX(ipt);
		spectra_K[ipt]=1.0E6*(A2/double(NMC))*spectra_K[ipt]/GetPTDPTDY_PHENIX(ipt);
	}
	ptbar_pi=ptbar_pi/double(Npi);
	ptbar_K=ptbar_K/double(NK);
	ptbar_p=ptbar_p/double(Np);
	
	//printf("Npi=%llu, NK=%llu, Np=%llu\n",Npi,NK,Np);
	//printf("ptbar_pi=%g, ptbar_K=%g, ptbar_p=%g\n",ptbar_pi,ptbar_K,ptbar_p);

	WriteSpectra(spectra_pi,"spectra/spectra_pi.txt");
	WriteSpectra(spectra_K,"spectra/spectra_K.txt");

	WriteSpectra(spectra_p,"spectra/spectra_p.txt");
	double chi2=CalcChiSquared_p(spectra_p, spectra_pi, spectra_K, error_p, error_pi, error_K);
	printf("%g\n",chi2);	
		
	return 0;
}