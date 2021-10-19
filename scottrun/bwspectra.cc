#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
#include "balhbt.h"

using namespace std;

void WriteSpectra(vector<double> &spectra,string filename,double normfactor){
	unsigned int ipt,Nspectra=spectra.size();
	double pt;
	FILE *fptr;
	fptr=fopen(filename.c_str(),"w");
	for(ipt=0;ipt<Nspectra;ipt++){
		if(ipt<20){
			pt=(ipt+0.5)*100.0;
		}
		else{
			pt=2000+(ipt-20+0.5)*200.0;
		}
		fprintf(fptr,"%4.0f %g\n",pt,normfactor*spectra[ipt]);
	}
	fclose(fptr);
}

int GetIPt_PHENIX(double pt){
	int ipt=-1;
	if(pt<2000){
		ipt=floorl(pt/100.0);
	}
	if(pt>=2000.0){
		ipt=20+floorl((pt-2000)/200.0);
	}
	return ipt;
}

double GetPTDPTDY_PHENIX(int ipt){
	double ptdptdy,ptplus,ptminus;
	if(ipt<20){
		ptminus=ipt*100.0;
		ptplus=ptminus+100.0;
		ptdptdy=PI*(ptplus*ptplus-ptminus*ptminus);
	}
	else{
		ptminus=2000.0+(ipt-20)*200.0;
		ptplus=ptminus+200.0;
		ptdptdy=PI*(ptplus*ptplus-ptminus*ptminus);
	}
	return ptdptdy;
}

double  CalcChiSquared(vector<double> &spectra,vector<double> &error,string exp_filename,double &normfactor){
	FILE *fptr=fopen(exp_filename.c_str(),"r");
	//char dummy[120];
	bool exists[30]={false};
	double sigma2[30]={0.0},spectra_phenix[30];
	double ptbar,ptlow,pthigh,stathigh,statlow,syshigh,syslow,chi2=0.0,dchi2;
	double norm_exp=0.0,norm=0.0,ptdptdy,sp;
	double meanpt=0.0,meanpt_exp=0.0;
	int ipt;
	do{
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf",&ptbar,&ptlow,&pthigh,&sp,&stathigh,&statlow,&syshigh,&syslow);
		if(!feof(fptr)){
			if(ptbar<2.0){
				ipt=GetIPt_PHENIX(ptbar*1000.0);
				spectra_phenix[ipt]=sp;
				exists[ipt]=true;
				sigma2[ipt]=0.25*(stathigh-statlow)*(stathigh-statlow)+0.25*(syshigh-syslow)*(syshigh-syslow);
				sigma2[ipt]+=error[ipt]*error[ipt];
				ptdptdy=GetPTDPTDY_PHENIX(ipt);
				norm_exp+=spectra_phenix[ipt]*ptdptdy;
				meanpt_exp+=spectra_phenix[ipt]*ptdptdy*ptbar;
				norm+=spectra[ipt]*ptdptdy;
				meanpt+=spectra[ipt]*ptdptdy*ptbar;
			}
		}
	}while(!feof(fptr));
	normfactor=norm_exp/norm;
	meanpt=meanpt/norm;
	meanpt_exp=meanpt_exp/norm_exp;
	//printf("normalization factor=%g, <pt>=%g, phenix <pt>=%g\n",norm_exp/norm,meanpt,meanpt_exp);
	for(ipt=0;ipt<30;ipt++){
		if(exists[ipt]){
			//printf("pt=%g, phenix_spectra=%g, bw_spectra=%g, sigma2=%g\n",ptbar,spectra_phenix[ipt],spectra[ipt],sigma2);
			dchi2=pow(spectra[ipt]*normfactor-spectra_phenix[ipt],2)/sigma2[ipt];
			chi2+=dchi2;
		}
	}
	fclose(fptr);
	return chi2;
}


int main(int argc,char *argv[]){
	int run_number=0;
	double A,BW_T,BW_UPERP,efficiency;
	if (argc != 3) {
		printf("Usage: bwspectra TF UPERP\n");
		exit(-1);
	}
	else{
		A=1.0;
		BW_T=atof(argv[1]);
		BW_UPERP=atof(argv[2]);
	}
	CBalHBT *balhbt=new CBalHBT(run_number,BW_T,BW_UPERP);
	balhbt->parmap.set("BF_SPECTRA_ONLY",true);
	balhbt->Init();
	double pt;
	double ptbar_pi=0.0,ptbar_K=0.0,ptbar_p=0.0;
	vector<vector<double>> bfnorm;
	unsigned int id,id1,i,iprod,ipt,Nspectra=30;
	long long unsigned int NK=0,Np=0,Npi=0;
	vector<double> spectra_pi,spectra_K,spectra_p;
	vector<double> error_pi,error_K,error_p;
	spectra_pi.resize(Nspectra);
	spectra_K.resize(Nspectra);
	spectra_p.resize(Nspectra);
	error_pi.resize(Nspectra);
	error_K.resize(Nspectra);
	error_p.resize(Nspectra);
	for(ipt=0;ipt<Nspectra;ipt++){
		spectra_pi[ipt]=spectra_K[ipt]=spectra_p[ipt]=0.0;
		error_pi[ipt]=error_K[ipt]=error_p[ipt]=0.0;
	}
	CHBTPart *part;
	long long int imc,pid;
	long long int NMC=balhbt->parmap.getI("BW_NMC_SPECTRA",100000);
	
	vector<CHBTPart *> partvec(1);
	for(id=0;id<1;id++){
		partvec[id]=new CHBTPart();
	}
	vector<vector<CHBTPart *>> productvec(1);
	
	for(imc=0;imc<NMC;imc++){
		balhbt->GetPart(balhbt->stablevec,id1);
		partvec[0]->resinfo=balhbt->stablevec[id1]->resinfo;
		if(balhbt->randy->ran()<0.5){
			partvec[0]->PartAntipart();
		}
		balhbt->bw->GetXP(partvec);
		//for(unsigned int i=0;i<partvec.size();i++){
		//	if(partvec[i]->resinfo->decay)
		//		printf("--------\npartvec[%d].y=%g\n",i,partvec[i]->y);
		//}
		balhbt->GetDecayProducts(partvec[0],productvec[0]);
		
		
		if(balhbt->randy->ran()<0.5){
			partvec[0]->PartAntipart();
			for(unsigned long int iprod=0;iprod<productvec[0].size();iprod++){
				productvec[0][iprod]->PartAntipart();
				//if(partvec[0]->resinfo->decay)
				//	printf("productvec[%ld].y=%g\n",iprod,productvec[iprod]->y);
			}
		}

		
		for(i=0;i<1;i++){
			for(iprod=0;iprod<productvec[i].size();iprod++){
				part=productvec[i][iprod];
				if(CBF::acceptancebal->acceptance_spectra(part,efficiency)){
					pid=abs(part->resinfo->code);
					if(pid==211 || pid==2212 || pid==321){
						pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
						ipt=GetIPt_PHENIX(pt);
						if(pid==211){
							ptbar_pi+=pt;
							Npi+=1;
							if(ipt<Nspectra)
								spectra_pi[ipt]+=1;
						}
						else if(pid==321){
							ptbar_K+=pt;
							NK+=1;
							if(ipt<Nspectra)
								spectra_K[ipt]+=1;
						}
						else if(pid==2212){
							ptbar_p+=pt;
							Np+=1;
							if(ipt<Nspectra)
								spectra_p[ipt]+=1;
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
	for(ipt=0;ipt<Nspectra;ipt++){
		error_pi[ipt]=sqrt(spectra_pi[ipt]);
		error_K[ipt]=sqrt(spectra_K[ipt]);
		error_p[ipt]=sqrt(spectra_p[ipt]);
		error_pi[ipt]=1.0E6*(A/double(NMC))*error_pi[ipt]/GetPTDPTDY_PHENIX(ipt);
		error_K[ipt]=1.0E6*(A/double(NMC))*error_K[ipt]/GetPTDPTDY_PHENIX(ipt);
		error_p[ipt]=1.0E6*(A/double(NMC))*error_p[ipt]/GetPTDPTDY_PHENIX(ipt);
		spectra_pi[ipt]=1.0E6*(A/double(NMC))*spectra_pi[ipt]/GetPTDPTDY_PHENIX(ipt);
		spectra_K[ipt]=1.0E6*(A/double(NMC))*spectra_K[ipt]/GetPTDPTDY_PHENIX(ipt);
		spectra_p[ipt]=1.0E6*(A/double(NMC))*spectra_p[ipt]/GetPTDPTDY_PHENIX(ipt);
	}
	ptbar_pi=ptbar_pi/double(Npi);
	ptbar_K=ptbar_K/double(NK);
	ptbar_p=ptbar_p/double(Np);
	
	//printf("Npi=%llu, NK=%llu, Np=%llu\n",Npi,NK,Np);
	//printf("ptbar_pi=%g, ptbar_K=%g, ptbar_p=%g\n",ptbar_pi,ptbar_K,ptbar_p);
	
	double chi2_pi,chi2_K,chi2_p,chi2,norm_pi,norm_K,norm_p;
	
	chi2_pi=CalcChiSquared(spectra_pi,error_pi,"../run/phenix_data/phenix_pion.txt",norm_pi);
	chi2_K=CalcChiSquared(spectra_K,error_K,"../run/phenix_data/phenix_kaon.txt",norm_K);
	chi2_p=CalcChiSquared(spectra_p,error_p,"../run/phenix_data/phenix_proton.txt",norm_p);
	chi2=chi2_pi+chi2_K+chi2_p;
	
	WriteSpectra(spectra_pi,"spectra/spectra_pi.txt",norm_pi);
	WriteSpectra(spectra_K,"spectra/spectra_K.txt",norm_K);
	WriteSpectra(spectra_p,"spectra/spectra_p.txt",norm_p);
	
	//printf("%g      %g %g %g\n",chi2,chi2_pi,chi2_K,chi2_p);	
	printf("%g\n",chi2);	
		
	return 0;
}