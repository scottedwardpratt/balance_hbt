#include "balhbt.h"
using namespace std;

void CBalHBT::GetDecayProducts(CHBTPart *part,vector<CHBTPart *> &products){
	int idaughter,ndaughters;
	for(int iproduct=0;iproduct<int(products.size());iproduct++){
		delete products[iproduct];
	}
	products.clear();
	
	array<CResInfo *,5> dresinfo;
	array<CHBTPart,5> daughter;
	GetDecayResInfo(part->resinfo,ndaughters,dresinfo);
	if(part->resinfo->decay){
		for(idaughter=0;idaughter<ndaughters;idaughter++){
			daughter[idaughter].resinfo=dresinfo[idaughter];
		}
		Decay(part,ndaughters,daughter);
	}
	else{
		ndaughters=1;
	}
	
	products.resize(ndaughters);
	if(ndaughters>1){
		for(idaughter=0;idaughter<ndaughters;idaughter++){
			products[idaughter]=new CHBTPart();
			products[idaughter]->Copy(&daughter[idaughter]);
		}
	}
	else{
		products[0]=new CHBTPart();
		products[0]->Copy(part);
	}
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		if(products[idaughter]->p[0]!=products[idaughter]->p[0]){
			printf("In GetDecayProducts, p0=NaN!\n");
		}
	}
}

void CBalHBT::GetDecayResInfo(CResInfo *resinfo0,int &ndaughters,array<CResInfo *,5> &daughter){
	//resinfo0->decay=false; //this turns off decays
	int idaughtertemp;
	int nmothers,imother,ndaughterstemp;
	array <CResInfo *,5> mother;
	array <CResInfo *,5> daughtertemp;
	/** Decay the i-particles */
	if(resinfo0->decay){
		mother[0]=resinfo0;
		nmothers=1;
		ndaughters=0;
	}
	else{
		ndaughters=1;
		nmothers=0;
		daughter[0]=resinfo0;
	}
	imother=0;
	while(imother<nmothers){
		mother[imother]->DecayGetResInfoPtr(ndaughterstemp,daughtertemp);
		for(idaughtertemp=0;idaughtertemp<ndaughterstemp;idaughtertemp++){
			if(daughtertemp[idaughtertemp]->decay){
				mother[nmothers]=daughtertemp[idaughtertemp];
				nmothers+=1;
			}
			else{
				daughter[ndaughters]=daughtertemp[idaughtertemp];
				ndaughters+=1;
			}
		}
		imother+=1;
	}
}

void CBalHBT::Decay(CHBTPart *mother,int &nbodies,array<CHBTPart,5> &daughter){
	int ibody,alpha;
	double mass[6],mtot;

	FourVector *p[6],kprime,qprime,pprime,u12,pp,u;
	double q,weight,wmax,sthet,cthet,phi;
	double p3mag,kprimemax,p3max,ppmax,kprimemax2,kprimemag2,qprimemax,qprimemax2,qprimemag2,ppmag;
	double e1prime,e2prime,e3prime,e4prime,e1max,e2max,e3max,e4max;

	mass[0]=mother->GetMass();
	p[0]=&mother->p;
	
	/* Create daughter objects */
	mtot=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		mass[ibody+1]=daughter[ibody].resinfo->mass;
	}
	for(ibody=0;ibody<nbodies;ibody++){
		mass[ibody+1]=daughter[ibody].resinfo->mass;
		mtot+=mass[ibody+1];
		p[ibody+1]=&daughter[ibody].p;
	}

	/* TWO-BODY DECAYS */
	if(nbodies==2){
		cthet=1.0-2.0*randy->ran();
		sthet=sqrt(1.0-cthet*cthet);
		phi=2.0*PI*randy->ran();
		q=sqrt(Misc::triangle(mass[0],mass[1],mass[2]));
		(*p[1])[3]=q*cthet;
		(*p[1])[1]=q*sthet*cos(phi);
		(*p[1])[2]=q*sthet*sin(phi);
		(*p[2])[3]=-(*p[1])[3];
		(*p[2])[2]=-(*p[1])[2];
		(*p[2])[1]=-(*p[1])[1];
		(*p[1])[0]=sqrt(mass[1]*mass[1]+(*p[1])[1]*(*p[1])[1]+(*p[1])[2]*(*p[1])[2]+(*p[1])[3]*(*p[1])[3]);
		(*p[2])[0]=sqrt(mass[2]*mass[2]+(*p[2])[1]*(*p[2])[1]+(*p[2])[2]*(*p[2])[2]+(*p[2])[3]*(*p[2])[3]);
	}
	/* THREE-BODY DECAYS */
	else if(nbodies==3){
		kprimemax2=Misc::triangle(mass[0]-mass[3],mass[1],mass[2]);
		kprimemax=sqrt(kprimemax2);
		p3max=sqrt(Misc::triangle(mass[0],mass[1]+mass[2],mass[3]));
		e1max=sqrt(pow(mass[1],2)+p3max*p3max);
		e2max=sqrt(pow(mass[2],2)+p3max*p3max);
		e3max=sqrt(pow(mass[3],2)+p3max*p3max);
		//wmax=p3max*(e1max*e2max/(mass[1]*mass[2]))*(mass[1]+mass[2])/(e1max+e2max);
		wmax=p3max*pow(e1max+e2max,2)*e3max/(mass[1]+mass[2]);
		do{
			TRY_AGAIN:
			do{
				kprime[1]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[2]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[3]=kprimemax*(2.0*randy->ran()-1.0);
				kprimemag2=kprime[1]*kprime[1]+
					kprime[2]*kprime[2]+kprime[3]*kprime[3];
			} while(kprimemag2>kprimemax2);
			e1prime=sqrt(kprimemag2+mass[1]*mass[1]);
			e2prime=sqrt(kprimemag2+mass[2]*mass[2]);
			if(e1prime+e2prime+mass[3]>mass[0]) goto TRY_AGAIN;
			p3mag=sqrt(Misc::triangle(mass[0],e1prime+e2prime,mass[3]));
			cthet=1.0-2.0*randy->ran();
			sthet=sqrt(1.0-cthet*cthet);
			phi=2.0*PI*randy->ran();
			(*p[3])[3]=p3mag*cthet;
			(*p[3])[1]=p3mag*sthet*cos(phi);
			(*p[3])[2]=p3mag*sthet*sin(phi);
			(*p[3])[0]=sqrt(p3mag*p3mag+mass[3]*mass[3]);
			for(alpha=1;alpha<4;alpha++)
				u12[alpha]=-(*p[3])[alpha]/(e1prime+e2prime);
			u12[0]=sqrt(1.0+u12[1]*u12[1]+u12[2]*u12[2]+u12[3]*u12[3]);
			kprime[0]=e1prime;
			Misc::lorentz(u12,kprime,*p[1]);
			kprime[0]=e2prime;
			for(alpha=1;alpha<=3;alpha++) kprime[alpha]=-kprime[alpha];
			Misc::lorentz(u12,kprime,(*p[2]));
			weight=p3mag*pow((*p[1])[0]+(*p[2])[0],2)*(*p[3])[0]/(e1prime+e2prime);
		} while(randy->ran()>weight/wmax);
	}
	/* FOUR-BODY DECAYS */
	else if(nbodies==4){
		kprimemax2=Misc::triangle(mass[0]-mass[3]-mass[4],mass[1],mass[2]);
		kprimemax=sqrt(kprimemax2);
		qprimemax2=Misc::triangle(mass[0]-mass[1]-mass[2],mass[3],mass[4]);
		qprimemax=sqrt(qprimemax2);
		
		ppmax=sqrt(Misc::triangle(mass[0],mass[1]+mass[2],mass[3]+mass[4]));
		e1max=sqrt(pow(mass[1],2)+ppmax*ppmax);
		e2max=sqrt(pow(mass[2],2)+ppmax*ppmax);
		e3max=sqrt(pow(mass[3],2)+ppmax*ppmax);
		e4max=sqrt(pow(mass[4],2)+ppmax*ppmax);
		wmax=ppmax*pow(e1max+e2max,2)*pow(e3max+e4max,2)/((mass[1]+mass[2])*(mass[3]+mass[4]));
		
		do{
			TRY_AGAIN_4:
			do{
				kprime[1]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[2]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[3]=kprimemax*(2.0*randy->ran()-1.0);
				kprimemag2=kprime[1]*kprime[1]+kprime[2]*kprime[2]+kprime[3]*kprime[3];
			} while(kprimemag2>kprimemax2);
			e1prime=sqrt(kprimemag2+mass[1]*mass[1]);
			e2prime=sqrt(kprimemag2+mass[2]*mass[2]);
			do{
				qprime[1]=qprimemax*(2.0*randy->ran()-1.0);
				qprime[2]=qprimemax*(2.0*randy->ran()-1.0);
				qprime[3]=qprimemax*(2.0*randy->ran()-1.0);
				qprimemag2=qprime[1]*qprime[1]+qprime[2]*qprime[2]+qprime[3]*qprime[3];
			} while(qprimemag2>qprimemax2);
			e3prime=sqrt(qprimemag2+mass[3]*mass[3]);
			e4prime=sqrt(qprimemag2+mass[4]*mass[4]);
			
			if(e1prime+e2prime+e3prime+e4prime>mass[0]) goto TRY_AGAIN_4;

			ppmag=Misc::triangle(mass[0],e1prime+e2prime,e3prime+e4prime);
			if(ppmag>0){
				ppmag=sqrt(ppmag);
				cthet=1.0-2.0*randy->ran();
				sthet=sqrt(1.0-cthet*cthet);
				phi=2.0*PI*randy->ran();
				pp[3]=ppmag*cthet;
				pp[1]=ppmag*sthet*cos(phi);
				pp[2]=ppmag*sthet*sin(phi);

				pp[0]=sqrt(ppmag*ppmag+pow(e1prime+e2prime,2));
				for(alpha=0;alpha<4;alpha++)
					u[alpha]=pp[alpha]/(e1prime+e2prime);
				kprime[0]=sqrt(mass[1]*mass[1]+kprimemag2);
				Misc::lorentz(u,kprime,*p[1]);
				kprime[0]=sqrt(mass[2]*mass[2]+kprimemag2);
				for(alpha=1;alpha<4;alpha++)
					kprime[alpha]=-kprime[alpha];
				Misc::lorentz(u,kprime,*p[2]);

				for(alpha=1;alpha<4;alpha++)
					pp[alpha]=-pp[alpha];
				pp[0]=sqrt(ppmag*ppmag+pow(e3prime+e4prime,2));
				for(alpha=0;alpha<4;alpha++)
					u[alpha]=pp[alpha]/(e3prime+e4prime);
				qprime[0]=sqrt(mass[3]*mass[3]+qprimemag2);
				Misc::lorentz(u,qprime,(*p[3]));
				qprime[0]=sqrt(mass[4]*mass[4]+qprimemag2);
				for(alpha=1;alpha<4;alpha++)
					qprime[alpha]=-qprime[alpha];
				Misc::lorentz(u,qprime,*p[4]);

				weight=ppmag*pow((*p[1])[0]+(*p[2])[0],2)*pow((*p[3])[0]+(*p[4])[0],2)/((e1prime+e2prime)*(e3prime+e4prime));
			}
			else weight=0.0;
		} while(randy->ran()>weight/wmax);
	}

	/* Boost the new particles */
	for(alpha=0;alpha<4;alpha++) u[alpha]=mother->p[alpha]/mass[0];
	for(ibody=0;ibody<nbodies;ibody++){
		Misc::lorentz(u,*p[ibody+1],pprime);
		for(alpha=0;alpha<4;alpha++){
			if(pprime[alpha]!=pprime[alpha]){
				printf("In Decay, pprime is NaN!\n");
				exit(1);
			}
			daughter[ibody].p[alpha]=pprime[alpha];
			daughter[ibody].x[alpha]=mother->x[alpha];
		}
	}
}
