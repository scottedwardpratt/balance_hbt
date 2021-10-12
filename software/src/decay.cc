#include "balhbt.h"
using namespace std;

void CBalHBT::GetDecayProducts(CHBTPart *part,vector<CHBTPart *> &products){
	int idaughter,ndaughters,alpha;
	double tau,mass;
	for(int iproduct=0;iproduct<int(products.size());iproduct++){
		delete products[iproduct];
	}
	products.clear();
	array<CResInfo *,5> dresinfo;
	vector<CHBTPart> daughter(5);
	vector<CHBTPart*> daughter_ptr(5);
	
	
	GetDecayResInfo(part->resinfo,ndaughters,dresinfo);
	if(ndaughters>5){
		printf("ndaughters=%d\n",ndaughters);
		part->resinfo->Print();
		for(idaughter=0;idaughter<ndaughters;idaughter++)
			printf("%d ",dresinfo[idaughter]->code);
		printf("\n");
		exit(1);
	}
	if(part->resinfo->decay){
		for(idaughter=0;idaughter<ndaughters;idaughter++){
			daughter[idaughter].resinfo=dresinfo[idaughter];
			daughter_ptr[idaughter]=&daughter[idaughter];
		}
		tau=HBARC/part->resinfo->width;
		mass=part->GetMass();
		for(alpha=0;alpha<4;alpha++){
			part->x[alpha]+=part->p[alpha]*tau/mass;
		}
		Decay(part,ndaughters,daughter_ptr);
	}
	else{
		ndaughters=1;
	}
	ndaughters=1;

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
		products[idaughter]->SetEtaYPt();
		if(products[idaughter]->p[0]!=products[idaughter]->p[0]){
			fprintf(logfile,"In GetDecayProducts, p0=NaN!\n");
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

void CBalHBT::Decay(CHBTPart *mother,int &nbodies,vector<CHBTPart *> &daughter){
	int ibody,alpha;
	vector<double> mass(6);
	vector<FourVector> p(5);
	double mtot;
	CHBTPart *dptr;
	FourVector u;
	mass[0]=mother->GetMass();
	
	// Make masses for daughters
	do{
		mtot=0.0;
		for(ibody=0;ibody<nbodies;ibody++){
			mass[ibody+1]=daughter[ibody]->resinfo->mass;
			mtot+=mass[ibody+1];
		}
		if(mtot>mass[0]){
			printf("mtot too big!!!\n");
			exit(1);
		}
	}while(mtot>mass[0]);

	
	if(nbodies==2){
		decay_nbody->SetMasses2(mass[0],mass[1],mass[2]);
		decay_nbody->GenerateMomenta2(p[0],p[1]);
	}
	else if(nbodies==3){
		decay_nbody->SetMasses3(mass[0],mass[1],mass[2],mass[3]);
		decay_nbody->GenerateMomenta3(p[0],p[1],p[1]);
	}
	else if(nbodies==4){
		decay_nbody->SetMasses4(mass[0],mass[1],mass[2],mass[3],mass[4]);
		decay_nbody->GenerateMomenta4(p[0],p[1],p[2],p[3]);
	}
	else{
		decay_nbody->SetMasses(nbodies,mass);
		decay_nbody->GenerateMomenta(p);
	}
	
	// Boost the new particles 
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=mother->p[alpha]/mass[0];
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		Misc::Boost(u,p[ibody]);
		for(alpha=0;alpha<4;alpha++)
			dptr->p[alpha]=p[ibody][alpha];
		dptr->eta=mother->eta;
		dptr->SetEtaYPt();
		for(alpha=0;alpha<4;alpha++){
			dptr->x[alpha]=mother->x[alpha];
		}
	}
}