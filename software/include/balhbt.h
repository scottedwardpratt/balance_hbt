#ifndef __BALHBT_H__
#define __BALHBT_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "parametermap.h"
#include "resonances.h"
#include "misc.h"
#include "randy.h"
#include "coral.h"
#include "constants.h"
#include "acceptancebal.h"
#include "decay_nbody.h"
#include <list>

using namespace std;

class CHBTPart;
class CStableInfo;
class CblastWave;
class CHBTCalc;
class CBF;

class CBalHBT{
public:
	int run_number;
	long long unsigned int NMC;
	double Tchem,taumax;
	bool STRANGEONLY,BARYONSONLY;
	vector<vector<double>> bfnorm;
	vector<CStableInfo *> stablevec;
	FILE *logfile;
	CBalHBT(int run_number);
	CBalHBT(int run_number,double BW_T,double BW_UPERP);
	void Init();
	~CBalHBT();
	CparameterMap parmap;
	CResList *reslist;
	CblastWave *bw;
	CRandy *randy;
	CDecay_NBody *decay_nbody;
	static CBF *bf;
	CHBTCalc *hbtcalc;
	void GetStableInfo(CResList *reslist,double taumax,vector<CStableInfo *> &stablevec,vector<vector<double>> &bfnormweight);
	void GetDecayProducts(CHBTPart *part,vector<CHBTPart *> &products);
	void GetDecayResInfo(CResInfo *resinfo0,int &ndaughters,array<CResInfo *,5> &daughter);
	void Decay(CHBTPart *mother,int &nbodies,vector<CHBTPart *> &daughter);
	void GetPart(vector<CStableInfo *> &stablevec,unsigned int &id);
	void InitHBT(vector<CStableInfo *> &stablevec,string hbtparsfilename);
	void CalcCFs();
	void CalcBFs();
	void CalcBFCFDenoms();
	void WriteResultsCF(string dirname);
	void WriteResultsBF(string dirname);
};

class CHBTPart{
public:
	FourVector p;
	FourVector x;
	double eta,y,pt;
	int motherid;
	CResInfo *resinfo;
	CHBTPart(){};
	~CHBTPart(){};
	CHBTPart(CResInfo *resinfoset){
		resinfo=resinfoset;
	}
	void BjBoost(double dely);
	void Copy(CHBTPart *);
	double GetMass();
	double GetRapidity();
	void Print();
	void PartAntipart();
	void SetEtaYPt();
	void Setp0();
	static CBalHBT *balhbt;
	double GetDCA();
};

class CStableInfo{
public:
	double density;
	CResInfo *resinfo;
	CStableInfo(CResInfo *resinfoset){
		resinfo=resinfoset;
		density=0.0;
	}
	static double denstot;
};
	
class CblastWave{
public:
	CparameterMap *parmap;
	CRandy *randy;
	CBalHBT *balhbt;
	double Tf,uperpmax,Rperp,etawidth,tau;
	double Ybeam,sigmaR;
	double sigma_eta;
	CResList *reslist;
	CblastWave(CparameterMap *parmapin,CRandy *randyset,CResList *reslistset);
	CblastWave(CBalHBT *balhbtset);
	void GenerateParts(vector<CResInfo *> &resinfovec,vector<CHBTPart *> &partvec);
	void GetXP(vector<CHBTPart *> &partvec);
	void GetXP(CHBTPart *part);
	void GetDecayMomenta(CHBTPart *mother,int &nbodies,vector<CHBTPart *> &daughterpartvec);
	void SetYbeam(double roots);
};

class CHBTCalc{
public:
	CHBTCalc(CparameterMap *parmapin);
	double GetPsiSquared(CHBTPart *part,CHBTPart *partprime,int id,int idprime);
	CparameterMap *parmap;
	vector<vector<CWaveFunction_generic *>> wf_same,wf_opp;
	CWaveFunction *wf_pp;
	CWaveFunction_classical *wf_classical;
};

class CCF_Arrays{
public:
	int NY,NPHI,NQ;
	double DELY,DELPHI,DELQ;
	CCF_Arrays(int NYset,double DELYset,int NPHIset,int NQset,double DELQset);
	vector<double> cf_inv,denom_inv,cf_out,denom_out,cf_side,denom_side,cf_long,denom_long;
	vector<double> cf_y,denom_y,cf_phi,denom_phi;
	void Increment(double dely,double delphi,double qinv,double qout,double qside,double qlong,double weight);
	void Zero();
	void Print();
	void WriteResultsCF(string dirname,long long int denom_count,int run_number);
	void WriteResultsBF(string dirname,long long int denom_count,int run_number);
	void WriteResultsBFCFDenoms(string dirname,long long int denom_count,int run_number);
};

class CBF{
public:
	CparameterMap *parmap;
	CBalHBT *balhbt;
	bool UseAllWFsForCF;
	int NYBINS,NPHIBINS,NQINVBINS;
	double DELPHI,DELY,DELQINV,QINVMAX;
	double NETWEIGHTsame,NETWEIGHTopp;
	
	CCF_Arrays *CF_pipluspiplus,*CF_pipluspiminus,*CF_KplusKplus,*CF_KplusKminus,*CF_pp,*CF_ppbar;
	CCF_Arrays *CF_piplusp,*CF_pipluspbar,*CF_piplusKplus,*CF_piplusKminus,*CF_Kplusp,*CF_Kpluspbar;
	
	void CBF_init(CparameterMap *parmapset);
	CBF(CBalHBT *balhbtset);
	static CAcceptanceBal *acceptancebal;
	void Zero();
	void Evaluate(vector<CHBTPart *> &partvec,vector<vector<CHBTPart *>> &productvec,
	vector<CHBTPart *> &partprimevec,vector<vector<CHBTPart *>> &productprimevec,int id0,int id1,int id0prime,int id1prime,double balweight,double balweightprime);
	void IncrementCF(CHBTPart *part,CHBTPart *partprime,double weight,double efficiency);
	void IncrementBFs(vector<CHBTPart *> partvec,vector<vector<CHBTPart *>> productvec,int id1,int id2,double balweight);
	void IncrementBFCFDenoms(vector<CHBTPart *> partvec,vector<vector<CHBTPart *>> productvec,int id1,int id2,double balweight);
	void WriteResultsCF(int run_number);
	void WriteResultsBF(int run_number);
	void WriteResultsBFCFDenoms(int run_number);
	double Getqinv(CHBTPart *part,CHBTPart *partprime);
	bool CHEAPPSISQUARED;
	double CheapPsiSquared(CHBTPart *part,CHBTPart *partprime);
	CHBTCalc *hbtcalc;
	CRandy *randy;
	long long int picount,Kcount,pcount,Ntest;
	
};

#endif

