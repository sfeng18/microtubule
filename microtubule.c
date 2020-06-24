#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <string.h>
#include <errno.h> 
#include "mpi.h"

#define MD_ON		0
#define HT_ON		0
#define FC_ON		1
#define FR_ON		1
#define MR_ON		1
#define CAP_ON		0	//+:GTP Cap -:GDP Cap
#define PR_ON		0
#define Di_Inertia	0
#define GROW_ON		0
#define BREAK_ON	0
#define HD_ON		0
#define RLW_ON		0
#define EOut_ON		0
#define MustGrow	0
#define GSteps		50000
#define G_Slow		0
#define GRStep		0
#define HD_Delay	0
#define RandLatP	0.5
#define LCompute	0
#define SlowTime	20
#define RS_OUT_ON	0

#define pi			3.1415926535897
#define cutoff		5
#define cutoff_s	2.7
#define cutoffsq	25
#define CellL		7
#define G_rate		7.0
#define EBase		0.71
#define Rate_In		0.6
#define Rate_Side	0.35
#define Rate_On		0.2
#define C_rate		1e-12
#define HD_rate		5.0e-7
#define S_rate		1e7
#define B_rate		5e-3
#define Afree		7.2e-11
#define dt_run		1e-2
#define FricCo		10.0
#define FricCo_o	10.0
#define GammaN		15
#define MaxFName	256
#define P_LT		5.0e3
#define W_In		1
#define W_Out		1
#define M_F			1.0
#define I_M			40.0
#define W_Lat_aa	0.1811
#define W_Lat_ab	0.1485
#define CTime		4000000
#define WallTime	400000
#define TipEnd		0
#define TipAll		0
#define TipEndPf	11

#define E_Conver	1e-10
#define T_Conver	200000
#define T_Hydrolyze	100000
#define Out_Step	10000

#define NLen		20000
#define DLen		0.00005
#define deltasqrt	DLen
#define deltalj		DLen
#define deltaexp	DLen

#define Detach_Curve 0
#define Temp_Kappa	0	//1:二分 2:渐增 
//#define Tub_R_Out 1
#define Test_Type	0
#define CAPBegin	221
#define CAPNo		-26	//+:GTP Cap -:GDP Cap
#define SecSteps	1
#define SecLen		0.1
#define BallR		10
#define Ballcut		400
#define InitFila	13
#define InitR		(0.98/sin(pi/FilaSet))//4.0
#define FilaMis		0
#define TipPf		FilaSet
#define TipPfs		1
#define TipLayer	7
#define OutputTime	1
#define para_m		1.0
#define para_l		1.5
#define para_s		2.0
#define RFType		0 //Distribution of Reinforced Units 0:Uniform 1:Gauss
#define RFLayers	0
#define RFUnits		39//26
#define RFRate		0.0
#define RFDensityOut 0
#define LongRFTest	0
#define CV_RFTest	0

#define Max_Nbr_Num		100
#define Max_Nbr_SqDis	27
#define Max_Add_Link	12

#define NInit		180
#define Times		20000
#define MTimes		800

#define My_Type		0
#define HELIX_ON	1
#define Unit_Type	1
#define Whole_Struc	1		//0:Sheet 1:Tube 2:Tube From File 
#define Lat_Length	1
#define TipType		1		//0:Sin 1:Mix 2:Power 3:Power Shift 4:Double Peak 5:Linear Mix 6:Decay
#define TipSinA		0
#define TipSinL		1.2
#define TipSinH		0
#define TipARatio	0
#define TipSinPhase	0
#define XI_0		5.0
#define BoxOut		0x00	//High: 1:PatchyPoint 2:nst Low: 1:fp Box 2:vmd Box 4:fp connect 8:vmd connect
#define DataOut		0
#define FixLink		1
#define PhaseFile	1		//0:Given by paras 1:Given by File 'PhaseData'	2:Given by File 'PhaseData07'

#define StepR		0.005
#define StepA		0.005

#define NCell			50
#define NCellX			40
#define NCellY			20
#define NCellZ			20
#define LCutSt			2.0
#define LCutEd			2.7
#define ACutSt			0//(pi/20)
#define ACutEd			(pi/3)
#define ELatAA			0.0906	//0.2286 0.0906 0.3628 0.1811
#define ELatBA			0.0743	//0.1874 0.0743 0.2975 0.1485
#define ELong			1.0
#define P_LJC			10
#define P_LJS			10//P_LJC
#define KLat			5
#define KLong			10
#define TLong			2.0e2//(30.0*TLat)
#define GLat			2e0
#define WeakAttr		0.0
#define TwistType		1
#define QDDelta			0.02
#define QDType			0x09	//High: chain Low: side 1:distribution 2: uniform 4: uni-axis 8: uni-Layer

#define Tst				3
#define Ted				50000
#define DeltaT			10
#define DeltaL			0e-2		//UniformForce: Strain Other: Distance
#define DeltaXi			0e-4
#define DeltaA			0e-3
#define BCType			0x10	//High:-end Low:+end 0:Free 1:Hinged 2:Slide 4:Fix
#define ForceXYZ		0x00	//Low: Direction:xyz(124) Rotate(8) High: 0:None 1:Uniform 2:+end 4:-end 8:ShearEdge

#if Whole_Struc && Test_Type!='S' && Test_Type!='K' && Test_Type!='V'
	#define FixEnd	1
#else
	#define FixEnd	0
#endif

typedef struct{
	int Type,Tub_No,Active,Tub_Id,HydroTime;
	double r[3],rinit[3];
#if MD_ON
	double v[3],a[3],apre[3],vomi[3],aomi[3],aomipre[3],mom[3],mombar[3],vomibar[3],rpre[3];
#endif
	double angle[5],xi;
	double n[3],s[3],t[3];
	double r1[3],r2[3],rbar1[3],rbar2[3];
	double chain1[3],chain2[3],side1[3],side2[3],side3[3],side4[3];
	double NS1[3],NS2[3],NS3[3],NS4[3];
	double MyPoten[12],QDRate[6];
	double poten,PotenB,PotenA,potent;
#if MD_ON
	double potenf,potenb,p_poten[3];
#endif
	double Phd,Puhd;
#if FixLink
	double PotenList[Max_Add_Link];
	int LinkList[Max_Add_Link],PosList[Max_Add_Link];
	char PPMark[Max_Add_Link];
	int AddNum; 
#endif
	int front,left,right,back,lseam,rseam;
	int ChainNo,CMark,LineNo,Tin,Ttop,Thd;//ChainNo:top->bottom,LineNo:bottom->top,Cmark:1~13
	int RFMark,WallMark;
}tubulin;

typedef struct{
	double r[3],Radius,BForce[3];
}ball;

int N=NInit,Nactive,Ntype,ttubno=13,BlockSize,TubSize,dtTime=(int)(dt_run/1e-2),BUGst=0,FilaSet=InitFila,MaxLayer=200;
int MsgTag,BUG=0,StopMsg=0,Ngrow=0,TGrow=GSteps,MpiSec,NStart=0,LStart=-1,DCCount=0,NA_A;
const int Nbox=NCellX*NCellY*NCellZ,NCsq=NCellX*NCellY,//BoxL=CellL*NCell/2,Box2L=CellL*NCell,
	BoxLX=CellL*NCellX/2,BoxLY=CellL*NCellY/2,BoxLZ=CellL*NCellZ/2,Box2LX=CellL*NCellX,Box2LY=CellL*NCellY,Box2LZ=CellL*NCellZ;
double KT=2e-3,TReal=0,KappaL=0,KappaH=20*pi/180,dt=dt_run;
double ansPrC[cutoff*NLen+2],ansPrS[cutoff*NLen+2],ansRDcyC[cutoff*NLen+2];
double ansPAC[2*NLen+2],ansADcyC[2*NLen+2];
double ansACos[2*NLen+2],ansSin[8*NLen+2],ansCos[8*NLen+2];
double PAngle1[2*NLen+2],PAngle2[2*NLen+2],PAngle3[2*NLen+2],PAngle4[2*NLen+2],RGamma[100];
double ForceR[cutoff*NLen+2],ForceAngle1[2*NLen+2],ForceAngle2[2*NLen+2],ForceAngle3[2*NLen+2],ForceAngle4[2*NLen+2];
#if !My_Type
	double alpha0=0*pi/180,omiga0=0*pi/180,xi0=XI_0*pi/180,xi1=12*pi/180,alpha1=3*pi/180,omiga1=3*pi/180;
	#if HELIX_ON
		double beta=103.45*pi/180,gamma0=13.8*pi/180;
	#else
		double beta=90*pi/180,gamma0=13.8*pi/180;
	#endif
#else
	double alpha0=0*pi/180,omiga0=0*pi/180,xi0=12*pi/180,xi1=XI_0*pi/180,alpha1=0*pi/180,omiga1=0*pi/180;
	#if HELIX_ON
		double beta=103.45*pi/180,gamma0=13.8*pi/180;
	#else
		double beta=90*pi/180,gamma0=13.8*pi/180;
	#endif
#endif
double SinBeta,SinGamma,CosBeta,CosGamma;
double cosxi_z[T_Hydrolyze+1],cosa_z[T_Hydrolyze+1],coso_z[T_Hydrolyze+1];
double sinxi_z[T_Hydrolyze+1],sina_z[T_Hydrolyze+1],sino_z[T_Hydrolyze+1];
double C0,C1,C2,Co0,Co1,Co2,sigmav,sigmar,Crv,SqCrv;
double *Rand_N;
double RAccept,at1,at2,*Tran;
double Distance[3],Dlength,MaxA=0,MinA=0,MyTipPhase[6]={0.0};
double *EPoten,*EVel,*ERot,*ETot,*FAvg,*FBall;
int TConver=0,TTot=0,MConver=0,Trec=0;
int *Head,*List,*CList,GNo,CNo,*Chain,**TGroup,*Top,*Tmark,*TubeEnd,*NbrList;
int icell,MaxEnd=0;
//double topr[4];
FILE *fout,*fmsg,*fbug,*fphase,*fshort;
tubulin *Tub,*ttub;
ball BigBall;
int point_info[10][4]={{0,-1,1,0},{1,-1,0,0},{2,-1,3,1},{2,-1,4,2},{3,-1,2,1},{3,-1,5,2},{4,-1,2,2},{4,-1,5,1},{5,-1,3,2},{5,-1,4,1}};
char fstring[MaxFName],tstring[MaxFName],*TubState;
/*
points chain:0,1 side:2,3,4,5
type All:-1 GTP:0 GDP:1
num same as points
potential 0,1,2,...
*/
//MPI
int myid,namelen;
int numprocs;
char processor_name[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;

inline double Dot(double *a,double *b){
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

inline double Dis2(double *a,double *b){
	return ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]));
}

inline int GetCellNo(double *r){
	return ((int)((r[0]+BoxLX)/CellL)+((int)((r[1]+BoxLY)/CellL))*NCellX+((int)((r[2]+BoxLZ)/CellL))*NCsq);
}

inline void PeriodBoundary(double *r){
	if(r[0]>BoxLX)r[0]-=Box2LX;
	else if(r[0]<-BoxLX)r[0]+=Box2LX;
	if(r[1]>BoxLY)r[1]-=Box2LY;
	else if(r[1]<-BoxLY)r[1]+=Box2LY;
	if(r[2]>BoxLZ)r[2]-=Box2LZ;
	else if(r[2]<-BoxLZ)r[2]+=Box2LZ;
}

double GammaRand(){	//生成均值为1 k等于GammaN 的gamma分布随机数 
	int i;
	double uk,gammak;
	gammak=0;
	for(i=0;i<GammaN;i++){
		uk=rand()/(RAND_MAX+0.1);
		gammak+=-log(uk)/GammaN;
	}
	return gammak;
}

void GaussRand(){	//用于生成12*N个符合标准正态分布的随机数 
	int i;
	double rand1,rand2,eta;
	free(Rand_N);
	Rand_N=(double *)malloc(12*N*sizeof(double));
	for(i=0;i<6*N;i++){
		do{
			rand1=2.0*rand()/(RAND_MAX+0.1)-1;
			rand2=2.0*rand()/(RAND_MAX+0.1)-1;
			eta=rand1*rand1+rand2*rand2;
		}while(eta>=1.0 || eta==0.0);
		eta=sqrt(-2*log(eta)/eta);
		Rand_N[2*i]=rand1*eta;
		Rand_N[2*i+1]=rand2*eta;
	}
}

double PGauss(double x){	//求符合标准正态分布的随机变量 X<x 的概率，精确到 1e-9 
	int i=0;
	double erf_i,erf_sum=0.0,x2;	//erf_i为第i项的值, erf_sum为总的值, x2为 x的平方, Fac_i为 i的阶乘(factorial) 
	x*=sqrt(0.5);
	erf_i=1.0*x;
	x2=1.0*x*x;
	while(fabs(erf_i)>1e-10){
		erf_sum+=erf_i;
		i++;
		erf_i*=-x2*(2*i-1)/(2*i+1)/i;
	}
	erf_sum=0.5+sqrt(1/pi)*erf_sum;
	return(erf_sum);
}

void RotMatrx(double *Mrot,double *Qin){	//依据四元数求旋转矩阵 
	double x,y,z,costheta,sintheta;
	x=Qin[0];
		y=Qin[1];
		z=Qin[2];

		costheta=cos(2*Qin[3]);
		sintheta=sin(2*Qin[3]);

		Mrot[0]=costheta+x*x*(1-costheta);
		Mrot[1]=sintheta*z+x*y*(1-costheta);
		Mrot[2]=-sintheta*y+x*z*(1-costheta);
		Mrot[3]=-sintheta*z+x*y*(1-costheta);
		Mrot[4]=costheta+y*y*(1-costheta);
		Mrot[5]=sintheta*x+y*z*(1-costheta);
		Mrot[6]=sintheta*y+x*z*(1-costheta);
		Mrot[7]=-sintheta*x+y*z*(1-costheta);
	    Mrot[8]=costheta+z*z*(1-costheta);//局部坐标
}

void CoodRot(double *Rout,double *Rin,double *Mrot){	//依据旋转矩阵转动坐标 
	double r[3];
	r[0]=Rin[0];
	r[1]=Rin[1];
	r[2]=Rin[2];
	Rout[0]=Mrot[0]*r[0]+Mrot[3]*r[1]+Mrot[6]*r[2];
	Rout[1]=Mrot[1]*r[0]+Mrot[4]*r[1]+Mrot[7]*r[2];
	Rout[2]=Mrot[2]*r[0]+Mrot[5]*r[1]+Mrot[8]*r[2];
}

void CrossL(double *a,double *b,double *c){ //c=a X b
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
}

void QuaTimes(double *a,double *b,double *c){	//四元数乘法  c=a * b
	int i;
	double c1,c2,c3,s1,s2,s3,s12,c1s2,c2s1;
	s1=a[3];
	c1=a[4];
	s2=b[3];
	c2=b[4];
	s12=s1*s2;
	c1s2=c1*s2;
	c2s1=c2*s1;
	c[0]=s12*(a[1]*b[2]-a[2]*b[1])+c1s2*b[0]+c2s1*a[0];
	c[1]=s12*(a[2]*b[0]-a[0]*b[2])+c1s2*b[1]+c2s1*a[1];
	c[2]=s12*(a[0]*b[1]-a[1]*b[0])+c1s2*b[2]+c2s1*a[2];
	s3=sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
	c[3]=s3;
	c[4]=c1*c2-s12*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
	if(s3<1e-100&&s3>-1e-100){
		c[0]=1;c[1]=0;c[2]=0;c[3]=s3;
	}
	else{
		s3=1/s3;
		c[0]=c[0]*s3;
		c[1]=c[1]*s3;
		c[2]=c[2]*s3;
	}
}

void Set_Temp(){
	double KB=1.380662e-23,Da=1.660538921e-27,Mdimer=1,dtReal,vScale,rScale,tScale;
	if(dtTime){
		dt=dt_run/dtTime;
		dtTime--;
		fprintf(stderr,"dt:%g\tdtT:%d\n",dt,dtTime);
	}
	fprintf(stdout,"Seting temperature parameters...\n");
//布朗动力学参数 	
	C0=exp(-dt*FricCo);
	C1=(1-C0)/(dt*FricCo);
	C2=(1-C1)/(dt*FricCo);
	Co0=exp(-dt*FricCo_o);
	Co1=(1-Co0)/(dt*FricCo_o);
	Co2=(1-Co1)/(dt*FricCo_o);
	Mdimer=2*55000*Da;
	KT=TReal*KB/Mdimer;
	rScale=1/2e-9;
	tScale=sqrt(18.5*300*KB/Mdimer)*rScale;
	vScale=rScale/tScale;
	dtReal=dt/tScale;
	sigmav=vScale*sqrt(KT*(1-exp(-2*dt*FricCo)));
	sigmar=rScale*dtReal*sqrt(KT/dt/FricCo*(2-(3-4*exp(-FricCo*dt)+exp(-2*dt*FricCo))/dt/FricCo));
	Crv=(1-C0)*(1-C0)/sqrt((2-(3-4*exp(-FricCo*dt)+exp(-2*dt*FricCo))/dt/FricCo)*(1-exp(-2*dt*FricCo))*dt*FricCo);
	SqCrv=sqrt(1-Crv*Crv);
	fprintf(stdout,"Done\n");
	fprintf(stdout,"Scale: r:%g t:%g v:%g\n",rScale,tScale,vScale);
	fprintf(stdout,"T:%g KT:%g sigv:%g sigr:%g Crv:%g\n",TReal,KT,sigmav,sigmar,Crv);
}

void Set_Tub(tubulin *Tubin,double *r,double *quaternion,double xi){//,double a,double b,double c){	//粒子状态设定  
	int j=0,MyHT;
	double x,y,z,w,sinxi,cosxi,costheta,sintheta;
	double ch1,ch2,si1,si2,si3,si4,tb1,tb2,sb1,sb2,tmod,smod,nmod;
	double cosbeta,sinbeta,cosgamma,singamma;
	double tlocal[3],slocal[3],nlocal[3],tbar1[3],tbar2[3],sbar1[3],sbar2[3];
	double CtS1[3],CtS2[3],CtS3[3],CtS4[3],NS1mod,NS2mod,NS3mod,NS4mod;
	
	ch1=ch2=si1=si2=si3=si4=tb1=tb2=sb1=sb2=0;
	tmod=smod=nmod=0;
	cosbeta=CosBeta;
	sinbeta=SinBeta;
	cosgamma=CosGamma;
	singamma=SinGamma;
	
		x=quaternion[0];
		y=quaternion[1];
		z=quaternion[2];
		w=quaternion[3];
		costheta=quaternion[4]*quaternion[4]-quaternion[3]*quaternion[3];
		sintheta=2*quaternion[3]*quaternion[4];
		tlocal[0]=costheta+x*x*(1-costheta);
		tlocal[1]=sintheta*z+x*y*(1-costheta);
		tlocal[2]=-sintheta*y+x*z*(1-costheta);
		slocal[0]=-sintheta*z+x*y*(1-costheta);
		slocal[1]=costheta+y*y*(1-costheta);
		slocal[2]=sintheta*x+y*z*(1-costheta);
		nlocal[0]=sintheta*y+x*z*(1-costheta);
		nlocal[1]=-sintheta*x+y*z*(1-costheta);
	    nlocal[2]=costheta+z*z*(1-costheta);//局部坐标
		j=1;
		MyHT=Tubin->HydroTime;
		if(Tubin->Type){
			#if HT_ON
			if(MyHT<T_Hydrolyze)MyHT++;
			#else
			MyHT=T_Hydrolyze;
			#endif
		}
		else MyHT=0;
		
			if(!MD_ON){
				xi=xi/2;
				sinxi=sin(xi);
				cosxi=cos(xi);
			}
			else{
				sinxi=sinxi_z[MyHT];
				cosxi=cosxi_z[MyHT];
			}
		for(j=0;j<3;j++){
			tbar1[j]=tlocal[j]*cosxi-slocal[j]*sinxi;
			tbar2[j]=-tlocal[j]*cosxi-slocal[j]*sinxi;
			sbar1[j]=tlocal[j]*sinxi+slocal[j]*cosxi;
			sbar2[j]=-tlocal[j]*sinxi+slocal[j]*cosxi;
		}
		for(j=0;j<3;j++){
			Tubin->r[j]=r[j];
			Tubin->angle[j]=quaternion[j];
			
			Tubin->rbar1[j]=tbar1[j];
			Tubin->rbar2[j]=tbar2[j];
			Tubin->r1[j]=r[j]+tbar1[j];//alpha（-）位置
			Tubin->r2[j]=r[j]+tbar2[j];//beta（+）位置
			
			Tubin->chain1[j]=tbar1[j]*cosa_z[MyHT]-sbar1[j]*sina_z[MyHT];
			
			Tubin->chain2[j]=tbar2[j]*coso_z[MyHT]-sbar2[j]*sino_z[MyHT];
				
			Tubin->side1[j]=(tbar1[j]*cosbeta-nlocal[j]*sinbeta)*cosgamma+sbar1[j]*singamma;//侧向结合位点方向
			Tubin->side2[j]=(-tbar1[j]*cosbeta+nlocal[j]*sinbeta)*cosgamma+sbar1[j]*singamma;
			Tubin->side3[j]=(-tbar2[j]*cosbeta-nlocal[j]*sinbeta)*cosgamma+sbar2[j]*singamma;//侧向结合位点方向
			Tubin->side4[j]=(tbar2[j]*cosbeta+nlocal[j]*sinbeta)*cosgamma+sbar2[j]*singamma;
			
			Tubin->t[j]=tlocal[j];
			Tubin->s[j]=slocal[j];
			Tubin->n[j]=nlocal[j];
		}
		CrossL(tlocal,Tubin->side1,CtS1);
		CrossL(tlocal,Tubin->side2,CtS2);
		CrossL(tlocal,Tubin->side3,CtS3);
		CrossL(tlocal,Tubin->side4,CtS4);
		NS1mod=1/sqrt(Dot(CtS1,CtS1));
		NS2mod=1/sqrt(Dot(CtS2,CtS2));
		NS3mod=1/sqrt(Dot(CtS3,CtS3));
		NS4mod=1/sqrt(Dot(CtS4,CtS4));
		for(j=0;j<3;j++){
			Tubin->NS1[j]=CtS1[j]*NS1mod;
			Tubin->NS2[j]=CtS2[j]*NS2mod;
			Tubin->NS3[j]=CtS3[j]*NS3mod;
			Tubin->NS4[j]=CtS4[j]*NS4mod;
		}
		Tubin->HydroTime=MyHT;
		Tubin->angle[3]=quaternion[3];
		Tubin->angle[4]=quaternion[4];
}

void Set_FileName(){
	int NameL,HTInt;
	char TypeChar,HelType,ttmp[MaxFName];
	time_t tnow;
//输出文件名字符串 
	fprintf(stdout,"Seting FileName...");
	if(My_Type){
		if(Unit_Type)TypeChar='M';
		else TypeChar='Y';
	}
	else if(Unit_Type)TypeChar='D';
	else TypeChar='T';
	if(HT_ON)HTInt=T_Hydrolyze/1000;
	else HTInt=0;
	if(HELIX_ON)HelType='T';
	else HelType='I';
	if(Whole_Struc==1)NameL=sprintf(fstring,"%c%dX%d_R%g-%gB%g-%gT%gG%g_L%gA%g_N%d%c%d",
		TypeChar,MTimes,Times,1.0*P_LJC,1.0*P_LJS,1.0*KLong,1.0*KLat,1.0*TLong,1.0*GLat,1.0*LCutEd,ACutEd*180/pi,N,HelType,LCompute);
	else if(Whole_Struc==2)NameL=sprintf(fstring,"F%dX%d_R%g-%gB%g-%gT%gG%g_L%gA%g_N%d%c%d",
		MTimes,Times,1.0*P_LJC,1.0*P_LJS,1.0*KLong,1.0*KLat,1.0*TLong,1.0*GLat,1.0*LCutEd,ACutEd*180/pi,N,HelType,LCompute);
	else NameL=sprintf(fstring,"%c%dX%d_R%g-%gB%g-%gT%gG%g_L%gA%g_N%d.%d",
		TypeChar,MTimes,Times,1.0*P_LJC,1.0*P_LJS,1.0*KLong,1.0*KLat,1.0*TLong,1.0*GLat,1.0*LCutEd,ACutEd*180/pi,N/Lat_Length,Lat_Length);
	if(TwistType==2)NameL+=sprintf(fstring+NameL,"_TT");
	NameL+=sprintf(fstring+NameL,"_EL%g",1.0*ELatAA);
	{
		char typel,typer;
		switch(BCType&0xF0){
			case 0x00:typel=0;break;
			case 0x10:typel='H';break;
			case 0x20:typel='S';break;
			case 0x40:typel='F';break;
			default:typel='N';
		}
		switch(BCType&0x0F){
			case 0x00:typer=0;break;
			case 0x01:typer='H';break;
			case 0x02:typer='S';break;
			case 0x04:typer='F';break;
			default:typer='N';
		}
		if(typel || typer){
			NameL+=sprintf(fstring+NameL,"_");
			if(typel)NameL+=sprintf(fstring+NameL,"l%c",typel);
			if(typer)NameL+=sprintf(fstring+NameL,"r%c",typer);
		}
	}
	if(ForceXYZ&0xF0){
		NameL+=sprintf(fstring+NameL,"_F");
		switch(ForceXYZ&0xF0){
			case 0x10:TypeChar='U';break;
			case 0x20:TypeChar='+';break;
			case 0x40:TypeChar='-';break;
			case 0x80:TypeChar='S';break;
			default:TypeChar='N';
		}
		NameL+=sprintf(fstring+NameL,"%c",TypeChar);
		switch(ForceXYZ&0x0F){
			case 0x01:TypeChar='x';break;
			case 0x02:TypeChar='y';break;
			case 0x04:TypeChar='z';break;
			case 0x08:TypeChar='t';break;
			default:TypeChar='N';
		}
		NameL+=sprintf(fstring+NameL,"%c",TypeChar);
	}
	if(DeltaL!=0.0 || DeltaXi!=0.0 || DeltaA!=0.0){
		NameL+=sprintf(fstring+NameL,"_Dt%ds%dd%d",DeltaT,Tst,Ted);
		if(DeltaL!=0.0)NameL+=sprintf(fstring+NameL,"L%g",1.0*DeltaL);
		if(DeltaXi!=0.0)NameL+=sprintf(fstring+NameL,"X%g",1.0*DeltaXi);
		if(DeltaA!=0.0)NameL+=sprintf(fstring+NameL,"A%g",1.0*DeltaA);
	}
	if(CAP_ON)NameL+=sprintf(fstring+NameL,"_Cap%d",CAP_ON);
	if(FR_ON||MR_ON)NameL+=sprintf(fstring+NameL,"_T%g",1.0*(MD_ON?TReal:KT));
	if(MD_ON){
		if(FR_ON)NameL+=sprintf(fstring+NameL,"_Fr");
		if(MR_ON)NameL+=sprintf(fstring+NameL,"_Mr");
	}
	else NameL+=sprintf(fstring+NameL,"_R%gA%g",1.0*StepR,1.0*StepA);
	if(W_Out!=1)NameL+=sprintf(fstring+NameL,"_Out%d",W_Out);
	if(W_Lat_aa!=0.1811)NameL+=sprintf(fstring+NameL,"_WL%g_%g",1.0*W_Lat_aa,1.0*W_Lat_ab);
	if(Test_Type)NameL+=sprintf(fstring+NameL,"_%c%dX%g",Test_Type,SecSteps,1.0*SecLen);
	if(Test_Type=='Y' || Test_Type=='K')NameL+=sprintf(fstring+NameL,".R%d",BallR);
	if(Test_Type=='T')NameL+=sprintf(fstring+NameL,".F%dL%d",CAPBegin,CAPNo);
	if(XI_0!=5.0)NameL+=sprintf(fstring+NameL,"_xi%g",1.0*XI_0);
	if(PR_ON)NameL+=sprintf(fstring+NameL,"_Pr%d",PR_ON);
	if(GROW_ON || BREAK_ON || HD_ON){
		NameL+=sprintf(fstring+NameL,"_");
		if(GROW_ON)NameL+=sprintf(fstring+NameL,"G");
		if(GROW_ON!=1)NameL+=sprintf(fstring+NameL,"%d",GROW_ON);
		if(BREAK_ON)NameL+=sprintf(fstring+NameL,"B");
#if GROW_ON==1 || GROW_ON==3
		NameL+=sprintf(fstring+NameL,"_G%g-%g.%dk",1.0*G_rate,1.0*EBase,GSteps/1000);
#elif GROW_ON==2
		NameL+=sprintf(fstring+NameL,"_G%g.%g.%g.%dk",1.0*Rate_In,1.0*Rate_Side,1.0*Rate_On,GSteps/1000);
#endif
		if(GROW_ON && G_Slow)NameL+=sprintf(fstring+NameL,"S%d",G_Slow);
		if(BREAK_ON)NameL+=sprintf(fstring+NameL,"C%g.%dW",1.0*C_rate,CTime/10000);
		if(BREAK_ON)NameL+=sprintf(fstring+NameL,"S%gB%g",1.0*S_rate,1.0*B_rate);
#if MustGrow
		NameL+=sprintf(fstring+NameL,"F");
#endif
		if(HD_ON==1)NameL+=sprintf(fstring+NameL,"HD%g",1.0*HD_rate);
		else if(HD_ON==2)NameL+=sprintf(fstring+NameL,"HDD%g",1.0*Afree);//dynamic
		if(HD_ON && HD_Delay)NameL+=sprintf(fstring+NameL,"D%dk",HD_Delay/1000);
	}
	else if(LongRFTest)NameL+=sprintf(fstring+NameL,"_HD%g",1.0*HD_rate);
	if(GRStep)NameL+=sprintf(fstring+NameL,"_Gr%d",GRStep);
	if(RLW_ON){
		NameL+=sprintf(fstring+NameL,"_RLW%gN%d",1.0*RandLatP,GammaN);
		if(RLW_ON<10)NameL+=sprintf(fstring+NameL,"N");
		if(RLW_ON>=10)NameL+=sprintf(fstring+NameL,"Z");
		if(RLW_ON>=15 && RLW_ON<20)NameL+=sprintf(fstring+NameL,"H1");
		else if(RLW_ON>=20 && RLW_ON<25)NameL+=sprintf(fstring+NameL,"H2");
		else if(RLW_ON>=25 && RLW_ON<30)NameL+=sprintf(fstring+NameL,"H3");
		else if(RLW_ON>=30 && RLW_ON<35)NameL+=sprintf(fstring+NameL,"H4");
	}
	if(EOut_ON)NameL+=sprintf(fstring+NameL,"_EO");
	if(RS_OUT_ON)NameL+=sprintf(fstring+NameL,"_RS");
	if(FilaSet!=13)NameL+=sprintf(fstring+NameL,"_IF%d",FilaSet);
	if(FilaMis!=0)NameL+=sprintf(fstring+NameL,"M%d",FilaMis);
	if(TipPf!=FilaSet)NameL+=sprintf(fstring+NameL,"_%dTPF%dL%d",TipPfs,TipPf,TipLayer);
	if(TipSinA>0){
		if(TipType==1)NameL+=sprintf(fstring+NameL,"_TM");
		else if(TipType==2)NameL+=sprintf(fstring+NameL,"_TP");
		else if(TipType==3)NameL+=sprintf(fstring+NameL,"_TPS");
		else if(TipType==4)NameL+=sprintf(fstring+NameL,"_TD");
		else if(TipType==5)NameL+=sprintf(fstring+NameL,"_TL");
		else if(TipType==6)NameL+=sprintf(fstring+NameL,"_TAD");
		else NameL+=sprintf(fstring+NameL,"_T");
		NameL+=sprintf(fstring+NameL,"A%dL%g",TipSinA,1.0*TipSinL);
		if(TipType==4&&TipARatio>0)NameL+=sprintf(fstring+NameL,"R%g",1.0*TipARatio);
		if((TipType==1||TipType==3||TipType==5||TipType==6)&&TipSinPhase!=0)NameL+=sprintf(fstring+NameL,"t%g",1.0*TipSinPhase);
		if(TipSinH>0)NameL+=sprintf(fstring+NameL,"H%d",TipSinH);
		if(PhaseFile>0)NameL+=sprintf(fstring+NameL,"F%d",PhaseFile);
	}
	if(RFLayers && RFUnits)NameL+=sprintf(fstring+NameL,"_%dRF%dL%d",RFType,RFUnits,RFLayers);
	if(RFRate)NameL+=sprintf(fstring+NameL,"_RFR%g",1.0*RFRate);
	if(TipEnd)NameL+=sprintf(fstring+NameL,"_T%dA%dP%d",TipEnd,TipAll,TipEndPf);
	if(LongRFTest)NameL+=sprintf(fstring+NameL,"_LRFT%d",LongRFTest);
	
	if(QDDelta)NameL+=sprintf(fstring+NameL,"_Q%gT%0.2X",1.0*QDDelta,QDType);
	if(FixLink){
		NameL+=sprintf(fstring+NameL,"_fl");
		if(cutoff_s<cutoff)NameL+=sprintf(fstring+NameL,"%g",cutoff_s);
	}
	if(OutputTime){
		tnow=time(0);
		strftime(tstring,sizeof(tstring),"%Y%m%d%H%M%S",localtime(&tnow));
		NameL+=sprintf(fstring+NameL,"_%s",tstring);
	}
	fprintf(stderr,"Filename length: %d\n",NameL);
	if(NameL>MaxFName){
		fprintf(stderr,"Filename too Long!\n");
		exit(0);
	}
	fprintf(stdout,"Done\n");
	fprintf(stdout,"File Name:\n%s\n",fstring);
}

void Set_Poten(){
	int i;
	double delta,para_a,para_b,para_c;
	char filename[MaxFName];
	FILE *ftest;
//常用三角函数值 	
	
	fprintf(stdout,"Generating Potential Lists...");
	CosBeta=cos(beta);
	CosGamma=cos(gamma0);	
	SinBeta=sin(beta);
	SinGamma=sin(gamma0);
	
	if(RLW_ON && myid==0){
		sprintf(filename,"Gamma_%s.txt",fstring);
		ftest=fopen(filename,"wt+");
	}
	for(i=0;i<100;i++){
		RGamma[i]=1+RandLatP*(GammaRand()-1);
		if(i==9 && RLW_ON>=15)RGamma[i]=1.5;
		if(RLW_ON>=20 && i==8)RGamma[i]=1.5;
		if(RLW_ON>=25 && i==10)RGamma[i]=1.5;
		if(RLW_ON>=30 && i==11)RGamma[i]=1.5;
		if(RLW_ON>=10 && myid==0)fprintf(ftest,"%4d\t%g\n",i,RGamma[i]);
	}
	if(RLW_ON && myid==0)fclose(ftest);	
	for(i=0;i<=T_Hydrolyze;i++){
		delta=1.0*i/T_Hydrolyze;
		cosxi_z[i]=cos((1-delta)*xi0/2+delta*xi1/2);
		cosa_z[i]=cos((1-delta)*alpha0+delta*alpha1);
		coso_z[i]=cos((1-delta)*omiga0+delta*omiga1);
		sinxi_z[i]=sin((1-delta)*xi0/2+delta*xi1/2);
		sina_z[i]=sin((1-delta)*alpha0+delta*alpha1);
		sino_z[i]=sin((1-delta)*omiga0+delta*omiga1);
	}	

//相关函数数组 

	for(i=0;i<=(int)(2*NLen);i++){
		ansACos[i]=acos(i*deltaexp-1);
	}
	ansACos[i]=ansACos[i-1];
	for(i=0;i<=(int)(8*NLen);i++){
		ansSin[i]=sin(i*deltaexp-4);
	}
	ansSin[i]=ansSin[i-1];
	for(i=0;i<=(int)(8*NLen);i++){
		ansCos[i]=cos(i*deltaexp-4);
	}
	ansCos[i]=ansCos[i-1];
	
	for(i=0;i<=(int)(2*NLen);i++){
		delta=acos(i*deltalj-1);
		if(delta<=ACutSt)
			ansADcyC[i]=1;
		else if(delta>=ACutEd)ansADcyC[i]=0;
		else {
			delta=sin(pi/(ACutEd-ACutSt)*(delta-ACutSt)/2.0);
			delta=delta*delta;
			ansADcyC[i]=1-delta*delta;
		}
	}
	ansADcyC[i]=ansADcyC[i-1];
	
	para_a=1.0/ACutEd/ACutEd;
	for(i=0;i<=(int)(2*NLen);i++){
		delta=acos(i*deltaexp-1);
		ansPAC[i]=delta*delta;
	}
	ansPAC[i]=ansPAC[i-1];
	for(i=0;i<=(int)(cutoff*NLen);i++){
		delta=i*deltalj;
		if(delta<=LCutSt)
			ansRDcyC[i]=1;
		else if(delta>=LCutEd)ansRDcyC[i]=0;
		else {
			delta=sin(pi/(LCutEd-LCutSt)*(delta-LCutSt)/2.0);
			delta=delta*delta;
			ansRDcyC[i]=1-delta*delta;
		}
	}
	ansRDcyC[i]=ansRDcyC[i-1];
	
	para_a=P_LJC;
	for(i=0;i<=(int)(cutoff/deltalj)+1;i++){
		delta=i*deltalj;
		ansPrC[i]=pow((1-exp(-para_a*(delta-2))),2)-1;
	}
	para_a=P_LJS;
	for(i=0;i<=(int)(cutoff/deltalj)+1;i++){
		delta=i*deltalj;
		ansPrS[i]=pow((1-exp(-para_a*(delta-2))),2)-1;
	}
	fprintf(stdout,"Done\n");
}

void Set_CellList(){
	int i,j,k,l,xup,xdown,yup,ydown,zup,zdown;
//网格邻居列表生成 & MPI数据传送表 	

	fprintf(stdout,"Generating CellList...");
	for(i=0;i<13;i++){
		Top[i]=-1;
		TubeEnd[i]=-1;
	}
	for(i=0,j=0;j<NCellZ;j++){
		if(j%NCellZ==0){
			zup=0;
			zdown=Nbox;
		}
		else if(j%NCellZ==NCellZ-1){
			zup=-Nbox;
			zdown=0;
		}
		else{
			zup=0;
			zdown=0;
		}
		for(k=0;k<NCellY;k++){
			if(k%NCellY==0){
				yup=0;
				ydown=NCsq;
			}
			else if(k%NCellY==NCellY-1){
				yup=-NCsq;
				ydown=0;
			}
			else{
				yup=0;
				ydown=0;
			}
			for(l=0;l<NCellX;l++,i++){
				if(l%NCellX==0){
					xup=0;
					xdown=NCellX;
				}
				else if(l%NCellX==NCellX-1){
					xup=-NCellX;
					xdown=0;
				}
				else{
					xup=0;
					xdown=0;
				}	
				Head[i]=-1;
				CList[26*i+0]=i+1+xup;
				CList[26*i+1]=i+NCellX-1+xdown+yup;
				CList[26*i+2]=i+NCellX+yup;
				CList[26*i+3]=i+NCellX+1+xup+yup;
				CList[26*i+4]=i+NCsq-NCellX-1+xdown+ydown+zup;
				CList[26*i+5]=i+NCsq-NCellX+ydown+zup;
				CList[26*i+6]=i+NCsq-NCellX+1+xup+ydown+zup;
				CList[26*i+7]=i+NCsq-1+xdown+zup;
				CList[26*i+8]=i+NCsq+zup;
				CList[26*i+9]=i+NCsq+1+xup+zup;
				CList[26*i+10]=i+NCsq+NCellX-1+xdown+yup+zup;
				CList[26*i+11]=i+NCsq+NCellX+yup+zup;
				CList[26*i+12]=i+NCsq+NCellX+1+xup+yup+zup;
				CList[26*i+13]=i-1+xdown;
				CList[26*i+14]=i-NCellX-1+xdown+ydown;
				CList[26*i+15]=i-NCellX+ydown;
				CList[26*i+16]=i-NCellX+1+xup+ydown;
				CList[26*i+17]=i-NCsq-NCellX-1+xdown+ydown+zdown;
				CList[26*i+18]=i-NCsq-NCellX+ydown+zdown;
				CList[26*i+19]=i-NCsq-NCellX+1+xup+ydown+zdown;
				CList[26*i+20]=i-NCsq-1+xdown+zdown;
				CList[26*i+21]=i-NCsq+zdown;
				CList[26*i+22]=i-NCsq+1+xup+zdown;
				CList[26*i+23]=i-NCsq+NCellX-1+xdown+yup+zdown;
				CList[26*i+24]=i-NCsq+NCellX+yup+zdown;
				CList[26*i+25]=i-NCsq+NCellX+1+xup+yup+zdown;

			}
		}
	}
	fprintf(stdout,"Done\n");
}

void Set_NeighbourList(int InitMark){
	int i,j,k,p,q;
	int icell,MyNbrNo,MaxNbrNum=0;
	double MyR1[3],MyR2[3],r[4][3],rmod[4];
	tubulin *Tubp,*Tub_j;
	
	//初始化 
	fprintf(stderr,"Generating NeighbourList (%d)...",InitMark);
	if(!InitMark)NbrList=(int *)malloc(Max_Nbr_Num*N*sizeof(int));
	for(i=0;i<N;i++){
		NbrList[Max_Nbr_Num*i]=0;
	}
	//填充邻居列表 
	for(i=0;i<N;i++){
		Tubp=&Tub[i];
		if(!Tubp->Active)continue;
		icell=GetCellNo(Tubp->r);
		for(k=0;k<3;k++){
			MyR1[k]=Tubp->r1[k];
			MyR2[k]=Tubp->r2[k];
		}
		MyNbrNo=Max_Nbr_Num*i;
		for(p=0;p<26;p++){
			q=CList[26*icell+p];
			for(j=Head[q];j!=-1;j=List[j]){
				if(j<=i)continue;
				Tub_j=&Tub[j];
				for(k=0;k<3;k++){
					r[0][k]=MyR1[k]-Tub_j->r1[k];
					r[1][k]=MyR1[k]-Tub_j->r2[k];
					r[2][k]=MyR2[k]-Tub_j->r1[k];
					r[3][k]=MyR2[k]-Tub_j->r2[k];
				}
				for(k=0;k<4;k++){
					PeriodBoundary(r[k]);
					rmod[k]=Dot(r[k],r[k]);
					if(rmod[k]<=Max_Nbr_SqDis)break;
					else continue;
				}
				if(k<4){
					NbrList[MyNbrNo+ ++NbrList[MyNbrNo]]=j;
					k=Max_Nbr_Num*j;
					NbrList[k+ ++NbrList[k]]=i;
				}
			}
		}
		for(j=Head[icell];j!=-1;j=List[j]){
			if(j<=i)continue;
			Tub_j=&Tub[j];
			for(k=0;k<3;k++){
				r[0][k]=MyR1[k]-Tub_j->r1[k];
				r[1][k]=MyR1[k]-Tub_j->r2[k];
				r[2][k]=MyR2[k]-Tub_j->r1[k];
				r[3][k]=MyR2[k]-Tub_j->r2[k];
			}
			for(k=0;k<4;k++){
				PeriodBoundary(r[k]);
				rmod[k]=Dot(r[k],r[k]);
				if(rmod[k]<=Max_Nbr_SqDis)break;
				else continue;
			}
			if(k<4){
				NbrList[MyNbrNo+ ++NbrList[MyNbrNo]]=j;
				k=Max_Nbr_Num*j;
				NbrList[k+ ++NbrList[k]]=i;
			}
		}
	}
	//验算 
	for(i=0;i<N;i++){
		k=Max_Nbr_Num*i;
		if(NbrList[k]>=Max_Nbr_Num){
			fprintf(stderr,"\nERROR: Too many neighbours for Unit %d (%d > Max_Nbr_Num-1)!\n",i,NbrList[k]);
			exit(-7);
		}
		if(MaxNbrNum<NbrList[k])MaxNbrNum=NbrList[k];
	}
	fprintf(stderr,"\tMaxNbrNum=%d\t",MaxNbrNum);
	fprintf(stderr,"Done\n");
}

void Set_ReinForce(void){
#if RFLayers && RFUnits
	int i,j,k=0,Nst,Ned,Hst,Ncount[RFLayers];
	char fname[MaxFName];
	double Ra,LRange,Pnow,PList[RFLayers*FilaSet],Psum=0.0;
	FILE *fout;
	
	fprintf(stdout,"Seting ReinForce Units...");
	Hst=N/FilaSet/2-RFLayers/2;
	Nst=13*Hst;
	Ned=Nst+RFLayers*FilaSet;
	if(RFUnits>=RFLayers*FilaSet){
		for(i=Nst;i<Ned;i++){
			Tub[i].Type=0;
			Tub[i].xi=xi0;
			Tub[i].RFMark=1;
		}
		return;
	}
	LRange=6.0/RFLayers;
	for(i=0;i<RFLayers;i++){
	#if RFType
		Pnow=(PGauss((-3+(i+1)*LRange))-PGauss((-3+i*LRange)))/FilaSet;
	#else
		Pnow=1.0/RFLayers/FilaSet;
	#endif
		for(j=0;j<FilaSet;j++){
			PList[i*FilaSet+j]=Pnow;
			Psum+=Pnow;
		}
	}
	for(i=0;i<RFLayers*FilaSet;i++)PList[i]/=Psum;
	
	for(i=0;i<RFLayers;i++)Ncount[i]=0;
	for(i=0;i<RFUnits;i++){
		Ra=rand()/(RAND_MAX+0.1);
		for(j=0;j<RFLayers*FilaSet;j++){
			if(Ra>0)Ra-=PList[j];
			else{
				j--;
				break;
			}
		}
		if(Tub[j+Nst].Type && j<RFLayers*FilaSet){
			Tub[j+Nst].Type=0;
			Tub[j+Nst].xi=xi0;
			Tub[j+Nst].RFMark=1;
			Ncount[j/FilaSet]++;
		}
		else i--;
		k++;
		if(k>1000)break;
	}
	fprintf(stdout,"Done\n");
#endif
	return;
}

void Set_QDRate(void){
	int i,j,k,p0,p1,p2;
	char fname[MaxFName];
	double QDList[10000+2],MyQD[10000+1],Dtmp,QDk,R1;
	tubulin *Tubp;
	FILE *fout;
	
	fprintf(stdout,"Generating QDRate...");
	QDList[0]=0;
	if(QDDelta){
		QDk=1.0/QDDelta;
		for(i=0;i<10001;i++){
			Dtmp=0.001*i;
			if(i)MyQD[i]=exp(QDk*log(QDk)+(QDk-1)*log(Dtmp)-QDk*Dtmp-lgamma(QDk));
			else MyQD[i]=0;
			if(i)QDList[i]=QDList[i-1]+0.001*(MyQD[i]+MyQD[i-1])/2;
		}
		QDList[10000+2]=1;
	}
	else{
		fprintf(stderr,"\nERROR: QDData should not be set when QDDelta==0.0 !!\n");
		exit(-1);
	}
	for(i=0;i<N;i++){
		Tubp=&Tub[i];
		if(!Tubp->Active)continue;
		if(QDType&0xF0){
			if((QDType&0x80) && (i%FilaSet)){
				k=i-i%FilaSet;
				for(j=0;j<2;j++)Tubp->QDRate[j]=Tub[k].QDRate[j];
			}
			else if((QDType&0x40) && (i/FilaSet)){
				k=i%FilaSet;
				for(j=0;j<2;j++)Tubp->QDRate[j]=Tub[k].QDRate[j];
			}
			else{
				for(j=0;j<2;j++){
					R1=1.0*rand()/(RAND_MAX+0.1);
					switch(QDType&0x30){
						case 0x20:Tubp->QDRate[j]==1.0+QDDelta*(2*R1-1);break;
						case 0x10:{
							p0=0;
							p2=10000+2;
							do{
								p1=(p0+p2)/2;
								if(R1>QDList[p1])p0=p1;
								if(R1<QDList[p1])p2=p1;
							}while(p2-p0>1);
							
							if(p0!=p2)Tubp->QDRate[j]=0.001*p0+0.001*(p2-p0)*(R1-QDList[p0])/(QDList[p2]-QDList[p0]);
							else Tubp->QDRate[j]=0.001*p1;
							break;
						}
					}
				}
				if(QDType&0x40)Tubp->QDRate[1]=Tubp->QDRate[0];
			}
		}
		if(QDType&0x0F){
			if((QDType&0x08) && (i%FilaSet)){
				k=i-i%FilaSet;
				for(j=2;j<6;j++)Tubp->QDRate[j]=Tub[k].QDRate[j];
			}
			else if((QDType&0x04) && (i/FilaSet)){
				k=i%FilaSet;
				for(j=2;j<6;j++)Tubp->QDRate[j]=Tub[k].QDRate[j];
			}
			else{
				for(j=2;j<6;j++){
					R1=1.0*rand()/(RAND_MAX+0.1);
					switch(QDType&0x03){
						case 0x02:Tubp->QDRate[j]==1.0+QDDelta*(2*R1-1);break;
						case 0x01:{
							p0=0;
							p2=10000+2;
							do{
								p1=(p0+p2)/2;
								if(R1>QDList[p1])p0=p1;
								if(R1<QDList[p1])p2=p1;
							}while(p2-p0>1);
							
							if(p0!=p2)Tubp->QDRate[j]=0.001*p0+0.001*(p2-p0)*(R1-QDList[p0])/(QDList[p2]-QDList[p0]);
							else Tubp->QDRate[j]=0.001*p1;
							break;
						}
					}
				}
				if(QDType&0x08){
					Tubp->QDRate[3]=Tubp->QDRate[2];
					Tubp->QDRate[5]=Tubp->QDRate[4];
				}
				else if(QDType&0x04){
					Tubp->QDRate[4]=Tubp->QDRate[2];
					Tubp->QDRate[5]=Tubp->QDRate[3];
				}
			}
		}
	}
	fprintf(stdout,"Done\n");
	return;
}

int Set_Position(int i,double *r,double *angle,int FilaMisSt,int FilaMisEd){
	int j,k,l;
	double Dtmp,MyHeight;
	if(Whole_Struc==1){
		if(TipSinA>0){
			if(HELIX_ON) r[0]=(i/FilaSet)*4.0+(i%FilaSet)*6.0/FilaSet-BoxLX+5;
			else r[0]=(i/FilaSet)*4.0-BoxLX+5;
			r[1]=-cos((i%FilaSet)*pi*2/FilaSet)*InitR;
			r[2]=sin((i%FilaSet)*pi*2/FilaSet)*InitR;
	
			angle[0]=1;
			angle[1]=0;//-(i%FilaSet)*pi*2/FilaSet;
			angle[2]=0;
			Dtmp=-(i%FilaSet)*pi*2/FilaSet/2;
			angle[3]=sin(Dtmp);
			angle[4]=cos(Dtmp);
			
			if(i<NInit)return 1;
			else{
				if(!TipType){
					MyHeight=TipSinA*(sin((i%FilaSet-0.5*FilaSet)*pi*2/TipSinL+pi/2)+1)/2+TipSinH;
					fprintf(stderr,"i: %d\tH: %g\tZ: %g\n",i,MyHeight,(i-NInit)/FilaSet+0.5);
					if((i-NInit)/FilaSet+0.5<MyHeight)return 1;
					else return 0;
				}
				else if(TipType==1){
					int L1,L2;
					L1=(int)TipSinL;
					L2=L1+1;
					Dtmp=(i%FilaSet)*pi*2/FilaSet;
					MyHeight=TipSinA*((sqrt(L2-TipSinL)*sin(L1*Dtmp)+sqrt(TipSinL-L1)*sin(L2*Dtmp+TipSinPhase*pi))-MinA)/(MaxA-MinA)+TipSinH;
					fprintf(stderr,"i: %d\tH: %g\tZ: %g\n",i,MyHeight,(i-NInit)/FilaSet+0.5);
					if((i-NInit)/FilaSet+0.5<MyHeight)return 1;
					else return 0;
				}
				else if(TipType==2){
					Dtmp=sin((i%FilaSet)*pi/FilaSet-pi/4);
					MyHeight=TipSinA*(1-pow(Dtmp*Dtmp,TipSinL))+TipSinH;
					fprintf(stderr,"i: %d\tH: %g\tZ: %g\n",i,MyHeight,(i-NInit)/FilaSet+0.5);
					if((i-NInit)/FilaSet+0.5<MyHeight)return 1;
					else return 0;
				}
				else if(TipType==3){
					Dtmp=sin((i%FilaSet)*pi/FilaSet+TipSinPhase*pi);
					MyHeight=TipSinA*(1-pow(Dtmp*Dtmp,TipSinL))+TipSinH;
					fprintf(stderr,"i: %d\tH: %g\tZ: %g\n",i,MyHeight,(i-NInit)/FilaSet+0.5);
					if((i-NInit)/FilaSet+0.5<MyHeight)return 1;
					else return 0;
				}
				else if(TipType==4){
					int MyPos=0,PeakNo,PeakX;
					double PeakL=0;
					PeakL=TipSinL*FilaSet;
					MyPos=i%FilaSet;
					if(MyPos<0.5*FilaSet)PeakNo=MyPos;
					else if(FilaSet-MyPos<0.5*FilaSet)PeakNo=FilaSet-MyPos;
					if(PeakNo<0.5*PeakL){
						PeakX=PeakNo;
						Dtmp=(cos(PeakX*2*pi/PeakL)+1)*0.5;
						MyHeight=TipSinA*Dtmp+TipSinH;
					}
					else {
						PeakX=0.5*FilaSet-PeakNo;
						Dtmp=(cos(PeakX*2*pi/(FilaSet-PeakL))+1)*0.5;
						MyHeight=TipSinA*TipARatio*Dtmp+TipSinH;
					}
					fprintf(stderr,"i: %d\tH: %g\tZ: %g\n",i,MyHeight,(i-NInit)/FilaSet+0.5);
					if((i-NInit)/FilaSet+0.5<MyHeight)return 1;
					else return 0;
				}
				else if(TipType==5){
					int L1,L2;
					L1=(int)TipSinL;
					L2=L1+1;
					Dtmp=(i%FilaSet)*pi*2/FilaSet;
					MyHeight=TipSinA*(((L2-TipSinL)*sin(L1*Dtmp)+(TipSinL-L1)*sin(L2*Dtmp+TipSinPhase*pi))-MinA)/(MaxA-MinA)+TipSinH;
					fprintf(stderr,"i: %d\tH: %g\tZ: %g\n",i,MyHeight,(i-NInit)/FilaSet+0.5);
					if((i-NInit)/FilaSet+0.5<MyHeight)return 1;
					else return 0;
				}
				else if(TipType==6){
					double A[6],Htmp;
					Dtmp=(i%FilaSet)*pi*2/FilaSet;
					Htmp=0;
					for(k=0;k<6;k++){
						A[k]=sqrt(pow(k+1,-TipSinL));
						Htmp+=A[k]*sin((k+1)*Dtmp+MyTipPhase[k]*pi);
					}
					MyHeight=TipSinA*(Htmp-MinA)/(MaxA-MinA)+TipSinH;
					fprintf(stderr,"i: %d\tH: %g\tZ: %g\n",i,MyHeight,(i-NInit)/FilaSet+0.5);
					if((i-NInit)/FilaSet+0.5<MyHeight)return 1;
					else return 0;
				}
			}
		}
		else if(FilaMisSt>0){
			if(i<FilaMisSt){
				if(HELIX_ON) r[0]=(i/13)*4.0+(i%13)*6.0/13-BoxLX+5;
				else r[0]=(i/13)*4.0-BoxLX+5;
				r[1]=-cos((i%13)*pi*2/13)*InitR;
				r[2]=sin((i%13)*pi*2/13)*InitR;
	
				angle[0]=1;
				angle[1]=0;//-(i%13)*pi*2/13;
				angle[2]=0;
				Dtmp=-(i%13)*pi*2/13/2;
				angle[3]=sin(Dtmp);
				angle[4]=cos(Dtmp);
			}
			else if(i>=FilaMisSt && i<FilaMisEd){
				k=(i-FilaMisSt)/FilaSet+FilaMisSt/13;
				l=(i-FilaMisSt)%FilaSet;
				if(HELIX_ON) r[0]=k*4.0+l*6.0/FilaSet-BoxLX+5;
				else r[0]=k*4.0-BoxLX+5;
				r[1]=-cos(l*pi*2/FilaSet)*InitR;
				r[2]=sin(l*pi*2/FilaSet)*InitR;
	
				angle[0]=1;
				angle[1]=0;//-(i%FilaSet)*pi*2/FilaSet;
				angle[2]=0;
				Dtmp=-l*pi*2/FilaSet/2;
				angle[3]=sin(Dtmp);
				angle[4]=cos(Dtmp);
			}
			else if(i>=FilaMisEd){
				k=(i-FilaMisEd)/13+FilaMisSt/13+FilaMis;
				l=(i-FilaMisEd)%13;
				if(HELIX_ON) r[0]=k*4.0+l*6.0/13-BoxLX+5;
				else r[0]=k*4.0-BoxLX+5;
				r[1]=-cos(l*pi*2/13)*InitR;
				r[2]=sin(l*pi*2/13)*InitR;
	
				angle[0]=1;
				angle[1]=0;//-(i%13)*pi*2/13;
				angle[2]=0;
				Dtmp=-l*pi*2/13/2;
				angle[3]=sin(Dtmp);
				angle[4]=cos(Dtmp);
			}
		}
		else{
			if(HELIX_ON) r[0]=(i/FilaSet)*4.0+(i%FilaSet)*6.0/FilaSet-BoxLX+5;
			else r[0]=(i/FilaSet)*4.0-BoxLX+5;
			r[1]=-cos((i%FilaSet)*pi*2/FilaSet)*InitR;
			r[2]=sin((i%FilaSet)*pi*2/FilaSet)*InitR;
	
			angle[0]=1;
			angle[1]=0;//-(i%FilaSet)*pi*2/FilaSet;
			angle[2]=0;
			Dtmp=-(i%FilaSet)*pi*2/FilaSet/2;
			angle[3]=sin(Dtmp);
			angle[4]=cos(Dtmp);
		}
	}
	else if(Whole_Struc==2){
		if(HELIX_ON) r[0]=(i/FilaSet)*4.0+(i%FilaSet)*6.0/FilaSet-BoxLX+5;
		else r[0]=(i/FilaSet)*4.0-BoxLX+5;
		r[1]=-cos((i%FilaSet)*pi*2/FilaSet)*InitR;
		r[2]=sin((i%FilaSet)*pi*2/FilaSet)*InitR;

		angle[0]=1;
		angle[1]=0;//-(i%FilaSet)*pi*2/FilaSet;
		angle[2]=0;
		Dtmp=-(i%FilaSet)*pi*2/FilaSet/2;
		angle[3]=sin(Dtmp);
		angle[4]=cos(Dtmp);
	}
	else{
		r[0]=(i/Lat_Length)*4.02-BoxLX+5;
		r[1]=0;
		r[2]=(i%Lat_Length)*2.02-BoxLX+5;

		angle[0]=1;
		angle[1]=0;//-(i%13)*pi*2/13;
		angle[2]=0;
		angle[3]=0;//-sin(pi/4);
		angle[4]=1;//cos(pi/4);
	}
	return 1;
} 

void ReadTub(){
	int i,j,k;
	FILE *fin;
	char Ctmp,*Strtmp,StrHead[256];
	
	fprintf(stdout,"Reading TubInput.txt ...");
	fin=fopen("TubInput.txt","rt");
	fgets(StrHead,256,fin);
	sscanf(StrHead,"%d%d",&FilaSet,&MaxLayer);
	do{
		TubState=(char *)malloc(FilaSet*MaxLayer*sizeof(char));
		Strtmp=(char *)malloc(MaxLayer*sizeof(char));
	}while(!TubState || !Strtmp);
	i=0;
	while(fgets(Strtmp,MaxLayer,fin)){
		fprintf(stderr,"i: %d\t%s",i,Strtmp);
		for(j=0;j<MaxLayer;j++){
			k=i+FilaSet*j;
			Ctmp=Strtmp[j];
			if(Ctmp<'A' || Ctmp>'z' || (Ctmp<'a' && Ctmp>'Z')){
				TubState[k]='N';
				break;
			}
			else TubState[k]=Ctmp;
		}
		if(j<MaxLayer){
			for(;j<MaxLayer;j++){
				k=i+FilaSet*j;
				TubState[k]='N';
			}
		}
		i++;
		if(i>=FilaSet)break;
	}
	fclose(fin);
	free(Strtmp);
	fprintf(stdout,"Done.\n");
}

void ReadPhase(){
	int i,PhaseNo;
	double KeyPara,Dtmp;
	char Strtmp[512],StrList[512]; 
	FILE *fp;
	if(PhaseFile==1){
		fp=fopen("PhaseData.txt","rt");
		PhaseNo=abs((int)TipSinPhase);
		while(fgets(Strtmp,512,fp)){
			sscanf(Strtmp,"%lf%[^\n]",&KeyPara,StrList);
			if(KeyPara==TipSinL){
				if(PhaseNo--)continue;
				fprintf(stderr,"TipPhase:");
				for(i=0;i<6;i++){
					strcpy(Strtmp,StrList);
					sscanf(Strtmp,"%lf%[^\n]",&Dtmp,StrList);
					MyTipPhase[i]=Dtmp;
					fprintf(stderr,"\t%g",Dtmp);
				}
				fprintf(stderr,"\n");
			}
			else continue;
		}
		fclose(fp);
	}
	else if(PhaseFile==2){
		int Phase2,Phase3;
		fp=fopen("PhaseData07.txt","rt");
		PhaseNo=abs((int)TipSinPhase);
		Phase2=PhaseNo/100;
		Phase3=PhaseNo%100;
		PhaseNo=Phase3*7+Phase2;
		while(fgets(Strtmp,512,fp)){
			sscanf(Strtmp,"%lf%[^\n]",&KeyPara,StrList);
			if(KeyPara==TipSinL){
				if(PhaseNo--)continue;
				fprintf(stderr,"TipPhase:");
				for(i=0;i<6;i++){
					strcpy(Strtmp,StrList);
					sscanf(Strtmp,"%lf%[^\n]",&Dtmp,StrList);
					MyTipPhase[i]=Dtmp;
					fprintf(stderr,"\t%g",Dtmp);
				}
				fprintf(stderr,"\n");
			}
			else continue;
		}
		fclose(fp);
	} 
}

void move(double *r,double *angle,double *xi){	//随机移动 
	int No;
	double rmod,anglemod;
	const double RandPara=2.0/(RAND_MAX);
	r[0]=rand()*RandPara-1;
	r[1]=rand()*RandPara-1;
	r[2]=rand()*RandPara-1;
	rmod=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
	rmod=sqrt(rmod);
	r[0]=0.001*(r[0]/rmod);
	r[1]=0.001*(r[1]/rmod);
	r[2]=0.001*(r[2]/rmod);
	*xi=0.0001*(rand()*RandPara-1);

	angle[0]=rand()*RandPara-1;
	angle[1]=rand()*RandPara-1;
	angle[2]=rand()*RandPara-1;
	anglemod=angle[0]*angle[0]+angle[1]*angle[1]+angle[2]*angle[2];
	anglemod=sqrt(anglemod);
	angle[0]=angle[0]/anglemod;
	angle[1]=angle[1]/anglemod;
	angle[2]=angle[2]/anglemod;
	angle[3]=0.0005*rand()*RandPara;

}

void move2(double *r,double *angle,double StepLen,int k){	//随机移动 
	int No1,No2,No3;
	double DeltaNo1,DeltaNo2,DeltaNo3;
	double s1,c1,s2,c2;
	double rmod,anglemod,RNormal,RU1,RU2;
	const double RandPara=2.0/(RAND_MAX);
	if(StepLen<0.01)StepLen=0.01;
	if(k!=1){
		DeltaNo1=(2-rand()*RandPara)*NLen;
		No1=(int)(DeltaNo1);
		DeltaNo1-=No1;
		RU1=(1-DeltaNo1)*ansACos[No1]+(DeltaNo1)*ansACos[No1+1];
		RU2=(rand()*RandPara-1)*pi;
		DeltaNo1=(RU1+4)*NLen;
		No1=(int)(DeltaNo1);
		DeltaNo1-=No1;
		DeltaNo2=(RU2+4)*NLen;
		No2=(int)(DeltaNo2);
		DeltaNo2-=No2;
		s1=(1-DeltaNo1)*ansSin[No1]+(DeltaNo1)*ansSin[No1+1];
		s2=(1-DeltaNo2)*ansSin[No2]+(DeltaNo2)*ansSin[No2+1];
		c1=(1-DeltaNo1)*ansCos[No1]+(DeltaNo1)*ansCos[No1+1];
		c2=(1-DeltaNo2)*ansCos[No2]+(DeltaNo2)*ansCos[No2+1];
		RNormal=StepLen*StepR;
		r[0]=RNormal*c2*s1;
		r[1]=RNormal*s2*s1;
		r[2]=RNormal*c1;
	}
	if(k){
		DeltaNo1=(2-rand()*RandPara)*NLen;
		No1=(int)(DeltaNo1);
		DeltaNo1-=No1;
		RU1=(1-DeltaNo1)*ansACos[No1]+(DeltaNo1)*ansACos[No1+1];
		RU2=(rand()*RandPara-1)*pi;
		DeltaNo1=(RU1+4)*NLen;
		No1=(int)(DeltaNo1);
		DeltaNo1-=No1;
		DeltaNo2=(RU2+4)*NLen;
		No2=(int)(DeltaNo2);
		DeltaNo2-=No2;
		s1=(1-DeltaNo1)*ansSin[No1]+(DeltaNo1)*ansSin[No1+1];
		s2=(1-DeltaNo2)*ansSin[No2]+(DeltaNo2)*ansSin[No2+1];
		c1=(1-DeltaNo1)*ansCos[No1]+(DeltaNo1)*ansCos[No1+1];
		c2=(1-DeltaNo2)*ansCos[No2]+(DeltaNo2)*ansCos[No2+1];
		angle[0]=c2*s1;
		angle[1]=s2*s1;
		angle[2]=c1;
		RNormal=StepLen*StepA;
		angle[3]=RNormal;
		angle[4]=sqrt(1-RNormal*RNormal);
	}
}

double potential(tubulin *Tub_i,tubulin *Tub_j,double *PotenAll){	//粒子间能量计算  
	int k,l;
	int No=0,DNo[6];
	double r[4][3],rmod[4],Pr[4],PrS[4],RDcy[4];
	double angle[6]={0.0},Poten[12]={0.0};
	double single[6],ADcy[6],Delta[6],Dtmp;
	double AD01,AD23,AD45;
	
	if(!Tub_i->Active || !Tub_j->Active)return 0;
    if(Tub_i->Tub_No==Tub_j->Tub_No)return 0;
	{
		for(l=0;l<3;l++){
			r[0][l]=Tub_i->r1[l]-Tub_j->r1[l];
			r[1][l]=Tub_i->r1[l]-Tub_j->r2[l];
			r[2][l]=Tub_i->r2[l]-Tub_j->r1[l];
			r[3][l]=Tub_i->r2[l]-Tub_j->r2[l];
		}
		for(k=0;k<4;k++){
			PeriodBoundary(r[k]);
			rmod[k]=Dot(r[k],r[k]);
			if(rmod[k]<=cutoffsq){
				rmod[k]=sqrt(rmod[k]);
				Delta[k]=1/rmod[k];
				for(l=0;l<3;l++)r[k][l]*=Delta[k];
				Delta[k]=rmod[k]*NLen;
				DNo[k]=(int)(Delta[k]);
				Delta[k]-=DNo[k];
				Pr[k]=(1-Delta[k])*ansPrC[DNo[k]]+(Delta[k])*ansPrC[DNo[k]+1];
				PrS[k]=(1-Delta[k])*ansPrS[DNo[k]]+(Delta[k])*ansPrS[DNo[k]+1];
				RDcy[k]=(1-Delta[k])*ansRDcyC[DNo[k]]+(Delta[k])*ansRDcyC[DNo[k]+1];
			}
			else{
				Pr[k]=0;
				PrS[k]=0;
				RDcy[k]=0;
			}
		}
		if(rmod[0]<=cutoff){
			
			angle[0]=-Dot(Tub_i->side1,r[0]);
			angle[1]=Dot(Tub_j->side2,r[0]);
			angle[2]=Dot(Tub_j->side1,r[0]);
			angle[3]=-Dot(Tub_i->side2,r[0]);
			for(k=0;k<4;k++){
				Delta[k]=(angle[k]+1)*NLen;
				DNo[k]=(int)(Delta[k]);
				Delta[k]-=DNo[k];
			}
			for(k=0;k<4;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
			AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[2]+Tub_j->QDRate[3])/2;
			AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[3]+Tub_j->QDRate[2])/2;
			if(RDcy[0]>0 && (AD01>0 || AD23>0)){
				double Dcni,Dcnj,Dtitj,Dnitj,Dnjti,MyAD;
				double *Pi,*Pj,*Ni,*Nj,*Cij,DTheta[6],Theta[6];
				if(AD01>0){
					Pi=Tub_i->side1;
					Ni=Tub_i->NS1;
					Pj=Tub_j->side2;
					Nj=Tub_j->NS2;
					Cij=r[0];
					MyAD=AD01;
				}
				else {
					Pi=Tub_i->side2;
					Ni=Tub_i->NS2;
					Pj=Tub_j->side1;
					Nj=Tub_j->NS1;
					Cij=r[0];
					MyAD=AD23;
				}
				Dcni=Dot(Cij,Ni);
				Dcnj=Dot(Cij,Nj);
				Dtitj=Dot(Tub_i->t,Tub_j->t);
				Dnitj=Dot(Tub_j->t,Ni);
				Dnjti=Dot(Tub_i->t,Nj);
				DTheta[0]=sqrt(1-Dcni*Dcni);			//Theta i bend
				DTheta[1]=-Dot(Cij,Pi)/DTheta[0];		//Theta i in plane
				DTheta[2]=Dtitj/sqrt(1-Dnitj*Dnitj);	//Theta titj in plane
				DTheta[3]=sqrt(1-Dcnj*Dcnj);			//Theta j bend
				DTheta[4]=Dot(Cij,Pj)/DTheta[3];		//Theta j in plane
				DTheta[5]=Dtitj/sqrt(1-Dnjti*Dnjti);	//Theta titj in plane
				for(k=0;k<6;k++){
					Delta[k]=(DTheta[k]+1)*NLen;
					DNo[k]=(int)(Delta[k]);
					Delta[k]-=DNo[k];
				}
				for(k=0;k<6;k++)Theta[k]=(1-Delta[k])*ansACos[DNo[k]]+(Delta[k])*ansACos[DNo[k]+1];
				if(Dot(Tub_j->t,Pi)<0)Theta[2]=-Theta[2];
				if(Dot(Tub_i->t,Pj)<0)Theta[5]=-Theta[5];
				single[0]=Theta[0]*Theta[0];
				single[1]=Theta[3]*Theta[3];
				single[2]=Theta[1]-0.5*Theta[2];
				single[2]*=single[2];
				single[3]=Theta[4]-0.5*Theta[5];
				single[3]*=single[3];
				Dtmp=RDcy[0]*MyAD;
				Poten[4]+=Dtmp*KLat*(single[0]+single[1]);
				Poten[5]+=Dtmp*GLat*(single[2]+single[3]);
			}
			if(PrS[0]<0){
				Poten[3]+=ELatAA*PrS[0]*(AD01+AD23);
			}
			else Poten[3]+=ELatAA*PrS[0];
			
			if(AD01>0 && PrS[0]*AD01<-0.35){
				Tub_i->right=Tub_j->Tub_No;
				Tub_j->left=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->right==Tub_j->Tub_No)Tub_i->right=-1;
				if(Tub_j->left==Tub_i->Tub_No)Tub_j->left=-1;
			}
			if(AD23>0 && PrS[0]*AD23<-0.35){
				Tub_i->left=Tub_j->Tub_No;
				Tub_j->right=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->left==Tub_j->Tub_No)Tub_i->left=-1;
				if(Tub_j->right==Tub_i->Tub_No)Tub_j->right=-1;
			}
		}
		if(rmod[1]<=cutoff){
			angle[0]=-Dot(Tub_i->side1,r[1]);
			angle[1]=Dot(Tub_j->side4,r[1]);
			angle[2]=-Dot(Tub_i->side2,r[1]);
			angle[3]=Dot(Tub_j->side3,r[1]);
			angle[4]=-Dot(Tub_i->chain1,r[1]);
			angle[5]=Dot(Tub_j->chain2,r[1]);
			for(k=0;k<6;k++){
				Delta[k]=(angle[k]+1)*NLen;
				DNo[k]=(int)(Delta[k]);
				Delta[k]-=DNo[k];
			}
			for(k=0;k<6;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
			
			AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[2]+Tub_j->QDRate[5])/2;
			AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[3]+Tub_j->QDRate[4])/2;
			AD45=ADcy[4]*ADcy[5]*(Tub_i->QDRate[0]+Tub_j->QDRate[1])/2;
			if(RDcy[1]>0 && AD45>0){
				double t0[3],Dn1n2,Dn1t,Dn2t,Dtw,Atw,DeltaAtw;
				int NoAtw;
				for(k=0;k<3;k++)t0[k]=r[1][k];
				Dn1n2=Dot(Tub_i->n,Tub_j->n);
				Dn1t=Dot(t0,Tub_i->n);
				Dn2t=Dot(t0,Tub_j->n);
				Dtw=(Dn1n2-Dn1t*Dn2t)/sqrt((1-Dn1t*Dn1t)*(1-Dn2t*Dn2t));
				DeltaAtw=(Dtw+1)*NLen;
				NoAtw=(int)(DeltaAtw);
				DeltaAtw-=NoAtw;
				Atw=(1-DeltaAtw)*ansPAC[NoAtw]+(DeltaAtw)*ansPAC[NoAtw+1];
				single[4]=(1-Delta[4])*ansPAC[DNo[4]]+(Delta[4])*ansPAC[DNo[4]+1];
				single[5]=(1-Delta[5])*ansPAC[DNo[5]]+(Delta[5])*ansPAC[DNo[5]+1];
				Dtmp=RDcy[1]*AD45;
				Poten[1]+=Dtmp*KLong*(single[4]+single[5]);
				Poten[2]+=Dtmp*TLong*Atw/(rmod[1]*rmod[1]);
			}
			if(RDcy[1]>0 && (AD01>0 || AD23>0)){
				double Dcni,Dcnj,Dtitj,Dnitj,Dnjti,MyAD;
				double *Pi,*Pj,*Ni,*Nj,*Cij,DTheta[6],Theta[6];
				if(AD01>0){
					Pi=Tub_i->side1;
					Ni=Tub_i->NS1;
					Pj=Tub_j->side4;
					Nj=Tub_j->NS4;
					Cij=r[1];
					MyAD=AD01;
				}
				else {
					Pi=Tub_i->side2;
					Ni=Tub_i->NS2;
					Pj=Tub_j->side3;
					Nj=Tub_j->NS3;
					Cij=r[1];
					MyAD=AD23;
				}
				Dcni=Dot(Cij,Ni);
				Dcnj=Dot(Cij,Nj);
				Dtitj=Dot(Tub_i->t,Tub_j->t);
				Dnitj=Dot(Tub_j->t,Ni);
				Dnjti=Dot(Tub_i->t,Nj);
				DTheta[0]=sqrt(1-Dcni*Dcni);			//Theta i bend
				DTheta[1]=-Dot(Cij,Pi)/DTheta[0];		//Theta i in plane
				DTheta[2]=Dtitj/sqrt(1-Dnitj*Dnitj);	//Theta titj in plane
				DTheta[3]=sqrt(1-Dcnj*Dcnj);			//Theta j bend
				DTheta[4]=Dot(Cij,Pj)/DTheta[3];		//Theta j in plane
				DTheta[5]=Dtitj/sqrt(1-Dnjti*Dnjti);	//Theta titj in plane
				for(k=0;k<6;k++){
					Delta[k]=(DTheta[k]+1)*NLen;
					DNo[k]=(int)(Delta[k]);
					Delta[k]-=DNo[k];
				}
				for(k=0;k<6;k++)Theta[k]=(1-Delta[k])*ansACos[DNo[k]]+(Delta[k])*ansACos[DNo[k]+1];
				if(Dot(Tub_j->t,Pi)<0)Theta[2]=-Theta[2];
				if(Dot(Tub_i->t,Pj)<0)Theta[5]=-Theta[5];
				single[0]=Theta[0]*Theta[0];
				single[1]=Theta[3]*Theta[3];
				single[2]=Theta[1]-0.5*Theta[2];
				single[2]*=single[2];
				single[3]=Theta[4]-0.5*Theta[5];
				single[3]*=single[3];
				Dtmp=RDcy[1]*MyAD;
				Poten[4]+=Dtmp*KLat*(single[0]+single[1]);
				Poten[5]+=Dtmp*GLat*(single[2]+single[3]);
			}
			if(Pr[1]<0)Poten[0]+=Pr[1]*ELong*AD45;
			else Poten[0]+=Pr[1]*ELong;
			if(PrS[1]<0)Poten[3]+=PrS[1]*ELatBA*(AD01+AD23);
			else Poten[3]+=PrS[1]*ELatBA;
			
			if(Pr[1]*AD45<-0.4){
				Tub_i->front=Tub_j->Tub_No;
				Tub_j->back=Tub_i->Tub_No;
			}
			if(AD23>0 && PrS[0]*AD23<-0.35){
				Tub_i->left=Tub_j->Tub_No;
				Tub_j->right=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->left==Tub_j->Tub_No)Tub_i->left=-1;
				if(Tub_j->right==Tub_i->Tub_No)Tub_j->right=-1;
			}
			#if HELIX_ON
			if(AD01>0 && PrS[0]*AD01<-0.35){
				Tub_i->rseam=Tub_j->Tub_No;
				Tub_j->lseam=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->rseam==Tub_j->Tub_No)Tub_i->rseam=-1;
				if(Tub_j->lseam==Tub_i->Tub_No)Tub_j->lseam=-1;
			}
			#endif
		}
		if(rmod[2]<=cutoff){
			angle[0]=Dot(Tub_j->side1,r[2]);
			angle[1]=-Dot(Tub_i->side4,r[2]);
			angle[2]=Dot(Tub_j->side2,r[2]);
			angle[3]=-Dot(Tub_i->side3,r[2]);
			angle[4]=-Dot(Tub_i->chain2,r[2]);
			angle[5]=Dot(Tub_j->chain1,r[2]);
			for(k=0;k<6;k++){
				Delta[k]=(angle[k]+1)*NLen;
				DNo[k]=(int)(Delta[k]);
				Delta[k]-=DNo[k];
			}
			for(k=0;k<6;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
			AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[5]+Tub_j->QDRate[2])/2;
			AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[4]+Tub_j->QDRate[3])/2;
			AD45=ADcy[4]*ADcy[5]*(Tub_i->QDRate[1]+Tub_j->QDRate[0])/2;
			if(RDcy[2]>0 && AD45>0){
				double t0[3],Dn1n2,Dn1t,Dn2t,Dtw,Atw,DeltaAtw;
				int NoAtw;
				for(k=0;k<3;k++)t0[k]=r[2][k];
				Dn1n2=Dot(Tub_i->n,Tub_j->n);
				Dn1t=Dot(t0,Tub_i->n);
				Dn2t=Dot(t0,Tub_j->n);
				Dtw=(Dn1n2-Dn1t*Dn2t)/sqrt((1-Dn1t*Dn1t)*(1-Dn2t*Dn2t));
				DeltaAtw=(Dtw+1)*NLen;
				NoAtw=(int)(DeltaAtw);
				DeltaAtw-=NoAtw;
				Atw=(1-DeltaAtw)*ansPAC[NoAtw]+(DeltaAtw)*ansPAC[NoAtw+1];
				single[4]=(1-Delta[4])*ansPAC[DNo[4]]+(Delta[4])*ansPAC[DNo[4]+1];
				single[5]=(1-Delta[5])*ansPAC[DNo[5]]+(Delta[5])*ansPAC[DNo[5]+1];
				Dtmp=RDcy[2]*AD45;
				Poten[7]+=Dtmp*KLong*(single[4]+single[5]);
				Poten[8]+=Dtmp*TLong*Atw/(rmod[2]*rmod[2]);
			}
			if(RDcy[2]>0 && (AD01>0 || AD23>0)){
				double Dcni,Dcnj,Dtitj,Dnitj,Dnjti,MyAD;
				double *Pi,*Pj,*Ni,*Nj,*Cij,DTheta[6],Theta[6];
				if(AD01>0){
					Pi=Tub_i->side4;
					Ni=Tub_i->NS4;
					Pj=Tub_j->side1;
					Nj=Tub_j->NS1;
					Cij=r[2];
					MyAD=AD01;
				}
				else {
					Pi=Tub_i->side3;
					Ni=Tub_i->NS3;
					Pj=Tub_j->side2;
					Nj=Tub_j->NS2;
					Cij=r[2];
					MyAD=AD23;
				}
				Dcni=Dot(Cij,Ni);
				Dcnj=Dot(Cij,Nj);
				Dtitj=Dot(Tub_i->t,Tub_j->t);
				Dnitj=Dot(Tub_j->t,Ni);
				Dnjti=Dot(Tub_i->t,Nj);
				DTheta[0]=sqrt(1-Dcni*Dcni);			//Theta i bend
				DTheta[1]=-Dot(Cij,Pi)/DTheta[0];		//Theta i in plane
				DTheta[2]=Dtitj/sqrt(1-Dnitj*Dnitj);	//Theta titj in plane
				DTheta[3]=sqrt(1-Dcnj*Dcnj);			//Theta j bend
				DTheta[4]=Dot(Cij,Pj)/DTheta[3];		//Theta j in plane
				DTheta[5]=Dtitj/sqrt(1-Dnjti*Dnjti);	//Theta titj in plane
				for(k=0;k<6;k++){
					Delta[k]=(DTheta[k]+1)*NLen;
					DNo[k]=(int)(Delta[k]);
					Delta[k]-=DNo[k];
				}
				for(k=0;k<6;k++)Theta[k]=(1-Delta[k])*ansACos[DNo[k]]+(Delta[k])*ansACos[DNo[k]+1];
				if(Dot(Tub_j->t,Pi)<0)Theta[2]=-Theta[2];
				if(Dot(Tub_i->t,Pj)<0)Theta[5]=-Theta[5];
				single[0]=Theta[0]*Theta[0];
				single[1]=Theta[3]*Theta[3];
				single[2]=Theta[1]-0.5*Theta[2];
				single[2]*=single[2];
				single[3]=Theta[4]-0.5*Theta[5];
				single[3]*=single[3];
				Dtmp=RDcy[2]*MyAD;
				Poten[10]+=Dtmp*KLat*(single[0]+single[1]);
				Poten[11]+=Dtmp*GLat*(single[2]+single[3]);
			}
			if(PrS[2]<0)Poten[9]+=PrS[2]*ELatBA*(AD01+AD23);
			else Poten[9]+=PrS[2]*ELatBA;
			if(Pr[2]<0)Poten[6]+=Pr[2]*ELong*AD45;
			else Poten[6]+=Pr[2]*ELong;
			
			if(Pr[2]*AD45<-0.4){
				Tub_i->back=Tub_j->Tub_No;
				Tub_j->front=Tub_i->Tub_No;
			}
			if(AD23>0 && PrS[2]*AD23<-0.35){
				Tub_i->right=Tub_j->Tub_No;
				Tub_j->left=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->right==Tub_j->Tub_No)Tub_i->right=-1;
				if(Tub_j->left==Tub_i->Tub_No)Tub_j->left=-1;
			}
			#if HELIX_ON
			if(AD01>0 && PrS[2]*AD01<-0.35){
				Tub_i->lseam=Tub_j->Tub_No;
				Tub_j->rseam=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->lseam==Tub_j->Tub_No)Tub_i->lseam=-1;
				if(Tub_j->rseam==Tub_i->Tub_No)Tub_j->rseam=-1;
			}
			#endif
		}
		if(rmod[3]<=cutoff){
			angle[0]=-Dot(Tub_i->side3,r[3]);
			angle[1]=Dot(Tub_j->side4,r[3]);
			angle[2]=Dot(Tub_j->side3,r[3]);
			angle[3]=-Dot(Tub_i->side4,r[3]);
			for(k=0;k<4;k++){
				Delta[k]=(angle[k]+1)*NLen;
				DNo[k]=(int)(Delta[k]);
				Delta[k]-=DNo[k];
			}
			for(k=0;k<4;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
			
			AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[4]+Tub_j->QDRate[5])/2;
			AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[5]+Tub_j->QDRate[4])/2;
			if(RDcy[3]>0 && (AD01>0 || AD23>0)){
				double Dcni,Dcnj,Dtitj,Dnitj,Dnjti,MyAD;
				double *Pi,*Pj,*Ni,*Nj,*Cij,DTheta[6],Theta[6];
				if(AD01>0){
					Pi=Tub_i->side3;
					Ni=Tub_i->NS3;
					Pj=Tub_j->side4;
					Nj=Tub_j->NS4;
					Cij=r[3];
					MyAD=AD01;
				}
				else {
					Pi=Tub_i->side4;
					Ni=Tub_i->NS4;
					Pj=Tub_j->side3;
					Nj=Tub_j->NS3;
					Cij=r[3];
					MyAD=AD23;
				}
				Dcni=Dot(Cij,Ni);
				Dcnj=Dot(Cij,Nj);
				Dtitj=Dot(Tub_i->t,Tub_j->t);
				Dnitj=Dot(Tub_j->t,Ni);
				Dnjti=Dot(Tub_i->t,Nj);
				DTheta[0]=sqrt(1-Dcni*Dcni);			//Theta i bend
				DTheta[1]=-Dot(Cij,Pi)/DTheta[0];		//Theta i in plane
				DTheta[2]=Dtitj/sqrt(1-Dnitj*Dnitj);	//Theta titj in plane
				DTheta[3]=sqrt(1-Dcnj*Dcnj);			//Theta j bend
				DTheta[4]=Dot(Cij,Pj)/DTheta[3];		//Theta j in plane
				DTheta[5]=Dtitj/sqrt(1-Dnjti*Dnjti);	//Theta titj in plane
				for(k=0;k<6;k++){
					Delta[k]=(DTheta[k]+1)*NLen;
					DNo[k]=(int)(Delta[k]);
					Delta[k]-=DNo[k];
				}
				for(k=0;k<6;k++)Theta[k]=(1-Delta[k])*ansACos[DNo[k]]+(Delta[k])*ansACos[DNo[k]+1];
				if(Dot(Tub_j->t,Pi)<0)Theta[2]=-Theta[2];
				if(Dot(Tub_i->t,Pj)<0)Theta[5]=-Theta[5];
				single[0]=Theta[0]*Theta[0];
				single[1]=Theta[3]*Theta[3];
				single[2]=Theta[1]-0.5*Theta[2];
				single[2]*=single[2];
				single[3]=Theta[4]-0.5*Theta[5];
				single[3]*=single[3];
				Dtmp=RDcy[3]*MyAD;
				Poten[10]+=Dtmp*KLat*(single[0]+single[1]);
				Poten[11]+=Dtmp*GLat*(single[2]+single[3]);
			}
			if(PrS[3]<0){
				Poten[9]+=ELatAA*PrS[3]*(AD01+AD23);
			}
			else Poten[9]+=ELatAA*PrS[3];
			
			if(AD01>0 && PrS[3]*AD01<-0.35){
				Tub_i->right=Tub_j->Tub_No;
				Tub_j->left=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->right==Tub_j->Tub_No)Tub_i->right=-1;
				if(Tub_j->left==Tub_i->Tub_No)Tub_j->left=-1;
			}
			if(AD23>0 && PrS[3]*AD23<-0.35){
				Tub_i->left=Tub_j->Tub_No;
				Tub_j->right=Tub_i->Tub_No;
			}
			else{
				if(Tub_i->left==Tub_j->Tub_No)Tub_i->left=-1;
				if(Tub_j->right==Tub_i->Tub_No)Tub_j->right=-1;
			}
		}
		for(k=0;k<12;k++)PotenAll[k]+=Poten[k];
		return 1;
	}
}
#if FixLink
void Set_LinkList(tubulin *Tub_i,tubulin *Tub_j){
	int i,j,k,l;
	int DNo[6],Link_i,Link_j;
	double r[4][3],rmod[4],RDcy[4];
	double angle[6]={0.0};
	double single[6],ADcy[6],Delta[6];
	double AD01,AD23,AD45;
	
	Link_i=Tub_i->AddNum;
	Link_j=Tub_j->AddNum;
	for(l=0;l<3;l++){
		r[0][l]=Tub_i->r1[l]-Tub_j->r1[l];
		r[1][l]=Tub_i->r1[l]-Tub_j->r2[l];
		r[2][l]=Tub_i->r2[l]-Tub_j->r1[l];
		r[3][l]=Tub_i->r2[l]-Tub_j->r2[l];
	}
	for(k=0;k<4;k++){
		PeriodBoundary(r[k]);
		rmod[k]=Dot(r[k],r[k]);
		if(rmod[k]<=cutoffsq){
			rmod[k]=sqrt(rmod[k]);
			Delta[k]=1/rmod[k];
			for(l=0;l<3;l++)r[k][l]*=Delta[k];
			
			Delta[k]=rmod[k]*NLen;
			DNo[k]=(int)(Delta[k]);
			Delta[k]-=DNo[k];
			RDcy[k]=(1-Delta[k])*ansRDcyC[DNo[k]]+(Delta[k])*ansRDcyC[DNo[k]+1];
		}
		else{
			RDcy[k]=0;
		}
	}
	if(rmod[0]<=cutoff_s){ 
		angle[0]=-Dot(Tub_i->side1,r[0]);
		angle[1]=Dot(Tub_j->side2,r[0]);
		angle[2]=Dot(Tub_j->side1,r[0]);
		angle[3]=-Dot(Tub_i->side2,r[0]);
		for(k=0;k<4;k++){
			Delta[k]=(angle[k]+1)*NLen;
			DNo[k]=(int)(Delta[k]);
			Delta[k]-=DNo[k];
		}
		for(k=0;k<4;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
		AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[2]+Tub_j->QDRate[3])/2;
		AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[3]+Tub_j->QDRate[2])/2;
		if((AD01>0 || AD23>0)){
			if(AD01>0){
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x34;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x43;
				Link_i++;
				Link_j++;
			}
			else {
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x43;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x34;
				Link_i++;
				Link_j++;
			}
		}
	}
	if(rmod[1]<=cutoff_s){
		angle[0]=-Dot(Tub_i->side1,r[1]);
		angle[1]=Dot(Tub_j->side4,r[1]);
		angle[2]=-Dot(Tub_i->side2,r[1]);
		angle[3]=Dot(Tub_j->side3,r[1]);
		angle[4]=-Dot(Tub_i->chain1,r[1]);
		angle[5]=Dot(Tub_j->chain2,r[1]);
		for(k=0;k<6;k++){
			Delta[k]=(angle[k]+1)*NLen;
			DNo[k]=(int)(Delta[k]);
			Delta[k]-=DNo[k];
		}
		for(k=0;k<6;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
		
		AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[2]+Tub_j->QDRate[5])/2;
		AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[3]+Tub_j->QDRate[4])/2;
		AD45=ADcy[4]*ADcy[5]*(Tub_i->QDRate[0]+Tub_j->QDRate[1])/2;
		if(AD45>0){
			Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
			Tub_i->PPMark[Link_i]=0x12|0x08;
			Tub_i->PosList[Link_i]=Link_j;
			Tub_j->PosList[Link_j]=Link_i;
			Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
			Tub_j->PPMark[Link_j]=0x21|0x08;
			Link_i++;
			Link_j++;
		}
		if((AD01>0 || AD23>0)){
			if(AD01>0){
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x36|0x08;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x63|0x08;
				Link_i++;
				Link_j++;
			}
			else {
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x45|0x08;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x54|0x08;
				Link_i++;
				Link_j++;
			}
		}
	}
	if(rmod[2]<=cutoff_s){
		angle[0]=Dot(Tub_j->side1,r[2]);
		angle[1]=-Dot(Tub_i->side4,r[2]);
		angle[2]=Dot(Tub_j->side2,r[2]);
		angle[3]=-Dot(Tub_i->side3,r[2]);
		angle[4]=-Dot(Tub_i->chain2,r[2]);
		angle[5]=Dot(Tub_j->chain1,r[2]);
		for(k=0;k<6;k++){
			Delta[k]=(angle[k]+1)*NLen;
			DNo[k]=(int)(Delta[k]);
			Delta[k]-=DNo[k];
		}
		for(k=0;k<6;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
		AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[5]+Tub_j->QDRate[2])/2;
		AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[4]+Tub_j->QDRate[3])/2;
		AD45=ADcy[4]*ADcy[5]*(Tub_i->QDRate[1]+Tub_j->QDRate[0])/2;
		if(AD45>0){
			Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
			Tub_i->PPMark[Link_i]=0x21|0x08;
			Tub_i->PosList[Link_i]=Link_j;
			Tub_j->PosList[Link_j]=Link_i;
			Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
			Tub_j->PPMark[Link_j]=0x12|0x08;
			Link_i++;
			Link_j++;
		}
		if((AD01>0 || AD23>0)){
			if(AD01>0){
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x63|0x08;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x36|0x08;
				Link_i++;
				Link_j++;
			}
			else {
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x54|0x08;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x45|0x08;
				Link_i++;
				Link_j++;
			}
		}
	}
	if(rmod[3]<=cutoff_s){
		angle[0]=-Dot(Tub_i->side3,r[3]);
		angle[1]=Dot(Tub_j->side4,r[3]);
		angle[2]=Dot(Tub_j->side3,r[3]);
		angle[3]=-Dot(Tub_i->side4,r[3]);
		for(k=0;k<4;k++){
			Delta[k]=(angle[k]+1)*NLen;
			DNo[k]=(int)(Delta[k]);
			Delta[k]-=DNo[k];
		}
		for(k=0;k<4;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
		
		AD01=ADcy[0]*ADcy[1]*(Tub_i->QDRate[4]+Tub_j->QDRate[5])/2;
		AD23=ADcy[2]*ADcy[3]*(Tub_i->QDRate[5]+Tub_j->QDRate[4])/2;
		if((AD01>0 || AD23>0)){
			if(AD01>0){
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x56;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x65;
				Link_i++;
				Link_j++;
			}
			else {
				Tub_i->LinkList[Link_i]=Tub_j->Tub_No;
				Tub_i->PPMark[Link_i]=0x65;
				Tub_i->PosList[Link_i]=Link_j;
				Tub_j->PosList[Link_j]=Link_i;
				Tub_j->LinkList[Link_j]=Tub_i->Tub_No;
				Tub_j->PPMark[Link_j]=0x56;
				Link_i++;
				Link_j++;
			}
		}
	}
	if(Link_i>Max_Add_Link){
		fprintf(stderr,"ERROR: Too Many Links (%d) of Tubulin %d!\n",Link_i,Tub_i->Tub_No);
		exit(-5);
	}
	if(Link_j>Max_Add_Link){
		fprintf(stderr,"ERROR: Too Many Links (%d) of Tubulin %d!\n",Link_j,Tub_j->Tub_No);
		exit(-5);
	}
	Tub_i->AddNum=Link_i;
	Tub_j->AddNum=Link_j;
}

double poten_fl(tubulin *Tub_i,double *PotenAll){	//粒子间能量计算 
	char MyMark;
	int i,j,k,l;
	int No=0,DNo[6];
	double rij[3],rmod,Pr,PrS,RDcy;
	double angle[6]={0.0};
	double single[6],ADcy[6],Delta[6],Dtmp;
	double AD01,AD23,AD45,MyAD;
	tubulin *Tub_j;
	
	if(!Tub_i->Active)return 0;
	for(j=0;j<Tub_i->AddNum;j++){
		MyMark=Tub_i->PPMark[j];
		PotenAll[j]=0;
		if(Tub_i->PPMark[j] & 0x80){
			PotenAll[j]=Tub_i->PotenList[j];
		}
		else{
			if(Tub_i->LinkList[j]<0)continue;
			Tub_j=&Tub[Tub_i->LinkList[j]];
			if(((MyMark>>4) & 0x07)<3){
				double *Pi,*Pj,MyQD;
				if((MyMark & 0x70)==0x10){
					for(l=0;l<3;l++)rij[l]=Tub_i->r1[l]-Tub_j->r2[l];
					Pi=Tub_i->chain1;
					Pj=Tub_j->chain2;
					MyQD=(Tub_i->QDRate[0]+Tub_j->QDRate[1])/2;
				}
				else if((MyMark & 0x70)==0x20){
					for(l=0;l<3;l++)rij[l]=Tub_i->r2[l]-Tub_j->r1[l];
					Pi=Tub_i->chain2;
					Pj=Tub_j->chain1;
					MyQD=(Tub_i->QDRate[1]+Tub_j->QDRate[0])/2;
				}
				PeriodBoundary(rij);
				rmod=Dot(rij,rij);
				if(rmod<=cutoffsq){
					rmod=sqrt(rmod);
					Dtmp=1/rmod;
					for(l=0;l<3;l++)rij[l]*=Dtmp;
					
					Dtmp=rmod*NLen;
					No=(int)(Dtmp);
					Dtmp-=No;
					Pr=(1-Dtmp)*ansPrC[No]+(Dtmp)*ansPrC[No+1];
					RDcy=(1-Dtmp)*ansRDcyC[No]+(Dtmp)*ansRDcyC[No+1];
					
					angle[0]=-Dot(Pi,rij);
					angle[1]=Dot(Pj,rij);
					for(k=0;k<2;k++){
						Delta[k]=(angle[k]+1)*NLen;
						DNo[k]=(int)(Delta[k]);
						Delta[k]-=DNo[k];
					}
					for(k=0;k<2;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
					MyAD=ADcy[0]*ADcy[1]*MyQD;
				}
				else{
					Pr=0;
					RDcy=0;
					MyAD=0;
				}
				if(RDcy>0){
					double Dn1n2,Dn1t,Dn2t,Dtw,Atw,DeltaAtw;
					int NoAtw;
					Dn1n2=Dot(Tub_i->n,Tub_j->n);
					Dn1t=Dot(rij,Tub_i->n);
					Dn2t=Dot(rij,Tub_j->n);
					Dtw=(Dn1n2-Dn1t*Dn2t)/sqrt((1-Dn1t*Dn1t)*(1-Dn2t*Dn2t));
					DeltaAtw=(Dtw+1)*NLen;
					NoAtw=(int)(DeltaAtw);
					DeltaAtw-=NoAtw;
					Atw=(1-DeltaAtw)*ansPAC[NoAtw]+(DeltaAtw)*ansPAC[NoAtw+1];
					single[0]=(1-Delta[0])*ansPAC[DNo[0]]+(Delta[0])*ansPAC[DNo[0]+1];
					single[1]=(1-Delta[1])*ansPAC[DNo[1]]+(Delta[1])*ansPAC[DNo[1]+1];
					Dtmp=RDcy*MyAD;
					PotenAll[j]+=Dtmp*(KLong*(single[0]+single[1])+TLong*Atw/(rmod*rmod));
				}
				if(Pr<0)PotenAll[j]+=Pr*ELong*MyAD;
				else PotenAll[j]+=Pr*ELong;
				if(Pr*MyAD<-0.4){
					if((MyMark & 0x70)==0x10){
						Tub_i->front=Tub_j->Tub_No;
						Tub_j->back=Tub_i->Tub_No;
					}
					else if((MyMark & 0x70)==0x20){
						Tub_i->back=Tub_j->Tub_No;
						Tub_j->front=Tub_i->Tub_No;
					}
				}
			}
			else{
				double *Pi,*Pj,*Ni,*Nj,*ri,*rj,QDi,QDj;
				switch(MyMark & 0x70){
					case 0x30:{
						Pi=Tub_i->side1;
						Ni=Tub_i->NS1;
						ri=Tub_i->r1;
						QDi=Tub_i->QDRate[2];
						break; 
					}
					case 0x40:{
						Pi=Tub_i->side2;
						Ni=Tub_i->NS2;
						ri=Tub_i->r1;
						QDi=Tub_i->QDRate[3];
						break; 
					}
					case 0x50:{
						Pi=Tub_i->side3;
						Ni=Tub_i->NS3;
						ri=Tub_i->r2;
						QDi=Tub_i->QDRate[4];
						break; 
					}
					case 0x60:{
						Pi=Tub_i->side4;
						Ni=Tub_i->NS4;
						ri=Tub_i->r2;
						QDi=Tub_i->QDRate[5];
						break; 
					}
				}
				switch(MyMark & 0x07){
					case 0x03:{
						Pj=Tub_j->side1;
						Nj=Tub_j->NS1;
						rj=Tub_j->r1;
						QDj=Tub_j->QDRate[2];
						break;
					}
					case 0x04:{
						Pj=Tub_j->side2;
						Nj=Tub_j->NS2;
						rj=Tub_j->r1;
						QDj=Tub_j->QDRate[3];
						break;
					}
					case 0x05:{
						Pj=Tub_j->side3;
						Nj=Tub_j->NS3;
						rj=Tub_j->r2;
						QDj=Tub_j->QDRate[4];
						break;
					}
					case 0x06:{
						Pj=Tub_j->side4;
						Nj=Tub_j->NS4;
						rj=Tub_j->r2;
						QDj=Tub_j->QDRate[5];
						break;
					}
				}
				
				for(l=0;l<3;l++)rij[l]=ri[l]-rj[l];
				PeriodBoundary(rij);
				rmod=Dot(rij,rij);
				if(rmod<=cutoffsq){
					rmod=sqrt(rmod);
					Dtmp=1/rmod;
					for(l=0;l<3;l++)rij[l]*=Dtmp;
					
					Dtmp=rmod*NLen;
					No=(int)(Dtmp);
					Dtmp-=No;
					PrS=(1-Dtmp)*ansPrS[No]+(Dtmp)*ansPrS[No+1];
					RDcy=(1-Dtmp)*ansRDcyC[No]+(Dtmp)*ansRDcyC[No+1];
					
					angle[0]=-Dot(Pi,rij);
					angle[1]=Dot(Pj,rij);
					for(k=0;k<2;k++){
						Delta[k]=(angle[k]+1)*NLen;
						DNo[k]=(int)(Delta[k]);
						Delta[k]-=DNo[k];
					}
					for(k=0;k<2;k++)ADcy[k]=(1-Delta[k])*ansADcyC[DNo[k]]+(Delta[k])*ansADcyC[DNo[k]+1];
					MyAD=ADcy[0]*ADcy[1]*(QDi+QDj)/2;
				}
				else{
					PrS=0;
					RDcy=0;
					MyAD=0;
				}
				
				if(RDcy>0){
					double Dcni,Dcnj,Dtitj,Dnitj,Dnjti;
					double DTheta[6],Theta[6];
					Dcni=Dot(rij,Ni);
					Dcnj=Dot(rij,Nj);
					Dtitj=Dot(Tub_i->t,Tub_j->t);
					Dnitj=Dot(Tub_j->t,Ni);
					Dnjti=Dot(Tub_i->t,Nj);
					DTheta[0]=sqrt(1-Dcni*Dcni);			//Theta i bend
					DTheta[1]=-Dot(rij,Pi)/DTheta[0];		//Theta i in plane
					DTheta[2]=Dtitj/sqrt(1-Dnitj*Dnitj);	//Theta titj in plane
					DTheta[3]=sqrt(1-Dcnj*Dcnj);			//Theta j bend
					DTheta[4]=Dot(rij,Pj)/DTheta[3];		//Theta j in plane
					DTheta[5]=Dtitj/sqrt(1-Dnjti*Dnjti);	//Theta titj in plane
					for(k=0;k<6;k++){
						Delta[k]=(DTheta[k]+1)*NLen;
						DNo[k]=(int)(Delta[k]);
						Delta[k]-=DNo[k];
					}
					for(k=0;k<6;k++)Theta[k]=(1-Delta[k])*ansACos[DNo[k]]+(Delta[k])*ansACos[DNo[k]+1];
					if(Dot(Tub_j->t,Pi)<0)Theta[2]=-Theta[2];
					if(Dot(Tub_i->t,Pj)<0)Theta[5]=-Theta[5];
					single[0]=Theta[0]*Theta[0];
					single[1]=Theta[3]*Theta[3];
					single[2]=Theta[1]-0.5*Theta[2];
					single[2]*=single[2];
					single[3]=Theta[4]-0.5*Theta[5];
					single[3]*=single[3];
					Dtmp=RDcy*MyAD;
					PotenAll[j]+=Dtmp*(KLat*(single[0]+single[1])+GLat*(single[2]+single[3]));
				}
				if(MyMark & 0x08){
					if(PrS<0)PotenAll[j]+=PrS*ELatBA*MyAD;
					else PotenAll[j]+=PrS*ELatBA;
				}
				else{
					if(PrS<0)PotenAll[j]+=PrS*ELatAA*MyAD;
					else PotenAll[j]+=PrS*ELatAA;
				}
			}
			Tub_i->PotenList[j]=PotenAll[j];
		}
	}
	return 1;
}
#endif
#if MD_ON
int Force(tubulin *Tubp_i,tubulin *Tubp_j){
	int k,l;
	long No=0;
	double r[5][3],rmod[5],lj[5],fr[4],acna,DeltaNo;
	double angle[20]={0.0},tangle=0.0,nangle=0.0,force_1[3],force_2[3],force_all[3],force_moment0[3],force_moment1[3],force_moment2[3],force_moment3[3],moment[3];
	double single[20],tsingle,nsingle,ftangle,fnangle,fangle[20],fpair[10][3],gpair[10],force_moment_t[3],force_moment_n[3];
	
	if(BUG)fprintf(stderr,"\nForce\ti: %d j: %d\n",Tubp_i->Tub_No,Tubp_j->Tub_No);
	if(Tubp_i->Active==0||Tubp_j->Active==0)return 0;
    if(Tubp_i->Tub_No==Tubp_j->Tub_No){
		return 0;
    }
    for(l=0;l<3;l++){
		force_1[l]=force_2[l]=force_all[l]=force_moment_n[l]=force_moment_t[l]=force_moment0[l]=force_moment1[l]=force_moment2[l]=force_moment3[l]=moment[l]=0;
	}

		for(l=0;l<3;l++){
			r[0][l]=Tubp_i->r1[l]-Tubp_j->r1[l];
			r[1][l]=Tubp_i->r1[l]-Tubp_j->r2[l];
			r[2][l]=Tubp_i->r2[l]-Tubp_j->r1[l];
			r[3][l]=Tubp_i->r2[l]-Tubp_j->r2[l];
			r[4][l]=Tubp_i->r[l]-Tubp_j->r[l];
		}
		for(k=0;k<5;k++){
			PeriodBoundary(r[k]);
			rmod[k]=Dot(r[k],r[k]);
			if(rmod[k]<cutoffsq){
				rmod[k]=sqrt(rmod[k]);
			}
		}

		if(rmod[0]>cutoff&&rmod[1]>cutoff&&rmod[2]>cutoff&&rmod[3]>cutoff){
			return 0;
    	}
    	lj[4]=1.0/(rmod[4]*rmod[4]);

		if(BUG)fprintf(stderr,"Rmod\n");
		for(k=0;k<4;k++){
			if(rmod[k]>=cutoff){
				lj[k]=0;
				fr[k]=0;//1e-9*(2.0*rand()/(RAND_MAX+0.1)-1);
			}
			else {
				No=(int)(rmod[k]/deltalj);
				DeltaNo=rmod[k]/deltalj;
				No=(int)(DeltaNo);
				DeltaNo=DeltaNo-No;
				lj[k]=(1-DeltaNo)*ansPrC[No]+DeltaNo*ansPrC[No+1];
				fr[k]=(1-DeltaNo)*ForceR[No]+DeltaNo*ForceR[No+1];
			}
			for(l=0;l<3;l++)r[k][l]=r[k][l]/rmod[k];
		}
		
		for(k=0;k<3;k++){
			angle[0]+=-Tubp_i->side1[k]*r[0][k];
			angle[1]+=Tubp_j->side2[k]*r[0][k];

			angle[2]+=Tubp_j->side1[k]*r[0][k];
			angle[3]+=-Tubp_i->side2[k]*r[0][k];

			angle[4]+=-Tubp_i->side3[k]*r[3][k];
			angle[5]+=Tubp_j->side4[k]*r[3][k];

			angle[6]+=Tubp_j->side3[k]*r[3][k];
			angle[7]+=-Tubp_i->side4[k]*r[3][k];

			angle[12]+=-Tubp_i->side1[k]*r[1][k];
			angle[13]+=Tubp_j->side4[k]*r[1][k];

			angle[14]+=Tubp_j->side1[k]*r[2][k];
			angle[15]+=-Tubp_i->side4[k]*r[2][k];

			angle[16]+=-Tubp_i->side2[k]*r[1][k];
			angle[17]+=Tubp_j->side3[k]*r[1][k];

			angle[18]+=Tubp_j->side2[k]*r[2][k];
			angle[19]+=-Tubp_i->side3[k]*r[2][k];
		}


		for(k=0;k<3;k++){

			angle[8]+=-Tubp_i->chain1[k]*r[1][k];
			angle[9]+=Tubp_j->chain2[k]*r[1][k];

			angle[10]+=-Tubp_i->chain2[k]*r[2][k];
			angle[11]+=Tubp_j->chain1[k]*r[2][k];

			tangle+=Tubp_i->t[k]*Tubp_j->t[k];
			nangle+=Tubp_i->n[k]*Tubp_j->n[k];
		}
		if(BUG)fprintf(stderr,"angle\n");
		for(k=0;k<8;k++){
			angle[k]=angle[k]+1;
			DeltaNo=angle[k]/deltaexp;
			No=(int)(DeltaNo);
			DeltaNo=DeltaNo-No;
			single[k]=(1-DeltaNo)*PAngle2[No]+DeltaNo*PAngle2[No+1];
			fangle[k]=(1-DeltaNo)*ForceAngle2[No]+DeltaNo*ForceAngle2[No+1];
			angle[k]=angle[k]-1;
			
		}
		for(;k<12;k++){
			angle[k]=angle[k]+1;
			DeltaNo=angle[k]/deltaexp;
			No=(int)(DeltaNo);
			DeltaNo=DeltaNo-No;
			single[k]=(1-DeltaNo)*PAngle1[No]+DeltaNo*PAngle1[No+1];
			fangle[k]=(1-DeltaNo)*ForceAngle1[No]+DeltaNo*ForceAngle1[No+1];
			angle[k]=angle[k]-1;
			
		}
		for(;k<20;k++){
			angle[k]=angle[k]+1;
			DeltaNo=angle[k]/deltaexp;
			No=(int)(DeltaNo);
			DeltaNo=DeltaNo-No;
			single[k]=(1-DeltaNo)*PAngle2[No]+DeltaNo*PAngle2[No+1];
			fangle[k]=(1-DeltaNo)*ForceAngle2[No]+DeltaNo*ForceAngle2[No+1];
			angle[k]=angle[k]-1;
			
		}
		
		tangle=tangle+1;
		nangle=nangle+1;
		DeltaNo=tangle/deltaexp;
		No=(int)(DeltaNo);
		DeltaNo=DeltaNo-No;
		tsingle=(1-DeltaNo)*PAngle3[No]+DeltaNo*PAngle3[No+1];
		ftangle=(1-DeltaNo)*ForceAngle3[No]+DeltaNo*ForceAngle3[No+1];
		DeltaNo=nangle/deltaexp;
		No=(int)(DeltaNo);
		DeltaNo=DeltaNo-No;
		nsingle=(1-DeltaNo)*PAngle4[No]+DeltaNo*PAngle4[No+1];
		fnangle=(1-DeltaNo)*ForceAngle4[No]+DeltaNo*ForceAngle4[No+1];

		tangle=tangle-1;
		nangle=nangle-1;
		
		if(fabs(nangle)<1)acna=acos(nangle);
		else acna=0;
	#if RLW_ON
		gpair[0]=Tubp_i->p_poten[1]*single[0]*single[1]*tsingle;
		gpair[1]=Tubp_i->p_poten[1]*single[2]*single[3]*tsingle;
		gpair[2]=Tubp_i->p_poten[1]*single[4]*single[5]*tsingle;
		gpair[3]=Tubp_i->p_poten[1]*single[6]*single[7]*tsingle;
		gpair[4]=Tubp_i->p_poten[0]*single[8]*single[9]*nsingle;
		gpair[5]=Tubp_i->p_poten[0]*single[10]*single[11]*nsingle;
		gpair[6]=Tubp_i->p_poten[2]*single[12]*single[13]*tsingle;
		gpair[7]=Tubp_i->p_poten[2]*single[14]*single[15]*tsingle;
		gpair[8]=Tubp_i->p_poten[2]*single[16]*single[17]*tsingle;
		gpair[9]=Tubp_i->p_poten[2]*single[18]*single[19]*tsingle;
	#else
		gpair[0]=single[0]*single[1]*tsingle;
		gpair[1]=single[2]*single[3]*tsingle;
		gpair[2]=single[4]*single[5]*tsingle;
		gpair[3]=single[6]*single[7]*tsingle;
		gpair[4]=single[8]*single[9]*nsingle;
		gpair[5]=single[10]*single[11]*nsingle;
		gpair[6]=single[12]*single[13]*tsingle;
		gpair[7]=single[14]*single[15]*tsingle;
		gpair[8]=single[16]*single[17]*tsingle;
		gpair[9]=single[18]*single[19]*tsingle;
	#endif
		if(lj[0]>0)gpair[0]=gpair[1]=1;
		if(lj[1]>0)gpair[4]=gpair[6]=gpair[8]=1;
		if(lj[2]>0)gpair[5]=gpair[7]=gpair[9]=1;
		if(lj[3]>0)gpair[2]=gpair[3]=1;
		for(k=0;k<3;k++){
			if(0);//lj[0]>0)fpair[0][k]=-fr[0]*r[0][k];
			else fpair[0][k]=-fr[0]*gpair[0]*r[0][k]+
				lj[0]*gpair[0]*(fangle[0]*(Tubp_i->side1[k]+angle[0]*r[0][k])+fangle[1]*(-Tubp_j->side2[k]+angle[1]*r[0][k]))/rmod[0];
			if(0);//lj[0]>0)fpair[1][k]=-fr[0]*r[0][k];
			else fpair[1][k]=-fr[0]*gpair[1]*r[0][k]+
				lj[0]*gpair[1]*(fangle[3]*(Tubp_i->side2[k]+angle[3]*r[0][k])+fangle[2]*(-Tubp_j->side1[k]+angle[2]*r[0][k]))/rmod[0];

			if(0);//lj[3]>0)fpair[2][k]=-fr[3]*r[3][k];
			else fpair[2][k]=-fr[3]*gpair[2]*r[3][k]+
				lj[3]*gpair[2]*(fangle[4]*(Tubp_i->side3[k]+angle[4]*r[3][k])+fangle[5]*(-Tubp_j->side4[k]+angle[5]*r[3][k]))/rmod[3];

			if(0);//lj[3]>0)fpair[3][k]=-fr[3]*r[3][k];
			else fpair[3][k]=-fr[3]*gpair[3]*r[3][k]+
				lj[3]*gpair[3]*(fangle[7]*(Tubp_i->side4[k]+angle[7]*r[3][k])+fangle[6]*(-Tubp_j->side3[k]+angle[6]*r[3][k]))/rmod[3];

			if(0);//lj[1]>0)fpair[6][k]=-fr[1]*r[1][k];
			else fpair[6][k]=-fr[1]*gpair[6]*r[1][k]+
				lj[1]*gpair[6]*(fangle[12]*(Tubp_i->side1[k]+angle[12]*r[1][k])+fangle[13]*(-Tubp_j->side4[k]+angle[13]*r[1][k]))/rmod[1];

			if(0);//lj[2]>0)fpair[7][k]=-fr[2]*r[2][k];
			else fpair[7][k]=-fr[2]*gpair[7]*r[2][k]+
				lj[2]*gpair[7]*(fangle[15]*(Tubp_i->side4[k]+angle[15]*r[2][k])+fangle[14]*(-Tubp_j->side1[k]+angle[14]*r[2][k]))/rmod[2];

			if(0);//lj[1]>0)fpair[8][k]=-fr[1]*r[1][k];
			else fpair[8][k]=-fr[1]*gpair[8]*r[1][k]+
				lj[1]*gpair[8]*(fangle[16]*(Tubp_i->side2[k]+angle[16]*r[1][k])+fangle[17]*(-Tubp_j->side3[k]+angle[17]*r[1][k]))/rmod[1];

			if(0);//lj[2]>0)fpair[9][k]=-fr[2]*r[2][k];
			else fpair[9][k]=-fr[2]*gpair[9]*r[2][k]+
				lj[2]*gpair[9]*(fangle[19]*(Tubp_i->side3[k]+angle[19]*r[2][k])+fangle[18]*(-Tubp_j->side2[k]+angle[18]*r[2][k]))/rmod[2];

			if(0);//lj[1]>0)fpair[4][k]=-fr[1]*r[1][k];
			else fpair[4][k]=-fr[1]*gpair[4]*r[1][k]+
				lj[1]*gpair[4]*(fangle[8]*(Tubp_i->chain1[k]+angle[8]*r[1][k])+fangle[9]*(-Tubp_j->chain2[k]+angle[9]*r[1][k]))/rmod[1];
				
			if(0);//lj[2]>0)fpair[5][k]=-fr[2]*r[2][k];
			else fpair[5][k]=-fr[2]*gpair[5]*r[2][k]+
				lj[2]*gpair[5]*(fangle[10]*(Tubp_i->chain2[k]+angle[10]*r[2][k])+fangle[11]*(-Tubp_j->chain1[k]+angle[11]*r[2][k]))/rmod[2];
			
//			if(lj[0]<0){
				force_moment0[k]=W_Lat_aa*lj[0]*(gpair[1]*fangle[3]*Tubp_i->side2[k]+gpair[0]*fangle[0]*Tubp_i->side1[k]);
				force_moment_t[k]+=W_Lat_aa*lj[0]*(gpair[1]+gpair[0])*ftangle*Tubp_i->t[k];
//			}
//			if(lj[1]<0){
				force_moment1[k]=lj[1]*(W_Lat_ab*(gpair[6]*fangle[12]*Tubp_i->side1[k]+gpair[8]*fangle[16]*Tubp_i->side2[k])+W_Out*gpair[4]*fangle[8]*Tubp_i->chain1[k]);
//				force_moment1[k]+=lj[1]*(gpair[4]*fangle[8]*Tubp_i->chain1[k]);
				force_moment_t[k]+=W_Lat_ab*lj[1]*(gpair[6]+gpair[8])*ftangle*Tubp_i->t[k];
				force_moment_n[k]+=W_Out*lj[1]*gpair[4]*fnangle*Tubp_i->n[k];
//			}
//			if(lj[2]<0){
				force_moment2[k]=lj[2]*(W_Lat_ab*(gpair[7]*fangle[15]*Tubp_i->side4[k]+gpair[9]*fangle[19]*Tubp_i->side3[k])+W_Out*gpair[5]*fangle[10]*Tubp_i->chain2[k]);
//				force_moment2[k]+=lj[2]*(gpair[5]*fangle[10]*Tubp_i->chain2[k]);
				force_moment_t[k]+=W_Lat_ab*lj[2]*(gpair[7]+gpair[9])*ftangle*Tubp_i->t[k];
				force_moment_n[k]+=W_Out*lj[2]*gpair[5]*fnangle*Tubp_i->n[k];
//			}
//			if(lj[3]<0){
				force_moment3[k]=W_Lat_aa*lj[3]*(gpair[2]*fangle[4]*Tubp_i->side3[k]+gpair[3]*fangle[7]*Tubp_i->side4[k]);
				force_moment_t[k]+=W_Lat_aa*lj[3]*(gpair[2]+gpair[3])*ftangle*Tubp_i->t[k];
//			}
			if(gpair[4]>0.5||gpair[5]>0.5){
				if(fabs(nangle)<1)force_moment_n[k]+=-2*P_LT*acna/sin(acna)*lj[4]*Tubp_i->n[k];
				else force_moment_n[k]+=-2*P_LT*lj[4]*Tubp_i->n[k];
			}
		}
		for(k=0;k<3;k++){
			force_1[k]=W_Lat_aa*(fpair[0][k]+fpair[1][k])+W_Out*fpair[4][k]+W_Lat_ab*(fpair[6][k]+fpair[8][k]);
			force_2[k]=W_Lat_aa*(fpair[2][k]+fpair[3][k])+W_Out*fpair[5][k]+W_Lat_ab*(fpair[7][k]+fpair[9][k]);
			force_all[k]=force_1[k]+force_2[k];
			if(gpair[4]>0.5||gpair[5]>0.5)force_all[k]+=2*P_LT*acna*lj[4]*acna*lj[4]*r[4][k];
		}
		moment[0]=	
			- (r[0][1]*force_moment0[2]-r[0][2]*force_moment0[1]) - (r[1][1]*force_moment1[2]-r[1][2]*force_moment1[1])
			- (r[2][1]*force_moment2[2]-r[2][2]*force_moment2[1]) - (r[3][1]*force_moment3[2]-r[3][2]*force_moment3[1])
			+ (Tubp_j->t[1]*force_moment_t[2]-Tubp_j->t[2]*force_moment_t[1]) + (Tubp_j->n[1]*force_moment_n[2]-Tubp_j->n[2]*force_moment_n[1])
			+ (Tubp_i->rbar1[1]*force_1[2]-Tubp_i->rbar1[2]*force_1[1]) + (Tubp_i->rbar2[1]*force_2[2]-Tubp_i->rbar2[2]*force_2[1]);
		moment[1]=	
			- (r[0][2]*force_moment0[0]-r[0][0]*force_moment0[2]) - (r[1][2]*force_moment1[0]-r[1][0]*force_moment1[2])
			- (r[2][2]*force_moment2[0]-r[2][0]*force_moment2[2]) - (r[3][2]*force_moment3[0]-r[3][0]*force_moment3[2])
			+ (Tubp_j->t[2]*force_moment_t[0]-Tubp_j->t[0]*force_moment_t[2]) + (Tubp_j->n[2]*force_moment_n[0]-Tubp_j->n[0]*force_moment_n[2])
			+ (Tubp_i->rbar1[2]*force_1[0]-Tubp_i->rbar1[0]*force_1[2]) + (Tubp_i->rbar2[2]*force_2[0]-Tubp_i->rbar2[0]*force_2[2]);
		moment[2]=	
			- (r[0][0]*force_moment0[1]-r[0][1]*force_moment0[0]) - (r[1][0]*force_moment1[1]-r[1][1]*force_moment1[0])
			- (r[2][0]*force_moment2[1]-r[2][1]*force_moment2[0]) - (r[3][0]*force_moment3[1]-r[3][1]*force_moment3[0])
			+ (Tubp_j->t[0]*force_moment_t[1]-Tubp_j->t[1]*force_moment_t[0]) + (Tubp_j->n[0]*force_moment_n[1]-Tubp_j->n[1]*force_moment_n[0])
			+ (Tubp_i->rbar1[0]*force_1[1]-Tubp_i->rbar1[1]*force_1[0]) + (Tubp_i->rbar2[0]*force_2[1]-Tubp_i->rbar2[1]*force_2[0]);

		for(k=0;k<3;k++){
			Tubp_i->a[k]+=force_all[k];
			Tubp_i->mom[k]+=moment[k];
		}

		Tubp_i->PotenB+=(W_Lat_aa*lj[0]*(gpair[0]+gpair[1])+W_Lat_ab*lj[1]*(gpair[6]+gpair[8])+W_Out*lj[1]*gpair[4])/2;
		Tubp_i->PotenA+=(W_Lat_aa*lj[3]*(gpair[2]+gpair[3])+W_Lat_ab*lj[2]*(gpair[7]+gpair[9])+W_Out*lj[2]*gpair[5])/2;
		if(gpair[4]>0.5||gpair[5]>0.5){
			Tubp_i->potent+=P_LT*acna*acna*lj[4]/2;
		}
		Tubp_i->poten=Tubp_i->PotenB+Tubp_i->PotenA;
		if(Tubp_i->Tub_No<Tubp_j->Tub_No && Tubp_i->Tin%100==0){
			if(lj[1]*gpair[4]<=-0.4){
				Tubp_i->front=Tubp_j->Tub_No;
				Tubp_j->back=Tubp_i->Tub_No;
			}
			if(lj[2]*gpair[5]<=-0.4){
				Tubp_i->back=Tubp_j->Tub_No;
				Tubp_j->front=Tubp_i->Tub_No;
			}
			if((lj[0]*gpair[1]+lj[3]*gpair[3])<=-0.35||lj[1]*gpair[8]<=-0.35){
					Tubp_i->left=Tubp_j->Tub_No;
					Tubp_j->right=Tubp_i->Tub_No;
			}
			else{
				if(Tubp_i->left==Tubp_j->Tub_No)Tubp_i->left=-1;
				if(Tubp_j->right==Tubp_i->Tub_No)Tubp_j->right=-1;
			}
			if((lj[0]*gpair[0]+lj[3]*gpair[2])<=-0.35||lj[2]*gpair[9]<=-0.35){
					Tubp_i->right=Tubp_j->Tub_No;
					Tubp_j->left=Tubp_i->Tub_No;
			}
			else {
				if(Tubp_i->right==Tubp_j->Tub_No)Tubp_i->right=-1;
				if(Tubp_j->left==Tubp_i->Tub_No)Tubp_j->left=-1;
			}
	#if HELIX_ON
			if(lj[2]*gpair[7]<=-0.35){
					Tubp_i->lseam=Tubp_j->Tub_No;
					Tubp_j->rseam=Tubp_i->Tub_No;
			}
			else{
				if(Tubp_i->lseam==Tubp_j->Tub_No)Tubp_i->lseam=-1;
				if(Tubp_j->rseam==Tubp_i->Tub_No)Tubp_j->rseam=-1;
			}
			if(lj[1]*gpair[6]<=-0.35){
					Tubp_i->rseam=Tubp_j->Tub_No;
					Tubp_j->lseam=Tubp_i->Tub_No;
			}
			else {
				if(Tubp_i->rseam==Tubp_j->Tub_No)Tubp_i->rseam=-1;
				if(Tubp_j->lseam==Tubp_i->Tub_No)Tubp_j->lseam=-1;
			}
	#endif
		}
		if(Tubp_i->front==Tubp_j->Tub_No){
			if((Tubp_i->chain1[0]+r[1][0])*Tubp_i->s[0]+(Tubp_i->chain1[1]+r[1][1])*Tubp_i->s[1]+(Tubp_i->chain1[2]+r[1][2])*Tubp_i->s[2]<0){
				Tubp_i->potenf=0.5+lj[1]*gpair[4]/2;
			}
			else Tubp_i->potenf=(-lj[1]*gpair[4]/2-0.5);
		}
		else if(Tubp_i->back==Tubp_j->Tub_No){
			if((Tubp_i->chain2[0]+r[2][0])*Tubp_i->s[0]+(Tubp_i->chain2[1]+r[2][1])*Tubp_i->s[1]+(Tubp_i->chain2[2]+r[2][2])*Tubp_i->s[2]>0){
				Tubp_i->potenb=0.5+lj[2]*gpair[5]/2;
			}
			else Tubp_i->potenb=(-lj[2]*gpair[5]/2-0.5);
		}
	
	return 0;
}
#endif

void init(){	//初始化  
	int i,j,k,l,icell,FilaMisSt=0,FilaMisEd=N;
	double r[3],angle[5],ra,Dtmp;
	tubulin *Tubp;
	
	if(PhaseFile)ReadPhase();
	if(Whole_Struc==1 && TipSinA>0){
		N=NInit+(TipSinA+TipSinH)*FilaSet;
		if(TipType==1 || TipType==5){
			double MySin;
			int L1,L2;
			L1=(int)TipSinL;
			L2=L1+1;
			for(i=0;i<1001;i++){
				Dtmp=2*pi*i*0.001;
//				MySin=sin(L1*Dtmp);
				if(TipType==1)MySin=sqrt(L2-TipSinL)*sin(L1*Dtmp)+sqrt(TipSinL-L1)*sin(L2*Dtmp+TipSinPhase*pi);
//				else if(TipType==2)MySin=sin(Dtmp-pi/4);
//				else if(TipType==3)MySin=sin(Dtmp+TipSinPhase*pi);
				else if(TipType==5)MySin=(L2-TipSinL)*sin(L1*Dtmp)+(TipSinL-L1)*sin(L2*Dtmp+TipSinPhase*pi);
				if(MySin>MaxA)MaxA=MySin;
				if(MySin<MinA)MinA=MySin;
			}
			fprintf(stderr,"MaxA: %g\tMinA: %g\n",MaxA,MinA);
		}
		else if(TipType==6){
			double MySin;
			for(i=0;i<1001;i++){
				Dtmp=2*pi*i*0.001;
				MySin=0;
				for(k=0;k<6;k++){
					MySin+=sqrt(pow(k+1,-TipSinL))*sin((k+1)*Dtmp+MyTipPhase[k]*pi);
				}
				if(MySin>MaxA)MaxA=MySin;
				if(MySin<MinA)MinA=MySin;
			}
			fprintf(stderr,"MaxA: %g\tMinA: %g\n",MaxA,MinA);
		}
	}
	if(Whole_Struc==2){
		ReadTub();
		N=FilaSet*MaxLayer;
	}
	Set_FileName();
	
//网格划分参数设定 	
	Nactive=Ntype=0;
	if(FilaSet!=13 && FilaMis>0){
		FilaMisSt=(N-FilaMis*FilaSet)/13/2*13;
		FilaMisEd=FilaMisSt+FilaMis*FilaSet;
		fprintf(stderr,"FilaMis %d From %d to %d.\n",FilaSet,FilaMisSt,FilaMisEd);
	}
	else if(FilaSet!=13 && FilaMis<0){
		FilaMisSt=(N+FilaMis*FilaSet)/13*13;
		fprintf(stderr,"FilaMis %d From %d to %d.\n",FilaSet,FilaMisSt,FilaMisEd);
	}
	if(N%numprocs)MpiSec=N/numprocs+1;
	else MpiSec=N/numprocs;
	if(MpiSec<numprocs){
		fprintf(stderr,"Bad Processes Number!");
		MPI_Abort(MPI_COMM_WORLD,1);
		exit(1);
	}
	TubSize=sizeof(tubulin);

	GNo=CNo=1;
	if(MD_ON)Set_Temp();
	Set_Poten();
//数组空间申请  
	fprintf(stderr,"PF\n");
	do{
		Tub=(tubulin *)malloc(N*(sizeof(tubulin)));
		ttub=(tubulin *)malloc(ttubno*(sizeof(tubulin)));
		Head=(int *)malloc(Nbox*(sizeof(int)));
		Top=(int *)malloc(13*(sizeof(int)));
		TubeEnd=(int *)malloc(13*(sizeof(int)));
		List=(int *)malloc(N*(sizeof(int)));
		CList=(int *)malloc(26*Nbox*(sizeof(int)));
		Chain=(int *)malloc(CNo*(sizeof(int)));
		TGroup=(int **)malloc(GNo*(sizeof(int*)));
		Tran=(double *)malloc(13*(sizeof(double)));
		Tmark=(int *)malloc(13*(sizeof(int)));
#if EOut_ON
		Trec=Times*MTimes/Out_Step;
		EPoten=(double *)malloc(Trec*(sizeof(double)));
		EVel=(double *)malloc(Trec*(sizeof(double)));
		ERot=(double *)malloc(Trec*(sizeof(double)));
		ETot=(double *)malloc(Trec*(sizeof(double)));
		FAvg=(double *)malloc(Trec*(sizeof(double)));
		FBall=(double *)malloc(Trec*(sizeof(double)));
#endif
	}while(!Tub || !Head || !List || !CList);
	srand((unsigned)time(NULL));
	if(numprocs!=1){
		for(i=0;i<myid;i++)srand((unsigned)rand());
	}
	fprintf(stderr,"rand\n");
#if EOut_ON
	for(i=0;i<Trec;i++){
		EPoten[i]=EVel[i]=ERot[i]=ETot[i]=FAvg[i]=FBall[i]=0;
	}
	Trec=0;
#endif
	fprintf(stderr,"init\n");

	Set_CellList();
	fprintf(stderr,"hl\n");
//初始粒子生成  
	if(myid==0){
		fprintf(stderr,"N: %d Ncell: %d Ncsq: %d Nbox: %d\nNumProcs: %d MpiSec: %d TubSize: %d\n",N,NCell,NCsq,Nbox,numprocs,MpiSec,TubSize);
		fprintf(fout,"N: %d Ncell: %d Ncsq: %d Nbox: %d\nNumProcs: %d MpiSec: %d TubSize: %d\n",N,NCell,NCsq,Nbox,numprocs,MpiSec,TubSize);
		for(i=0;i<N;i++){
			
			Tubp=&Tub[i];
			Tubp->Active=Set_Position(i,r,angle,FilaMisSt,FilaMisEd);
			if(Unit_Type){
				Tubp->xi=xi1;
				Tubp->Type=1;
			}
			else{
				Tubp->xi=xi0;
				Tubp->Type=0;
			}
			if(TipSinA>0 && i>=NInit){
				Tubp->RFMark=1;
				Tubp->xi=xi0;
				Tubp->Type=0;
			}
			if(Whole_Struc==2){
				switch(TubState[i]){
					case 'd':
					case 'D':{
						Tubp->Active=1;
						Tubp->RFMark=0;
						Tubp->xi=xi1;
						Tubp->Type=1;
						break;
					}
					case 't':
					case 'T':{
						Tubp->Active=1;
						Tubp->RFMark=1;
						Tubp->xi=xi0;
						Tubp->Type=0;
						break;
					}
					case 'n':
					case 'N':
					default:Tubp->Active=0;
				}
			}
#if CAP_ON
			if(CAP_ON>0 && i>=N-CAP_ON){
				Tubp->RFMark=1;
				Tubp->xi=xi0;
				Tubp->Type=0;
			}
			else if(CAP_ON<0 && i>=N+CAP_ON){
				Tubp->RFMark=0;
				Tubp->xi=xi1;
				Tubp->Type=1;
			}
#endif
#if Test_Type=='T'
			if(CAPNo>0 && i>=CAPBegin && i<CAPBegin+CAPNo){
				Tubp->RFMark=1;
				Tubp->xi=xi0;
				Tubp->Type=0;
			}
			else if(CAPNo<0 && i>=CAPBegin && i<CAPBegin-CAPNo){
				Tubp->RFMark=0;
				Tubp->xi=xi1;
				Tubp->Type=1;
			}
#endif
			PeriodBoundary(r);
			Set_Tub(Tubp,r,angle,Tubp->xi);
			for(j=0;j<6;j++)Tubp->QDRate[j]=1.0;
			for(j=0;j<3;j++){
				Tubp->rinit[j]=r[j];
			}
			Tubp->HydroTime=0;
			Tubp->RFMark=0;
		#if MD_ON
			for(j=0;j<3;j++){
				Tubp->v[j]=0;
				Tubp->a[j]=0;
				Tubp->apre[j]=0;
				Tubp->vomi[j]=0;
				Tubp->aomi[j]=0;
				Tubp->aomipre[j]=0;
				Tubp->mom[j]=0;
				Tubp->p_poten[j]=1;
			}
		#endif
			
			fprintf(fout,"\ni: %d t0: %g t1: %g t2: %g r0: %g r1: %g r2: %g ",i,Tubp->t[0],Tubp->t[1],Tubp->t[2],Tubp->r[0],Tubp->r[1],Tubp->r[2]);
			Tubp->poten=0;
			Tubp->PotenB=0;
			Tubp->PotenA=0;
			List[i]=-1;
			icell=GetCellNo(Tubp->r);
			List[i]=Head[icell];
			Head[icell]=i;
			Tubp->front=Tubp->back=Tubp->left=Tubp->right=Tubp->lseam=Tubp->rseam=-1;
			Tubp->WallMark=0;
			Tubp->Tub_No=i;
			Tubp->CMark=-1;
			Tubp->ChainNo=-1;
			Tubp->LineNo=-1;
			Tubp->Tub_Id=i/MpiSec;
			ra=1.0*rand()/(RAND_MAX+0.1);
			Tubp->Thd=-log(1-ra)/HD_rate+HD_Delay;
			Tubp->Phd=1.0*rand()/(RAND_MAX+0.1);
			Tubp->Puhd=1.0;
			if(GROW_ON && HD_ON)Tubp->Tin=(N-i)*GSteps/10;
			else Tubp->Tin=0;
			if(RFRate){
				ra=1.0*rand()/(RAND_MAX+0.1);
				if(ra<RFRate){
					Tubp->RFMark=1;
					Tubp->xi=xi0;
					Tubp->Type=0;
				}
			}
#if CV_RFTest==1
			if(i>N-100){
				Tubp->RFMark=0;
				Tubp->xi=xi1;
				Tubp->Type=1;
			}
			else if(i>=26 && !Tubp->RFMark){
				Tubp->xi=xi1;
				Tubp->Type=1;
			}
#elif CV_RFTest==2
			if(i<26){
				Tubp->RFMark=0;
				Tubp->xi=xi1;
				Tubp->Type=1;
			}
			else if(i>=N-13){
				Tubp->RFMark=1;
				Tubp->xi=xi0;
				Tubp->Type=0;
			}
#endif
#if LongRFTest==1 && HD_ON!=2
			Tubp->Tin=(N-i)*GSteps;
#elif LongRFTest==2
			Tubp->Tin=-i*GSteps;
#endif
#if Test_Type=='X'
			if(i>=N-13*6){
				j=i%13;
				if(j>6)j=13-j;
				if(i>=N-13*(6-j))Tubp->Active=0;
			}
#endif
		}
	}
#if RFUnits && RFLayers
	Set_ReinForce();
#endif
	if(QDDelta)Set_QDRate();
	Set_NeighbourList(0);
	if(QDDelta){
		FILE *ftmp;
		char Stmp[MaxFName];
		sprintf(Stmp,"QD_Q%gT%0.2X_%s.txt",1.0*QDDelta,QDType,tstring);
		ftmp=fopen(Stmp,"wt");
		for(i=0;i<N;i++){
			Tubp=&Tub[i];
			if(!Tubp->Active)continue;
			fprintf(stderr,"QDRate\ti:%d",i);
			fprintf(ftmp,"%d",i);
			for(j=0;j<6;j++){
				fprintf(stderr,"\t%g",Tubp->QDRate[j]);
				fprintf(ftmp,"\t%g",Tubp->QDRate[j]);
			}
			fprintf(stderr,"\n");
			fprintf(ftmp,"\n");
		}
		fclose(ftmp);
	}
	
	BigBall.r[0]=Tub[N/2].r[0];
	BigBall.r[1]=0;
	BigBall.r[2]=5+BallR;
//广播初始数据  
	
	if(numprocs!=1){
		MPI_Bcast(Tub,N*TubSize,MPI_CHAR,0,MPI_COMM_WORLD); 
		MPI_Bcast(Head,Nbox,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(List,N,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&BigBall,sizeof(ball),MPI_CHAR,0,MPI_COMM_WORLD);
	}
	fprintf(stderr,"bcast\n");
//初态能量确定 & 构型列表生成  

	for(i=0;i<13;i++)Top[i]=-1;
#if MD_ON
		for(i=0;i<N;i++){
			Tubp=&Tub[i];
			icell=GetCellNo(Tubp->r);
			for(k=0;k<26;k++){
				l=CList[26*icell+k];
				for(j=Head[l];j!=-1;j=List[j]){
					Force(Tubp,&Tub[j]);
				}
			}

			for(j=Head[icell];j!=-1;j=List[j]){
				Force(Tubp,&Tub[j]);
			}
			for(j=0;j<3;j++){		
				Tubp->v[j]=0;
				Tubp->a[j]=0;
				Tubp->apre[j]=0;
				Tubp->vomi[j]=0;
				Tubp->aomi[j]=0;
				Tubp->aomipre[j]=0;
				Tubp->mom[j]=0;
			}
		}
	fprintf(stderr,"force\n");
#elif FixLink
	{
		double PotenTmp[Max_Add_Link];
		for(i=0;i<N;i++){
			Tubp=&Tub[i];
			if(!Tubp->Active)continue;
			for(j=0;j<Max_Add_Link;j++){
				Tubp->LinkList[j]=-1;
				Tubp->PosList[j]=-1;
				Tubp->PotenList[j]=0;
				Tubp->PPMark[j]=0;
			}
			Tubp->AddNum=0;
		}
		for(i=0;i<N;i++){
			Tubp=&Tub[i];
			if(!Tubp->Active)continue;
			icell=GetCellNo(Tubp->r);
			for(k=0;k<26;k++){
				l=CList[26*icell+k];
				for(j=Head[l];j!=-1;j=List[j]){
					if(!Tub[j].Active)continue;
					if(j<=i)continue;
					Set_LinkList(Tubp,&Tub[j]);
				}
			}
			for(j=Head[icell];j!=-1;j=List[j]){
				if(!Tub[j].Active)continue;
				if(j<=i)continue;
				Set_LinkList(Tubp,&Tub[j]);
			}
			for(j=0;j<Tubp->AddNum;j++)PotenTmp[j]=0;
			poten_fl(Tubp,PotenTmp);
		}
		for(i=0;i<N;i++){
			Tubp=&Tub[i];
			if(!Tubp->Active)continue;
			fprintf(stderr,"i: %d\tAddN: %d\n",i,Tubp->AddNum);
			fprintf(stderr,"LinkTub:");
			for(j=0;j<Tubp->AddNum;j++)fprintf(stderr,"\t%d",Tubp->LinkList[j]);
			fprintf(stderr,"\nLinkMark:");
			for(j=0;j<Tubp->AddNum;j++)fprintf(stderr,"\t%0.2X",Tubp->PPMark[j]);
			fprintf(stderr,"\nLinkPoten:");
			for(j=0;j<Tubp->AddNum;j++)fprintf(stderr,"\t%g",Tubp->PotenList[j]);
			fprintf(stderr,"\n");
		}
	}
#else
	{
		int p,q;
		for(i=0;i<N;i++){
			Tubp=&Tub[i];
			if(!Tubp->Active)continue;
			for(j=0;j<12;j++)Tubp->MyPoten[j]=0;
			q=Max_Nbr_Num*i;
			for(p=1;p<=NbrList[q];p++){
				j=NbrList[q+p];
				potential(Tubp,&Tub[j],Tubp->MyPoten);
			}
		}
	}
#endif
	for(i=0;i<13;i++){
		for(j=i,k=0;j!=-1;j=Tub[j].front,k++){
			Tub[j].CMark=i;
			Tub[j].LineNo=k;
			if(Tub[j].front==-1){
				Tub[j].ChainNo=1;
				Top[i]=j;
			}
		}
	}
	for(i=0;i<13;i++){
		for(j=Top[i],k=1;j!=-1;j=Tub[j].back,k++){
			Tub[j].ChainNo=k;
		}
	}
	for(i=0;i<N;i++){
		Tubp=&Tub[i];
		if(!Tubp->Active)continue;
		fprintf(stderr,"i: %d\tf: %d\tb: %d\tl: %d\tr: %d\n",Tubp->Tub_No,Tubp->front,Tubp->back,Tubp->left,Tubp->right);
	}
	fprintf(stderr,"done\n");
}

int ReInit(char *SourceFile){
	int i;
	FILE *fsource;

	fprintf(stderr,"Restart from %s!\n",SourceFile);	//初始化  
	fsource=fopen(SourceFile,"r");
	if(!fsource){
		fprintf(stderr,"Cannot open file %s!\n",SourceFile);
		exit(1);
	}
	init();
	
	fread(&N,sizeof(int),1,fsource);
	Tub=(tubulin *)realloc(Tub,N*(sizeof(tubulin)));
	fread(Tub,sizeof(tubulin),N,fsource);
	fclose(fsource);
	return 0;
}

void AddTub(){
	int i,j,k,l,p,q,icell,TGno,topno=13;
	double r[3],angle[4],dxi,dangle[4],ra,length,dr[3],xi,SumProb=0;
	tubulin *Tubp;
	double TmpProb[13];
	
	if(myid==0){
		fprintf(stderr,"N: %d Ncell: %d Ncsq: %d Nbox: %d\nNumProcs: %d BlockSize: %d\n",N,NCell,NCsq,Nbox,numprocs,BlockSize);
		for(l=0;l<13;l++){
			if(Top[l]==-1){
				fprintf(stderr,"%d (%d)\n",Top[l],l);
				for(i=0;i<13;i++){
					if(Tub[i].back>=0)continue;
					k=Tub[i].CMark;
					if(k==-1)continue;
					fprintf(stderr,"i:%d\tk:%d\t",i,k);
					for(j=i,p=1;j>=0;j=Tub[j].front,p++){
						Tub[j].LineNo=p;
						if(Tub[j].front<0)Top[k]=j;
						fprintf(stderr,"%d ",j);
					}
					fprintf(stderr,"top %d\n",Top[k]);
					for(j=Top[k],p=1;j>=0;j=Tub[j].back,p++){
						Tub[j].ChainNo=p;
					}
				}
			}
		}
		
		fprintf(stderr,"\nTop:\t");
		fprintf(fout,"\nTop:\t");
		for(j=0;Top[j]!=-1&&j<13;j++){
			fprintf(stderr,"%d (%d)\t",Top[j],j);
			fprintf(fout,"%d (%d)\t",Top[j],j);
		}
			
		for(p=N-1;!Tub[p].Active;p--);
		fprintf(stderr,"\ngrow\n");
		fprintf(fout,"\ngrow\n");
		if(Tub[p].back>0){
			TGno=0;
			for(i=0;i<ttubno;i++){
				Tubp=&Tub[Top[i]];

				r[0]=Tubp->r[0]+4.0*Tubp->t[0];
				r[1]=Tubp->r[1]+4.0*Tubp->t[1];
				r[2]=Tubp->r[2]+4.0*Tubp->t[2];

				move2(dr,dangle,10.0,2);
				QuaTimes(Tubp->angle,dangle,angle);
				xi=xi0;
				
				Tubp=&ttub[i];
				PeriodBoundary(r);
				for(j=0;j<3;j++){
					Tubp->rinit[j]=r[j];
				}
				Tubp->HydroTime=0;
				Tubp->xi=xi0;
				Tubp->Type=0;
			#if MD_ON
				for(j=0;j<3;j++){
					Tubp->v[j]=0;
					Tubp->a[j]=0;
					Tubp->apre[j]=0;
					Tubp->vomi[j]=0;
					Tubp->aomi[j]=0;
					Tubp->aomipre[j]=0;
					Tubp->mom[j]=0;
					Tubp->p_poten[j]=1;
				}
			#endif
				
				Set_Tub(Tubp,r,angle,xi);
				Tubp->poten=0;
				Tubp->PotenB=0;
				Tubp->PotenA=0;
				Tubp->front=Tubp->back=Tubp->left=Tubp->right=Tubp->lseam=Tubp->rseam=-1;
				Tubp->Active=1;
				Tubp->WallMark=0;
				Tubp->CMark=i;
				Tubp->ChainNo=-1;
				Tubp->Tub_No=-2-i;
				Tubp->Tin=0;
				Tubp->LineNo=Tub[Top[i]].LineNo+1;
				ra=1.0*rand()/(RAND_MAX+0.1);
				Tubp->Thd=-log(1-ra)/HD_rate+HD_Delay;
				Tubp->Phd=1.0*rand()/(RAND_MAX+0.1);
				Tubp->Puhd=1.0;
				Tubp->RFMark=0;
				for(j=0;j<6;j++)Tubp->QDRate[j]=1.0;
				if(RFRate){
					ra=1.0*rand()/(RAND_MAX+0.1);
					if(ra<RFRate){
						Tubp->RFMark=1;
						Tubp->xi=xi0;
						Tubp->Type=0;
					}
				}
			#if MD_ON
				icell=GetCellNo(Tubp->r);
				for(p=0;p<26;p++){
					q=CList[26*icell+p];
					for(j=Head[q];j!=-1;j=List[j]){
						Force(Tubp,&Tub[j]);
					}
				}
				for(j=Head[icell];j!=-1;j=List[j]){
					Force(Tubp,&Tub[j]);
				}
				for(j=0;j<3;j++){
					Tubp->v[j]=0;
					Tubp->a[j]=0;
					Tubp->apre[j]=0;
					Tubp->vomi[j]=0;
					Tubp->aomi[j]=0;
					Tubp->aomipre[j]=0;
					Tubp->mom[j]=0;
				}
			#else
				{
					for(j=0;j<12;j++)Tubp->MyPoten[j]=0;
					icell=GetCellNo(Tubp->r);
					for(k=0;k<26;k++){
						l=CList[26*icell+k];
						for(j=Head[l];j!=-1;j=List[j]){
							if(!Tub[j].Active)continue;
							potential(Tubp,&Tub[j],Tubp->MyPoten);
						}
					}
					for(j=Head[icell];j!=-1;j=List[j]){
						if(!Tub[j].Active)continue;
						potential(Tubp,&Tub[j],Tubp->MyPoten);
					}
					Tubp->poten=0;
					for(j=0;j<4;j++)Tubp->poten+=Tubp->MyPoten[3*j];
				}
			#endif
				ra=1.0*rand()/(RAND_MAX+0.1);
#if GROW_ON==1
				if(exp(-G_rate*(Tubp->poten+EBase))>ra){
					Tmark[i]=1;
					TGno++;
				}
				else Tmark[i]=0;
#elif GROW_ON==2
				if(Tubp->poten<-0.63 && Rate_In>ra){
					Tmark[i]=1;
					TGno++;
				}
				else if(Tubp->poten<-0.53 && Rate_Side>ra){
					Tmark[i]=1;
					TGno++;
				}
				else if(Tubp->poten<-0.40 && Rate_On>ra){
					Tmark[i]=1;
					TGno++;
				}
				else Tmark[i]=0;
#elif GROW_ON==3
				Tmark[i]=0;
				TmpProb[i]=4;
				if(Tubp->left>=0)TmpProb[i]-=1;
				if(Tubp->right>=0)TmpProb[i]-=1;
				TmpProb[i]*=exp(-G_rate*(Tubp->poten+EBase));
				Tmark[i]=0;
				SumProb+=TmpProb[i];
#endif
				fprintf(stderr,"j %d\tr %g\t %g\t %g\tpoten %g\t%d\tProb %g\n",i,Tubp->r[0],Tubp->r[1],Tubp->r[2],Tubp->poten,Tmark[i],TmpProb[i]);
			}
#if GROW_ON==3
			if(!MustGrow && SumProb<6.0){
				SumProb=6.0;
			}
			ra=SumProb*rand()/(RAND_MAX+0.1);
			if(ra>0){
				for(i=0;i<13;i++){
					ra-=TmpProb[i];
					if(ra<0)break;
				}
				if(i<13){
					Tmark[i]=1;
					TGno=1;
				}
			}
			else TGno=0;
#endif
#if MustGrow
			if(MaxEnd<15 && TGno==0){
				double Dtmp=0;
				k=-1;
				for(i=0;i<ttubno;i++){
					Tubp=&ttub[i];
					if(Dtmp>Tubp->poten){
						Dtmp=Tubp->poten;
						k=i;
					}
				}
				if(k>0){
					i=k;
					TGno=1;
					Tmark[i]=1;
				}
			}
#endif
			if(TGno==0){
				fprintf(stderr,"Bad grow!(%g)\n",ra);
				fprintf(fout,"Bad grow!(%g)\n",ra);
				return;
			}
			else{
				l=rand()%TGno+1;
				k=0;
				for(j=0;j<topno;j++){
					if(Tmark[j]){
						k++;
						if(k!=l)Tmark[j]=0;
					}
				}
				fprintf(stderr,"grow 1 in %d: %d\n",TGno,l);
				fprintf(fout,"grow 1 in %d: %d\n",TGno,l);
				TGno=1;
			}
			Tub=(tubulin *)realloc(Tub,(N+TGno)*(sizeof(tubulin)));
			List=(int *)realloc(List,(N+TGno)*(sizeof(int)));
			for(j=0;j<topno;j++){
				if(Tmark[j]){

					Tub[N]=ttub[j];
					Tubp=&Tub[N];
					Top[j]=N;
					fprintf(stderr,"New j:%d r %g\t%g\t%g\t poten %g\tsize %d\n",j,Tubp->r[0],Tubp->r[1],Tubp->r[2],Tubp->poten,TubSize);
					icell=GetCellNo(Tubp->r);
					fprintf(stderr,"icell %d\tN %d\tj %d\tTGno %d ra %g\n",icell,N,j,TGno,ra);
					fprintf(fout,"icell %d\tN %d\tj %d\tTGno %d ra %g\n",icell,N,j,TGno,ra);
					List[N]=Head[icell];
					Head[icell]=N;
					Tubp->Tub_No=N;
					Tubp->CMark=j;
					Tubp->Tub_Id=N/MpiSec;
					for(q=N,p=1;q>=0;q=Tub[q].back,p++){
						Tub[q].ChainNo=p;
					}
					fprintf(stderr,"No: %d\tRFMark: %d\tThd: %d\tPhd: %g\n",N,Tubp->RFMark,Tubp->Thd,Tubp->Phd);
					#if FixLink
						for(i=0;i<Max_Add_Link;i++){
							Tubp->LinkList[i]=-1;
							Tubp->PosList[i]=-1;
							Tubp->PotenList[i]=0;
							Tubp->PPMark[i]=0;
						}
						Tubp->AddNum=0;
						
						icell=GetCellNo(Tubp->r);
						for(k=0;k<26;k++){
							l=CList[26*icell+k];
							for(q=Head[l];q!=-1;q=List[q]){
								if(!Tub[q].Active)continue;
								Set_LinkList(Tubp,&Tub[q]);
							}
						}
						for(q=Head[icell];q!=-1;q=List[q]){
							if(!Tub[q].Active)continue;
							if(q==N)continue;
							Set_LinkList(Tubp,&Tub[q]);
						}
						
						fprintf(stderr,"j: %d\tAddN: %d\n",j,Tubp->AddNum);
						if(!Tubp->AddNum){
							Tubp->Active=0;
							Top[Tubp->CMark]=Tubp->back;
							for(q=Tubp->back,p=1;q>=0;q=Tub[q].back,p++){
								Tub[q].ChainNo=p;
							}
							fprintf(stderr,"Fails to grow.\n");
						}
						else{
							fprintf(stderr,"LinkTub:");
							for(q=0;q<Tubp->AddNum;q++)fprintf(stderr,"\t%d",Tubp->LinkList[q]);
							fprintf(stderr,"\nLinkMark:");
							for(q=0;q<Tubp->AddNum;q++)fprintf(stderr,"\t%0.2X",Tubp->PPMark[q]);
							fprintf(stderr,"\nLinkPoten:");
							for(q=0;q<Tubp->AddNum;q++)fprintf(stderr,"\t%g",Tubp->PotenList[q]);
							fprintf(stderr,"\n");
						}
					#endif
					N+=1;
				}
			}
		#if !MD_ON && !FixLink
			NbrList=(int *)realloc(NbrList,Max_Nbr_Num*N*sizeof(int));
			Set_NeighbourList(1);
		#endif
			fprintf(stderr,"Grow done.\n");
			fprintf(fout,"Grow done.\n");
			if(N%numprocs)MpiSec=N/numprocs+1;
			else MpiSec=N/numprocs;
		}
		else{
			fprintf(stderr,"Not Grow!\n");
			fprintf(fout,"Not Grow!\n");
			Tubp=&Tub[p];
			length=10;
			for(j=0;Top[j]!=-1&&j<topno;j++){
				length=(Tubp->r[0]-Tub[Top[j]].r[0])*(Tubp->r[0]-Tub[Top[j]].r[0])+(Tubp->r[1]-Tub[Top[j]].r[1])*(Tubp->r[1]-Tub[Top[j]].r[1])+(Tubp->r[2]-Tub[Top[j]].r[2])*(Tubp->r[2]-Tub[Top[j]].r[2]);
				fprintf(stderr,"\nTest:\tNo:%d\tj:%d\tlength:%g",p,Top[j],length);
				fprintf(fout,"\nTest:\tNo:%d\tj:%d\tlength:%g",p,Top[j],length);
				if(length<(cutoff+2)*(cutoff+2)&&p!=Top[j]){
					length=-1;break;
				}
			}
			if(length>0){
				fprintf(stderr,"\nDestroy:\n%d\tlen: %g",p,length);
				fprintf(fout,"\nDestroy:\n%dlen: %g",p,length);
				Tubp->Active=0;
				Top[Tubp->CMark]=Tubp->back;
				for(j=Tubp->back,p=1;j>=0;j=Tub[j].back,p++){
					Tub[j].ChainNo=p;
				}
				if(Tubp->front>=0){Tub[Tubp->front].back=-1;Tubp->front=-1;}
				if(Tubp->back>=0){Tub[Tubp->back].front=-1;Tubp->back=-1;}
				if(Tubp->left>=0){Tub[Tubp->left].right=-1;Tubp->left=-1;}
				if(Tubp->right>=0){Tub[Tubp->right].left=-1;Tubp->right=-1;}
			}
			else if(FixLink && !Tubp->AddNum){
				for(i=0;i<Max_Add_Link;i++){
					Tubp->LinkList[i]=-1;
					Tubp->PosList[i]=-1;
					Tubp->PotenList[i]=0;
					Tubp->PPMark[i]=0;
				}
				Tubp->AddNum=0;
				
				icell=GetCellNo(Tubp->r);
				for(k=0;k<26;k++){
					l=CList[26*icell+k];
					for(q=Head[l];q!=-1;q=List[q]){
						if(!Tub[q].Active)continue;
						Set_LinkList(Tubp,&Tub[q]);
					}
				}
				for(q=Head[icell];q!=-1;q=List[q]){
					if(!Tub[q].Active)continue;
					if(q==N)continue;
					Set_LinkList(Tubp,&Tub[q]);
				}
				
				fprintf(stderr,"j: %d\tAddN: %d\n",j,Tubp->AddNum);
				fprintf(stderr,"LinkTub:");
				for(q=0;q<Tubp->AddNum;q++)fprintf(stderr,"\t%d",Tubp->LinkList[q]);
				fprintf(stderr,"\nLinkMark:");
				for(q=0;q<Tubp->AddNum;q++)fprintf(stderr,"\t%0.2X",Tubp->PPMark[q]);
				fprintf(stderr,"\nLinkPoten:");
				for(q=0;q<Tubp->AddNum;q++)fprintf(stderr,"\t%g",Tubp->PotenList[q]);
				fprintf(stderr,"\n");
			}
		}
	}
}

#if FixLink
void Monte_Carlo(int M){
	int i,j,k,l,p,q,icell,jcell,ijeq=0;
	double poten,PotenB,PotenA,dr[3],dangle[5],r[3],angle[5],dv,ra,xi,ra1,Vxi,anglemod;//,a,b,c;
	int d,topno=13,minno=-1,TGno;
	double PotenTmp[Max_Add_Link];
	const double RandPara=1.0/(RAND_MAX+0.1);
	tubulin *Tubp,tubtmp;

	RAccept=1.0*Times*N;
	at1=at2=0.0;
	NA_A=0;

fprintf(stderr,"Monte Carlo\t%d\n",M);
	for(k=0;k<Times;k++){
		for(l=0;l<N;l++){
			#if Whole_Struc
	//			if(M%50==0)
					i=(int)(rand()*RandPara*N);
	//			else i=13+(int)(rand()/(RAND_MAX+1.0)*(N-13));
			#else
				i=(int)(rand()*RandPara*N);
			#endif

			tubtmp=Tub[i];
			Tubp=&tubtmp;
			if(!Tubp->Active)continue;

			if(!Tubp->Type)xi=xi0;
			else xi=xi1;
			
			for(j=0;j<Tubp->AddNum;j++)PotenTmp[j]=0;
			
			poten_fl(Tubp,PotenTmp);
			Tubp->poten=0;
			for(j=0;j<Tubp->AddNum;j++)Tubp->poten+=PotenTmp[j];

			if(KT>1)move2(dr,dangle,3.0,l%2);
			else move2(dr,dangle,1.0,l%2);
			
			#if BCType
			if(M || k%5){
				#if BCType & 0xF0
					#if Whole_Struc
					if(i<13)
					#else
					if(i==0 || i==Lat_Length-1)
					#endif
					{
						#if BCType & 0x10
						for(j=0;j<3;j++)dr[j]=0;
						#elif BCType & 0x20
							dr[1]=0;
							#if ForceXYZ & 0x01
							dr[0]=0;
							#endif
							#if ForceXYZ & 0x02
							dr[2]=0;
							#endif
							#if ForceXYZ & 0x04
							dr[2]=0;
							#endif
						#elif BCType & 0x40
						for(j=0;j<3;j++)dr[j]=0;
						dangle[3]=0;
						dangle[4]=1;
						#endif
					}
				#endif
				#if BCType & 0x0F
					#if Whole_Struc
					if(i>=N-13)
					#else
					if(i==N-1 || i==N-Lat_Length)
					#endif
					{
						#if BCType & 0x01
						for(j=0;j<3;j++)dr[j]=0;
						#elif BCType & 0x02
							dr[1]=0;
							#if ForceXYZ & 0x01
							dr[0]=0;
							#endif
							#if ForceXYZ & 0x02
							dr[2]=0;
							#endif
							#if ForceXYZ & 0x04
							dr[2]=0;
							#endif
						#elif BCType & 0x04
						for(j=0;j<3;j++)dr[j]=0;
						dangle[3]=0;
						dangle[4]=1;
						#endif
					}
				#endif
			}
			#endif

			if(l%2){
				QuaTimes(Tubp->angle,dangle,angle);
				Set_Tub(Tubp,Tubp->r,angle,xi);
			}
			else{
				for(j=0;j<3;j++){
					r[j]=Tubp->r[j]+dr[j];
				}
				PeriodBoundary(r);
				Set_Tub(Tubp,r,Tubp->angle,xi);
			}
			for(j=0;j<Tubp->AddNum;j++){
				PotenTmp[j]=0;
				Tubp->PPMark[j]=Tubp->PPMark[j]&0x7F;
			}
			poten_fl(Tubp,PotenTmp);
			poten=0;
			for(j=0;j<Tubp->AddNum;j++)poten+=PotenTmp[j];
			
			if(poten>Tubp->poten)
			{				
				at1+=1.0;
				dv=-(poten-Tubp->poten)/KT;
				ra=rand()*RandPara;
				if(exp(dv)<ra){
					if(icell!=jcell)ijeq--;
					RAccept-=1.0;
					at2+=1.0;
					if(l%2)NA_A++;
					
					Tub[i].poten=Tubp->poten/2;
					continue;
				}
			}
			for(j=0;j<Tubp->AddNum;j++){
				Tubp->PPMark[j]=Tubp->PPMark[j]|0x80;
				p=Tubp->LinkList[j];
				q=Tubp->PosList[j];
				Tub[p].PotenList[q]=Tubp->PotenList[j];
				Tub[p].PPMark[q]=Tub[p].PPMark[q]|0x80;
			}
			Tubp->poten=poten/2;
			Tub[i]=tubtmp;
		}
		if(k%20==10){
			
		#if HD_ON
			for(i=0;i<N;i++){
				Tubp=&Tub[i];
				if(Tubp->Type!=1 && i>25 && Tubp->Tin>0){
			#if HD_ON==1
					if(Tubp->front>=0 && Tubp->Thd<Tubp->Tin)Tubp->Type=1;
			#elif HD_ON==2
					double ActiveE,Cnow;
					ActiveE=0;
					if(Tubp->front>=0)ActiveE+=0.161+Tubp->potenf;
					if(Tubp->back>=0)ActiveE+=0.161+Tubp->potenb;
					Cnow=7.2e11*exp(18.5*ActiveE);
					Cnow=exp(-Cnow*SlowTime);
					Tubp->Phd-=Tubp->Puhd*(1-Cnow);
					if(Tubp->Phd<0)Tubp->Type=1;
					else Tubp->Puhd*=Cnow;
			#endif
					if(Tubp->RFMark)Tubp->Type=0;
			#if Test_Type=='W'
					if(i<130)Tubp->Type=0;
			#endif
					if(Tubp->Type)Tubp->xi=xi1;
					else Tubp->xi=xi0;
				}
				Tubp->Tin+=20;
			}
		#endif
		#if GROW_ON
			if(Ngrow>=TGrow){
				Ngrow=0;
#if Test_Type=='W'
				if(MaxEnd<40)AddTub();
#else
				AddTub();
#endif
			#if GRStep
				TGrow=TGrow+(int)(GRStep*(2.0*rand()/(RAND_MAX+0.1)-1));
				if(TGrow<1000)TGrow=1000;
			#endif
			}
			Ngrow+=20;
		#endif
		
			for(i=0;i<Nbox;i++)Head[i]=-1;
			for(i=0;i<N;i++){
				Tubp=&Tub[i];
				icell=GetCellNo(Tubp->r);
				List[i]=Head[icell];
				Head[icell]=i;
				anglemod=Tubp->angle[0]*Tubp->angle[0]+Tubp->angle[1]*Tubp->angle[1]+Tubp->angle[2]*Tubp->angle[2];
				anglemod=1/sqrt(anglemod);
				Tubp->angle[0]*=anglemod;
				Tubp->angle[1]*=anglemod;
				Tubp->angle[2]*=anglemod;
				anglemod=1/sqrt(Tubp->angle[3]*Tubp->angle[3]+Tubp->angle[4]*Tubp->angle[4]);
				Tubp->angle[3]*=anglemod;
				Tubp->angle[4]*=anglemod;
				Set_Tub(Tubp,Tubp->r,Tubp->angle,Tubp->xi);
			}
		}
		#if GROW_ON
		if(k%500==10){
			for(i=0;i<N;i++){
				Tubp=&Tub[i];
				if(!Tubp->Active)continue;
				for(j=0;j<Max_Add_Link;j++){
					Tubp->LinkList[j]=-1;
					Tubp->PosList[j]=-1;
					Tubp->PotenList[j]=0;
					Tubp->PPMark[j]=0;
				}
				Tubp->AddNum=0;
			}
			for(i=0;i<N;i++){
				Tubp=&Tub[i];
				if(!Tubp->Active)continue;
				icell=GetCellNo(Tubp->r);
				for(p=0;p<26;p++){
					l=CList[26*icell+p];
					for(j=Head[l];j!=-1;j=List[j]){
						if(!Tub[j].Active)continue;
						if(j<=i)continue;
						Set_LinkList(Tubp,&Tub[j]);
					}
				}
				for(j=Head[icell];j!=-1;j=List[j]){
					if(!Tub[j].Active)continue;
					if(j<=i)continue;
					Set_LinkList(Tubp,&Tub[j]);
				}
			}
		}
		#endif
	}
}
#else
void Monte_Carlo(int M){
	int i,j,k,l,p,q,icell,jcell,ijeq=0;
	double poten,PotenB,PotenA,dr[3],dangle[5],r[3],angle[5],dv,ra,xi,ra1,Vxi,anglemod;//,a,b,c;
	int d,topno=13,minno=-1,TGno;
	double PotenTmp[12];
	const double RandPara=1.0/(RAND_MAX+0.1);
	tubulin *Tubp,tubtmp;

	RAccept=1.0*Times*N;
	at1=at2=0.0;
	NA_A=0;

	for(k=0;k<Times;k++){
		for(l=0;l<N;l++){
			#if Whole_Struc
	//			if(M%50==0)
					i=(int)(rand()*RandPara*N);
	//			else i=13+(int)(rand()/(RAND_MAX+1.0)*(N-13));
			#else
				i=(int)(rand()*RandPara*N);
			#endif

			tubtmp=Tub[i];
			Tubp=&tubtmp;
			if(!Tubp->Active)continue;

			if(!Tubp->Type)xi=xi0;
			else xi=xi1;
			for(j=0;j<12;j++)Tubp->MyPoten[j]=0;
			
			q=Max_Nbr_Num*i;
			for(p=1;p<=NbrList[q];p++){
				j=NbrList[q+p];
				potential(Tubp,&Tub[j],Tubp->MyPoten);
			}
			Tubp->PotenB=0;
			Tubp->PotenA=0;
			for(j=0;j<6;j++){
				Tubp->PotenB+=Tubp->MyPoten[j];
				Tubp->PotenA+=Tubp->MyPoten[6+j];
			}
			Tubp->poten=Tubp->PotenB+Tubp->PotenA;

			if(KT>1)move2(dr,dangle,3.0,l%2);
			else move2(dr,dangle,1.0,l%2);
			
			#if BCType
			if(M || k%5){
				#if BCType & 0xF0
					#if Whole_Struc
					if(i<13)
					#else
					if(i==0 || i==Lat_Length-1)
					#endif
					{
						#if BCType & 0x10
						for(j=0;j<3;j++)dr[j]=0;
						#elif BCType & 0x20
							dr[1]=0;
							#if ForceXYZ & 0x01
							dr[0]=0;
							#endif
							#if ForceXYZ & 0x02
							dr[2]=0;
							#endif
							#if ForceXYZ & 0x04
							dr[2]=0;
							#endif
						#elif BCType & 0x40
						for(j=0;j<3;j++)dr[j]=0;
						dangle[3]=0;
						dangle[4]=1;
						#endif
					}
				#endif
				#if BCType & 0x0F
					#if Whole_Struc
					if(i>=N-13)
					#else
					if(i==N-1 || i==N-Lat_Length)
					#endif
					{
						#if BCType & 0x01
						for(j=0;j<3;j++)dr[j]=0;
						#elif BCType & 0x02
							dr[1]=0;
							#if ForceXYZ & 0x01
							dr[0]=0;
							#endif
							#if ForceXYZ & 0x02
							dr[2]=0;
							#endif
							#if ForceXYZ & 0x04
							dr[2]=0;
							#endif
						#elif BCType & 0x04
						for(j=0;j<3;j++)dr[j]=0;
						dangle[3]=0;
						dangle[4]=1;
						#endif
					}
				#endif
			}
			#endif

			if(l%2){
				QuaTimes(Tubp->angle,dangle,angle);
				Set_Tub(Tubp,Tubp->r,angle,xi);
			}
			else{
				for(j=0;j<3;j++){
					r[j]=Tubp->r[j]+dr[j];
				}
				PeriodBoundary(r);
				Set_Tub(Tubp,r,Tubp->angle,xi);
			}
			for(j=0;j<12;j++)PotenTmp[j]=0;
			q=Max_Nbr_Num*i;
			for(p=1;p<=NbrList[q];p++){
				j=NbrList[q+p];
				potential(Tubp,&Tub[j],PotenTmp);
			}
			PotenB=0;
			PotenA=0;
			for(j=0;j<6;j++){
				PotenB+=PotenTmp[j];
				PotenA+=PotenTmp[6+j];
			}
			poten=PotenB+PotenA;
			
			if(poten>Tubp->poten)
			{				
				at1+=1.0;
				dv=-(poten-Tubp->poten)/KT;
				ra=rand()*RandPara;
				if(exp(dv)<ra){
					if(icell!=jcell)ijeq--;
					RAccept-=1.0;
					at2+=1.0;
					if(l%2)NA_A++;
					
					Tub[i].poten=Tubp->poten/2;
					Tub[i].PotenB=Tubp->PotenB/2;
					Tub[i].PotenA=Tubp->PotenA/2;
					for(j=0;j<12;j++)Tub[i].MyPoten[j]=Tubp->MyPoten[j]/2;

					continue;
				}
			}
			Tubp->poten=poten/2;
			Tubp->PotenB=PotenB/2;
			Tubp->PotenA=PotenA/2;
			for(j=0;j<12;j++)Tubp->MyPoten[j]/=2;
			Tub[i]=tubtmp;
		}

		if(k%20==10){
			
		#if HD_ON
			for(i=0;i<N;i++){
				Tubp=&Tub[i];
				if(Tubp->Type!=1 && i>25 && Tubp->Tin>0){
			#if HD_ON==1
					if(Tubp->front>=0 && Tubp->Thd<Tubp->Tin)Tubp->Type=1;
			#elif HD_ON==2
					double ActiveE,Cnow;
					ActiveE=0;
					if(Tubp->front>=0)ActiveE+=0.161+Tubp->potenf;
					if(Tubp->back>=0)ActiveE+=0.161+Tubp->potenb;
					Cnow=7.2e11*exp(18.5*ActiveE);
					Cnow=exp(-Cnow*SlowTime);
					Tubp->Phd-=Tubp->Puhd*(1-Cnow);
					if(Tubp->Phd<0)Tubp->Type=1;
					else Tubp->Puhd*=Cnow;
			#endif
					if(Tubp->RFMark)Tubp->Type=0;
	//							if(exp(-C_rate*(Tubp->Tin-CTime))<ra)Tubp->Type=1;
			#if Test_Type=='W'
					if(i<130)Tubp->Type=0;
			#endif
					if(Tubp->Type)Tubp->xi=xi1;
					else Tubp->xi=xi0;
				}
				Tubp->Tin+=20;
			}
		#endif
		#if GROW_ON
			if(Ngrow>=TGrow){
				Ngrow=0;
#if Test_Type=='W'
				if(MaxEnd<40)AddTub();
#else
				AddTub();
#endif
			#if GRStep
				TGrow=TGrow+(int)(GRStep*(2.0*rand()/(RAND_MAX+0.1)-1));
				if(TGrow<1000)TGrow=1000;
			#endif
			}
			Ngrow+=20;
		#endif
		
			for(i=0;i<Nbox;i++)Head[i]=-1;
			for(i=0;i<N;i++){
				Tubp=&Tub[i];
				icell=GetCellNo(Tubp->r);
				List[i]=Head[icell];
				Head[icell]=i;
				anglemod=Tubp->angle[0]*Tubp->angle[0]+Tubp->angle[1]*Tubp->angle[1]+Tubp->angle[2]*Tubp->angle[2];
				anglemod=1/sqrt(anglemod);
				Tubp->angle[0]*=anglemod;
				Tubp->angle[1]*=anglemod;
				Tubp->angle[2]*=anglemod;
				anglemod=1/sqrt(Tubp->angle[3]*Tubp->angle[3]+Tubp->angle[4]*Tubp->angle[4]);
				Tubp->angle[3]*=anglemod;
				Tubp->angle[4]*=anglemod;
				Set_Tub(Tubp,Tubp->r,Tubp->angle,Tubp->xi);
			}
			Set_NeighbourList(1);
		}
	}
}
#endif

#if MD_ON
void StdEQuit(tubulin *Tubin){
	int i,j,k,l;
	tubulin *Tubp;
	fprintf(stderr,"No:%d \nTy %d Act %d CNo %d CMark %d LNo %d Tin %d WallM %d\n"
		,Tubin->Tub_No,Tubin->Type,Tubin->Active,Tubin->ChainNo,Tubin->CMark,Tubin->LineNo,Tubin->Tin,Tubin->WallMark);
	fprintf(stderr,"P %g P1 %g P2 %g fr %d bk %d lf %d rt %d\n"
		,Tubin->poten,Tubin->PotenB,Tubin->PotenA,Tubin->front,Tubin->back,Tubin->left,Tubin->right);
	for(i=0;i<3;i++){
		fprintf(stderr,"i:%d r %g rpre %g t %g s %g n %g ch1 %g ch2 %g si1 %g si2 %g si3 %g si4 %g\n"
			,i,Tubin->r[i],Tubin->rpre[i],Tubin->t[i],Tubin->s[i],Tubin->n[i],Tubin->chain1[i],Tubin->chain2[i],Tubin->side1[i],Tubin->side2[i],Tubin->side3[i],Tubin->side4[i]);
		fprintf(stderr,"a %g v %g apre %g aomi %g vomi %g aomipre %g mom %g \n"
			,Tubin->a[i],Tubin->v[i],Tubin->apre[i],Tubin->aomi[i],Tubin->vomi[i],Tubin->aomipre[i],Tubin->mom[i]);
		Tubin->r[i]=Tubin->rpre[i];
	}
	Set_Tub(Tubin,Tubin->r,Tubin->angle,Tubin->xi);
	BUG=1;
	{
		Tubp=Tubin;
		Tubp->poten=Tubp->PotenB=Tubp->PotenA=0;
		for(k=0;k<26;k++){
			l=CList[26*icell+k];
			for(j=Head[l];j!=-1;j=List[j]){
				if(!Tub[j].Active)continue;
				Force(Tubp,&Tub[j]);	
			}
		}
		for(j=Head[icell];j!=-1;j=List[j]){
			if(!Tub[j].Active)continue;
			Force(Tubp,&Tub[j]);
		}
	}
	BUG=0;
	for(i=0;i<Nbox;i++){
		if(Head[i]!=-1)fprintf(stderr,"\ncell:%d\t",i);
		for(j=Head[i];j>=0;j=List[j]){
			fprintf(stderr,"%d ",j);
		}
	}
}

void MoleDynamic(int MoleDynamicNo){
	int i,j,k,l,p,t,m,BcastSize,Npre,Maxline,mincom,VirtualL;
	double dangle[4],angle[4],xi,anglemod,Tnow,ra,ActiveE,Cnow;
	double EFTrans[9];
	tubulin *Tubp,tubtmp;
#if Di_Inertia
	double x,y,z,costheta,sintheta,Tra[9],TraInv[9];
#endif
	fprintf(stderr,"MoleDynamic id:%d\n",myid);
	for(t=0;t<Times;t++){
		TTot++;
		if(BUG)fprintf(stderr,"\nt:%d",t);
		for(m=NStart+myid*MpiSec;m<NStart+(myid+1)*MpiSec && m<N;m++){
#if FixEnd
			if(t && N>60 && m<26)continue;
#endif
			tubtmp=Tub[m];
			Tubp=&tubtmp;
			if(!Tubp->Active)continue;
			if(Whole_Struc && Tubp->LineNo<=LStart)continue;
#if 0//GROW_ON || BREAK_ON
			if(Tubp->Type && Tubp->ChainNo>=25 && (Tubp->left!=-1&&Tubp->right!=-1) && Tubp->Tin-CTime>WallTime){
				if(Tubp->WallMark<50){
					Tubp->WallMark++;
					continue;
				}
				else Tubp->WallMark=1;
			}
			else Tubp->WallMark=0;
#endif
#if Whole_Struc && BREAK_ON
			if(Tubp->left<0 || Tubp->right<0){
				if(Tubp->WallMark<0){
					Tubp->WallMark--;
				}
				else Tubp->WallMark=-1;
			}
			else if(Tubp->WallMark<0){
				Tubp->WallMark=1;
			}
#endif
			Tubp->poten=Tubp->PotenB=Tubp->PotenA=Tubp->potent=0;
			icell=GetCellNo(Tubp->r);
			for(k=0;k<26;k++){
				l=CList[26*icell+k];
				for(j=Head[l];j!=-1;j=List[j]){
					if(!Tub[j].Active)continue;
					Force(Tubp,&Tub[j]);

				}
			}
			for(j=Head[icell];j!=-1;j=List[j]){
				if(!Tub[j].Active)continue;
				if(m==j)continue;
				Force(Tubp,&Tub[j]);
			}
			Tub[m]=tubtmp;
		}

		if(BUG)fprintf(stderr,"t:%d half: %d\n",t,myid);
		GaussRand();
//		if(1)
			{
			Tnow=0;
				for(m=NStart+myid*MpiSec;m<NStart+(myid+1)*MpiSec	&& m<N;m++){
#if FixEnd
					if(t && N>60 && m<26)continue;
#endif
					tubtmp=Tub[m];
					Tubp=&tubtmp;
					if(!Tubp->Active)continue;
					if(Whole_Struc && Tubp->LineNo<=LStart)continue;
#if Di_Inertia
					x=Tubp->angle[0];
					y=Tubp->angle[1];
					z=Tubp->angle[2];
					costheta=cos(2*Tubp->angle[3]);
					sintheta=sin(-2*Tubp->angle[3]);

					Tra[0]=costheta+x*x*(1-costheta);
					Tra[1]=sintheta*z+x*y*(1-costheta);
					Tra[2]=-sintheta*y+x*z*(1-costheta);
					Tra[3]=-sintheta*z+x*y*(1-costheta);
					Tra[4]=costheta+y*y*(1-costheta);
					Tra[5]=sintheta*x+y*z*(1-costheta);
					Tra[6]=sintheta*y+x*z*(1-costheta);
					Tra[7]=-sintheta*x+y*z*(1-costheta);
				    Tra[8]=costheta+z*z*(1-costheta);
			    
			    	Tubp->mombar[0]=(Tra[0]*Tubp->mom[0]+Tra[3]*Tubp->mom[1]+Tra[6]*Tubp->mom[2])/I_M;
				    Tubp->mombar[1]=(Tra[1]*Tubp->mom[0]+Tra[4]*Tubp->mom[1]+Tra[7]*Tubp->mom[2])/I_M/3.5;
				    Tubp->mombar[2]=(Tra[2]*Tubp->mom[0]+Tra[5]*Tubp->mom[1]+Tra[8]*Tubp->mom[2])/I_M/3.5;
				    
				    sintheta=sin(2*Tubp->angle[3]);

					TraInv[0]=costheta+x*x*(1-costheta);
					TraInv[1]=sintheta*z+x*y*(1-costheta);
					TraInv[2]=-sintheta*y+x*z*(1-costheta);
					TraInv[3]=-sintheta*z+x*y*(1-costheta);
					TraInv[4]=costheta+y*y*(1-costheta);
					TraInv[5]=sintheta*x+y*z*(1-costheta);
					TraInv[6]=sintheta*y+x*z*(1-costheta);
					TraInv[7]=-sintheta*x+y*z*(1-costheta);
				    TraInv[8]=costheta+z*z*(1-costheta);
			    
			    	Tubp->aomi[0]=(TraInv[0]*Tubp->mombar[0]+TraInv[3]*Tubp->mombar[1]+TraInv[6]*Tubp->mombar[2]);
				    Tubp->aomi[1]=(TraInv[1]*Tubp->mombar[0]+TraInv[4]*Tubp->mombar[1]+TraInv[7]*Tubp->mombar[2]);
				    Tubp->aomi[2]=(TraInv[2]*Tubp->mombar[0]+TraInv[5]*Tubp->mombar[1]+TraInv[8]*Tubp->mombar[2]);
#else
				    Tubp->aomi[0]=Tubp->mom[0]/I_M;
				    Tubp->aomi[1]=Tubp->mom[1]/I_M;
				    Tubp->aomi[2]=Tubp->mom[2]/I_M;
#endif 
					Tubp->a[0]=Tubp->a[0]/M_F;
					Tubp->a[1]=Tubp->a[1]/M_F;
					Tubp->a[2]=Tubp->a[2]/M_F;
	    
					dangle[3]=0;

					for(p=0;p<3;p++){
#if !Whole_Struc
	#if Test_Type && Test_Type!='p'
						if(m==0)Tubp->aomi[p]=Tubp->vomi[p]=Tubp->a[p]=Tubp->v[p]=0;
	#endif
	#if Test_Type=='E' || Test_Type=='C'
						if(m==N-1)Tubp->aomi[p]=Tubp->vomi[p]=Tubp->a[p]=Tubp->v[p]=0;
	#elif Test_Type=='B' || Test_Type=='L'
						if(m==N-1 && p!=0)Tubp->a[p]=Tubp->v[p]=0;
	#elif Test_Type=='F'
						if(m==N-1 && p!=2)Tubp->a[p]=Tubp->v[p]=0;
	#endif
#else
	#if Test_Type && Test_Type!='K'  && Test_Type!='V'
						if(m<13)Tubp->aomi[p]=Tubp->vomi[p]=Tubp->a[p]=Tubp->v[p]=0;
	#endif
	#if Test_Type=='E'
						if(m>=N-13)Tubp->aomi[p]=Tubp->vomi[p]=Tubp->a[p]=Tubp->v[p]=0;
	#elif Test_Type=='B' || Test_Type=='L'
						if(m>=N-13 && p!=0)Tubp->a[p]=Tubp->v[p]=0;
	#elif Test_Type=='F'
						if(m>=N-13 && p!=2)Tubp->a[p]=Tubp->v[p]=0;
	#elif Test_Type=='Y'
						if(m>=N-13)Tubp->a[p]=Tubp->v[p]=0;
	#elif Test_Type=='K'  || Test_Type=='V'
						if(m<13 && p!=0)Tubp->aomi[p]=Tubp->vomi[p]=Tubp->a[p]=Tubp->v[p]=0;
						if(m>=N-13 && p!=0)Tubp->a[p]=Tubp->v[p]=0;
	#endif
#endif
#if TipEnd
//						if(MoleDynamicNo<10 && m>=(N/FilaSet-TipEnd)*FilaSet)Tubp->a[p]=Tubp->v[p]=0;
//						if(MoleDynamicNo<2)Tubp->a[p]=Tubp->v[p]=0;
#endif
						Tubp->rpre[p]=Tubp->r[p];
#if FC_ON
	#if FR_ON
						Tubp->v[p]=C0*Tubp->v[p]+(C2*Tubp->a[p]+(C1-C2)*Tubp->apre[p])*dt+sigmav*Rand_N[4*p+1+12*Tubp->Tub_No]*(Crv*Rand_N[4*p+12*Tubp->Tub_No]+SqCrv)/M_F;
						Tubp->r[p]=Tubp->r[p]+C1*Tubp->v[p]*dt+C2*Tubp->a[p]*dt*dt+sigmar*Rand_N[4*p+12*Tubp->Tub_No]/M_F;
	#else
						Tubp->v[p]=C0*Tubp->v[p]+(C2*Tubp->a[p]+(C1-C2)*Tubp->apre[p])*dt;
						Tubp->r[p]=Tubp->r[p]+C1*Tubp->v[p]*dt+C2*Tubp->a[p]*dt*dt;
	#endif
#else											
						Tubp->v[p]=Tubp->v[p]+(Tubp->a[p]+Tubp->apre[p])*dt/2;	
						Tubp->r[p]=Tubp->r[p]+Tubp->v[p]*dt+Tubp->a[p]*dt*dt/2;
#endif
						Tubp->apre[p]=Tubp->a[p];
						Tubp->a[p]=0;
#if FC_ON						
	#if MR_ON
						Tubp->vomi[p]=Co0*Tubp->vomi[p]+(Co2*Tubp->aomi[p]+(Co1-Co2)*Tubp->aomipre[p])*dt+sigmav*Rand_N[4*p+3+12*Tubp->Tub_No]*(Crv*Rand_N[4*p+2+12*Tubp->Tub_No]+SqCrv)/I_M;
						dangle[p]=Co1*Tubp->vomi[p]*dt+Co2*Tubp->aomi[p]*dt*dt+sigmar*Rand_N[4*p+2+12*Tubp->Tub_No]/I_M;
	#else
						Tubp->vomi[p]=Co0*Tubp->vomi[p]+(Co2*Tubp->aomi[p]+(Co1-Co2)*Tubp->aomipre[p])*dt;
						dangle[p]=Co1*Tubp->vomi[p]*dt+Co2*Tubp->aomi[p]*dt*dt;
	#endif
#else						
						Tubp->vomi[p]=Tubp->vomi[p]+(Tubp->aomi[p]+Tubp->aomipre[p])*dt/2;
						dangle[p]=Tubp->vomi[p]*dt+Tubp->aomi[p]*dt*dt/2;
#endif
						dangle[3]+=dangle[p]*dangle[p]; 
						Tubp->aomipre[p]=Tubp->aomi[p];
						Tubp->aomi[p]=Tubp->mom[p]=0;

					}
					PeriodBoundary(Tubp->r);
					dangle[3]=sqrt(dangle[3]);
					if(dangle[3]<1e-100&&dangle[3]>-1e-100){
						dangle[1]=dangle[2]=0;
						dangle[0]=1;
					}
					else{
						for(p=0;p<3;p++)dangle[p]=dangle[p]/dangle[3];
					}
					dangle[3]=dangle[3]/2;
					QuaTimes(dangle,Tubp->angle,angle);

					if(Tubp->Type)xi=xi1;
					else xi=xi0;
					Set_Tub(Tubp,Tubp->r,angle,xi);
#if EOut_ON
	#if Di_Inertia
					x=Tubp->angle[0];
					y=Tubp->angle[1];
					z=Tubp->angle[2];
					costheta=cos(2*Tubp->angle[3]);
					sintheta=sin(-2*Tubp->angle[3]);

					Tra[0]=costheta+x*x*(1-costheta);
					Tra[1]=sintheta*z+x*y*(1-costheta);
					Tra[2]=-sintheta*y+x*z*(1-costheta);
					Tra[3]=-sintheta*z+x*y*(1-costheta);
					Tra[4]=costheta+y*y*(1-costheta);
					Tra[5]=sintheta*x+y*z*(1-costheta);
					Tra[6]=sintheta*y+x*z*(1-costheta);
					Tra[7]=-sintheta*x+y*z*(1-costheta);
				    Tra[8]=costheta+z*z*(1-costheta);
					Tubp->vomibar[0]=(Tra[0]*Tubp->vomi[0]+Tra[3]*Tubp->vomi[1]+Tra[6]*Tubp->vomi[2]);
				    Tubp->vomibar[1]=(Tra[1]*Tubp->vomi[0]+Tra[4]*Tubp->vomi[1]+Tra[7]*Tubp->vomi[2]);
				    Tubp->vomibar[2]=(Tra[2]*Tubp->vomi[0]+Tra[5]*Tubp->vomi[1]+Tra[8]*Tubp->vomi[2]);
				    ERot[Trec]+=Tubp->vomibar[0]*Tubp->vomibar[0]+3.5*Tubp->vomibar[1]*Tubp->vomibar[1]+3.5*Tubp->vomibar[2]*Tubp->vomibar[2];
	#else
				    ERot[Trec]+=Tubp->vomi[0]*Tubp->vomi[0]+Tubp->vomi[1]*Tubp->vomi[1]+Tubp->vomi[2]*Tubp->vomi[2];
	#endif
				    EVel[Trec]+=Tubp->v[0]*Tubp->v[0]+Tubp->v[1]*Tubp->v[1]+Tubp->v[2]*Tubp->v[2];	
					FAvg[Trec]+=Tubp->apre[0]*Tubp->apre[0]+Tubp->apre[1]*Tubp->apre[1]+Tubp->apre[2]*Tubp->apre[2];
					EPoten[Trec]+=Tubp->poten+Tubp->potent;
#endif
					if((Tubp->apre[0]<-1e6||Tubp->apre[0]>1e6) || (Tubp->apre[1]<-1e6||Tubp->apre[1]>1e6) || (Tubp->apre[2]<-1e6||Tubp->apre[2]>1e6)){
						fprintf(stderr,"Abort (apre)\n");
						StdEQuit(Tubp);
						MPI_Abort(MPI_COMM_WORLD,1);
						exit(1);
					}
					if((Tubp->angle[0]*Tubp->angle[0]+Tubp->angle[1]*Tubp->angle[1]+Tubp->angle[2]*Tubp->angle[2]>50) || (Tubp->angle[3]*Tubp->angle[3]>50)){
						fprintf(stderr,"Abort (anglemod2)\n");
						StdEQuit(Tubp);
						MPI_Abort(MPI_COMM_WORLD,1);
						exit(1);
					}
					Tub[m]=tubtmp;
				}
//			}

			if(numprocs!=1){
				MPI_Barrier(MPI_COMM_WORLD);
				if(myid!=0){
					if(myid!=numprocs-1)BcastSize=MpiSec;
					else BcastSize=N-NStart-myid*MpiSec;
					MPI_Send(&Tub[NStart+myid*MpiSec],BcastSize*TubSize,MPI_CHAR,0,myid,MPI_COMM_WORLD);
#if EOut_ON
					if(!(TTot%Out_Step)){
						EFTrans[0]=EPoten[Trec];
						EFTrans[1]=EVel[Trec];
						EFTrans[2]=ERot[Trec];
						EFTrans[3]=ETot[Trec];
						EFTrans[4]=FAvg[Trec];
						EFTrans[5]=BigBall.BForce[0];
						EFTrans[6]=BigBall.BForce[1];
						EFTrans[7]=BigBall.BForce[2];
						MPI_Send(EFTrans,9,MPI_DOUBLE,0,myid,MPI_COMM_WORLD);
					}
#endif
				}
				else{
					for(i=1;i<numprocs;i++){
						if(i!=numprocs-1)BcastSize=MpiSec;
						else BcastSize=N-NStart-i*MpiSec;
						MPI_Recv(&Tub[NStart+i*MpiSec],BcastSize*TubSize,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
#if EOut_ON
						if(!(TTot%Out_Step)){
							MPI_Recv(EFTrans,9,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
							EPoten[Trec]+=EFTrans[0];
							EVel[Trec]+=EFTrans[0];
							ERot[Trec]+=EFTrans[0];
							ETot[Trec]+=EFTrans[0];
							FAvg[Trec]+=EFTrans[0];
							BigBall.BForce[0]+=EFTrans[0];
							BigBall.BForce[1]+=EFTrans[0];
							BigBall.BForce[2]+=EFTrans[0];
						}
#endif
					}
				}
			}	
#if EOut_ON
			if(myid==0 && !(TTot%Out_Step)){
	#if Test_Type=='Y'
				FBall[Trec]=sqrt(BigBall.BForce[0]*BigBall.BForce[0]+BigBall.BForce[1]*BigBall.BForce[1]+BigBall.BForce[2]*BigBall.BForce[2]);
	#endif
				EPoten[Trec]=EPoten[Trec]/N/Out_Step;
				EVel[Trec]=M_F*EVel[Trec]/2/N/Out_Step;
				ERot[Trec]=I_M*ERot[Trec]/2/N/Out_Step;
				ETot[Trec]=EPoten[Trec]+EVel[Trec]+ERot[Trec];
				FAvg[Trec]=sqrt(FAvg[Trec]/N/Out_Step);
	#if Test_Type=='Y'
				fprintf(fphase,"%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n",Trec*Out_Step,EPoten[Trec],EVel[Trec],ERot[Trec],ETot[Trec],FAvg[Trec],FBall[Trec]);
	#else
				fprintf(fphase,"%d,%.16g,%.16g,%.16g,%.16g,%.16g\n",Trec*Out_Step,EPoten[Trec],EVel[Trec],ERot[Trec],ETot[Trec],FAvg[Trec]);
	#endif
				if(!MConver && fabs(ETot[Trec-1]-ETot[Trec])<E_Conver)TConver++;
				else if(!MConver)TConver=0;
				if(!MConver && TConver>=T_Conver){
					TConver=TTot;
					MConver=1;
					StopMsg=1;
				}
					Trec++;
			}
#endif
			BigBall.BForce[0]=BigBall.BForce[1]=BigBall.BForce[2]=0;		
			if(t%SlowTime==0){
				for(i=0;i<Nbox;i++)Head[i]=-1;
				for(i=0;i<N;i++){
					Tubp=&Tub[i];
					if(!Tubp->Active)continue;
					icell=GetCellNo(Tubp->r);
					if(icell<0 || icell>Nbox){
						fprintf(stderr,"%d Tub %d icell= %d! (NBox %d)\n",myid,i,icell,Nbox);
						fprintf(stderr,"r0: %g r1: %g r2: %g\n",Tubp->r[0],Tubp->r[1],Tubp->r[2]);
					}
					List[i]=Head[icell];
					Head[icell]=i;
					Tubp->Tub_Id=i/MpiSec;
					if(myid==0){
#if Whole_Struc && BREAK_ON
						if(Tubp->Type==1)
							{
							double Ymax=0,Ymin=0,Zmax=0,Zmin=0,YZ;
							ra=(rand()/(RAND_MAX+0.1));
							if(Test_Type=='W' && MaxEnd<17)VirtualL=17-MaxEnd;
							else if(!Tubp->Type)VirtualL=(Tubp->ChainNo<=5)?-Tubp->ChainNo:-5;
							else VirtualL=0;
							if(Tubp->Tub_No>25 && Tubp->WallMark<-SlowTime && (Tubp->left<0||Tubp->right<0) && exp(-B_rate*(Tubp->ChainNo+VirtualL)*(Tubp->ChainNo+VirtualL)*(Tubp->poten+1.184)/225)<ra){
								fprintf(stderr,"\nBreak:\n%d",i);
								for(j=i,p=0;j>=0;j=Tub[j].front,p++){
									fprintf(stderr,"\t%d",j);
									if(!p){
										Ymax=Ymax=Tubp->r[1];
										Zmax=Zmax=Tubp->r[2];
									}
									else{
										if(Tub[j].r[1]>Ymax)Ymax=Tub[j].r[1];
										else if(Tub[j].r[1]<Ymin)Ymin=Tub[j].r[1];
										if(Tub[j].r[2]>Zmax)Zmax=Tub[j].r[2];
										else if(Tub[j].r[2]<Zmin)Zmin=Tub[j].r[2];
									}
								}
								YZ=sqrt((Ymax-Ymin)*(Ymax-Ymin)+(Zmax-Zmin)*(Zmax-Zmin));
								if(YZ/Tubp->ChainNo>0.3){
									fprintf(fout,"\nBreak:\n%d",i);
									Top[Tubp->CMark]=Tubp->back;
									j=Tubp->back;
									if(j!=-1)Tub[j].front=-1;
									for(p=1;j>=0;j=Tub[j].back,p++){
										Tub[j].ChainNo=p;
									}
									for(j=i;j>=0;j=Tub[j].front){
										fprintf(stderr,"\t%d",j);
										fprintf(fout,"\t%d",j);
										Tub[j].Active=0;
										if(Tub[j].back>=0){
											Tub[Tub[j].back].front=-1;
											Tub[j].back=-1;
										}
										if(Tub[j].left>=0){
											Tub[Tub[j].left].right=-1;
											Tub[j].left=-1;
										}
										if(Tub[j].right>=0){
											Tub[Tub[j].right].left=-1;
											Tub[j].right=-1;
										}
									}
									fprintf(stderr,"\tBreak done.\n");
								}
								else fprintf(stderr,"\tBreak cancel.\n");
							}
						}
#endif
#if Whole_Struc && (GROW_ON || BREAK_ON)
	#if HD_ON
						if(Tubp->Type!=1 && i>25 && Tubp->Tin>0){
		#if HD_ON==1
							if(Tubp->Thd<Tubp->Tin)Tubp->Type=1;
		#elif HD_ON==2
							ActiveE=0;
							if(Tubp->front>=0)ActiveE+=0.161+Tubp->potenf;
							if(Tubp->back>=0)ActiveE+=0.161+Tubp->potenb;
							Cnow=7.2e11*exp(18.5*ActiveE);
							Cnow=exp(-Cnow*SlowTime);
							Tubp->Phd-=Tubp->Puhd*(1-Cnow);
							if(Tubp->Phd<0)Tubp->Type=1;
							else Tubp->Puhd*=Cnow;
		#endif
							if(Tubp->RFMark)Tubp->Type=0;
//							if(exp(-C_rate*(Tubp->Tin-CTime))<ra)Tubp->Type=1;
		#if Test_Type=='W'
							if(i<130)Tubp->Type=0;
		#endif
						}
	#endif
						else if(Tubp->Type==1 && Top[Tubp->CMark]==i){
							ra=rand()/(RAND_MAX+0.1);
							if(exp(-C_rate*(Tubp->Ttop-CTime))<ra){
								fprintf(fout,"\nBreak Single:\n%d",i);
								Top[Tubp->CMark]=Tubp->back;
								j=Tubp->back;
								if(j!=-1)Tub[j].front=-1;
								for(p=1;j>=0;j=Tub[j].back,p++){
									Tub[j].ChainNo=p;
								}
								for(j=i;j>=0;j=Tub[j].front){
									fprintf(stderr,"\t%d",j);
									fprintf(fout,"\t%d",j);
									Tub[j].Active=0;
									if(Tub[j].back>=0){
										Tub[Tub[j].back].front=-1;
										Tub[j].back=-1;
									}
									if(Tub[j].left>=0){
										Tub[Tub[j].left].right=-1;
										Tub[j].left=-1;
									}
									if(Tub[j].right>=0){
										Tub[Tub[j].right].left=-1;
										Tub[j].right=-1;
									}
								}
								fprintf(stderr,"\tBreak done.\n");
							}
						}
#elif LongRFTest
	#if HD_ON
						if(Tubp->Type!=1 && i>25 && Tubp->Tin>0){
							ra=rand()/(RAND_MAX+0.1);
		#if HD_ON==1
							if(Tubp->Thd<Tubp->Tin)Tubp->Type=1;
		#elif HD_ON==2
							ActiveE=0;
							if(Tubp->front>=0)ActiveE+=0.161+Tubp->potenf;
							if(Tubp->back>=0)ActiveE+=0.161+Tubp->potenb;
							Cnow=Afree*exp(18.5*ActiveE);
							Cnow=exp(-Cnow*SlowTime);
							Tubp->Phd-=Tubp->Puhd*(1-Cnow);
							if(Tubp->Phd<0)Tubp->Type=1;
							else Tubp->Puhd*=Cnow;
		#endif
							if(Tubp->RFMark)Tubp->Type=0;
						}
	#endif
#endif
					}
					Tubp->Tin+=SlowTime;
				}
				for(i=0;i<13;i++){
					Tub[Top[i]].Ttop+=SlowTime;
				}
			}
			if(myid==0){
					#if Detach_Curve==1
					if(t%(500*SlowTime)==0)DC1_Out();
					#endif
				if(t%SlowTime==0){
					for(i=0;i<N;i++){
						Tubp=&Tub[i];
						if(!Tubp->Active)continue;
						anglemod=Tubp->angle[0]*Tubp->angle[0]+Tubp->angle[1]*Tubp->angle[1]+Tubp->angle[2]*Tubp->angle[2];
						anglemod=sqrt(anglemod);
						Tubp->angle[0]=Tubp->angle[0]/anglemod;
						Tubp->angle[1]=Tubp->angle[1]/anglemod;
						Tubp->angle[2]=Tubp->angle[2]/anglemod;
						if(Tubp->angle[3]>pi)Tubp->angle[3]=Tubp->angle[3]-2*pi*((int)(Tubp->angle[3]/pi/2));
						else if(Tubp->angle[3]<-pi)Tubp->angle[3]=Tubp->angle[3]+2*pi*((int)(Tubp->angle[3]/pi/2));
						Set_Tub(Tubp,Tubp->r,Tubp->angle,Tubp->xi);
					}
				}
#if GROW_ON
				if(Ngrow>=TGrow && t%SlowTime==0){
					Ngrow=0;
					fprintf(stderr,"t: %d \n",t);
	#if Test_Type=='W'
					if(MaxEnd<40)AddTub();
	#else
					AddTub();
	#endif
	#if GRStep
					TGrow=TGrow+(int)(GRStep*(2.0*rand()/(RAND_MAX+0.1)-1));
					if(TGrow<1000)TGrow=1000;
	#endif
				}
				Ngrow++;
#endif
			}
			if(numprocs!=1){
				if(t%SlowTime==0)
				MPI_Barrier(MPI_COMM_WORLD);
					{
					Npre=N;
					MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
					if(N>Npre){
						Tub=(tubulin *)realloc(Tub,N*(sizeof(tubulin)));
						List=(int *)realloc(List,N*(sizeof(int)));
					}
				MPI_Bcast(Tub,N*TubSize,MPI_CHAR,0,MPI_COMM_WORLD);
				}
				if(t%SlowTime==0 && N>Npre){
					icell=GetCellNo(Tub[N-1].r);
					List[N-1]=Head[icell];
					Head[icell]=N-1;
#if Whole_Struc	&& LCompute
					for(i=0;i<13;i++){
						if(Maxline<Tub[Top[i]].LineNo)Maxline=Tub[Top[i]].LineNo;
					}
					if(Maxline>LCompute){
						LStart=Maxline-LCompute;
						NStart=LStart*13;
						mincom=N-NStart;
					}
					else {
						mincom=N;
						NStart=0;
					}
					if(mincom%numprocs)MpiSec=mincom/numprocs+1;
					else MpiSec=mincom/numprocs;
#else
					if(N%numprocs)MpiSec=N/numprocs+1;
					else MpiSec=N/numprocs;
#endif
				}
			}
		}
	}
}
#endif

void output(int M){
	int i,j,k,p,Nst1=1,Nst2=1;
	char filename[MaxFName],Stmp[MaxFName],Ctype,Stype[6];
 	time_t tnow;
	double zero=0.0,BoxLout,ttheta=0,ntheta=0;
	double BoxLoutX,BoxLoutY,BoxLoutZ;
	tubulin *Tubp;
	FILE *fdata,*fp,*fvmd,*frs,*frf;

	BoxLoutX=BoxLX;
	BoxLoutY=BoxLY;
	BoxLoutZ=BoxLZ;
	if(myid==0){
		if(N>=10000){
			fprintf(stderr,"Abort (N>9999)\n");
			MPI_Abort(MPI_COMM_WORLD,1);
			exit(1);
		}

		sprintf(filename,"%.4d.pdb",M+1);
		fp=fopen(filename,"wt+");
		if(fp==NULL){
			fprintf(stderr,"Cannot open PDB file!");
			fprintf(stderr, "open file :%s\n", strerror(errno));
			fprintf(stderr,"Trying again...");
			fp=fopen(filename,"wt+");
			if(fp==NULL){
				fprintf(stderr,"Cannot open PDB file!");
				fprintf(stderr, "open file :%s\n", strerror(errno));
				exit(EXIT_FAILURE);
			}
		}
		sprintf(filename,"%.4d.data",M+1);
		fdata=fopen(filename,"wt+");
		if(fdata==NULL){
			fprintf(stderr,"Cannot open DATA file!");
			fprintf(stderr, "open file :%s\n", strerror(errno));
			fprintf(stderr,"Trying again...");
			fdata=fopen(filename,"wt+");
			if(fdata==NULL){
				fprintf(stderr,"Cannot open DATA file!");
				fprintf(stderr, "open file :%s\n", strerror(errno));
				exit(EXIT_FAILURE);
			}
		}
		
		sprintf(filename,"vMT_%s.pdb",fstring);
		if(M==0)fvmd=fopen(filename,"wt+");
		else fvmd=fopen(filename,"at+");
		if(fvmd==NULL){
			fprintf(stderr,"Cannot open VMD file!");
			fprintf(stderr, "open file :%s\n", strerror(errno));
			fprintf(stderr,"Trying again...");
			if(M==0)fvmd=fopen(filename,"wt+");
			else fvmd=fopen(filename,"at+");
			if(fvmd==NULL){
				fprintf(stderr,"Cannot open VMD file!");
				fprintf(stderr, "open file :%s\n", strerror(errno));
				exit(EXIT_FAILURE);
			}
		}
		
#if RFDensityOut
		sprintf(filename,"RFDensity_%s.csv",fstring);
		if(M==0)frf=fopen(filename,"wt+");
		else frf=fopen(filename,"at+");
		if(frf==NULL){
			fprintf(stderr,"Cannot open RFDensity file!");
			fprintf(stderr, "open file :%s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
#endif
		Nactive=0;
		Ntype=0;
		fprintf(fdata,"   No  Active  Type     r0      r1      r2    Angle0  Angle1  Angle2  Angle3  Angle4    xi     front   back    left    right   lseam    rseam   CMark  ChainNo  Tin  LineNo  Potent\n");

#if BoxOut & 0x01
		fprintf(fp,"ATOM  %5d  K   LEU N   1    %8.3f%8.3f%8.3f  1.00 %6.3f\n",1,-BoxLoutX,-BoxLoutY,-BoxLoutZ,zero);
		fprintf(fp,"ATOM  %5d  K   LEU N   2    %8.3f%8.3f%8.3f  1.00 %6.3f\n",2,-BoxLoutX,BoxLoutY,-BoxLoutZ,zero);
		fprintf(fp,"ATOM  %5d  K   LEU N   3    %8.3f%8.3f%8.3f  1.00 %6.3f\n",3,BoxLoutX,-BoxLoutY,-BoxLoutZ,zero);
		fprintf(fp,"ATOM  %5d  K   LEU N   4    %8.3f%8.3f%8.3f  1.00 %6.3f\n",4,BoxLoutX,BoxLoutY,-BoxLoutZ,zero);
		fprintf(fp,"ATOM  %5d  K   LEU N   5    %8.3f%8.3f%8.3f  1.00 %6.3f\n",5,-BoxLoutX,-BoxLoutY,BoxLoutZ,zero);
		fprintf(fp,"ATOM  %5d  K   LEU N   6    %8.3f%8.3f%8.3f  1.00 %6.3f\n",6,-BoxLoutX,BoxLoutY,BoxLoutZ,zero);
		fprintf(fp,"ATOM  %5d  K   LEU N   7    %8.3f%8.3f%8.3f  1.00 %6.3f\n",7,BoxLoutX,-BoxLoutY,BoxLoutZ,zero);
		fprintf(fp,"ATOM  %5d  K   LEU N   8    %8.3f%8.3f%8.3f  1.00 %6.3f\n",8,BoxLoutX,BoxLoutY,BoxLoutZ,zero);
		Nst1=9;
#endif
#if BoxOut & 0x02
		fprintf(fvmd,"ATOM  %5d  K   LEU N   1    %8.3f%8.3f%8.3f  1.00 %6.3f\n",1,-BoxLoutX,-BoxLoutY,-BoxLoutZ,zero);
		fprintf(fvmd,"ATOM  %5d  K   LEU N   2    %8.3f%8.3f%8.3f  1.00 %6.3f\n",2,-BoxLoutX,BoxLoutY,-BoxLoutZ,zero);
		fprintf(fvmd,"ATOM  %5d  K   LEU N   3    %8.3f%8.3f%8.3f  1.00 %6.3f\n",3,BoxLoutX,-BoxLoutY,-BoxLoutZ,zero);
		fprintf(fvmd,"ATOM  %5d  K   LEU N   4    %8.3f%8.3f%8.3f  1.00 %6.3f\n",4,BoxLoutX,BoxLoutY,-BoxLoutZ,zero);
		fprintf(fvmd,"ATOM  %5d  K   LEU N   5    %8.3f%8.3f%8.3f  1.00 %6.3f\n",5,-BoxLoutX,-BoxLoutY,BoxLoutZ,zero);
		fprintf(fvmd,"ATOM  %5d  K   LEU N   6    %8.3f%8.3f%8.3f  1.00 %6.3f\n",6,-BoxLoutX,BoxLoutY,BoxLoutZ,zero);
		fprintf(fvmd,"ATOM  %5d  K   LEU N   7    %8.3f%8.3f%8.3f  1.00 %6.3f\n",7,BoxLoutX,-BoxLoutY,BoxLoutZ,zero);
		fprintf(fvmd,"ATOM  %5d  K   LEU N   8    %8.3f%8.3f%8.3f  1.00 %6.3f\n",8,BoxLoutX,BoxLoutY,BoxLoutZ,zero);
		Nst2=9;
#endif
		strcpy(Stype,"TUB");
		for(j=0;j<N;j++){
			Tubp=&Tub[j];
			if(Tubp->Type)Ntype++;
			if(Tubp->Active){
				Nactive++;
				#if FixLink
				Tubp->PotenA=0;
				Tubp->PotenB=0;
				for(k=0;k<Tubp->AddNum;k++){
					switch(Tubp->PPMark[k] & 0x70){
						case 0x10:
						case 0x30:
						case 0x40:Tubp->PotenB+=Tubp->PotenList[k];break;
						case 0x20:
						case 0x50:
						case 0x60:Tubp->PotenA+=Tubp->PotenList[k];break;
					}
				}
				Tubp->PotenA*=0.5;
				Tubp->PotenB*=0.5;
				#endif
				switch(Tubp->Type){
					case 0:Ctype='N';strcpy(Stype,"LEU");break;
					case 1:Ctype='B';strcpy(Stype,"GLY");break;
					case 2:Ctype='C';break;
					case 3:Ctype='O';break;
					default:Ctype='P';
				}
				fprintf(fp,"ATOM  %5d  %c   %3s %c%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					2*j+Nst1,Ctype,Stype,Tubp->Type?'D':'T',j,(Tubp->r1[0]),(Tubp->r1[1]),(Tubp->r1[2]),Tubp->PotenB);
				fprintf(fp,"ATOM  %5d  %c   %3s %c%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					2*j+Nst1+1,Ctype+1,"PHE",'A',j,(Tubp->r2[0]),(Tubp->r2[1]),(Tubp->r2[2]),Tubp->PotenA);
				fprintf(fvmd,"ATOM  %5d  %c   %3s %c%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,Ctype,Stype,Tubp->Type?'D':'T',j,(Tubp->r1[0]),(Tubp->r1[1]),(Tubp->r1[2]),Tubp->PotenB);
				fprintf(fvmd,"ATOM  %5d  %c   %3s %c%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,Ctype+1,"PHE",'A',j,(Tubp->r2[0]),(Tubp->r2[1]),(Tubp->r2[2]),Tubp->PotenA);
			#if BoxOut & 0x10
					//Lateral
				fprintf(fvmd,"ATOM  %5d  S1  %3s S%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r1[0]+Tubp->side1[0]),(Tubp->r1[1]+Tubp->side1[1]),(Tubp->r1[2]+Tubp->side1[2]),Tubp->PotenB);
				fprintf(fvmd,"ATOM  %5d  S2  %3s S%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r1[0]+Tubp->side2[0]),(Tubp->r1[1]+Tubp->side2[1]),(Tubp->r1[2]+Tubp->side2[2]),Tubp->PotenB);
				fprintf(fvmd,"ATOM  %5d  S3  %3s S%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r2[0]+Tubp->side3[0]),(Tubp->r2[1]+Tubp->side3[1]),(Tubp->r2[2]+Tubp->side3[2]),Tubp->PotenA);
				fprintf(fvmd,"ATOM  %5d  S4  %3s S%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r2[0]+Tubp->side4[0]),(Tubp->r2[1]+Tubp->side4[1]),(Tubp->r2[2]+Tubp->side4[2]),Tubp->PotenA);
					//Longitudinal
				fprintf(fvmd,"ATOM  %5d  P1  %3s C%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r1[0]+Tubp->chain1[0]),(Tubp->r1[1]+Tubp->chain1[1]),(Tubp->r1[2]+Tubp->chain1[2]),Tubp->PotenB);
				fprintf(fvmd,"ATOM  %5d  P2  %3s C%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r2[0]+Tubp->chain2[0]),(Tubp->r2[1]+Tubp->chain2[1]),(Tubp->r2[2]+Tubp->chain2[2]),Tubp->PotenA);
			#endif
			#if BoxOut & 0x20
					//n
				fprintf(fvmd,"ATOM  %5d  H1  %3s N%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r[0]+0.6*Tubp->n[0]),(Tubp->r[1]+0.6*Tubp->n[1]),(Tubp->r[2]+0.6*Tubp->n[2]),Tubp->PotenB);
					//s
				fprintf(fvmd,"ATOM  %5d  H2  %3s B%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r[0]+0.6*Tubp->s[0]),(Tubp->r[1]+0.6*Tubp->s[1]),(Tubp->r[2]+0.6*Tubp->s[2]),Tubp->PotenB);
					//t
				fprintf(fvmd,"ATOM  %5d  H3  %3s T%4d    %8.3f%8.3f%8.3f 1.00 %7.4f\n",
					Nst2++,"TUB",j,(Tubp->r[0]+0.6*Tubp->t[0]),(Tubp->r[1]+0.6*Tubp->t[1]),(Tubp->r[2]+0.6*Tubp->t[2]),Tubp->PotenB);
			#endif
			}
			fprintf(fdata,"%4d%8d%6d  %8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8.4f\n",
				j,Tubp->Active,Tubp->Type,Tubp->r[0],Tubp->r[1],Tubp->r[2],Tubp->angle[0],Tubp->angle[1],Tubp->angle[2],Tubp->angle[3],Tubp->angle[4],
				Tubp->xi,Tubp->front,Tubp->back,Tubp->left,Tubp->right,Tubp->lseam,Tubp->rseam,Tubp->CMark,Tubp->ChainNo,Tubp->Tin,Tubp->LineNo,Tubp->potent);
			if(j%100==99)
				fprintf(fdata,"\n   No  Active  Type     r0      r1      r2    Angle0  Angle1  Angle2  Angle3  Angle4    xi     front   back    left    right   lseam    rseam   CMark  ChainNo  Tin  LineNo   Potent\n");
		}

#if BoxOut & 0x01
	fprintf(fp,"CONECT %4d%5d\n",1,2);
	fprintf(fp,"CONECT %4d%5d\n",1,3);
	fprintf(fp,"CONECT %4d%5d\n",2,4);
	fprintf(fp,"CONECT %4d%5d\n",3,4);
	fprintf(fp,"CONECT %4d%5d\n",5,6);
	fprintf(fp,"CONECT %4d%5d\n",5,7);
	fprintf(fp,"CONECT %4d%5d\n",6,8);
	fprintf(fp,"CONECT %4d%5d\n",7,8);
	fprintf(fp,"CONECT %4d%5d\n",1,5);
	fprintf(fp,"CONECT %4d%5d\n",2,6);
	fprintf(fp,"CONECT %4d%5d\n",3,7);
	fprintf(fp,"CONECT %4d%5d\n",4,8);
#endif
	Nst2=1;
#if BoxOut & 0x02
	fprintf(fvmd,"CONECT %4d%5d\n",1,2);
	fprintf(fvmd,"CONECT %4d%5d\n",1,3);
	fprintf(fvmd,"CONECT %4d%5d\n",2,4);
	fprintf(fvmd,"CONECT %4d%5d\n",3,4);
	fprintf(fvmd,"CONECT %4d%5d\n",5,6);
	fprintf(fvmd,"CONECT %4d%5d\n",5,7);
	fprintf(fvmd,"CONECT %4d%5d\n",6,8);
	fprintf(fvmd,"CONECT %4d%5d\n",7,8);
	fprintf(fvmd,"CONECT %4d%5d\n",1,5);
	fprintf(fvmd,"CONECT %4d%5d\n",2,6);
	fprintf(fvmd,"CONECT %4d%5d\n",3,7);
	fprintf(fvmd,"CONECT %4d%5d\n",4,8);
	Nst2=9;
#endif
	for(j=0;j<N;j++){
		if(Tub[j].Active){
		#if BoxOut & 0x04
			fprintf(fp,"CONECT %4d%5d\n",2*j+Nst1,2*j+Nst1+1);
		#endif
		#if BoxOut & 0x08
			fprintf(fvmd,"CONECT %4d%5d\n",Nst2,Nst2+1);
		#endif
			Nst2+=2;
			if(BoxOut & 0x10)Nst2+=6;
			if(BoxOut & 0x20)Nst2+=3;
		}
	}
//	ATOM      1   C  LEU Y 129     -887.19  -32.24  -48.15  1.00  0.5184778284491873
		fprintf(fp,"END\n");
		fclose(fp);
		fprintf(fvmd,"END\n");
		fclose(fvmd);
		fprintf(fdata,"END\n");
		fclose(fdata);
		
		if(RS_OUT_ON && M%10==0){
			tnow=time(0);
			strftime(Stmp,sizeof(Stmp),"%Y%m%d%H%M%S",localtime(&tnow));
			sprintf(filename,"RestrartSourse%s.tubulin",Stmp);
			frs=fopen(filename,"w");
			fwrite(&N,sizeof(int),1,frs);
			fwrite(Tub,sizeof(tubulin),N,frs);
			fclose(frs);
		}
	}
#if DataOut && !(BoxOut & 0xF0)
	{
		FILE *fpp;
		double BoxPoten;
		sprintf(filename,"PotenData_%s.txt",tstring);
		if(M==0)fpp=fopen(filename,"wt+");
		else fpp=fopen(filename,"at+");
		if(fpp==NULL){
			fprintf(stderr,"Cannot open PotenData file!");
			fprintf(stderr, "open file :%s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		if(BoxOut & 0x02){
			BoxPoten=Tub[0].PotenB;
			for(j=1;j<=8;j++)fprintf(fpp,"%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n",M,j,BoxPoten,BoxPoten,BoxPoten,BoxPoten,BoxPoten,BoxPoten);
			Nst2=9;
		}
		else Nst2=1;
		for(j=0;j<N;j++){
			Tubp=&Tub[j];
			if(Tubp->Active){
				fprintf(fpp,"%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n",M,Nst2++,Tubp->MyPoten[0],Tubp->MyPoten[1],Tubp->MyPoten[2],Tubp->MyPoten[3],Tubp->MyPoten[4],Tubp->MyPoten[5]);
				fprintf(fpp,"%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n",M,Nst2++,Tubp->MyPoten[6],Tubp->MyPoten[7],Tubp->MyPoten[8],Tubp->MyPoten[9],Tubp->MyPoten[10],Tubp->MyPoten[11]);
			}
		}
		fclose(fpp);
	}
#endif
	if(DeltaL!=0.0 || DeltaXi!=0.0 || DeltaA!=0.0){
		int Nall=0,N1=0,N13=0;
		double PotenAll=0,Poten1=0,Poten13=0,Dtmp,Ptest=0;
		sprintf(filename,"MyAvgPoten_%s.txt",tstring);
		if(!M)fp=fopen(filename,"wt");
		else fp=fopen(filename,"at");
		for(j=0;j<N;j++){
			if(Tub[j].Active){
				Dtmp=Tub[j].poten;
				Nall++;
				PotenAll+=Dtmp;
				if(j && j<N-1 && j!=Lat_Length-1 && j!=N-Lat_Length){
					N1++;
					Poten1+=Dtmp;
					#if Whole_Struc
					if(j>12 && j<N-13)
					#else
					if(j%Lat_Length && (j+1)%Lat_Length && j/Lat_Length && (N-1-j)/Lat_Length)
					#endif
					{
						N13++;
						Poten13+=Dtmp;
					}
				}
			}
		}
		fprintf(fp,"%d\t%g\t%g\t%g\t%g\n",M,PotenAll/Nall,Poten1/N1,Poten13/N13,Ptest/N1);
		//Averange Energy Per: Dimer | Dimer except 1&N | Dimer between 13 and N-13 | twist energy
		fclose(fp);
	}
}

int main(int argc,char **argv){	
	int i,j,k;
	double Dtmp;
	double start,end,totaltime;
	char fname[MaxFName],SRestart[6]="-r";
	tubulin *Tubp;
//MPI初始化  
	MPI_Init(&argc,&argv);

	start=clock()*1.0;
//	fprintf(stderr,"MPIinit done.\n");
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);/*得到总的进程数*/
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);/*得到当前进程号*/
	MPI_Get_processor_name(processor_name,&namelen);/*得到机器名*/
//	fprintf(stderr,"Hello World! Process %d of %d on %s\n",	myid, numprocs, processor_name);
	if(myid==0){
		fout=(FILE *)malloc(sizeof(FILE*));
		fout=fopen("printout.txt","wt+");
		fprintf(stderr,"Initializing...\n");
		fprintf(fout,"Initializing...\n");
	}
	if(argc>1 && !strcmp(SRestart,argv[1]))ReInit(argv[2]);
//	MPI_Barrier(MPI_COMM_WORLD);
	else init();
	if(myid==0){
		fprintf(stderr,"\nOK: %d processes.\n",numprocs);
		fprintf(fout,"\nOK: %d processes.\n",numprocs);
	}
	if(myid==0)output(0);
#if MD_ON
		for(i=0;i<MTimes;i++){
			if(!GROW_ON && StopMsg && myid==0)break;
			if(myid==0){
				fprintf(stderr,"%4.2f%%\n",i*100.0/MTimes);
				fprintf(fout,"%4.2f%%\n",i*100.0/MTimes);
			}
			if(i<10)BUG=0;
			else BUG=0;
				MoleDynamic(i);
			if(dtTime)Set_Temp();
			if(myid==0){
			
				output(i+1);
				if(i%10==9){
					srand(rand());
					fclose(fout);
					fout=fopen("printout.txt","at+");
				}
			}
			#if Temp_Kappa==1
				if(KappaH-KappaL<1e-3){
					StopMsg=2;
					break;
				}
			#endif
		}
		end=clock()*1.0;
		totaltime=(end-start)/CLOCKS_PER_SEC;
		if(myid==0){
			fprintf(stderr,"\nTime: %g s ( %d h %d m %g s).\n",
			totaltime,(int)(totaltime/3600),(int)(totaltime/60)-60*(int)(totaltime/3600),totaltime-60*(int)(totaltime/60));
			fprintf(fout,"\nTime: %g s ( %d h %d m %g s).\n",
			totaltime,(int)(totaltime/3600),(int)(totaltime/60)-60*(int)(totaltime/3600),totaltime-60*(int)(totaltime/60));
		}
#else
		if(numprocs>1){
			fprintf(stderr,"No more than 1 process in MC.");
			fprintf(fout,"No more than 1 process in MC.");
			fclose(fout);
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		fprintf(stderr,"RAND_MAX: %d\t%g\n",RAND_MAX,(double)RAND_MAX);
		fprintf(fout,"RAND_MAX: %d\t%g\n",RAND_MAX,(double)RAND_MAX);
		
		fprintf(stderr,"\tAccept\t\tTempAcpt\tAAcpt\tMoleNum(T/A/N)");
		fprintf(fout,"\tAccept\t\tTempAcpt\tAAcpt\tMoleNum(T/A/N)");
		output(0);
		for(i=0;i<MTimes;i++){
			if(i>Tst && i<Ted && i%DeltaT==0){
				if(DeltaXi!=0.0){
					xi0+=DeltaXi;
					Tubp=&Tub[j];
					for(j=0;j<N;j++){
						Set_Tub(Tubp,Tubp->r,Tubp->angle,xi0);
					}
				}
				if(DeltaA!=0.0){
					#if !Whole_Struc
						{
						double DAngle[5],MyAngle,NewAngle[5];
						DAngle[0]=1;
						DAngle[1]=0;
						DAngle[2]=0;
						#if ForceXYZ & 0x10
						for(j=0;j<N;j++){
							Tubp=&Tub[j];
							MyAngle=DeltaA*j/N/2;
							DAngle[3]=sin(MyAngle);
							DAngle[4]=cos(MyAngle);
							QuaTimes(Tubp->angle,DAngle,NewAngle);
							Set_Tub(Tubp,Tubp->r,NewAngle,Tubp->xi);
						}
						#else
							#if ForceXYZ & 0x20
							j=N-1;
							#elif ForceXYZ & 0x40
							j=0;
							#endif
							{
								Tubp=&Tub[j];
								MyAngle=DeltaA/2;
								DAngle[3]=sin(MyAngle);
								DAngle[4]=cos(MyAngle);
								QuaTimes(Tubp->angle,DAngle,NewAngle);
								Set_Tub(Tubp,Tubp->r,NewAngle,Tubp->xi);
							}
						#endif
					}
					#elif ForceXYZ & 0x08
						{
						double y,z,R,DAngle[5],MyAngle,NewAngle[5];
						DAngle[0]=1;
						DAngle[1]=0;
						DAngle[2]=0;
						for(j=0;j<N;j++)
						#if ForceXYZ & 0x20
							if(j>=N-13)
						#elif ForceXYZ & 0x40
							if(j<13)
						#else
							if(0)
						#endif
						{
							Tubp=&Tub[j];
							y=Tubp->r[1];
							z=Tubp->r[2];
							R=sqrt(y*y+z*z);
							MyAngle=acos(y/R);
							if(z<0)MyAngle=-MyAngle;
							MyAngle+=DeltaA;
							Tubp->r[1]=R*cos(MyAngle);
							Tubp->r[2]=R*sin(MyAngle);
							DAngle[3]=sin(DeltaA);
							DAngle[4]=cos(DeltaA);
							QuaTimes(Tubp->angle,DAngle,NewAngle);
							Set_Tub(Tubp,Tubp->r,NewAngle,Tubp->xi);
						}
					}
					#else
					;
					#endif
				}
				if(DeltaL!=0.0){
					#if ForceXYZ & 0x10
					for(j=0;j<N;j++)
					{
						Tubp=&Tub[j];
						#if ForceXYZ & 0x01
						Tubp->r[0]+=DeltaL*Tubp->rinit[0];
						#endif
						#if ForceXYZ & 0x02
						Tubp->r[1]+=DeltaL*Tubp->rinit[1];
						#endif
						#if ForceXYZ & 0x04
						Tubp->r[2]+=DeltaL*Tubp->rinit[2];
						#endif
						Set_Tub(Tubp,Tubp->r,Tubp->angle,Tubp->xi);
					}
					#else
						for(j=0;j<N;j++)
						#if ForceXYZ & 0x20
							#if Whole_Struc
							if(j>=N-13)
							#else
						//	if(j==N-1 || j==N-Lat_Length)
							if(j/Lat_Length==(N-1)/Lat_Length)
							#endif
						#elif ForceXYZ & 0x40
							#if Whole_Struc
							if(j<13)
							#else
						//	if(j==0 || j==Lat_Length-1)
							if(!(j/Lat_Length))
							#endif
						#elif ForceXYZ & 0x80
							#if Whole_Struc
							if(0)
							#else
							if(((ForceXYZ & 0x04) && !(j/Lat_Length)) || ((ForceXYZ & 0x01) && !(j%Lat_Length)))
							#endif
						#endif
						{
							Tubp=&Tub[j];
							#if ForceXYZ & 0x01
							Tubp->r[0]+=DeltaL;
							#endif
							#if ForceXYZ & 0x02
							Tubp->r[1]+=DeltaL;
							#endif
							#if ForceXYZ & 0x04
							Tubp->r[2]+=DeltaL;
							#endif
							Set_Tub(Tubp,Tubp->r,Tubp->angle,Tubp->xi);
						}
					#endif
				}
			}
			Monte_Carlo(i);
			output(i+1);
			fprintf(stderr,"\n%d\t%g  \t%g  \t%g  \t%d/%d/%d",i,RAccept/N/Times,1-at2/at1,1-2.0*NA_A/Times/N,Ntype,Nactive,N);
			fprintf(fout,"\n%d\t%g  \t%g  \t%g  \t%d/%d/%d",i,RAccept/N/Times,1-at2/at1,1-2.0*NA_A/Times/N,Ntype,Nactive,N);
			if(i%10==9){
				fclose(fout);
				fout=fopen("printout.txt","at+");
			}
			#if GROW_ON && G_Slow
			if(TGrow<Times)TGrow+=G_Slow;
			else if(TGrow<50*Times)TGrow+=GSteps;
			else TGrow+=Times;
			#endif
		}
	
		end=clock()*1.0;
		totaltime=(end-start)/CLOCKS_PER_SEC;

		fprintf(stderr,"\nTime: %g s ( %d d %d h %d m %g s).\n",
		totaltime,(int)(totaltime/3600/24),(int)(totaltime/3600)-24*(int)(totaltime/3600/24),(int)(totaltime/60)-60*(int)(totaltime/3600),totaltime-60*(int)(totaltime/60));
		fprintf(fout,"\nTime: %g s ( %d d %d h %d m %g s).\n",
		totaltime,(int)(totaltime/3600/24),(int)(totaltime/3600)-24*(int)(totaltime/3600/24),(int)(totaltime/60)-60*(int)(totaltime/3600),totaltime-60*(int)(totaltime/60));
		fprintf(stderr,"\nData:\nG_rate:\t%g\nC_rate:\t%g\nS_rate:\t%g\nB_rate:\t%g\n",G_rate,C_rate,S_rate,B_rate);
		fprintf(fout,"\nData:\nG_rate:\t%g\nC_rate:\t%g\nS_rate:\t%g\nB_rate:\t%g\n",G_rate,C_rate,S_rate,B_rate);
		fprintf(stderr,"No:\nAll:\t%d\nGTP:\t%d\nGDP:\t%d\nNon-active:\t%d\n",N,N-Ntype,Ntype,N-Nactive);
		fprintf(fout,"No:\nAll:\t%d\nGTP:\t%d\nGDP:\t%d\nNon-active:\t%d\n",N,N-Ntype,Ntype,N-Nactive);
#endif
	if(myid==0){
		fclose(fout);
#if EOut_ON
		fprintf(fshort,"\nConvergence Condition: E:%g T:%d\nConvergence Time:%8d\nEND",E_Conver,T_Conver,TConver);
		fclose(fshort);
		fclose(fphase);
#endif
	}
	{
	free(Tub);
	if(ttub)free(ttub);
	if(Head)free(Head);
	if(Top)free(Top);
	if(List)free(List);
	if(CList)free(CList);
	if(Chain)free(Chain);
	if(TGroup)free(TGroup);
	if(Tran)free(Tran);
	if(Tmark)free(Tmark);
	if(TubState)free(TubState);
#if EOut_ON
		free(EPoten);
		free(EVel);
		free(ERot);
		free(ETot);
		free(FAvg);
		free(FBall);
#endif
	}
	if(!myid)switch(StopMsg){
		case 1:fprintf(stderr,"Finished because of convergence.\n");break;
		case 2:fprintf(stderr,"Finished at the critical point.\n");break;
		default:break;
	}
	if(numprocs!=1){
		if(StopMsg){
			MPI_Abort(MPI_COMM_WORLD,1);
			exit(0);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();/*结束*/
	return 1;
}


