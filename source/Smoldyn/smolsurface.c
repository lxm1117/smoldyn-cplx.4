/* Steven Andrews, started 10/22/2001.
 This is a library of functions for the Smoldyn program.
 See documentation called Smoldyn_doc1.pdf and Smoldyn_doc2.pdf, and the Smoldyn
 website, which is at www.smoldyn.org.
 Copyright 2003-2013 by Steven Andrews.  This work is distributed under the terms
 of the Gnu Lesser General Public License (LGPL). */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "Geometry.h"
#include "math2.h"
#include "random2.h"
#include "Rn.h"
#include "RnSort.h"
#include "SurfaceParam.h"
#include "smoldyn.h"
#include "smoldynfuncs.h"
#include "Sphere.h"
#include "Zn.h"
#include <sstream>
#include <string>

#include "smoldynconfigure.h"


/******************************************************************************/
/********************************** Surfaces **********************************/
/******************************************************************************/


/******************************************************************************/
/****************************** Local declarations ****************************/
/******************************************************************************/

// enumerated types
enum SrfAction surfstring2act(char *string);
char *surfact2string(enum SrfAction act,char *string);

// low level utilities
panelptr readpanelname(simptr sim,surfaceptr srf,const char *str);
int panelpoints(enum PanelShape ps,int dim);
int surfpanelparams(enum PanelShape ps,int dim);
void panelmiddle(panelptr pnl,double *middle,int dim,int onpanel);
int srfsamestate(enum MolecState ms1,enum PanelFace face1,enum MolecState ms2,enum MolecState *ms3ptr);
void srfreverseaction(enum MolecState ms1,enum PanelFace face1,enum MolecState ms2,enum MolecState *ms3ptr,enum PanelFace *face2ptr,enum MolecState *ms4ptr);
void srftristate2index(enum MolecState ms,enum MolecState ms1,enum MolecState ms2,enum MolecState *ms3ptr,enum PanelFace *faceptr,enum MolecState *ms4ptr);
void srfindex2tristate(enum MolecState ms3,enum PanelFace face,enum MolecState ms4,enum MolecState *msptr,enum MolecState *ms1ptr,enum MolecState *ms2ptr);
int withincmptcheck(simptr sim, moleculeptr mptr);

// memory management
surfactionptr surfaceactionalloc(int species);
void surfaceactionfree(surfactionptr actdetails);
int panelsalloc(surfaceptr srf,int dim,int maxpanel,int maxspecies,enum PanelShape ps);
void panelfree(panelptr pnl);
int emittersalloc(surfaceptr srf,enum PanelFace face,int oldmaxspecies,int maxspecies);
surfaceptr surfacealloc(surfaceptr srf,int oldmaxspecies,int maxspecies,int dim);
void surfacefree(surfaceptr srf,int maxspecies);
surfacessptr surfacessalloc(surfacessptr srfss,int maxsurface,int maxspecies,int dim);

// data structure output

// structure set up
int surfsetmaxpanel(surfaceptr srf,int dim,enum PanelShape ps,int maxpanel);
int surfsetemitterabsorption(simptr sim);
double srfcalcrate(simptr sim,surfaceptr srf,int i,enum MolecState ms1,enum PanelFace face,enum MolecState ms2);
double srfcalcprob(simptr sim,surfaceptr srf,int i,enum MolecState ms1,enum PanelFace face,enum MolecState ms2);
int surfupdateparams(simptr sim);
int surfupdatelists(simptr sim);

// core simulation functions
void panelnormal(panelptr pnl,double *pos,enum PanelFace face,int dim,double *norm);
void movept2panel(double *pt,panelptr pnl,int dim,double margin);
int closestpanelpt(panelptr pnl,int dim,double *testpt,double *pnlpt);
void surfacereflect(moleculeptr mptr,panelptr pnl,double *crsspt,int dim,enum PanelFace face);
int surfacejump(moleculeptr mptr,panelptr pnl,double *crsspt,enum PanelFace face,int dim);
int rotating_jump(moleculeptr mptr,panelptr pnl, double *crsspt,enum PanelFace face, int dim);
int dosurfinteract(simptr sim,moleculeptr mptr,int ll,int m,panelptr pnl,enum PanelFace face,double *crsspt, int *actptr);
int checksurfaces_cplx(simptr sim, moleculeptr mptr, int m, int ll, int reborn);

/******************************************************************************/
/****************************** enumerated types ******************************/
/******************************************************************************/

/* surfstring2face */
enum PanelFace surfstring2face(char *string) {
	enum PanelFace ans;

	if(strbegin(string,"front",0)) ans=PFfront;
	else if(strbegin(string,"back",0)) ans=PFback;
	else if(strbegin(string,"all",0) || strbegin(string,"both",0)) ans=PFboth;
	else ans=PFnone;
	return ans; }


/* surfface2string */
char *surfface2string(enum PanelFace face,char *string) {
	if(face==PFfront) strcpy(string,"front");
	else if(face==PFback) strcpy(string,"back");
	else if(face==PFboth) strcpy(string,"both");
	else strcpy(string,"none");
	return string; }


/* surfstring2act */
enum SrfAction surfstring2act(char *string) {
	enum SrfAction ans;

	if(strbegin(string,"reflect",0)) ans=SAreflect;
	else if(strbegin(string,"transmit",0)) ans=SAtrans;
	else if(strbegin(string,"absorb",0)) ans=SAabsorb;
	else if(strbegin(string,"jump",0)) ans=SAjump;	
	else if(!strcmp(string,"periodic")) ans=SAjump;
	else if(strbegin(string, "rotate_jump",0)) ans=SArotate_jump;		// cplx
	else if(!strcmp(string,"rotate_periodic")) ans=SArotate_jump;
	else if(!strcmp(string,"port")) ans=SAport;
	else if(strbegin(string,"multiple",0)) ans=SAmult;
	else if(strbegin(string,"no",0)) ans=SAno;
	else if(strbegin(string,"adsorb",0)) ans=SAadsorb;
	else if(strbegin(string,"revdes",0)) ans=SArevdes;
	else if(strbegin(string,"irrevdes",0)) ans=SAirrevdes;
	else if(strbegin(string,"flip",0)) ans=SAflip;
	else ans=SAnone;
	return ans; }


/* surfact2string */
char *surfact2string(enum SrfAction act,char *string) {
	if(act==SAreflect) strcpy(string,"reflect");
	else if(act==SAtrans) strcpy(string,"transmit");
	else if(act==SAabsorb) strcpy(string,"absorb");
	else if(act==SAjump) strcpy(string,"jump");
	else if(act==SAport) strcpy(string,"port");
	else if(act==SAmult) strcpy(string,"multiple");
	else if(act==SAno) strcpy(string,"no");
	else if(act==SAadsorb) strcpy(string,"adsorb");
	else if(act==SArevdes) strcpy(string,"revdes");
	else if(act==SAirrevdes) strcpy(string,"irrevdes");
	else if(act==SAflip) strcpy(string,"flip");
	else if(act==SArotate_jump) strcpy(string,"rotate_jump");
	else strcpy(string,"none");
	return string; }


/* surfstring2ps */
enum PanelShape surfstring2ps(char *string) {
	enum PanelShape ans;

	if(strbegin(string,"rectangle",0)) ans=PSrect;
	else if(strbegin(string,"triangle",0)) ans=PStri;
	else if(strbegin(string,"sphere",0)) ans=PSsph;
	else if(strbegin(string,"cylinder",0)) ans=PScyl;
	else if(strbegin(string,"hemisphere",0)) ans=PShemi;
	else if(strbegin(string,"disk",0)) ans=PSdisk;
	else if(strbegin(string,"all",0)) ans=PSall;
	else ans=PSnone;
	return ans; }


/* surfps2string */
char *surfps2string(enum PanelShape ps,char *string) {
	if(ps==PSrect) strcpy(string,"rect");
	else if(ps==PStri) strcpy(string,"tri");
	else if(ps==PSsph) strcpy(string,"sph");
	else if(ps==PScyl) strcpy(string,"cyl");
	else if(ps==PShemi) strcpy(string,"hemi");
	else if(ps==PSdisk) strcpy(string,"disk");
	else if(ps==PSall) strcpy(string,"all");
	else strcpy(string,"none");
	return string; }


/* surfstring2dm */
enum DrawMode surfstring2dm(char *string) {
	enum DrawMode ans;

	if(strbegin(string,"none",0)) ans=DMno;
	else if(!strcmp(string,"ve") || !strcmp(string,"ev")) ans=DMve;
	else if(!strcmp(string,"vf") || !strcmp(string,"fv")) ans=DMvf;
	else if(!strcmp(string,"vef") || !strcmp(string,"vfe") || !strcmp(string,"evf") || !strcmp(string,"efv") || !strcmp(string,"fev") || !strcmp(string,"fve")) ans=DMvef;
	else if(strbegin(string,"vertex",0)) ans=DMvert;
	else if(strbegin(string,"edge",0)) ans=DMedge;
	else if(strbegin(string,"face",0)) ans=DMface;
	else ans=DMnone;
	return ans; }


/* surfdm2string */
char *surfdm2string(enum DrawMode dm,char *string) {
	if(dm==DMno) strcpy(string,"no");
	else if(dm==DMvert) strcpy(string,"vert");
	else if(dm==DMedge) strcpy(string,"edge");
	else if(dm==DMve) strcpy(string,"ve");
	else if(dm==DMface) strcpy(string,"face");
	else if(dm==DMvf) strcpy(string,"vf");
	else if(dm==DMef) strcpy(string,"ef");
	else if(dm==DMvef) strcpy(string,"vef");
	else strcpy(string,"none");
	return string; }


/******************************************************************************/
/***************************** low level utilities ****************************/
/******************************************************************************/


/* readsurfacename. */
int readsurfacename(simptr sim,const char *str,enum PanelShape *psptr,int *pptr) {
	char snm[STRCHAR],pnm[STRCHAR],*colon;
	int itct,s,p;
	enum PanelShape ps;

	if(!str) return -1;
	if(!sim->srfss || !sim->srfss->nsrf) return -2;

	itct=sscanf(str,"%s",snm);
	if(itct!=1) return -3;														// cannot read name
	colon=strchr(snm,':');
	if(colon) {
		strcpy(pnm,colon+1);
		*colon='\0'; }
	else pnm[0]='\0';
	ps=PSnone;
	p=-1;

	if(!strcmp(snm,"all")) {													// all surfaces
		s=-5;
		if(pnm[0]=='\0');																// all surface, no panel
		else if(!strcmp(pnm,"all")) {										// all:all
			ps=PSall;
			p=-5; }
		else																						// all:panel
			p=-2; }
	else {																						// specific surface
		s=stringfind(sim->srfss->snames,sim->srfss->nsrf,snm);
		if(s==-1)																				// unknown surface
			s=-4;
		else {
			if(pnm[0]=='\0');															// no panel name
			else if(!strcmp(pnm,"all")) {									// surface:all
				ps=PSall;
				p=-5; }
			else if(VCellDefined && strstr(pnm,"tri_")==pnm) {		// surface:tri_#_#_#	(VCELL)
				ps=PStri;
				int memIndex;
				int globalIndex;
				sscanf(pnm,"tri_%d_%d_%d",&p,&globalIndex,&memIndex); }
			else {																				// surface:panel
				for(ps=(PanelShape)0;p==-1 && ps<PSMAX;ps=(PanelShape)(ps+1))
					p=stringfind(sim->srfss->srflist[s]->pname[ps],sim->srfss->srflist[s]->npanel[ps],pnm);
				if(p==-1) {
					ps=PSnone;
					p=-3; }
				else ps=(PanelShape)(ps-1); }}}

	if(psptr) *psptr=ps;
	if(pptr) *pptr=p;
	return s; }


/* readpanelname */
panelptr readpanelname(simptr sim,surfaceptr srf,const char *str) {
	enum PanelShape ps;
	char name[STRCHAR];
	int s,p;

	if(strchr(str,':')) strcpy(name,str);
	else if(srf) sprintf(name,"%s:%s",srf->sname,str);
	else return NULL;
	s=readsurfacename(sim,name,&ps,&p);
	if(s<0 || p<0) return NULL;
	return sim->srfss->srflist[s]->panels[ps][p]; }


/* panelpoints */
int panelpoints(enum PanelShape ps,int dim) {
	int npts;

	if(ps==PSrect) {
		if(dim==1) npts=1;
		else if(dim==2) npts=4;
		else npts=8; }
	else if(ps==PStri) {
		if(dim==1) npts=1;
		else if(dim==2) npts=4;
		else npts=6; }
	else if(ps==PSsph) npts=2;
	else if(ps==PScyl && dim>1) npts=5;
	else if(ps==PShemi && dim>1) npts=3;
	else if(ps==PSdisk && dim>1) npts=2;
	else npts=0;
	return npts; }


/* surfpanelparams */
int surfpanelparams(enum PanelShape ps,int dim) {
	int n;

	if(ps==PSrect) n=2*dim-1;
	else if(ps==PStri) n=dim*dim;
	else if(ps==PSsph) n=2*dim;
	else if(ps==PScyl && dim>1) n=(dim==2)?5:9;
	else if(ps==PShemi && dim>1) n=3*dim;
	else if(ps==PSdisk && dim>1) n=(dim==2)?5:8;
	else n=0;
	return n; }


/* panelmiddle */
void panelmiddle(panelptr pnl,double *middle,int dim,int onpanel) {
	enum PanelShape ps;
	double **point,norm[DIMMAX];
	int d;

	ps=pnl->ps;
	point=pnl->point;
	if(ps==PSrect) {
		Geo_RectCenter(point,middle,dim); }
	else if(ps==PStri) {
		Geo_TriCenter(point,middle,dim); }
	else if(ps==PSsph) {
		for(d=0;d<dim;d++)
			middle[d]=point[0][d];
		if(onpanel)
			middle[0]+=point[1][0]; }
	else if(ps==PScyl) {
		Geo_LineCenter(point,middle,dim);
		if(onpanel) {
			if(dim==2) Geo_LineNormal2D(point[0],point[1],middle,norm);
			else Geo_LineNormal3D(point[0],point[1],middle,norm);
			for(d=0;d<dim;d++)
				middle[d]+=norm[d]*point[2][0]; }}
	else if(ps==PShemi) {
		for(d=0;d<dim;d++)
			middle[d]=point[0][d];
		if(onpanel)
			for(d=0;d<dim;d++)
				middle[d]-=point[2][d]*point[1][0]; }
	else if(ps==PSdisk) {
		for(d=0;d<dim;d++)
			middle[d]=point[0][d]; }
	return; }


/* panelarea */
double panelarea(panelptr pnl,int dim) {
	enum PanelShape ps;
	double **point,*front,area;
	int d;

	ps=pnl->ps;
	point=pnl->point;
	front=pnl->front;

	if(dim==1) {
		if(ps==PSrect || ps==PStri) area=1;
		else if(ps==PSsph) area=2;
		else area=0; }
	else if(dim==2) {
		if(ps==PSrect) area=fabs(point[1][(int)front[2]]-point[0][(int)front[2]]);
		else if(ps==PStri) area=sqrt((point[1][0]-point[0][0])*(point[1][0]-point[0][0])+(point[1][1]-point[0][1])*(point[1][1]-point[0][1]));
		else if(ps==PSsph) area=2.0*PI*point[1][0];
		else if(ps==PScyl) area=2.0*sqrt((point[1][0]-point[0][0])*(point[1][0]-point[0][0])+(point[1][1]-point[0][1])*(point[1][1]-point[0][1]));
		else if(ps==PShemi) area=PI*point[1][0];
		else if(ps==PSdisk) area=2.0*point[1][0];
		else area=0; }
	else if(dim==3) {
		if(ps==PSrect) {
			d=0;
			while(d==(int)front[1] || d==(int)front[2]) d++;
			area=fabs((point[2][(int)front[2]]-point[0][(int)front[2]])*(point[2][d]-point[0][d])); }
		else if(ps==PStri) area=Geo_TriArea3(point[0],point[1],point[2],front);
		else if(ps==PSsph) area=4.0*PI*point[1][0]*point[1][0];
		else if(ps==PScyl) area=2.0*PI*point[2][0]*sqrt((point[1][0]-point[0][0])*(point[1][0]-point[0][0])+(point[1][1]-point[0][1])*(point[1][1]-point[0][1])+(point[1][2]-point[0][2])*(point[1][2]-point[0][2]));
		else if(ps==PShemi) area=2.0*PI*point[1][0]*point[1][0];
		else if(ps==PSdisk) area=PI*point[1][0]*point[1][0];
		else area=0; }
	else area=0;

	return area; }


/* surfacearea */
double surfacearea(surfaceptr srf,int dim,int *totpanelptr) {
	int ps,p,totpanel;
	double area;

	totpanel=0;
	area=0;
	for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
		for(p=0;p<srf->npanel[ps];p++) {
			totpanel++;
			area+=panelarea(srf->panels[ps][p],dim); }
	if(totpanelptr) *totpanelptr=totpanel;
	return area; }


/* surfacearea2 */
double surfacearea2(simptr sim,int surface,enum PanelShape ps,char *pname,int *totpanelptr) {
	int panel,dim,totpanel;
	double area;
	surfaceptr srf;
	int slo,shi,pslo,pshi,plo,phi,s,p;

	dim=sim->dim;
	if(ps==PSnone) {area=0;totpanel=0;}											// ps is none

	else if(surface>=0 && ps!=PSall && pname && strcmp(pname,"all")) {		// specific panel
		srf=sim->srfss->srflist[surface];
		panel=stringfind(srf->pname[ps],srf->npanel[ps],pname);
		if(panel<0) {
			area=0;
			totpanel=0; }
		else {
			area=panelarea(srf->panels[ps][panel],sim->dim);
			totpanel=1; }}

	else {																									// multiple panels
		slo=(surface>=0)?surface:0;
		shi=(surface>=0)?surface+1:sim->srfss->nsrf;
		pslo=(ps!=PSall)?ps:0;
		pshi=(ps!=PSall)?ps+1:PSMAX;

		area=0;
		totpanel=0;
		for(s=slo;s<shi;s++)
			for(ps=(PanelShape)pslo;ps<pshi;ps=(PanelShape)(ps+1)) {
				srf=sim->srfss->srflist[s];
				if(!pname || !strcmp(pname,"all")) {plo=0;phi=srf->npanel[ps];}
				else if((panel=stringfind(srf->pname[ps],srf->npanel[ps],pname))<0) plo=phi=0;
				else {plo=panel;phi=panel+1;}
				for(p=plo;p<phi;p++) {
					area+=surfacearea2(sim,s,ps,srf->pname[ps][p],NULL);
					totpanel++; }}}

	if(totpanelptr) *totpanelptr=totpanel;
	return area; }


/* panelrandpos */
void panelrandpos(panelptr pnl,double *pos,int dim) {
	int d;
	double **point,*front,flt1,flt2,v1[DIMMAX],dcm[DIMMAX*DIMMAX];
	enum PanelShape ps;

	point=pnl->point;
	front=pnl->front;
	ps=pnl->ps;

	if(dim==1) {
		if(ps==PSrect || ps==PStri)
			pos[0]=point[0][0];
		else if(ps==PSsph)
			pos[0]=point[0][0]+(coinrandD(0.5)?point[1][0]:-point[1][0]);
		else
			pos[0]=0; }

	else if(dim==2) {
		if(ps==PSrect) {
			pos[(int)front[1]]=point[0][(int)front[1]];
			pos[(int)front[2]]=unirandCCD(point[0][(int)front[2]],point[1][(int)front[2]]); }
		else if(ps==PStri) {
			flt1=randCCD();
			pos[0]=point[0][0]+flt1*(point[1][0]-point[0][0]);
			pos[1]=point[0][1]+flt1*(point[1][1]-point[0][1]); }
		else if(ps==PSsph) {
			flt1=unirandCOD(0,2.0*PI);
			pos[0]=point[0][0]+point[1][0]*cos(flt1);
			pos[1]=point[0][1]+point[1][0]*sin(flt1); }
		else if(ps==PScyl) {
			flt1=randCCD();
			flt2=coinrandD(0.5)?1:-1;
			pos[0]=point[0][0]+flt1*(point[1][0]-point[0][0])+flt2*point[2][0]*front[0];
			pos[1]=point[0][1]+flt1*(point[1][1]-point[0][1])+flt2*point[2][0]*front[1]; }
		else if(ps==PShemi) {
			flt1=unirandCCD(0,PI)+atan2(point[2][1],point[2][0])+PI/2.0;
			pos[0]=point[0][0]+point[1][0]*cos(flt1);
			pos[1]=point[0][1]+point[1][0]*sin(flt1); }
		else if(ps==PSdisk) {
			flt1=unirandCCD(-1,1);
			pos[0]=point[0][0]-flt1*point[1][0]*front[1];
			pos[1]=point[0][1]+flt1*point[1][0]*front[0]; }
		else
			pos[0]=pos[1]=0; }

	else if(dim==3) {
		if(ps==PSrect) {
			pos[0]=unirandCCD(point[0][0],point[2][0]);
			pos[1]=unirandCCD(point[0][1],point[2][1]);
			pos[2]=unirandCCD(point[0][2],point[2][2]);
			pos[(int)front[1]]=point[0][(int)front[1]]; }
		else if(ps==PStri) {
			trianglerandCD(point[0],point[1],point[2],dim,pos); }
		else if(ps==PSsph) {
			sphererandCCD(pos,point[1][0],point[1][0]);
			for(d=0;d<dim;d++) pos[d]+=point[0][d]; }
		else if(ps==PScyl) {
			flt1=randCCD();
			flt2=unirandCOD(0,2.0*PI);
			v1[0]=point[1][0]-point[0][0];
			v1[1]=point[1][1]-point[0][1];
			v1[2]=point[1][2]-point[0][2];
			Sph_Newz2Dcm(v1,0,dcm);
			v1[0]=point[2][0]*cos(flt2);
			v1[1]=point[2][0]*sin(flt2);
			v1[2]=0;
			dotMVD(dcm,v1,pos,3,3);
			pos[0]+=point[0][0]+flt1*(point[1][0]-point[0][0]);
			pos[1]+=point[0][1]+flt1*(point[1][1]-point[0][1]);
			pos[2]+=point[0][2]+flt1*(point[1][2]-point[0][2]); }
		else if(ps==PShemi) {
			flt1=thetarandCCD();
			if(flt1<PI/2.0) flt1=PI-flt1;
			flt2=unirandCOD(0,2.0*PI);
			v1[0]=point[1][0]*sin(flt1)*cos(flt2);
			v1[1]=point[1][0]*sin(flt1)*sin(flt2);
			v1[2]=point[1][0]*cos(flt1);
			Sph_Newz2Dcm(point[2],0,dcm);
			dotMVD(dcm,v1,pos,3,3);
			pos[0]+=point[0][0];
			pos[1]+=point[0][1];
			pos[2]+=point[0][2]; }
		else if(ps==PSdisk) {
			flt1=unirandCOD(0,2.0*PI);
			flt2=radrandcircCCD(point[1][0]);
			v1[0]=flt2*cos(flt1);
			v1[1]=flt2*sin(flt1);
			v1[2]=0;
			Sph_Newz2Dcm(front,0,dcm);
			dotMVD(dcm,v1,pos,3,3);
			pos[0]+=point[0][0];
			pos[1]+=point[0][1];
			pos[2]+=point[0][2]; }
		else
			pos[0]=pos[1]=0; }

	return; }


/* surfrandpos */
panelptr surfrandpos(surfaceptr srf,double *pos,int dim) {
	panelptr pnl;

	if(!srf->totpanel) return NULL;
	pnl=srf->paneltable[intrandpD(srf->totpanel,srf->areatable)];
	panelrandpos(pnl,pos,dim);
	return pnl; }


/* issurfprod */
int issurfprod(simptr sim,int i,enum MolecState ms) {
	surfacessptr srfss;
	int s,i1;
	surfaceptr srf;
	enum MolecState ms1;
	enum PanelFace face;
	surfactionptr actdetails;

	srfss=sim->srfss;
	if(!srfss) return 0;

	for(s=0;s<srfss->nsrf;s++) {
		srf=srfss->srflist[s];
		i1=i;															// first try no species change at surface
		for(ms1=(MolecState)0;ms1<MSMAX;ms1=(MolecState)(ms1+1))
			for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1)) {
				actdetails=srf->actdetails[i1][ms1][face];
				if(actdetails)
					if(actdetails->srfrate[ms]>0 || actdetails->srfprob[ms]>0)
						if(actdetails->srfnewspec[ms]==i)
							return 1; }
		for(i1=0;i1<srfss->maxspecies;i1++)	// failed, so try with species change at surface
			for(ms1=(MolecState)0;ms1<MSMAX;ms1=(MolecState)(ms1+1))
				for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1)) {
					actdetails=srf->actdetails[i1][ms1][face];
					if(actdetails)
						if(actdetails->srfrate[ms]>0 || actdetails->srfprob[ms]>0)
							if(actdetails->srfnewspec[ms]==i)
								return 1; }}
		return 0; }


/* srfsamestate */
int srfsamestate(enum MolecState ms1,enum PanelFace face1,enum MolecState ms2,enum MolecState *ms3ptr) {
	int same;
	enum MolecState ms3;

	if(face1==PFfront && ms2==MSsoln) same=1;
	else if(face1==PFback && ms2==MSbsoln) same=1;
	else if(face1==PFnone && ms1==ms2) same=1;
	else same=0;

	if(ms3ptr) {
		if(ms1==MSsoln) {
			if(face1==PFfront) ms3=MSsoln;
			else if(face1==PFback) ms3=MSbsoln;
			else ms3=MSnone; }
		else {
			if(face1==PFfront) ms3=MSsoln;
			else if(face1==PFback) ms3=MSbsoln;
			else ms3=ms1; }
		*ms3ptr=ms3; }
	return same; }


/* srfreverseaction */
void srfreverseaction(enum MolecState ms1,enum PanelFace face1,enum MolecState ms2,enum MolecState *ms3ptr,enum PanelFace *face2ptr,enum MolecState *ms4ptr) {
	enum MolecState ms3,ms4;
	enum PanelFace face2;

	if(ms1==MSsoln && face1==PFnone) {			// soln but neither face, which is impossible
		ms3=ms4=MSnone;
		face2=PFnone; }

	else if(ms1==MSsoln) {									// soln and either face, for simple collision
		if(ms2==MSsoln || ms2==MSbsoln) {
			ms3=MSsoln;
			face2=(ms2==MSsoln)?PFfront:PFback; }
		else {
			ms3=ms2;
			face2=PFnone; }
		ms4=(face1==PFfront)?MSsoln:MSbsoln; }

	else if(face1==PFnone) {								// bound and no face, for desorb etc.
		if(ms2==MSsoln || ms2==MSbsoln) {
			ms3=MSsoln;
			face2=(ms2==MSsoln)?PFfront:PFback; }
		else {
			ms3=ms2;
			face2=PFnone; }
		ms4=ms1; }

	else {																	// bound and face, for surface-bound collision
		if(ms2==MSsoln || ms2==MSbsoln) {
			ms3=ms1;
			face2=(ms2==MSsoln)?PFfront:PFback;
			ms4=(face1==PFfront)?MSsoln:MSbsoln; }
		else {
			ms3=ms2;
			face2=PFboth;
			ms4=ms1; }}

	if(ms3ptr) *ms3ptr=ms3;
	if(face2ptr) *face2ptr=face2;
	if(ms4ptr) *ms4ptr=ms4;
	return; }


/* srftristate2index */
void srftristate2index(enum MolecState ms,enum MolecState ms1,enum MolecState ms2,enum MolecState *ms3ptr,enum PanelFace *faceptr,enum MolecState *ms4ptr) {
	enum MolecState ms3,ms4;
	enum PanelFace face;

	if(ms==MSnone) ms=MSsoln;

	if(ms==MSsoln && (ms1==MSsoln || ms1==MSbsoln)) ms3=MSsoln;
	else if(ms==MSsoln) ms3=ms1;
	else ms3=ms;

	if(ms1==MSsoln) face=PFfront;
	else if(ms1==MSbsoln) face=PFback;
	else face=PFnone;

	ms4=ms2;

	if(ms!=MSsoln && ms1!=MSsoln && ms1!=MSbsoln && ms1!=ms) {
		ms3=MSnone;
		face=PFnone;
		ms4=MSnone; }

	if(ms3ptr) *ms3ptr=ms3;
	if(faceptr) *faceptr=face;
	if(ms4ptr) *ms4ptr=ms4;
	return; }


/* srfindex2tristate */
void srfindex2tristate(enum MolecState ms3,enum PanelFace face,enum MolecState ms4,enum MolecState *msptr,enum MolecState *ms1ptr,enum MolecState *ms2ptr) {
	enum MolecState ms,ms1,ms2;
	ms=ms3;

	if(face==PFfront) ms1=MSsoln;
	else if(face==PFback) ms1=MSbsoln;
	else ms1=ms3;

	ms2=ms4;

	if(msptr) *msptr=ms;
	if(ms1ptr) *ms1ptr=ms1;
	if(ms2ptr) *ms2ptr=ms2;
	return; }


/******************************************************************************/
/****************************** memory management *****************************/
/******************************************************************************/


/* surfaceactionalloc */
surfactionptr surfaceactionalloc(int species) {
	surfactionptr actdetails;
	enum MolecState ms;

	actdetails=NULL;
	actdetails=(surfactionptr) malloc(sizeof(struct surfactionstruct));
	if(!actdetails) return NULL;
	actdetails->srfnewspec=NULL;
	actdetails->srfrate=NULL;
#if OPTION_VCELL
	actdetails->srfRateValueProvider=NULL;
#endif
	actdetails->srfprob=NULL;
	actdetails->srfcumprob=NULL;
	actdetails->srfdatasrc=NULL;
	actdetails->srfrevprob=NULL;

	CHECKMEM(actdetails->srfnewspec=(int*) calloc(MSMAX1,sizeof(int)));
	for(ms=(MolecState)0;ms<MSMAX1;ms=(MolecState)(ms+1)) actdetails->srfnewspec[ms]=species;

	CHECKMEM(actdetails->srfrate=(double*) calloc(MSMAX1,sizeof(double)));
	for(ms=(MolecState)0;ms<MSMAX1;ms=(MolecState)(ms+1)) actdetails->srfrate[ms]=0;

#if OPTION_VCELL
	CHECKMEM(actdetails->srfRateValueProvider=(ValueProvider**) calloc(MSMAX1,sizeof(ValueProvider*)));
	for(ms=(MolecState)0;ms<MSMAX1;ms=(MolecState)(ms+1)) actdetails->srfRateValueProvider[ms]=0;
#endif

	CHECKMEM(actdetails->srfprob=(double*) calloc(MSMAX1,sizeof(double)));
	for(ms=(MolecState)0;ms<MSMAX1;ms=(MolecState)(ms+1)) actdetails->srfprob[ms]=0;

	CHECKMEM(actdetails->srfcumprob=(double*) calloc(MSMAX1,sizeof(double)));
	for(ms=(MolecState)0;ms<MSMAX1;ms=(MolecState)(ms+1)) actdetails->srfcumprob[ms]=0;

	CHECKMEM(actdetails->srfdatasrc=(int*) calloc(MSMAX1,sizeof(int)));
	for(ms=(MolecState)0;ms<MSMAX1;ms=(MolecState)(ms+1)) actdetails->srfdatasrc[ms]=0;

	CHECKMEM(actdetails->srfrevprob=(double*) calloc(MSMAX1,sizeof(double)));
	for(ms=(MolecState)0;ms<MSMAX1;ms=(MolecState)(ms+1)) actdetails->srfrevprob[ms]=0;

	return actdetails;

 failure:
	surfaceactionfree(actdetails);
	simLog(NULL,10,"Unable to allocate memory in surfaceactionalloc");
	return NULL; }


/* surfaceactionfree */
void surfaceactionfree(surfactionptr actdetails) {
	if(!actdetails) return;
	free(actdetails->srfrevprob);
	free(actdetails->srfdatasrc);
	free(actdetails->srfcumprob);
	free(actdetails->srfprob);
	free(actdetails->srfrate);
	free(actdetails->srfnewspec);
	free(actdetails);
	return; }


/* panelsalloc. */
int panelsalloc(surfaceptr srf,int dim,int maxpanel,int maxspecies,enum PanelShape ps) {
	panelptr *newpnls,pnl;
	int p,npts,pt,d,oldmaxpanel;
	char **newpname,string[STRCHAR];

	npts=panelpoints(ps,dim);
	newpname=NULL;
	newpnls=NULL;
	CHECKBUG(srf,"missing surface parameter in panelsalloc");
	oldmaxpanel=srf->maxpanel[ps];

	if(maxpanel<=0 || maxpanel<oldmaxpanel) return 0;
	else if(maxpanel==oldmaxpanel) return 1;

	CHECKMEM(newpname=(char**) calloc(maxpanel,sizeof(char*)));
	for(p=0;p<maxpanel;p++) newpname[p]=NULL;
	for(p=0;p<oldmaxpanel;p++)
		newpname[p]=srf->pname[ps][p];
	for(;p<maxpanel;p++) {
		CHECKMEM(newpname[p]=EmptyString());
		sprintf(newpname[p],"%s%i",surfps2string(ps,string),p); }

	CHECKMEM(newpnls=(panelptr*) calloc(maxpanel,sizeof(panelptr)));
	for(p=0;p<maxpanel;p++) newpnls[p]=NULL;
	for(p=0;p<oldmaxpanel;p++)
		newpnls[p]=srf->panels[ps][p];
	for(;p<maxpanel;p++) {
		CHECKMEM(newpnls[p]=(panelptr) malloc(sizeof(struct panelstruct)));
		pnl=newpnls[p];
		pnl->pname=newpname[p];
		pnl->ps=ps;
		pnl->srf=srf;
		pnl->npts=npts;
		pnl->point=NULL;
		pnl->maxneigh=0;
		pnl->nneigh=0;
		pnl->neigh=NULL;
		pnl->emitterabsorb[PFfront]=NULL;
		pnl->emitterabsorb[PFback]=NULL;

		CHECKMEM(pnl->point=(double**) calloc(npts,sizeof(double*)));
		for(pt=0;pt<npts;pt++) pnl->point[pt]=NULL;
		for(pt=0;pt<npts;pt++) {
			CHECKMEM(pnl->point[pt]=(double*) calloc(dim,sizeof(double)));
			for(d=0;d<dim;d++) pnl->point[pt][d]=0; }
		pnl->front[0]=pnl->front[1]=pnl->front[2]=0;
		pnl->jumpp[0]=pnl->jumpp[1]=NULL;
		pnl->jumpf[0]=pnl->jumpf[1]=PFnone; }

	srf->maxpanel[ps]=maxpanel;
	free(srf->pname[ps]);
	srf->pname[ps]=newpname;
	free(srf->panels[ps]);
	srf->panels[ps]=newpnls;

	if(srf->maxemitter[PFfront]) {		// add emitter stuff to new panels
		CHECK(emittersalloc(srf,PFfront,maxspecies,maxspecies)==0); }
	if(srf->maxemitter[PFback]) {
		CHECK(emittersalloc(srf,PFback,maxspecies,maxspecies)==0); }
		
	return 1;

 failure:
	if(ErrorType!=1) simLog(NULL,10,"Unable to allocate memory in panelsalloc");
	return 0; }


/* panelfree */
void panelfree(panelptr pnl) {
	int pt;

	if(!pnl) return;
	free(pnl->emitterabsorb[PFback]);
	free(pnl->emitterabsorb[PFfront]);
	free(pnl->neigh);
	if(pnl->npts && pnl->point) {
		for(pt=0;pt<pnl->npts;pt++)
			if(pnl->point[pt]) free(pnl->point[pt]);
		free(pnl->point); }
	free(pnl);
	return; }


/* emittersalloc */
int emittersalloc(surfaceptr srf,enum PanelFace face,int oldmaxspecies,int maxspecies) {
	int i1,p;
	enum PanelShape ps;
	panelptr pnl;
	int *newmaxemitter,*newnemitter;
	double **newemitteramount,***newemitterpos,*newemitterabsorb;

	newmaxemitter=NULL;
	newnemitter=NULL;
	newemitteramount=NULL;
	newemitterpos=NULL;
	newemitterabsorb=NULL;

	if(!srf->maxemitter[face])				// creating emitters from scratch
		oldmaxspecies=0;

	if(maxspecies>oldmaxspecies) {		// surface structure data, for either from scratch or larger maxspecies
		CHECKMEM(newmaxemitter=(int*) calloc(maxspecies,sizeof(int)));
		for(i1=0;i1<oldmaxspecies;i1++) newmaxemitter[i1]=srf->maxemitter[face][i1];
		for(;i1<maxspecies;i1++) newmaxemitter[i1]=0;

		CHECKMEM(newnemitter=(int*) calloc(maxspecies,sizeof(int)));
		for(i1=0;i1<oldmaxspecies;i1++) newnemitter[i1]=srf->nemitter[face][i1];
		for(;i1<maxspecies;i1++) newnemitter[i1]=0;
		
		CHECKMEM(newemitteramount=(double**) calloc(maxspecies,sizeof(double*)));
		for(i1=0;i1<oldmaxspecies;i1++) newemitteramount[i1]=srf->emitteramount[face][i1];
		for(;i1<maxspecies;i1++) newemitteramount[i1]=NULL;

		CHECKMEM(newemitterpos=(double***) calloc(maxspecies,sizeof(double**)));
		for(i1=0;i1<oldmaxspecies;i1++) newemitterpos[i1]=srf->emitterpos[face][i1];
		for(;i1<maxspecies;i1++) newemitterpos[i1]=NULL;

		free(srf->maxemitter[face]);
		free(srf->nemitter[face]);
		free(srf->emitteramount[face]);
		free(srf->emitterpos[face]);
		srf->maxemitter[face]=newmaxemitter;
		srf->nemitter[face]=newnemitter;
		srf->emitteramount[face]=newemitteramount;
		srf->emitterpos[face]=newemitterpos; }

	for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))						// panel structure data
		for(p=0;p<srf->maxpanel[ps];p++) {
			pnl=srf->panels[ps][p];
			if(!pnl->emitterabsorb[face] || maxspecies>oldmaxspecies) {
				CHECKMEM(newemitterabsorb=(double*) calloc(maxspecies,sizeof(double)));
				i1=0;
				if(maxspecies>oldmaxspecies)
					for(i1=0;i1<oldmaxspecies;i1++) newemitterabsorb[i1]=pnl->emitterabsorb[face][i1];
				for(;i1<maxspecies;i1++) newemitterabsorb[i1]=0;
				free(pnl->emitterabsorb[face]);
				pnl->emitterabsorb[face]=newemitterabsorb; }}

	return 0;

 failure:
	simLog(NULL,10,"Unable to allocate memory in emittersalloc");
	return 1; }


/* surfacealloc */
surfaceptr surfacealloc(surfaceptr srf,int oldmaxspecies,int maxspecies,int dim,int maxsitecode) {
	int i,freesrf;
	enum PanelShape ps;
	enum MolecState ms;
	enum SrfAction ***newaction;
	surfactionptr ***newactdetails;

	if(srf && oldmaxspecies==maxspecies) return srf;

	newaction=NULL;
	newactdetails=NULL;
	freesrf=0;
	
	if(!srf) {
		srf=(surfaceptr) malloc(sizeof(struct surfacestruct));
		if(!srf) return NULL;
		freesrf=1;
		srf->sname=NULL;
		srf->srfss=NULL;
		srf->selfindex=-1;
		srf->action=NULL;
		srf->actdetails=NULL;
		srf->neighhop=0;
		srf->fcolor[0]=srf->fcolor[1]=srf->fcolor[2]=0;
		srf->bcolor[0]=srf->bcolor[1]=srf->bcolor[2]=0;
		srf->fcolor[3]=srf->bcolor[3]=1;
		srf->edgepts=1;
		srf->edgestipple[0]=1;
		srf->edgestipple[1]=0xFFFF;
		srf->fdrawmode=(dim==3?DMface:DMedge);
		srf->bdrawmode=(dim==3?DMface:DMedge);
		srf->fshiny=0;
		srf->bshiny=0;
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) srf->maxpanel[ps]=0;
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) srf->npanel[ps]=0;
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) srf->pname[ps]=NULL;
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) srf->panels[ps]=NULL;
		srf->port[PFfront]=NULL;
		srf->port[PFback]=NULL;
		srf->totarea=0;
		srf->totpanel=0;
		srf->areatable=NULL;
		srf->paneltable=NULL;
		srf->maxemitter[PFfront]=srf->maxemitter[PFback]=NULL;
		srf->nemitter[PFfront]=srf->nemitter[PFback]=NULL;
		srf->emitteramount[PFfront]=srf->emitteramount[PFback]=NULL;
		srf->emitterpos[PFfront]=srf->emitterpos[PFback]=NULL; }
	
	if(maxspecies) {
		CHECKMEM(newaction=(enum SrfAction***) calloc(maxspecies,sizeof(enum SrfAction**)));
		for(i=0;i<maxspecies;i++) newaction[i]=NULL;
		for(i=0;i<oldmaxspecies;i++)
			newaction[i]=srf->action[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newaction[i]=(enum SrfAction**) calloc(MSMAX,sizeof(enum SrfAction*)));
			for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1))
				newaction[i][ms]=NULL;
			for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1)) {
				CHECKMEM(newaction[i][ms]=(enum SrfAction*) calloc(3,sizeof(enum SrfAction)));
				newaction[i][ms][PFfront]=newaction[i][ms][PFback]=SAtrans;
				newaction[i][ms][PFnone]=SAno;
		}}
		
		CHECKMEM(newactdetails=(surfactionptr***) calloc(maxspecies,sizeof(surfactionptr**)));
		for(i=0;i<maxspecies;i++) newactdetails[i]=NULL;
		for(i=0;i<oldmaxspecies;i++)
			newactdetails[i]=srf->actdetails[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newactdetails[i]=(surfactionptr**) calloc(MSMAX,sizeof(surfactionptr*)));
			for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1))
				newactdetails[i][ms]=NULL;
			for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1)) {
				CHECKMEM(newactdetails[i][ms]=(surfactionptr*) calloc(3,sizeof(surfactionptr)));
				newactdetails[i][ms][PFfront]=NULL;
				newactdetails[i][ms][PFback]=NULL;
				newactdetails[i][ms][PFnone]=NULL; }}

		if(srf->maxemitter[PFfront]) {	// these calls are for larger maxspecies
			CHECK(emittersalloc(srf,PFfront,oldmaxspecies,maxspecies)==0); }
		if(srf->maxemitter[PFback]) {
			CHECK(emittersalloc(srf,PFback,oldmaxspecies,maxspecies)==0); }
		
		free(srf->action);
		srf->action=newaction;
		free(srf->actdetails);
		srf->actdetails=newactdetails; }
	
	return srf;
	
failure:
	if(freesrf) surfacefree(srf,maxspecies);
	simLog(NULL,10,"Unable to allocate memory in surfacealloc");
	return NULL; }


/* surfacefree */
void surfacefree(surfaceptr srf,int maxspecies) {
	int p,i,emit;
	enum PanelFace face;
	enum PanelShape ps;
	enum MolecState ms;

	if(!srf) return;

	for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
		if(srf->emitterpos[face]) {
			for(i=0;i<maxspecies;i++)
				if(srf->emitterpos[face][i]) {
					for(emit=0;emit<srf->maxemitter[face][i];emit++)
						free(srf->emitterpos[face][i][emit]);
					free(srf->emitterpos[face][i]); }
				free(srf->emitterpos[face]); }
		if(srf->emitteramount[face]) {
			for(i=0;i<maxspecies;i++)
				free(srf->emitteramount[face][i]);
			free(srf->emitteramount[face]); }
		free(srf->nemitter[face]);
		free(srf->maxemitter[face]); }
	
	free(srf->paneltable);
	free(srf->areatable);
	
	for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) {
		for(p=0;p<srf->maxpanel[ps];p++) {
			if(srf->panels[ps]) panelfree(srf->panels[ps][p]);
			if(srf->pname[ps]) free(srf->pname[ps][p]); }
		free(srf->pname[ps]);
		free(srf->panels[ps]); }
	
	for(i=0;i<maxspecies;i++)
		if(srf->actdetails[i]) {
			for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1))
				if(srf->actdetails[i][ms]) {
					for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1))
						surfaceactionfree(srf->actdetails[i][ms][face]);
					free(srf->actdetails[i][ms]); }
			free(srf->actdetails[i]); }
	free(srf->actdetails);
	
	for(i=0;i<maxspecies;i++) {
		if(srf->action[i]) {
			for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1))
				free(srf->action[i][ms]);
			free(srf->action[i]); 
		}
	}
	free(srf->action);
	
	free(srf);
	return; }


/* surfacessalloc */
surfacessptr surfacessalloc(surfacessptr srfss,int maxsurface,int maxspecies,int dim,int maxsitecode) {
	int s,newsrfss,oldmaxsrf;
	char **newnames;
	surfaceptr *newsrflist;

	if(maxsurface<1 || maxspecies<0) return NULL;

	newsrfss=0;
	newnames=NULL;
	newsrflist=NULL;

	if(!srfss) {												// new allocation
		srfss=(surfacessptr) malloc(sizeof(struct surfacesuperstruct));
		if(!srfss) return NULL;
		newsrfss=1;
		srfss->condition=SCinit;
		srfss->sim=NULL;
		srfss->maxspecies=maxspecies;
		srfss->maxsrf=0;
		srfss->nsrf=0;
		srfss->epsilon=100*DBL_EPSILON;
		srfss->margin=100*DBL_EPSILON;
		srfss->neighdist=1000*DBL_EPSILON;
		srfss->snames=NULL;
		srfss->srflist=NULL;
		srfss->maxmollist=0;
		srfss->nmollist=0;
		srfss->srfmollist=NULL; }
	else {																// checks, and update maxspecies if reallocation
		if(maxsurface<srfss->maxsrf) return NULL;
		if(maxspecies<srfss->maxspecies) return NULL;
		if(maxspecies>srfss->maxspecies) {
			for(s=0;s<srfss->maxsrf;s++) {
				CHECK(surfacealloc(srfss->srflist[s],srfss->maxspecies,maxspecies,dim,maxsitecode)!=NULL);
				srfss->srflist[s]->selfindex=s; }}
		srfss->maxspecies=maxspecies; }

	if(maxsurface>srfss->maxsrf) {			// allocate any new surface names and surfaces
		oldmaxsrf=srfss->maxsrf;
		CHECKMEM(newnames=(char**) calloc(maxsurface,sizeof(char*)));		// surface names
		for(s=0;s<maxsurface;s++) newnames[s]=NULL;
		for(s=0;s<oldmaxsrf;s++) newnames[s]=srfss->snames[s];
		for(;s<maxsurface;s++)
			CHECKMEM(newnames[s]=EmptyString());

		CHECKMEM(newsrflist=(surfaceptr*) calloc(maxsurface,sizeof(surfaceptr)));	// surface list
		for(s=0;s<maxsurface;s++) newsrflist[s]=NULL;
		for(s=0;s<oldmaxsrf;s++)
			newsrflist[s]=srfss->srflist[s];
		for(;s<maxsurface;s++) {
			CHECK(newsrflist[s]=surfacealloc(NULL,0,maxspecies,dim,maxsitecode));
			newsrflist[s]->srfss=srfss;
			newsrflist[s]->sname=newnames[s];
			newsrflist[s]->selfindex=s; }

		srfss->maxsrf=maxsurface;
		free(srfss->snames);
		srfss->snames=newnames;
		free(srfss->srflist);
		srfss->srflist=newsrflist;

		if(srfss->sim && srfss->sim->mols && srfss->sim->mols->surfdrift) {
			CHECK(molexpandsurfdrift(srfss->sim,srfss->sim->mols->maxspecies,oldmaxsrf)==0); }}
	
	return srfss;

 failure:
 	if(newsrfss) surfacessfree(srfss);
	if(ErrorType!=1) simLog(NULL,10,"Unable to allocate memory in surfacessalloc");
 	return NULL; }


/* surfacessfree */
void surfacessfree(surfacessptr srfss) {
	int s;

	if(!srfss) return;

	free(srfss->srfmollist);
	if(srfss->srflist) {
		for(s=0;s<srfss->maxsrf;s++)
			surfacefree(srfss->srflist[s],srfss->maxspecies);
		free(srfss->srflist); }

	if(srfss->snames) {
		for(s=0;s<srfss->maxsrf;s++)
			free(srfss->snames[s]);
		free(srfss->snames); }

	free(srfss);
	return; }


/******************************************************************************/
/**************************** data structure output ***************************/
/******************************************************************************/

/* surfaceoutput */
void surfaceoutput(simptr sim) {
	int s,p,i,dim,nspecies,jumpfrnt,jumpback,ll,p2,emit,d,same,none;
	surfacessptr srfss;
	surfaceptr srf;
	double **point,*front,prob;
	char **pname,string[STRCHAR];
	panelptr pnl;
	enum SrfAction act,***action;
	enum MolecState ms,ms2,ms3,ms4,ms5;
	enum PanelFace face;
	surfactionptr actdetails;

	simLog(sim,2,"SURFACE PARAMETERS\n");
	srfss=sim->srfss;
	dim=sim->dim;
	nspecies=sim->mols?sim->mols->nspecies:0;
	if(!srfss) {
		simLog(sim,2," No internal surfaces\n\n");
		return; }

	simLog(sim,1," Allocated for %i species\n",srfss->maxspecies-1);

	simLog(sim,2," Surface epsilon, margin, and neighbor distances: %g %g %g\n",srfss->epsilon,srfss->margin,srfss->neighdist);

	if(sim->mols) {
		simLog(sim,2," Molecule lists checked after diffusion:");
		for(ll=0;ll<srfss->nmollist;ll++)
			if(srfss->srfmollist[ll]&SMLdiffuse) simLog(sim,2," %s",sim->mols->listname[ll]);
		simLog(sim,2,"\n");
		simLog(sim,2," Molecule lists checked after reactions:");
		for(ll=0;ll<srfss->nmollist;ll++)
			if(srfss->srfmollist[ll]&SMLreact) simLog(sim,2," %s",sim->mols->listname[ll]);
		simLog(sim,2,"\n");
		simLog(sim,2," Molecule lists checked for surface-bound molecules:");
		for(ll=0;ll<srfss->nmollist;ll++)
			if(srfss->srfmollist[ll]&SMLsrfbound) simLog(sim,2," %s",sim->mols->listname[ll]);
		simLog(sim,2,"\n"); }

	simLog(sim,2," Surfaces allocated: %i, surfaces defined: %i\n\n",srfss->maxsrf,srfss->nsrf);
	for(s=0;s<srfss->nsrf;s++) {
		srf=srfss->srflist[s];
		simLog(sim,2," Surface: %s\n",srfss->snames[s]);
		if(srf->port[PFfront]) simLog(sim,2,"  The front of this surface is part of port %s\n",srf->port[PFfront]->portname);
		if(srf->port[PFback]) simLog(sim,2,"  The back of this surface is part of port %s\n",srf->port[PFback]->portname);

		if(sim->mols && srf->action) {
			simLog(sim,2,"  actions for molecules:\n");
			action=srf->action;
			for(i=1;i<nspecies;i++) {
				for(face=(PanelFace)0;face<2;face=(PanelFace)(face+1)) {
					same=1;
					act=action[i][MSsoln][face];
					for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1))
						if(action[i][ms][face]!=act) same=0;
						if(same) {
							if(act!=SAmult)
								simLog(sim,2,"   %s(all) at %s: %s\n",sim->mols->spname[i],face==PFfront?"front":"back",surfact2string(act,string)); }
						else {
							for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1)) {
								act=action[i][ms][face];
								if(sim->mols->exist[i][ms] && act!=SAmult) {
									simLog(sim,2,"   %s(%s)",sim->mols->spname[i],molms2string(ms,string));
									simLog(sim,2," at %s: %s\n",face==PFfront?"front":"back",surfact2string(act,string)); }}}}}

			simLog(sim,2,"  rates for molecules:");
			none=1;
			for(i=1;i<nspecies;i++) {
				for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1))
					for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1))
						if(action[i][ms][face]==SAmult && srf->actdetails) {
							actdetails=srf->actdetails[i][ms][face];
							for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1)) {
								prob=actdetails->srfprob[ms2];
								if(prob>0 && !srfsamestate(ms,face,ms2,NULL)) {
									if(none) {
										none=0;
										simLog(sim,2,"\n"); }
									srfindex2tristate(ms,face,ms2,&ms3,&ms4,&ms5);
									simLog(sim,2,"    %s(%s)",sim->mols->spname[i],molms2string(ms3,string));
									simLog(sim,2," %s ->",molms2string(ms4,string));
									simLog(sim,2," %s",molms2string(ms5,string));
									if(actdetails->srfnewspec[ms2]!=i) simLog(sim,2," (convert to %s)",sim->mols->spname[actdetails->srfnewspec[ms2]]);
									if(actdetails->srfdatasrc[ms2]==1) simLog(sim,2,", requested rate=%g",actdetails->srfrate[ms2]);
									else if(actdetails->srfdatasrc[ms2]==2) simLog(sim,2,", requested prob=%g",actdetails->srfprob[ms2]);
									simLog(sim,2,", actual rate=%g, prob=%g\n",srfcalcrate(sim,srf,i,ms,face,ms2),prob); }}}}
		if(none) simLog(sim,2," no stochastic interactions\n"); }

		if(sim->mols) {
			simLog(sim,2,"  surface-bound molecules %s hop between neighboring panels upon collision\n",srf->neighhop?"do":"do not"); }

		simLog(sim,2,"  front color: %g %g %g %g\n",srf->fcolor[0],srf->fcolor[1],srf->fcolor[2],srf->fcolor[3]);
		simLog(sim,2,"  back color: %g %g %g %g\n",srf->bcolor[0],srf->bcolor[1],srf->bcolor[2],srf->bcolor[3]);
		simLog(sim,2,"  edge points: %g, polygon modes: %s %s\n",srf->edgepts,surfdm2string(srf->fdrawmode,string),surfdm2string(srf->bdrawmode,string));
		if(srf->edgestipple[1]!=0xFFFF) simLog(sim,2,"   edge stippling: %ui %X\n",srf->edgestipple[0],srf->edgestipple[1]);
		if(srf->fshiny!=0) simLog(sim,2,"  front shininess: %g\n",srf->fshiny);
		if(srf->bshiny!=0) simLog(sim,2,"  back shininess: %g\n",srf->bshiny);
		
		for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
			if(srf->nemitter[face])
				for(i=0;i<nspecies;i++)
					if(srf->nemitter[face][i]) {
						simLog(sim,2,"  %i %s-side emitters defined, of %i allocated:\n",srf->nemitter[face][i],surfface2string(face,string),srf->maxemitter[face][i]);
						for(emit=0;emit<srf->nemitter[face][i];emit++) {
							simLog(sim,2,"   %g at (%g",srf->emitteramount[face][i][emit],srf->emitterpos[face][i][emit][0]);
							for(d=1;d<dim;d++)
								simLog(sim,2,",%g",srf->emitterpos[face][i][emit][d]);
							simLog(sim,2,")\n"); }}}
		
		jumpfrnt=jumpback=0;
		if(srf->action) {
			for(i=1;i<nspecies;i++) {
				if(srf->action[i][MSsoln][PFfront]==SAjump) jumpfrnt=1;
				if(srf->action[i][MSsoln][PFback]==SAjump) jumpback=1; 
		}}
		
		if(srf->maxpanel[PSrect]) {
			pname=srf->pname[PSrect];
			simLog(sim,2,"  rectangle panels allocated: %i, defined: %i\n",srf->maxpanel[PSrect],srf->npanel[PSrect]);
			for(p=0;p<srf->npanel[PSrect] && p<21;p++) {
				if(p==20) {
					simLog(sim,2,"   ...\n");
					p=srf->npanel[PSrect]-1; }
				pnl=srf->panels[PSrect][p];
				point=pnl->point;
				front=pnl->front;
				if(dim==1) simLog(sim,2,"   %s: %g, facing %c0",pname[p],point[0][0],front[0]==1?'+':'-');
				else if(dim==2) simLog(sim,2,"   %s: (%g,%g), (%g,%g), facing: %c%1.0f, length: %g",pname[p],point[0][0],point[0][1],point[1][0],point[1][1],front[0]==1?'+':'-',front[1],Geo_LineLength(point[0],point[1],2));
				else simLog(sim,2,"   %s: (%g,%g,%g), (%g,%g,%g), (%g,%g,%g), (%g,%g,%g), facing: %c%1.0f, area: %g",pname[p],point[0][0],point[0][1],point[0][2],point[1][0],point[1][1],point[1][2],point[2][0],point[2][1],point[2][2],point[3][0],point[3][1],point[3][2],front[0]==1?'+':'-',front[1],Geo_QuadArea(point[0],point[1],point[2],point[3],3));
				if(jumpfrnt && pnl->jumpp[PFfront]) simLog(sim,2,"; front jump: %s, %s",pnl->jumpp[PFfront]->pname,surfface2string(pnl->jumpf[PFfront],string));
				else if(jumpfrnt) simLog(sim,2,"; front jump: NO PANEL");
				if(jumpback && pnl->jumpp[PFback]) simLog(sim,2,"; back jump: %s, %s",pnl->jumpp[PFback]->pname,surfface2string(pnl->jumpf[PFback],string));
				else if(jumpback) simLog(sim,2,"; back jump: NO PANEL");
				for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
					if(pnl->emitterabsorb[face]) {
						simLog(sim,2,"; %s absorb probs.:",surfface2string(face,string));
						for(i=1;i<nspecies;i++)
							simLog(sim,2," %g",pnl->emitterabsorb[face][i]); }}
				simLog(sim,2,"\n");
				if(pnl->maxneigh) {
					simLog(sim,1,"    neighbors allocated: %i, defined: %i",pnl->maxneigh,pnl->nneigh);
					for(p2=0;p2<pnl->nneigh;p2++)
						simLog(sim,1,", %s",pnl->neigh[p2]->pname);
					simLog(sim,1,"\n"); }}}
		if(srf->maxpanel[PStri]) {
			pname=srf->pname[PStri];
			simLog(sim,2,"  triangle panels allocated: %i, defined: %i\n",srf->maxpanel[PStri],srf->npanel[PStri]);
			for(p=0;p<srf->npanel[PStri] && p<21;p++) {
				if(p==20) {
					simLog(sim,2,"   ...\n");
					p=srf->npanel[PStri]-1; }
				pnl=srf->panels[PStri][p];
				point=pnl->point;
				front=pnl->front;
				if(dim==1) simLog(sim,2,"   %s: %g, facing %c0",pname[p],point[0][0],front[0]==1?'+':'-');
				else if(dim==2) simLog(sim,2,"   %s: (%g,%g), (%g,%g), facing: (%g,%g), length: %g",pname[p],point[0][0],point[0][1],point[1][0],point[1][1],front[0],front[1],Geo_LineLength(point[0],point[1],2));
				else simLog(sim,2,"   %s: (%g,%g,%g), (%g,%g,%g), (%g,%g,%g), facing: (%g,%g,%g), area: %g",pname[p],point[0][0],point[0][1],point[0][2],point[1][0],point[1][1],point[1][2],point[2][0],point[2][1],point[2][2],front[0],front[1],front[2],Geo_TriArea3(point[0],point[1],point[2],front));
				if(jumpfrnt && pnl->jumpp[PFfront]) simLog(sim,2,"; front jump: %s, %s",pnl->jumpp[PFfront]->pname,surfface2string(pnl->jumpf[PFfront],string));
				else if(jumpfrnt) simLog(sim,2,"; front jump: NO PANEL");
				if(jumpback && pnl->jumpp[PFback]) simLog(sim,2,"; back jump: %s, %s",pnl->jumpp[PFback]->pname,surfface2string(pnl->jumpf[PFback],string));
				else if(jumpback) simLog(sim,2,"; back jump: NO PANEL");
				for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
					if(pnl->emitterabsorb[face]) {
						simLog(sim,2,"; %s absorb probs.:",surfface2string(face,string));
						for(i=1;i<nspecies;i++)
							simLog(sim,2," %g",pnl->emitterabsorb[face][i]); }}
				simLog(sim,2,"\n");
				if(pnl->maxneigh) {
					simLog(sim,1,"    neighbors allocated: %i, defined: %i",pnl->maxneigh,pnl->nneigh);
					for(p2=0;p2<pnl->nneigh;p2++)
						simLog(sim,1,", %s",pnl->neigh[p2]->pname);
					simLog(sim,1,"\n"); }}}
		if(srf->maxpanel[PSsph]) {
			pname=srf->pname[PSsph];
			simLog(sim,2,"  sphere panels allocated: %i, defined: %i\n",srf->maxpanel[PSsph],srf->npanel[PSsph]);
			for(p=0;p<srf->npanel[PSsph] && p<21;p++) {
				if(p==20) {
					simLog(sim,2,"   ...\n");
					p=srf->npanel[PSsph]-1; }
				pnl=srf->panels[PSsph][p];
				point=pnl->point;
				front=pnl->front;
				if(dim==1) simLog(sim,2,"   %s: %g, R=%g, facing: %s",pname[p],point[0][0],point[1][0],front[0]==1?"out":"in");
				else if(dim==2) simLog(sim,2,"   %s: (%g,%g), R=%g, facing: %s, draw: %g, length: %g",pname[p],point[0][0],point[0][1],point[1][0],front[0]==1?"out":"in",point[1][1],2.0*PI*point[1][0]);
				else simLog(sim,2,"   %s: (%g,%g,%g), R=%g, facing: %s, draw: %g %g, area: %g",pname[p],point[0][0],point[0][1],point[0][2],point[1][0],front[0]==1?"out":"in",point[1][1],point[1][2],4.0*PI*point[1][0]*point[1][0]);
				if(jumpfrnt && pnl->jumpp[PFfront]) simLog(sim,2,"; front jump: %s, %s",pnl->jumpp[PFfront]->pname,surfface2string(pnl->jumpf[PFfront],string));
				else if(jumpfrnt) simLog(sim,2,"; front jump: NO PANEL");
				if(jumpback && pnl->jumpp[PFback]) simLog(sim,2,"; back jump: %s, %s",pnl->jumpp[PFback]->pname,surfface2string(pnl->jumpf[PFback],string));
				else if(jumpback) simLog(sim,2,"; back jump: NO PANEL");
				for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
					if(pnl->emitterabsorb[face]) {
						simLog(sim,2,"; %s absorb probs.:",surfface2string(face,string));
						for(i=1;i<nspecies;i++)
							simLog(sim,2," %g",pnl->emitterabsorb[face][i]); }}
				simLog(sim,2,"\n");
				if(pnl->maxneigh) {
					simLog(sim,1,"    neighbors allocated: %i, defined: %i",pnl->maxneigh,pnl->nneigh);
					for(p2=0;p2<pnl->nneigh;p2++)
						simLog(sim,1,", %s",pnl->neigh[p2]->pname);
					simLog(sim,1,"\n"); }}}
		if(srf->maxpanel[PScyl]) {
			pname=srf->pname[PScyl];
			simLog(sim,2,"  cylinder panels allocated: %i, defined: %i\n",srf->maxpanel[PScyl],srf->npanel[PScyl]);
			for(p=0;p<srf->npanel[PScyl] && p<21;p++) {
				if(p==20) {
					simLog(sim,2,"   ...\n");
					p=srf->npanel[PScyl]-1; }
				pnl=srf->panels[PScyl][p];
				point=pnl->point;
				front=pnl->front;
				if(dim==1) simLog(sim,10,"   error, cylinders are not implemented in 1-D");
				else if(dim==2) simLog(sim,2,"   %s: (%g,%g) to (%g,%g), R=%g, facing: %s, length: %g",pname[p],point[0][0],point[0][1],point[1][0],point[1][1],point[2][0],front[2]==1?"out":"in",2.0*Geo_LineLength(point[0],point[1],2));
				else simLog(sim,2,"   %s: (%g,%g,%g) to (%g,%g,%g), R=%g, facing: %s, draw: %g %g, area: %g",pname[p],point[0][0],point[0][1],point[0][2],point[1][0],point[1][1],point[1][2],point[2][0],front[2]==1?"out":"in",point[2][1],point[2][2],2.0*PI*point[2][0]*Geo_LineLength(point[0],point[1],3));
				if(jumpfrnt && pnl->jumpp[PFfront]) simLog(sim,2,"; front jump: %s, %s",pnl->jumpp[PFfront]->pname,surfface2string(pnl->jumpf[PFfront],string));
				else if(jumpfrnt) simLog(sim,2,"; front jump: NO PANEL");
				if(jumpback && pnl->jumpp[PFback]) simLog(sim,2,"; back jump: %s, %s",pnl->jumpp[PFback]->pname,surfface2string(pnl->jumpf[PFback],string));
				else if(jumpback) simLog(sim,2,"; back jump: NO PANEL");
				for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
					if(pnl->emitterabsorb[face]) {
						simLog(sim,2,"; %s absorb probs.:",surfface2string(face,string));
						for(i=1;i<nspecies;i++)
							simLog(sim,2," %g",pnl->emitterabsorb[face][i]); }}
				simLog(sim,2,"\n");
				if(pnl->maxneigh) {
					simLog(sim,1,"    neighbors allocated: %i, defined: %i",pnl->maxneigh,pnl->nneigh);
					for(p2=0;p2<pnl->nneigh;p2++)
						simLog(sim,1,", %s",pnl->neigh[p2]->pname);
					simLog(sim,1,"\n"); }}}
		if(srf->maxpanel[PShemi]) {
			pname=srf->pname[PShemi];
			simLog(sim,2,"  hemisphere panels allocated: %i, defined: %i\n",srf->maxpanel[PShemi],srf->npanel[PShemi]);
			for(p=0;p<srf->npanel[PShemi] && p<21;p++) {
				if(p==20) {
					simLog(sim,2,"   ...\n");
					p=srf->npanel[PShemi]-1; }
				pnl=srf->panels[PShemi][p];
				point=pnl->point;
				front=pnl->front;
				if(dim==1) simLog(sim,2,"   error, hemispheres are not implemented in 1-D");
				else if(dim==2) simLog(sim,2,"   %s: (%g,%g), R=%g, facing: %s, opening: (%g,%g), draw: %g, length: %g",pname[p],point[0][0],point[0][1],point[1][0],front[0]==1?"out":"in",point[2][0],point[2][1],point[1][1],PI*point[1][0]);
				else simLog(sim,2,"   %s: (%g,%g,%g), R=%g, facing: %s, opening: (%g,%g,%g), draw: %g %g, area: %g",pname[p],point[0][0],point[0][1],point[0][2],point[1][0],front[0]==1?"out":"in",point[2][0],point[2][1],point[2][2],point[1][1],point[1][2],2.0*PI*point[1][0]*point[1][0]);
				if(jumpfrnt && pnl->jumpp[PFfront]) simLog(sim,2,"; front jump: %s, %s",pnl->jumpp[PFfront]->pname,surfface2string(pnl->jumpf[PFfront],string));
				else if(jumpfrnt) simLog(sim,2,"; front jump: NO PANEL");
				if(jumpback && pnl->jumpp[PFback]) simLog(sim,2,"; back jump: %s, %s",pnl->jumpp[PFback]->pname,surfface2string(pnl->jumpf[PFback],string));
				else if(jumpback) simLog(sim,2,"; back jump: NO PANEL");
				for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
					if(pnl->emitterabsorb[face]) {
						simLog(sim,2,"; %s absorb probs.:",surfface2string(face,string));
						for(i=1;i<nspecies;i++)
							simLog(sim,2," %g",pnl->emitterabsorb[face][i]); }}
				simLog(sim,2,"\n");
				if(pnl->maxneigh) {
					simLog(sim,1,"    neighbors allocated: %i, defined: %i",pnl->maxneigh,pnl->nneigh);
					for(p2=0;p2<pnl->nneigh;p2++)
						simLog(sim,1,", %s",pnl->neigh[p2]->pname);
					simLog(sim,1,"\n"); }}}
		if(srf->maxpanel[PSdisk]) {
			pname=srf->pname[PSdisk];
			simLog(sim,2,"  disk panels allocated: %i, defined: %i\n",srf->maxpanel[PSdisk],srf->npanel[PSdisk]);
			for(p=0;p<srf->npanel[PSdisk] && p<21;p++) {
				if(p==20) {
					simLog(sim,2,"   ...\n");
					p=srf->npanel[PSdisk]-1; }
				pnl=srf->panels[PSdisk][p];
				point=pnl->point;
				front=pnl->front;
				if(dim==1) simLog(sim,10,"   error, disks are not implemented in 1-D");
				else if(dim==2) simLog(sim,2,"   %s: (%g,%g), R=%g, facing: (%g,%g), length: %g",pname[p],point[0][0],point[0][1],point[1][0],front[0],front[1],2.0*point[1][0]);
				else simLog(sim,2,"   %s: (%g,%g,%g), R=%g, facing: (%g,%g,%g), draw: %g, area: %g",pname[p],point[0][0],point[0][1],point[0][2],point[1][0],front[0],front[1],front[2],point[1][1],PI*point[1][0]*point[1][0]);
				if(jumpfrnt && pnl->jumpp[PFfront]) simLog(sim,2,"; front jump: %s, %s",pnl->jumpp[PFfront]->pname,surfface2string(pnl->jumpf[PFfront],string));
				else if(jumpfrnt) simLog(sim,2,"; front jump: NO PANEL");
				if(jumpback && pnl->jumpp[PFback]) simLog(sim,2,"; back jump: %s, %s",pnl->jumpp[PFback]->pname,surfface2string(pnl->jumpf[PFback],string));
				else if(jumpback) simLog(sim,2,"; back jump: NO PANEL");
				for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
					if(pnl->emitterabsorb[face]) {
						simLog(sim,2,"; %s absorb probs.:",surfface2string(face,string));
						for(i=1;i<nspecies;i++)
							simLog(sim,2," %g",pnl->emitterabsorb[face][i]); }}
				simLog(sim,2,"\n");
				if(pnl->maxneigh) {
					simLog(sim,1,"    neighbors allocated: %i, defined: %i",pnl->maxneigh,pnl->nneigh);
					for(p2=0;p2<pnl->nneigh;p2++)
						simLog(sim,1,", %s",pnl->neigh[p2]->pname);
					simLog(sim,1,"\n"); }}}

		simLog(sim,2,"\n"); }
	return; }


/* writesurfaces */
void writesurfaces(simptr sim,FILE *fptr) {
	int s,i,d,d2,j,p,dim,nspecies;
	surfacessptr srfss;
	surfaceptr srf;
	enum MolecState ms1,ms2;
	enum PanelFace face;
	char string[STRCHAR];
	enum PanelShape ps;
	panelptr pnl;
	double **point,*front;
	surfactionptr actdetails;
	enum SrfAction act;

	dim=sim->dim;
	nspecies=sim->mols?sim->mols->nspecies:0;
	fprintf(fptr,"# Surface parameters\n");
	if(!sim->srfss) return;
	srfss=sim->srfss;
	fprintf(fptr,"max_surface %i\n",srfss->maxsrf);
	fprintf(fptr,"epsilon %g\n",srfss->epsilon);
	fprintf(fptr,"margin %g\n",srfss->margin);
	fprintf(fptr,"neighbor_dist %g\n",srfss->neighdist);
	fprintf(fptr,"\n");
	for(s=0;s<srfss->nsrf;s++) {
		srf=srfss->srflist[s];
		fprintf(fptr,"start_surface %s\n",srfss->snames[s]);
		if(srf->action) {
			for(i=1;i<nspecies;i++) {
				for(ms1=(MolecState)0;ms1<MSMAX;ms1=(MolecState)(ms1+1))
					for(face=(PanelFace)0;face<2;face=(PanelFace)(face+1)) {
						act=srf->action[i][ms1][face];
						if(act!=SAtrans && act!=SAmult) {
							fprintf(fptr,"action %s(%s)",sim->mols->spname[i],molms2string(ms1,string));
							fprintf(fptr," %s %s\n",face==PFfront?"front":"back",surfact2string(act,string)); }}}}

		if(srf->action && srf->actdetails) {
			for(i=1;i<nspecies;i++) {
				for(ms1=(MolecState)0;ms1<MSMAX;ms1=(MolecState)(ms1+1))
					for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1))
						if(srf->action[i][ms1][face]==SAmult && srf->actdetails[i][ms1][face]) {
							actdetails=srf->actdetails[i][ms1][face];
							for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1)) {
								if(actdetails->srfrate[ms2]>0) {
									fprintf(fptr,"rate %s(%s)",sim->mols->spname[i],molms2string(ms1,string));
									fprintf(fptr," %s",(face==PFnone)?molms2string(ms1,string):(face==PFfront?"fsoln":"bsoln"));
									fprintf(fptr," %s %g",molms2string(ms2,string),actdetails->srfrate[ms2]);
									if(actdetails->srfnewspec[ms2]!=i)
										fprintf(fptr," %s",sim->mols->spname[actdetails->srfnewspec[ms2]]);
										fprintf(fptr,"\n"); }
								else if(actdetails->srfrate[ms2]<0 && actdetails->srfprob[ms2]>0) {
									fprintf(fptr,"rate_internal %s(%s)",sim->mols->spname[i],molms2string(ms1,string));
									fprintf(fptr," %s",(face==PFnone)?molms2string(ms1,string):(face==PFfront?"fsoln":"bsoln"));
									fprintf(fptr," %s %g",molms2string(ms2,string),actdetails->srfprob[ms2]);
									if(actdetails->srfnewspec[ms2]!=i)
										fprintf(fptr," %s",sim->mols->spname[actdetails->srfnewspec[ms2]]);
										fprintf(fptr,"\n"); }}}}}

		if(srf->neighhop) {
			fprintf(fptr,"neighbor_action hop\n"); }
		
		fprintf(fptr,"color front %g %g %g %g\n",srf->fcolor[0],srf->fcolor[1],srf->fcolor[2],srf->fcolor[3]);
		fprintf(fptr,"color back %g %g %g %g\n",srf->bcolor[0],srf->bcolor[1],srf->bcolor[2],srf->bcolor[3]);
		if(srf->fshiny!=0) fprintf(fptr,"shininess front %g\n",srf->fshiny);
		if(srf->bshiny!=0) fprintf(fptr,"shininess back %g\n",srf->bshiny);
		fprintf(fptr,"thickness %g\n",srf->edgepts);
		if(srf->edgestipple[1]!=0xFFFF)
			printf("   stipple %ui 0x%X\n",srf->edgestipple[0],srf->edgestipple[1]);
		fprintf(fptr,"polygon front %s\n",surfdm2string(srf->fdrawmode,string));
		fprintf(fptr,"polygon back %s\n",surfdm2string(srf->bdrawmode,string));

		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) {
			if(srf->maxpanel[ps]) fprintf(fptr,"max_panels %s %i\n",surfps2string(ps,string),srf->maxpanel[ps]);
			for(p=0;p<srf->npanel[ps];p++) {
				pnl=srf->panels[ps][p];
				point=pnl->point;
				front=pnl->front;
				fprintf(fptr,"panel %s",surfps2string(ps,string));
				if(ps==PSrect) {
					fprintf(fptr," %c%i",pnl->front[0]>0?'+':'-',(int)(front[1]));
					for(d=0;d<dim;d++) fprintf(fptr," %g",point[0][d]);
					for(d=0;d<dim;d++)
						if(d!=pnl->front[1]) fprintf(fptr," %g",point[dim==2?1:2][d]-point[0][d]); }
				else if(ps==PStri) {
					for(d2=0;d2<dim;d2++)
						for(d=0;d<dim;d++)
							fprintf(fptr," %g",point[d2][d]); }
				else if(ps==PSsph) {
					for(d=0;d<dim;d++) fprintf(fptr," %g",point[0][d]);
					fprintf(fptr," %c%g",front[0]>0?'+':'-',point[1][0]);
					if(dim>1) fprintf(fptr," %g",point[1][1]);
					if(dim>2) fprintf(fptr," %g",point[1][2]); }
				else if(ps==PScyl) {
					for(d=0;d<dim;d++) fprintf(fptr," %g",point[0][d]);
					for(d=0;d<dim;d++) fprintf(fptr," %g",point[1][d]);
					fprintf(fptr," %c%g",front[dim==2?2:0]>0?'+':'-',point[2][0]);
					if(dim>2) fprintf(fptr," %g %g",point[2][1],point[2][2]); }
				else if(ps==PShemi) {
					for(d=0;d<dim;d++) fprintf(fptr," %g",point[0][d]);
					fprintf(fptr," %c%g",front[0]>0?'+':'-',point[1][0]);
					for(d=0;d<dim;d++) fprintf(fptr," %g",point[2][d]);
					if(dim>=2) fprintf(fptr," %g",point[1][1]);
					if(dim>2) fprintf(fptr," %g",point[1][2]); }
				else if(ps==PSdisk) {
					for(d=0;d<dim;d++) fprintf(fptr," %g",point[0][d]);
					fprintf(fptr," %g",point[1][0]);
					for(d=0;d<dim;d++) fprintf(fptr," %g",front[d]);
					if(dim>2) fprintf(fptr," %g",point[1][1]); }
				fprintf(fptr," %s\n",pnl->pname); }}

		if(srf->action) {
			for(i=1;i<nspecies;i++)
				for(face=(PanelFace)0;face<2;face=(PanelFace)(face+1))
					if(srf->action[i][MSsoln][face]==SAjump)
						for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
							for(p=0;p<srf->npanel[ps];p++) {
								pnl=srf->panels[ps][p];
								if(pnl->jumpp[face] && (pnl->jumpf[face]==PFfront || pnl->jumpf[face]==PFback))
									fprintf(fptr,"jump %s %s -> %s %s\n",pnl->pname,surfface2string(face,string),pnl->jumpp[face]->pname,surfface2string(pnl->jumpf[face],string)); }}

		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
			for(p=0;p<srf->npanel[ps];p++) {
				pnl=srf->panels[ps][p];
				if(pnl->nneigh) {
					fprintf(fptr,"neighbors %s",pnl->pname);
					for(j=0;j<pnl->nneigh;j++) {
						if(pnl->neigh[j]->srf==srf)
							fprintf(fptr," %s",pnl->neigh[j]->pname);
						else
							fprintf(fptr," %s:%s",pnl->neigh[j]->srf->sname,pnl->neigh[j]->pname); }
					fprintf(fptr,"\n"); }}
		fprintf(fptr,"end_surface\n\n"); }
	return; }


/* checksurfaceparams */
int checksurfaceparams(simptr sim,int *warnptr) {
	int error,warn;
	int d,c,s,i,dim,nspecies,p,pt,fjump,bjump,fport,bport;
	surfacessptr srfss;
	surfaceptr srf;
	enum MolecState ms,ms2;
	double syslen,**point,*front,lowwall[3],highwall[3],norm[3],prob;
	enum PanelShape ps;
	enum PanelFace face;
	panelptr pnl;
	char string[STRCHAR],string2[STRCHAR];
	surfactionptr actdetails;

	error=warn=0;

	dim=sim->dim;
	nspecies=sim->mols?sim->mols->nspecies:0;
	systemcorners(sim,lowwall,highwall);
	syslen=0;
	for(d=0;d<dim;d++) syslen+=(highwall[d]-lowwall[d])*(highwall[d]-lowwall[d]);
	syslen=sqrt(syslen);

	srfss=sim->srfss;
	if(!srfss) {
		if(warnptr) *warnptr=warn;
		return 0; }

	if(srfss->condition!=SCok) {											// condition
		warn++;
		simLog(sim,7," WARNING: surface structure %s\n",simsc2string(srfss->condition,string)); }

	if(sim->mols)																			// maxspecies
		if(sim->mols->maxspecies!=sim->srfss->maxspecies) {error++;simLog(sim,10," SMOLDYN BUG: number of molecule species differ between mols and srfss\n");}

	if(sim->mols)																			// nmollist
		if(srfss->nmollist!=sim->mols->nlist) {error++;simLog(sim,10," SMOLDYN BUG: mismatch between number of molecule lists in surface and in molecule structures");}
																										// possible superstructure warnings, errors
	if(srfss->nsrf==0) {warn++;simLog(sim,5," WARNING: Surface superstructure is defined, but no surfaces are defined\n");}
	if(srfss->nsrf>srfss->maxsrf) {error++;simLog(sim,10," SMOLDYN BUG: More surfaces are defined than allocated\n");}
	if(srfss->epsilon>0.01*syslen) {warn++;simLog(sim,5," WARNING: surface epsilon value is large compared to system size\n");}
	if(srfss->margin>0.01*syslen) {warn++;simLog(sim,5," WARNING: surface margin value is large compared to system size\n");}
	if(srfss->neighdist>0.01*syslen) {warn++;simLog(sim,5," WARNING: surface neighbor distance value is large compared to system size\n");}

	for(s=0;s<srfss->nsrf;s++) {											// surface missing and incorrect elements, not panels
		if(strlen(srfss->snames[s])==0) {error++;simLog(sim,9," ERROR: surface %i is unnamed\n",s);}
		if(!srfss->srflist[s]) {error++;simLog(sim,10," SMOLDYN BUG: pointer to surface %i is NULL",s);}
		srf=srfss->srflist[s];
		if(srf->sname!=srfss->snames[s]) {error++;simLog(sim,10," SMOLDYN BUG: surface %i name pointer does not match surface superstructure pointer\n",s);}
		for(c=0;c<4;c++) {
			if(srf->fcolor[c]<0 || srf->fcolor[c]>1) {error++;simLog(sim,10," SMOLDYN BUG: surface %i front color %i value is out of bounds\n",s,c);}
			if(srf->bcolor[c]<0 || srf->bcolor[c]>1) {error++;simLog(sim,10," SMOLDYN BUG: surface %i front color %i value is out of bounds\n",s,c);} }
		if(srf->edgepts<0) {error++;simLog(sim,10," SMOLDYN BUG: surface %i drawing thickness is negative\n",s);}
		if(srf->fdrawmode<DMno || srf->fdrawmode>DMvef) {error++;simLog(sim,10," SMOLDYN BUG: surface %i front drawing mode is out of bounds\n",s);}
		if(srf->bdrawmode<DMno || srf->bdrawmode>DMvef) {error++;simLog(sim,10," SMOLDYN BUG: surface %i back drawing mode is out of bounds\n",s);}
		if(srf->fshiny<0 || srf->fshiny>128) {error++;simLog(sim,10," SMOLDYN BUG: surface %i front shininess is out of bounds\n",s);}
		if(srf->bshiny<0 || srf->bshiny>128) {error++;simLog(sim,10," SMOLDYN BUG: surface %i back shininess is out of bounds\n",s);} }

	for(s=0;s<srfss->nsrf;s++) {											// surface missing and incorrect elements, about panels
		srf=srfss->srflist[s];
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) {
			if(srf->npanel[ps]>srf->maxpanel[ps]) {error++;simLog(sim,10," SMOLDYN BUG: surface %i has more panels of shape %i defined than allocated\n",s,ps);}
			if(srf->npanel[ps]>0) {
				if(!srf->pname || !srf->pname[ps]) {error++;simLog(sim,10," SMOLDYN BUG: surface %i has no name strings allocated for shape %i\n",s,ps);}
				if(!srf->panels || !srf->panels[ps]) {error++;simLog(sim,10," SMOLDYN BUG: surface %i has panels listed but not allocated for shape %i\n",s,ps);} }
			for(p=0;p<srf->npanel[ps];p++) {
				if(strlen(srf->pname[ps][p])==0) {error++;simLog(sim,10," ERROR: surface %i, panel shape %i, panel %i is unnamed\n",s,ps,p);}
				if(!srf->panels[ps][p]) {error++;simLog(sim,10," SMOLDYN BUG: surface %i, panel shape %i, panel %i pointer is NULL\n",s,ps,p);} }}}
	
	for(s=0;s<srfss->nsrf;s++) {											// panel missing and incorrect elements
		srf=srfss->srflist[s];
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
			for(p=0;p<srf->npanel[ps];p++) {
				pnl=srf->panels[ps][p];
				if(pnl->pname!=srf->pname[ps][p]) {error++;simLog(sim,10," SMOLDYN BUG: surface %i, panel shape %i, panel %i name pointer does not match surface pointer\n",s,ps,p);}
				if(pnl->ps!=ps) {error++;simLog(sim,10," SMOLDYN BUG: surface %i, shape %i, panel %i listing of shape does not match actual shape\n",s,ps,p);}
				if(pnl->srf!=srf) {error++;simLog(sim,10," SMOLDYN BUG: surface %i, shape %i, panel %i listing of surface does not match actual surface\n",s,ps,p);}
				if(pnl->npts!=panelpoints(ps,dim)) {error++;simLog(sim,10," SMOLDYN BUG: surface %i, shape %i, panel %i number of points (npts) is incorrect\n",s,ps,p);}
				for(pt=0;pt<pnl->npts;pt++)
					if(!pnl->point || !pnl->point[pt]) {error++;simLog(sim,10," SMOLDYN BUG: surface %i, shape %i, panel %i points element is not fully allocated\n",s,ps,p);} }}

	for(s=0;s<srfss->nsrf;s++) {											// panel points and front
		srf=srfss->srflist[s];
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
			for(p=0;p<srf->npanel[ps];p++) {
				pnl=srf->panels[ps][p];
				point=pnl->point;
				front=pnl->front;
				if(ps==PSrect) {
					if(dim==1) {
						if(point[0][0]<lowwall[0]) {warn++;simLog(sim,5," WARNING: surface %s, panel %s is entirely outside the system volume\n",srf->sname,pnl->pname);}
						if(point[0][0]>highwall[0]) {warn++;simLog(sim,5," WARNING: surface %s, panel %s is entirely outside the system volume\n",srf->sname,pnl->pname);}
						if(!(front[0]==-1 || front[0]==1)) {error++;simLog(sim,10," SMOLDYN BUG: surface %s, panel %s front vector is set incorrectly\n",srf->sname,pnl->pname);}
						if(front[1]!=0) {error++;simLog(sim,10," SMOLDYN BUG: surface %s, panel %s front vector is set incorrectly\n",srf->sname,pnl->pname);} }
					else if(dim==2) {
						Geo_LineNormal(point[0],point[1],norm);
						if(!Geo_LineXaabb2(point[0],point[1],norm,lowwall,highwall)) {warn++;simLog(sim,5," WARNING: surface %s, panel %s is entirely outside the system volume\n",srf->sname,pnl->pname);}
						if(!(front[0]==-1 || front[0]==1)) {error++;simLog(sim,10," SMOLDYN BUG: surface %s, panel %s front vector is set incorrectly\n",srf->sname,pnl->pname);} }}}}
				// should also check points and front for other shapes and dimensions

	for(s=0;s<srfss->nsrf;s++) {											// jump surfaces and panels
		srf=srfss->srflist[s];
		if(srf->action && nspecies>1) {
			fjump=bjump=0;
			for(i=1;i<nspecies;i++) {
				if(srf->action[i][MSsoln][PFfront]==SAjump) fjump=1;
				if(srf->action[i][MSsoln][PFback]==SAjump) bjump=1; 
			}
			if(fjump)
				for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
					for(p=0;p<srf->npanel[ps];p++) {
						pnl=srf->panels[ps][p];
						if(!pnl->jumpp[PFfront]) {warn++;simLog(sim,5," WARNING: front of surface %s has jump action but panel %s has no jump destination\n",srf->sname,pnl->pname);}
						else if(pnl->jumpp[PFfront]==pnl) {warn++;simLog(sim,5," WARNING: surface %s, panel %s front is a jump type panel that is its own destination\n",srf->sname,pnl->pname);}
						if(pnl->jumpp[PFfront] && !(pnl->jumpf[PFfront]==PFfront || pnl->jumpf[PFfront]==PFback)) {error++;simLog(sim,10," SMOLDYN BUG: surface %s, panel %s jumps to an undefined panel face\n",srf->sname,pnl->pname);} }
			if(bjump)
				for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
					for(p=0;p<srf->npanel[ps];p++) {
						pnl=srf->panels[ps][p];
						if(!pnl->jumpp[PFback]) {warn++;simLog(sim,5," WARNING: back of surface %s has jump action but panel %s has no jump destination\n",srf->sname,pnl->pname);}
						else if(pnl->jumpp[PFback]==pnl) {warn++;simLog(sim,5," WARNING: surface %s, panel %s back is a jump type panel that is its own destination\n",srf->sname,pnl->pname);}
						if(pnl->jumpp[PFback] && !(pnl->jumpf[PFback]==PFfront || pnl->jumpf[PFback]==PFback)) {error++;simLog(sim,10," SMOLDYN BUG: surface %s, panel %s jumps to an undefined panel face\n",srf->sname,pnl->pname);} }}}

	for(s=0;s<srfss->nsrf;s++) {											// make sure panel panel neighbors are stored correctly
		srf=srfss->srflist[s];
		for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
			for(p=0;p<srf->npanel[ps];p++) {
				pnl=srf->panels[ps][p];
				if(pnl->nneigh>0 && !pnl->neigh) {error++;simLog(sim,10," SMOLDYN BUG: surface %s, panel %s has %i neighbors, but none listed\n",srf->sname,pnl->pname,pnl->nneigh);} }}

	for(s=0;s<srfss->nsrf;s++) {											// make sure porting actions are linked to ports
		srf=srfss->srflist[s];
		if(srf->action) {
			fport=bport=0;
			for(i=0;i<nspecies;i++) {
				if(srf->action[i][MSsoln][PFfront]==SAport) fport=1;
				if(srf->action[i][MSsoln][PFback]==SAport) bport=1; 
			}
			if(fport && !srf->port[PFfront]) {
				error++;
				simLog(sim,9," ERROR: surface %s front face has porting action, but is not linked to a port\n",srf->sname);
				for(i=0;i<nspecies;i++) 
					if(srf->action[i][MSsoln][PFfront]==SAport) 
						srf->action[i][MSsoln][PFfront]=SAreflect; 
			}
			if(bport && !srf->port[PFback]) {
				error++;
				simLog(sim,9," ERROR: surface %s back face has porting action, but is not linked to a port\n",srf->sname);
				for(i=0;i<nspecies;i++) 
						if(srf->action[i][MSsoln][PFback]==SAport) 
							srf->action[i][MSsoln][PFback]=SAreflect; 

			}}}

	if(sim->mols)														// check that requested surface rates can be achieved
		for(s=0;s<srfss->nsrf;s++) {
			srf=srfss->srflist[s];
			if(srf->action && srf->actdetails) {
				for(i=1;i<sim->mols->nspecies;i++)
					for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1))
						for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1))
							if(srf->action[i][ms][face]==SAmult) {
								actdetails=srf->actdetails[i][ms][face];
								if(!actdetails) {
									simLog(sim,10," SMOLDYN BUG: action is SAmult but action details not allocated");error++;}
								else {
									for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1)) {
										prob=actdetails->srfprob[ms2];
										if(prob<0 || prob>1) {
											error++;
											simLog(sim,10," SMOLDYN BUG: surface interaction probability is <0 or >1\n"); }
										else {
											if(actdetails->srfdatasrc[ms2]==1)
												if(fabs((actdetails->srfrate[ms2]-srfcalcrate(sim,srf,i,ms,face,ms2))/actdetails->srfrate[ms2])>0.01) {
													warn++;
													simLog(sim,5," WARNING: requested surface '%s' rate for %s(%s) -> %s(%s) cannot be achieved\n",srf->sname,sim->mols->spname[i],molms2string(ms,string),sim->mols->spname[i],molms2string(ms2,string2)); }}}}}}}
	
	if(sim->mols)																			// check that don't have conflicting surface actions and emitters
		for(s=0;s<srfss->nsrf;s++) {
			srf=srfss->srflist[s];
			if(srf->action) {
				for(face=(PanelFace)0;face<2;face=(PanelFace)(face+1)) {
					for(i=1;i<sim->mols->nspecies;i++)
						if(srf->nemitter[face] && srf->nemitter[face][i])
							if(srf->action[i][MSsoln][face]!=SAreflect) {
								warn++;
								simLog(sim,5," WARNING: surface '%s' %s-side action for %s is listed as %s, but molecules will absorb instead because of unbounded emitters\n",srf->sname,surfface2string(face,string),sim->mols->spname[i],surfact2string(srf->action[i][MSsoln][face],string2)); }}}}

	if(sim->mols)																			// check that emitter absorptions aren't pegged to 1
		for(s=0;s<srfss->nsrf;s++) {
			srf=srfss->srflist[s];
			for(face=PFfront;face<=PFback;face=(PanelFace)(face+1))
				for(i=1;i<sim->mols->nspecies;i++)
					if(srf->nemitter[face] && srf->nemitter[face][i])
						for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
							for(p=0;p<srf->npanel[ps];p++)
								if(srf->panels[ps][p]->emitterabsorb[face][i]==1.0) {
									warn++;
									simLog(sim,5," WARNING: surface '%s', panel %s, %s-side absorption probability for %s cannot be made large enough for accurate effective unbounded diffusion\n",srf->sname,srf->panels[ps][p]->pname,surfface2string(face,string),sim->mols->spname[i]); }}

// check panels that normals are normal and normalized, with proper signs (error)
// check that panels are at least partially inside system volume (warning)
// check winding direction, if possible
// check for asymmetric jumping and that jumping is enabled for jumpable panels
// check that epsilon, margin, and neighbor_dist values are reasonable

	if(warnptr) *warnptr=warn;
	return error; }


/******************************************************************************/
/****************************** structure set up ******************************/
/******************************************************************************/


/* surfenablesurfaces */
int surfenablesurfaces(simptr sim,int maxsurf) {
	surfacessptr srfss;

	if(sim->srfss)									// check for redundant function call
		if(maxsurf==-1 || sim->srfss->maxsrf==maxsurf)
			if((sim->mols && sim->srfss->maxspecies==sim->mols->maxspecies) || (!sim->mols && sim->srfss->maxspecies==0))
				return 0;
	srfss=surfacessalloc(sim->srfss,maxsurf<0?5:maxsurf,sim->mols?sim->mols->maxspecies:0,sim->dim,sim->mols?sim->mols->maxsitecode:0);
	if(!srfss) return 1;
	sim->srfss=srfss;
	srfss->sim=sim;
	boxsetcondition(sim->boxs,SCparams,0);
	surfsetcondition(sim->srfss,SClists,0);
	return 0; }


/* surfexpandmaxspecies */
int surfexpandmaxspecies(simptr sim,int maxspecies) {
	surfacessptr srfss;

	if(!sim->srfss) return 0;
	srfss=sim->srfss;
	if(srfss->maxspecies>=maxspecies) return 0;
	srfss=surfacessalloc(srfss,srfss->maxsrf,maxspecies,sim->dim,sim->mols->maxsitecode);
	if(!srfss) return 1;
	return 0; }


/* surfaddsurface */
surfaceptr surfaddsurface(simptr sim,const char *surfname) {
	int er,s;
	surfacessptr srfss;
	surfaceptr srf;

	if(!sim->srfss) {
		er=surfenablesurfaces(sim,-1);
		if(er) return NULL; }
	srfss=sim->srfss;

	s=stringfind(srfss->snames,srfss->nsrf,surfname);
	if(s<0) {
		if(srfss->nsrf==srfss->maxsrf) {
			er=surfenablesurfaces(sim,srfss->nsrf*2+1);
			if(er) return NULL; }
		s=srfss->nsrf++;
		strncpy(srfss->snames[s],surfname,STRCHAR-1);
		srfss->snames[s][STRCHAR-1]='\0';
		srf=srfss->srflist[s];
		surfsetcondition(srfss,SClists,0); }
	else
		srf=srfss->srflist[s];

	surfsetcondition(sim->srfss,SClists,0);
	return srf; }


/* surfsetcondition. */
void surfsetcondition(surfacessptr surfss,enum StructCond cond,int upgrade) {
	if(!surfss) return;
	if(upgrade==0 && surfss->condition>cond) surfss->condition=cond;
	else if(upgrade==1 && surfss->condition<cond) surfss->condition=cond;
	else if(upgrade==2) surfss->condition=cond;
	if(surfss->sim && surfss->condition<surfss->sim->condition) {
		cond=surfss->condition;
		simsetcondition(surfss->sim,cond==SCinit?SClists:cond,0); }
	return; }


/* surfsetepsilon */
int surfsetepsilon(simptr sim,double epsilon) {
	int er;
	
	if(!sim->srfss) {
		er=surfenablesurfaces(sim,-1);
		if(er) return 2; }
	if(epsilon<=0) return 3;
	sim->srfss->epsilon=epsilon;
	return 0; }


/* surfsetmargin */
int surfsetmargin(simptr sim,double margin) {
	int er;
	
	if(!sim->srfss) {
		er=surfenablesurfaces(sim,-1);
		if(er) return 2; }
	if(margin<0) return 3;
	sim->srfss->margin=margin;
	return 0; }


/* surfsetneighdist */
int surfsetneighdist(simptr sim,double neighdist) {
	int er;
	
	if(!sim->srfss) {
		er=surfenablesurfaces(sim,-1);
		if(er) return 2; }
	if(neighdist<=0) return 3;
	sim->srfss->neighdist=neighdist;
	return 0; }


/* surfsetneighhop */
int surfsetneighhop(surfaceptr srf,int neighhop) {
	if(!srf) return 1;
	srf->neighhop=neighhop;
	return 0; }


/* surfsetcolor */
int surfsetcolor(surfaceptr srf,enum PanelFace face,double *rgba) {
	int col;

	if(!srf) return 1;
	for(col=0;col<4;col++)
		if(rgba[col]<0 || rgba[col]>1) return 2;

	if(face==PFfront || face==PFboth)
		for(col=0;col<4;col++) srf->fcolor[col]=rgba[col];
	if(face==PFback || face==PFboth)
		for(col=0;col<4;col++) srf->bcolor[col]=rgba[col];

	return 0; }


/* surfsetedgepts */
int surfsetedgepts(surfaceptr srf,double value) {
	if(!srf) return 1;
	if(value<0) return 2;
	srf->edgepts=value;
	return 0; }


/* surfsetstipple */
int surfsetstipple(surfaceptr srf,int factor,int pattern) {
	if(!srf) return 1;
	if(factor>=0) {
		if(factor==0) return 2;
		srf->edgestipple[0]=(unsigned int) factor; }
	if(pattern>=0) {
		if(pattern>0xFFFF) return 2;
		srf->edgestipple[1]=(unsigned int) pattern; }
	return 0; }


/* surfsetdrawmode */
int surfsetdrawmode(surfaceptr srf,enum PanelFace face,enum DrawMode dm) {
	if(!srf) return 1;
	if(dm==DMnone) return 2;
	if(face==PFfront || face==PFboth) srf->fdrawmode=dm;
	if(face==PFback || face==PFboth) srf->bdrawmode=dm;
	return 0; }


/* surfsetshiny */
int surfsetshiny(surfaceptr srf,enum PanelFace face,double shiny) {
	if(!srf) return 1;
	if(shiny<0 || shiny>128) return 2;
	if(face==PFfront || face==PFboth) srf->fshiny=shiny;
	if(face==PFback || face==PFboth) srf->bshiny=shiny;
	return 0; }



/* surfsetaction */
int surfsetaction(surfaceptr srf,int ident,enum MolecState ms,enum PanelFace face,enum SrfAction act) {
	int i;
	enum MolecState ms1,mslo,mshi;
	simptr sim;

	sim=srf->srfss->sim;
	if(ms==MSbsoln || ms==MSnone) return 2;
	else if(ms==MSall) {mslo=MolecState(0);mshi=MolecState(MSMAX-1);}
	else mslo=mshi=ms;

	if(face==PFfront || face==PFback || face==PFboth) {		// face is front, back, or both
		if(!(act==SAreflect || act==SAtrans || act==SAabsorb || act==SAjump || act==SAport || act==SAmult || act==SArotate_jump))			// cplx	
			return 3;

		i=0;
		while(1) {
			if(ident>0) {						// single species
				if(i==ident) break;
				i=ident; }
			else if(ident==-5) {		// all species
				if(++i==sim->mols->nspecies) break; }
			else if(ident==-6) {		// wildcard selected species
				if((i=molwildcardname(sim->mols,NULL,0,0))<0) break; }
			
			for(ms1=mslo;ms1<=mshi;ms1=MolecState(ms1 + 1)) {
				if(face==PFfront || face==PFboth) {
					srf->action[i][ms1][PFfront]=act;
				}
				if(face==PFback || face==PFboth) {
					srf->action[i][ms1][PFback]=act;
				} 
			}}}
	else {																// face is PFnone
		if(!(act==SAmult || act==SAno)) return 3;										

		i=0;
		while(1) {
			if(ident>0) {						// single species
				if(i==ident) break;
				i=ident; }
			else if(ident==-5) {		// all species
				if(++i==sim->mols->nspecies) break; }
			else if(ident==-6) {		// wildcard selected species
				if((i=molwildcardname(sim->mols,NULL,0,0))<0) break; }
			
			for(ms1=mslo;ms1<=mshi;ms1=MolecState(ms1 + 1)) {
				srf->action[i][ms1][face]=act; 
			}}}														// cplx originially there.

	surfsetcondition(srf->srfss,SCparams,0);
	return 0; }


/* surfsetrate */
int surfsetrate(surfaceptr srf,int ident,enum MolecState ms,enum MolecState ms1,enum MolecState ms2,int newident,double value,int which) {
	int i;
	enum MolecState ms3,ms4;
	enum PanelFace face;
	surfactionptr actdetails;
	simptr sim;

	sim=srf->srfss->sim;
	if(ms==MSbsoln || ms==MSall) return 2;

	if(ms1>MSbsoln) return 3;
	else if(ms!=MSsoln && ms1!=MSsoln && ms1!=MSbsoln && ms1!=ms) return 3;

	if(ms2>MSbsoln) return 4;
	else if(ms2==ms1) return 4;

	if((newident!=-5 && newident<0) || newident>=srf->srfss->maxspecies) return 5;

	if(value<0) return 6;
	else if(which==2 && value>1) return 6;

	srftristate2index(ms,ms1,ms2,&ms3,&face,&ms4);

	i=0;
	while(1) {
		if(ident>0) {						// single species
			if(i==ident) break;
			i=ident; }
		else if(ident==-5) {		// all species
			if(++i==sim->mols->nspecies) break; }
		else if(ident==-6) {		// wildcard selected species
			if((i=molwildcardname(sim->mols,NULL,0,0))<0) break; }
		
		if(!srf->actdetails[i][ms3][face]) {
			actdetails=surfaceactionalloc(i);				// allocate action details
			if(!actdetails) return -1;
			srf->actdetails[i][ms3][face]=actdetails; }
		actdetails=srf->actdetails[i][ms3][face];

		srf->action[i][ms3][face]=SAmult;						// set action to mult
		if(which==1) {															// and set up action details
			actdetails->srfrate[ms4]=value;
			actdetails->srfdatasrc[ms4]=1; }
		else if(which==2) {
			actdetails->srfprob[ms4]=value;
			actdetails->srfdatasrc[ms4]=2; }
		actdetails->srfnewspec[ms4]=(newident==-5)?i:newident; }

	surfsetcondition(srf->srfss,SCparams,0);
	return 0; }


/* surfsetmaxpanel */
int surfsetmaxpanel(surfaceptr srf,int dim,enum PanelShape ps,int maxpanel) {
	int ok;

	if(!srf) return 1;
	if(ps>=PSMAX) return 2;
	if(dim==1 && ps>PSsph) return 2;
	if(maxpanel==srf->maxpanel[ps]) return 0;
	if(maxpanel<srf->maxpanel[ps]) return 3;
	ok=panelsalloc(srf,dim,maxpanel,srf->srfss->maxspecies,ps);
	if(!ok) return -1;
	return 0; }


/* surfaddpanel */
int surfaddpanel(surfaceptr srf,int dim,enum PanelShape ps,const char *string,double *params,const char *name) {
	int p,ok,pdim,pt,d,axis01,axis12;
	panelptr pnl;
	double point[8][3],front[3],length;
	char ch;
	enum PanelShape ps2;
	
	if(!srf) return 1;
	if(ps>=PSMAX) return 2;
	if(dim==1 && ps>PSsph) return 2;

	if(ps==PSrect) {									// rectangle
		if(!string) return 3;
		ch=string[0];
		if(!(ch=='+' || ch=='-')) return 3;
		if(string[1]=='0' || string[1]=='x') pdim=0;			// pdim is the perpendicular axis
		else if(string[1]=='1' || string[1]=='y') pdim=1;
		else if(string[1]=='2' || string[1]=='z') pdim=2;
		else pdim=4;
		if(pdim>=dim) return 3;

		if(dim==1) {
			point[0][0]=params[0];
			front[0]=(double)(ch=='+'?1:-1);
			front[1]=0; }
		else if(dim==2) {									// params[0] and [1] are the given starting point, [2] is the length
			pt=1;
			if(ch=='-') pt*=-1;
			if(pdim==0) pt*=-1;
			if(params[2]<0) pt*=-1;
			if(pt<0) pt=0;									// pt is the point index of the given starting point
			point[pt][0]=params[0];
			point[pt][1]=params[1];
			point[!pt][0]=params[0]+(pdim==0?0:params[2]);
			point[!pt][1]=params[1]+(pdim==1?0:params[2]);
			front[0]=(double)(ch=='+'?1:-1);
			front[1]=(double)pdim;
			front[2]=(double)(!pdim);
			point[2][0]=front[0]*front[1];
			point[2][1]=-front[0]*front[2];
			point[3][0]=-front[0]*front[1];
			point[3][1]=front[0]*front[2]; }
		else if(dim==3) {									// params[0], [1], and [2] are the given starting point
			pt=1;														// params[3] and [4] are the dimensions on the other axes
			if(ch=='-') pt*=-1;
			if(params[3]*params[4]<0) pt*=-1;
			if(pt<0) pt=3;
			point[0][0]=params[0];					// point[0] is the given starting point
			point[0][1]=params[1];
			point[0][2]=params[2];
			point[pt][0]=params[0]+(pdim==2?params[3]:0);
			point[pt][1]=params[1]+(pdim==0?params[3]:0);
			point[pt][2]=params[2]+(pdim==1?params[4]:0);
			point[2][0]=params[0]+(pdim==0?0:params[3]);
			point[2][1]=params[1]+(pdim==1?0:(pdim==0?params[3]:params[4]));
			point[2][2]=params[2]+(pdim==2?0:params[4]);
			point[4-pt][0]=params[0]+(pdim==1?params[3]:0);
			point[4-pt][1]=params[1]+(pdim==2?params[4]:0);
			point[4-pt][2]=params[2]+(pdim==0?params[4]:0);
			front[0]=(double)(ch=='+'?1:-1);
			front[1]=(double)pdim;
			if(pdim==0) axis01=(pt==1?1:2);				// axis01 is the axis parallel to edge from point 0 to point 1
			else if(pdim==1) axis01=(pt==1?2:0);
			else axis01=(pt==1?0:1);
			axis12=(axis01+1)%3;									// axis12 is the axis parallel to edge from point 1 to point 2
			if(axis12==pdim) axis12=(axis12+1)%3;
			front[2]=(double)axis01;
			for(d=0;d<3;d++)
				point[4][d]=point[5][d]=point[6][d]=point[7][d]=0;
			point[4][axis12]=(point[2][axis12]>point[1][axis12]?-1:1);
			point[5][axis01]=(point[1][axis01]>point[0][axis01]?1:-1);
			point[6][axis12]=(point[2][axis12]>point[1][axis12]?1:-1);
			point[7][axis01]=(point[1][axis01]>point[0][axis01]?-1:1); }}

	else if(ps==PStri) {							// triangle
		if(dim==1) {
			point[0][0]=params[0];
			front[0]=(double)1; }
		else if(dim==2) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[1][0]=params[2];
			point[1][1]=params[3];
			Geo_LineNormal(point[0],point[1],front);
			point[2][0]=front[1];
			point[2][1]=-front[0];
			point[3][0]=-front[1];
			point[3][1]=front[0]; }
		else if(dim==3) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[0][2]=params[2];
			point[1][0]=params[3];
			point[1][1]=params[4];
			point[1][2]=params[5];
			point[2][0]=params[6];
			point[2][1]=params[7];
			point[2][2]=params[8];
			Geo_TriNormal(point[0],point[1],point[2],front);
			Geo_UnitCross(point[0],point[1],NULL,front,point[3]);
			Geo_UnitCross(point[1],point[2],NULL,front,point[4]);
			Geo_UnitCross(point[2],point[0],NULL,front,point[5]); }}
	
	else if(ps==PSsph) {							// sphere
		if(dim==1) {
			point[0][0]=params[0];
			point[1][0]=fabs(params[1]);
			front[0]=(double)sign(params[1]); }
		else if(dim==2) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[1][0]=fabs(params[2]);
			point[1][1]=params[3];
			if(point[1][1]<=0) return 4;
			front[0]=(double)sign(params[2]); }
		else if(dim==3) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[0][2]=params[2];
			point[1][0]=fabs(params[3]);
			point[1][1]=params[4];
			if(point[1][1]<=0) return 4;
			point[1][2]=params[5];
			if(point[1][2]<=0) return 4;
			front[0]=(double)sign(params[3]); }}
	
	else if(ps==PScyl) {								// cylinder
		if(dim==2) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[1][0]=params[2];
			point[1][1]=params[3];
			if(point[0][0]==point[1][0] && point[0][1]==point[1][1]) return 5;
			point[2][0]=fabs(params[4]);
			Geo_LineNormal(point[0],point[1],front);
			front[2]=(double)sign(params[4]);
			point[3][0]=front[1];
			point[3][1]=-front[0];
			point[4][0]=-front[1];
			point[4][1]=front[0]; }
		else if(dim==3) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[0][2]=params[2];
			point[1][0]=params[3];
			point[1][1]=params[4];
			point[1][2]=params[5];
			if(point[0][0]==point[1][0] && point[0][1]==point[1][1] && point[0][2]==point[1][2]) return 5;
			point[2][0]=fabs(params[6]);
			point[2][1]=params[7];
			if(point[2][1]<=0) return 4;
			point[2][2]=params[8];
			if(point[2][2]<=0) return 4;
			front[2]=(double)sign(params[6]);
			length=Geo_LineLength(point[0],point[1],3);
			point[3][0]=-(point[1][0]-point[0][0])/length;
			point[3][1]=-(point[1][1]-point[0][1])/length;
			point[3][2]=-(point[1][2]-point[0][2])/length;
			point[4][0]=(point[1][0]-point[0][0])/length;
			point[4][1]=(point[1][1]-point[0][1])/length;
			point[4][2]=(point[1][2]-point[0][2])/length; }}
	
	else if(ps==PShemi) {								// hemisphere
		if(dim==2) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[1][0]=fabs(params[2]);
			point[1][1]=params[5];
			if(point[1][1]<=0) return 4;
			point[2][0]=params[3];
			point[2][1]=params[4];
			if(normalizeVD(point[2],2)<=0) return 6;
			front[0]=(double)sign(params[2]); }
		else if(dim==3) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[0][2]=params[2];
			point[1][0]=fabs(params[3]);
			point[1][1]=params[7];
			if(point[1][1]<=0) return 4;
			point[1][2]=params[8];
			if(point[1][2]<=0) return 4;
			point[2][0]=params[4];
			point[2][1]=params[5];
			point[2][2]=params[6];
			if(normalizeVD(point[2],3)<=0) return 6;
			front[0]=(double)sign(params[3]); }}
	
	else if(ps==PSdisk) {								// disk
		if(dim==2) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[1][0]=params[2];
			if(point[1][0]<=0) return 7;
			front[0]=params[3];
			front[1]=params[4];
			if(normalizeVD(front,2)<=0) return 8; }
		else if(dim==3) {
			point[0][0]=params[0];
			point[0][1]=params[1];
			point[0][2]=params[2];
			point[1][0]=params[3];
			if(point[1][0]<=0) return 7;
			point[1][1]=params[7];
			if(point[1][1]<=0) return 4;
			front[0]=params[4];
			front[1]=params[5];
			front[2]=params[6];
			if(normalizeVD(front,3)<=0) return 8; }}
	
	p=-1;
	if(name && name[0]!='\0') {
		for(ps2=(PanelShape)0;ps2<PSMAX && p==-1;ps2=(PanelShape)(ps2+1))
			p=stringfind(srf->pname[ps2],srf->npanel[ps2],name);
		if(p!=-1 && (ps2-1)!=ps) return 9; }      // panel name was used before for a different shape

	if(p==-1) {																	// new panel
		if(srf->npanel[ps]==srf->maxpanel[ps]) {
			ok=panelsalloc(srf,dim,srf->maxpanel[ps]*2+1,srf->srfss->maxspecies,ps);
			if(!ok) return -1; }
		p=srf->npanel[ps]++; }

	pnl=srf->panels[ps][p];
	for(pt=0;pt<pnl->npts;pt++)
		for(d=0;d<dim;d++)
			pnl->point[pt][d]=point[pt][d];
	for(d=0;d<3;d++)
		pnl->front[d]=front[d];
	if(name && name[0]!='\0') strcpy(srf->pname[ps][p],name);

	surfsetcondition(srf->srfss,SClists,0);
	boxsetcondition(srf->srfss->sim->boxs,SCparams,0);
	compartsetcondition(srf->srfss->sim->cmptss,SCparams,0);
	
	return 0; }


/* surfsetemitterabsorption */
int surfsetemitterabsorption(simptr sim) {
	surfacessptr srfss;
	surfaceptr srf;
	panelptr pnl;
	int dim,s,i,nspecies,emit,er,p;
	enum PanelFace face;
	enum PanelShape ps;
	double difc,middle[DIMMAX],normal[DIMMAX],numer,denom,amount,*pos,dist,vdiff[DIMMAX],kappa,prob;

	srfss=sim->srfss;
	dim=sim->dim;
	nspecies=sim->mols->nspecies;
	er=0;
	for(s=0;s<srfss->nsrf;s++)
		for(face=PFfront;face<=PFback;face=(PanelFace)(face+1)) {
			srf=srfss->srflist[s];
			if(srf->nemitter[face])
				for(i=1;i<nspecies;i++)
					if(srf->nemitter[face][i]) {
						difc=sim->mols->difc[i][MSsoln];
						for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1))
							for(p=0;p<srf->npanel[ps];p++) {
								pnl=srf->panels[ps][p];
								panelmiddle(pnl,middle,dim,1);						// middle position
								panelnormal(pnl,middle,face==PFfront?PFback:PFfront,dim,normal);	// normal vector
								numer=0;
								denom=0;
								for(emit=0;emit<srf->nemitter[face][i];emit++) {
									amount=srf->emitteramount[face][i][emit];
									pos=srf->emitterpos[face][i][emit];
									dist=distanceVVD(middle,pos,dim);
									if(!(dist>0)) er=1;
									denom+=amount/dist;
									sumVD(1.0,middle,-1.0,pos,vdiff,dim);
									numer+=amount*dotVVD(vdiff,normal,dim)/(dist*dist*dist); }
								kappa=difc*numer/denom;
								prob=surfaceprob(kappa,0,sim->dt,difc,NULL,SPAirrAds);
								pnl->emitterabsorb[face][i]=prob; }}}

	if(er) simLog(sim,5,"WARNING: an unbounded emitter is at a surface panel which will cause inaccurate operation");
	return er; }


/* surfsetjumppanel */
int surfsetjumppanel(surfaceptr srf,panelptr pnl1,enum PanelFace face1,int bidirect,panelptr pnl2,enum PanelFace face2) {
	if(!srf) return 1;
	if(!pnl1) return 2;
	if(face1!=PFfront && face1!=PFback) return 3;
	if(bidirect!=0 && bidirect!=1) return 4;
	if(!pnl2 || pnl2==pnl1 || pnl2->ps!=pnl1->ps) return 5;
	if(face2!=PFfront && face2!=PFback) return 6;

	pnl1->jumpp[face1]=pnl2;
	pnl1->jumpf[face1]=face2;
	if(bidirect) {
		pnl2->jumpp[face2]=pnl1;
		pnl2->jumpf[face2]=face1; }
	return 0; }


/* srfcalcrate */
double srfcalcrate(simptr sim,surfaceptr srf,int i,enum MolecState ms1,enum PanelFace face,enum MolecState ms2) {
	double rate,prob,probrev,sum,dt,difc;
	enum MolecState ms3,ms4,ms5;
	surfactionptr actdetails;
	enum PanelFace face2;
	
	if(ms1==MSsoln && face==PFnone) return -1;		// impossible situation
	if(srf->action[i][ms1][face]!=SAmult) return -1; // could be computed but isn't
	actdetails=srf->actdetails[i][ms1][face];
	if(!actdetails) return -1;							// no data available
	prob=actdetails->srfprob[ms2];
	if(prob<0 || prob>1) return -2;					// impossible probability
	if(prob==0) return 0;															// rate is 0

	srfreverseaction(ms1,face,ms2,&ms3,&face2,&ms4);	// find probability of reverse action
	if(face2==PFboth) probrev=0;
	else if(srf->actdetails[i][ms3][face2])
		probrev=srf->actdetails[i][ms3][face2]->srfprob[ms4];
	else
		probrev=0;
	if(probrev<0) probrev=0;

	dt=sim->dt;
	difc=sim->mols->difc[i][MSsoln];

	sum=0;
	if(ms1!=MSsoln && face==PFnone) {							// if surface-bound, find sum of probabilities
		for(ms5=(MolecState)0;ms5<MSMAX1;ms5=(MolecState)(ms5+1)) {
			if(ms5!=ms1 && actdetails->srfprob[ms5]>=0)
				sum+=actdetails->srfprob[ms5]; }}

	if(ms1==MSsoln) {															// find rates
		if((face==PFfront && ms2==MSsoln) || (face==PFback && ms2==MSbsoln))
			rate=-3;																			// solution to solution (reflect), can't compute
		else if(ms2==MSsoln || ms2==MSbsoln) {					// solution to solution (transmit)
			if(probrev<=0) rate=surfacerate(prob,0,dt,difc,NULL,SPAirrTrans);
			else rate=surfacerate(prob,probrev,dt,difc,NULL,SPArevTrans); }
		else {																					// solution to surface (adsorb)
			if(probrev<=0) rate=surfacerate(prob,0,dt,difc,NULL,SPAirrAds);
			else rate=surfacerate(prob,probrev,dt,difc,NULL,SPArevAds); }}

	else {
		if(face==PFnone) {
			if(ms2==MSsoln || ms2==MSbsoln) {						// surface to solution (desorb)
				if(probrev<=0) rate=surfacerate(prob,sum,dt,difc,NULL,SPAirrDes);
				else rate=surfacerate(prob,probrev,dt,difc,NULL,SPArevDes); }
			else {																			// surface to surface (flip)
				if(ms1==ms2) rate=-3;											// same state
				else if(probrev<=0) rate=surfacerate(prob,sum,dt,difc,NULL,SPAirrFlip);
				else rate=surfacerate(prob,probrev,dt,difc,NULL,SPArevFlip); }}
		else {
			if((face==PFfront && ms2==MSsoln) || (face==PFback && ms2==MSbsoln))
				rate=-3;																	// surface-bound reflect, can't compute
			else if(ms2==MSsoln || ms2==MSbsoln) {			// surface-bound transmit
				if(probrev<=0) rate=surfacerate(prob,0,dt,difc,NULL,SPAirrTrans);
				else rate=surfacerate(prob,probrev,dt,difc,NULL,SPArevTrans); }
			else {																			// surface-bound hop to new surface
				rate=surfacerate(prob,0,dt,difc,NULL,SPAirrAds); }}}
	
	return rate;
	return 0; }


/* srfcalcprob */
double srfcalcprob(simptr sim,surfaceptr srf,int i,enum MolecState ms1,enum PanelFace face,enum MolecState ms2) {
	double rate,prob,raterev,sum,dt,difc;
	enum MolecState ms3,ms4,ms5;
	surfactionptr actdetails;
	enum PanelFace face2;

	if(ms1==MSsoln && face==PFnone) return 0;			// impossible situation
	if(srf->action[i][ms1][face]!=SAmult) return -1;	// could be computed but isn't
	actdetails=srf->actdetails[i][ms1][face];
	if(!actdetails) return -1;							// no data available
	if(actdetails->srfdatasrc[ms2]==2) return actdetails->srfprob[ms2];
	rate=actdetails->srfrate[ms2];
	if(rate<0) return -2;															// impossible rate
	if(rate==0) return 0;															// prob is 0

	srfreverseaction(ms1,face,ms2,&ms3,&face2,&ms4);	// find rate of reverse action
	if(face2==PFboth) raterev=0;
	else if(srf->actdetails[i][ms3][face2])
		raterev=srf->actdetails[i][ms3][face2]->srfrate[ms4];
	else
		raterev=0;
	if(raterev<0) raterev=0;

	dt=sim->dt;
	difc=sim->mols->difc[i][MSsoln];

	sum=0;
	if(ms1!=MSsoln && face==PFnone) {					// if surface-bound, find sum of rates
		for(ms5=(MolecState)0;ms5<MSMAX1;ms5=(MolecState)(ms5+1)) {
			if(ms5!=ms1 && actdetails->srfrate[ms5]>=0)
				sum+=actdetails->srfrate[ms5]; }}

	if(ms1==MSsoln) {																	// find probs
		if((face==PFfront && ms2==MSsoln) || (face==PFback && ms2==MSbsoln))
			prob=0;																				// solution to solution (reflect), can't compute
		else if(ms2==MSsoln || ms2==MSbsoln) {					// solution to solution (transmit)
			if(raterev<=0) prob=surfaceprob(rate,0,dt,difc,NULL,SPAirrTrans);
			else prob=surfaceprob(rate,raterev,dt,difc,NULL,SPArevTrans); }
		else {																					// solution to surface (adsorb)
			if(raterev<=0) prob=surfaceprob(rate,0,dt,difc,NULL,SPAirrAds);
			else prob=surfaceprob(rate,raterev,dt,difc,NULL,SPArevAds); }}

	else {
		if(face==PFnone) {
			if(ms2==MSsoln || ms2==MSbsoln) {						// surface to solution (desorb)
				if(raterev<=0) prob=surfaceprob(rate,sum,dt,difc,NULL,SPAirrDes);
				else prob=surfaceprob(rate,raterev,dt,difc,NULL,SPArevDes); }
			else {																			// surface to surface (flip)
				if(ms1==ms2) prob=0;											// same state
				else if(raterev<=0) prob=surfaceprob(rate,sum,dt,difc,NULL,SPAirrFlip);
				else prob=surfaceprob(rate,raterev,dt,difc,NULL,SPArevFlip); }}
		else {
			if((face==PFfront && ms2==MSsoln) || (face==PFback && ms2==MSbsoln))
				prob=0;																	// surface-bound reflect, can't compute
			else if(ms2==MSsoln || ms2==MSbsoln) {			// surface-bound transmit
				if(raterev<=0) prob=surfaceprob(rate,0,dt,difc,NULL,SPAirrTrans);
				else prob=surfaceprob(rate,raterev,dt,difc,NULL,SPArevTrans); }
			else {																			// surface-bound hop to new surface
				prob=surfaceprob(rate,0,dt,difc,NULL,SPAirrAds); }}}
	
	return prob; }


/* surfsetneighbors */
int surfsetneighbors(panelptr pnl,panelptr *neighlist,int nneigh,int add) {
	int newmaxneigh,p,p2;
	panelptr *newneigh;

	if(add) {
		if(pnl->nneigh+nneigh>pnl->maxneigh) { // allocate more space
			newmaxneigh=pnl->nneigh+nneigh;
			newneigh=(panelptr*) calloc(newmaxneigh,sizeof(panelptr));
			if(!newneigh) return 1;
			for(p2=0;p2<pnl->nneigh;p2++) newneigh[p2]=pnl->neigh[p2];
			for(;p2<newmaxneigh;p2++) newneigh[p2]=NULL;
			free(pnl->neigh);
			pnl->maxneigh=newmaxneigh;
			pnl->neigh=newneigh; }
		for(p=0;p<nneigh;p++) {
			for(p2=0;p2<pnl->nneigh && pnl->neigh[p2]!=neighlist[p];p2++);
			if(p2==pnl->nneigh)
				pnl->neigh[pnl->nneigh++]=neighlist[p]; }}
	else if(!neighlist) {
		pnl->nneigh=0; }
	else {
		for(p=0;p<nneigh;p++) {
			for(p2=0;p2<pnl->nneigh && pnl->neigh[p2]!=neighlist[p];p2++);
			if(p2<pnl->nneigh)
				pnl->neigh[p2]=pnl->neigh[--pnl->nneigh]; }}

	return 0; }


/* surfaddemitter */
int surfaddemitter(surfaceptr srf,enum PanelFace face,int i,double amount,double *pos,int dim) {
	int er,oldmax,newmax,emit,d;
	double *newamount,**newpos;

	newamount=NULL;
	newpos=NULL;
	if(!srf->maxemitter[face]) {
		er=emittersalloc(srf,face,srf->srfss->maxspecies,srf->srfss->maxspecies);
		CHECK(!er); }

	if(srf->nemitter[face][i]==srf->maxemitter[face][i]) {	// allocate emitter space if needed
		oldmax=srf->maxemitter[face][i];
		newmax=oldmax*2+1;

		CHECKMEM(newamount=(double*) calloc(newmax,sizeof(double)));	// emitteramount
		for(emit=0;emit<oldmax;emit++)
			newamount[emit]=srf->emitteramount[face][i][emit];
		for(;emit<newmax;emit++) newamount[emit]=0;

		CHECKMEM(newpos=(double**) calloc(newmax,sizeof(double*)));	// emitterpos
		for(emit=0;emit<oldmax;emit++)
			newpos[emit]=srf->emitterpos[face][i][emit];
		for(;emit<newmax;emit++)
			newpos[emit]=NULL;
		for(emit=oldmax;emit<newmax;emit++) {
			CHECKMEM(newpos[emit]=(double*) calloc(dim,sizeof(double)));
			for(d=0;d<dim;d++) newpos[emit][d]=0; }

		free(srf->emitteramount[face][i]);						// replace old with new
		srf->emitteramount[face][i]=newamount;
		free(srf->emitterpos[face][i]);
		srf->emitterpos[face][i]=newpos;
		srf->maxemitter[face][i]=newmax; }

	emit=srf->nemitter[face][i]++;									// add new emitter
	srf->emitteramount[face][i][emit]=amount;
	for(d=0;d<dim;d++)
		srf->emitterpos[face][i][emit][d]=pos[d];

	surfsetcondition(srf->srfss,SCparams,0);

	return 0;

 failure:
	free(newamount);
	free(newpos);
	if(ErrorType!=1) simLog(NULL,10,"Unable to allocate memory in surfaddemitter");
	return 1; }



/* surfreadstring */
surfaceptr surfreadstring(simptr sim,ParseFilePtr pfp,surfaceptr srf,const char *word,char *line2) {
	char nm[STRCHAR],nm1[STRCHAR],nm2[STRCHAR],facenm[STRCHAR],shapenm[STRCHAR],actnm[STRCHAR];
	int dim,i,p,p2,i1,i2,i3,itct,er,s;
	const int maxpnllist=128;
	panelptr pnl,pnllist[128];
	double fltv1[9],f1;
	surfacessptr srfss;
	enum PanelShape ps;
	enum PanelFace face,face2;
	enum SrfAction act;
	enum MolecState ms,ms1,ms2;
	enum DrawMode dm;
	int *ident,*sites;

	dim=sim->dim;
	srfss=sim->srfss;

	if(!strcmp(word,"name")) {								// name
		itct=sscanf(line2,"%s",nm);
		CHECKS(itct==1,"error reading surface name");
		srf=surfaddsurface(sim,nm);
		CHECKS(srf,"failed to add surface");
		CHECKS(!strnword(line2,2),"unexpected text following name"); }

	else if(!strcmp(word,"action")) {							// action
		CHECKS(srf,"need to enter surface name before action");
		CHECKS(srfss->maxspecies,"need to enter molecules before action");
		itct=sscanf(line2,"%s %s %s",nm1,nm2,actnm);
		//printf("surf.c line 2693, line2=%s, actnm=%s\n",line2,actnm);
		CHECKS(itct==3,"action format: species[(state)] face action");
		i=readmolname(sim,nm1,&ms1,0);
		face=surfstring2face(nm2);
		if(face==PFnone || (i<=0 && i>-5)) {	// try old format
			face=surfstring2face(nm1);
			i=readmolname(sim,nm2,&ms1,0);
			// i=readmolname_cplx(sim,nm2,&ident,&sites);
			// printf("surfstring, nm2:%s, i=%d\n",nm2,i);
			if(face==PFnone && (i<=0 && i>-5)) {	// back to new format
				i=readmolname(sim,nm1,&ms1,0);
				face=surfstring2face(nm2);
				CHECKS(face!=PFnone,"in action, face name needs to be 'front', 'back', or 'both'");
				CHECKS(i>0 || i==-5 || i==-6,"in action, molecule name not recognized"); }
			CHECKS(face!=PFnone,"in action, face name needs to be 'front', 'back', or 'both'");
			CHECKS(i>0 || i==-5 || i==-6,"in action, molecule name not recognized"); }
		act=surfstring2act(actnm);
		// CHECKS(act<=SAmult || act==SAadsorb,"in action statement, action not recognized or not permitted");
		er=surfsetaction(srf,i,ms1,face,act);
		CHECKS(!er,"BUG: error with surfsetaction statement");
		CHECKS(!strnword(line2,4),"unexpected text following action"); }

	else if(!strcmp(word,"rate") || !strcmp(word,"rate_internal")) { // rate, rate_internal
		CHECKS(srf,"need to enter surface name first");
		CHECKS(srfss->maxspecies,"need to enter molecules first");
		
#if OPTION_VCELL
		bool constRate = true;
		stringstream ss(line2);
		ss >> nm >> nm1 >> nm2;
		string rawStr;
		getline(ss, rawStr);
		size_t found = rawStr.find(";");
		ValueProvider* valueProvider = NULL;
		if(found!=string::npos)
		{
			string rateExpStr = rawStr.substr(0, found);
			valueProvider = sim->valueProviderFactory->createValueProvider(rateExpStr);
			try {
				double rate = valueProvider->getConstantValue();
				f1=rate;
				constRate = true;			
				delete valueProvider;
			} catch (...) {
				constRate = false;
			}
		}
		else
#endif
		{
			itct=sscanf(line2,"%s %s %s %lg",nm,nm1,nm2,&f1);
			CHECKS(itct==4,"format: species[(state)] state1 state2 value [new_species]");
		}
		
		i=readmolname(sim,nm,&ms,0);
		CHECKS(i>0 || i==-5,"molecule name not recognized");
		CHECKS(ms<MSbsoln,"state is not recognized or not permitted");
		ms1=molstring2ms(nm1);
		CHECKS(ms1<=MSbsoln,"state1 is not recognized or not permitted");
		ms2=molstring2ms(nm2);
		CHECKS(ms2<=MSbsoln,"state2 is not recognized or not permitted");
		if(ms==MSsoln) {
			CHECKS(ms1!=ms2,"it is not permitted for state1 to equal state2"); }
		else {
			CHECKS(ms1==ms || ms1==MSsoln || ms1==MSbsoln,"state1 does not make sense"); }
#if OPTION_VCELL
		if(found!=string::npos)
		{
			if(constRate == true)
				CHECKS(f1>=0,"negative surface rate values are not permitted");
			string name = rawStr.substr(found+2); //after the ";" denoting the end of rate, the found move one more position(the space) to get to the end of the line, which would be the species name
			char * tempLine = new char[name.size() + 1];
			strcpy(tempLine, name.c_str());
			line2 = tempLine;
		}
        else
#endif
		{
			CHECKS(f1>=0,"negative surface rate values are not permitted");
			line2=strnword(line2,5);
		}

		i3=i;
		if(line2) {
			itct=sscanf(line2,"%s",nm);
			CHECKS(itct==1,"cannot read new species name");
			i3=stringfind(sim->mols->spname,sim->mols->nspecies,nm);
			CHECKS(i3!=-1,"new species name not recognized");
			line2=strnword(line2,2); }
#if OPTION_VCELL
		if(found!=string::npos)
		{
			if(constRate)
			{
				if(!strcmp(word,"rate"))
					er=surfsetrate(srf,i,ms,ms1,ms2,i3,f1,1);
				else {
					CHECKS(f1<=1,"surface interaction probabilities cannot be greater than 1");
					er=surfsetrate(srf,i,ms,ms1,ms2,i3,f1,2); }
				CHECKS(er!=-1,"out of memory");
				CHECKS(!er,"BUG: error in surfsetrate");
			}
			else //rate is a function
			{
				if(!strcmp(word,"rate"))
					er=surfSetRateExp(srf, i, ms, ms1, ms2, i3, valueProvider,1);
				else {
					CHECKS(f1<=1,"surface interaction probabilities cannot be greater than 1");
					er=surfSetRateExp(srf, i, ms, ms1, ms2, i3, valueProvider, 2); }
				CHECKS(er!=-1,"out of memory");
				CHECKS(!er,"BUG: error in surfsetrate");
			}
		}
		else
#endif
		{
			if(!strcmp(word,"rate")){
				er=surfsetrate(srf,i,ms,ms1,ms2,i3,f1,1);
			}
			else {
				CHECKS(f1<=1,"surface interaction probabilities cannot be greater than 1");
				er=surfsetrate(srf,i,ms,ms1,ms2,i3,f1,2);
			 }
			CHECKS(er!=-1,"out of memory");
			CHECKS(!er,"BUG: error in surfsetrate");
		}
		CHECKS(!line2,"unexpected text at end of line"); }

	else if(!strcmp(word,"color") || !strcmp(word,"colour")) {		// color
		CHECKS(srf,"need to enter surface name before color");
		itct=sscanf(line2,"%s",facenm);
		CHECKS(itct==1,"color format: face color");
		face=surfstring2face(facenm);
		CHECKS(face!=PFnone,"in color, face name needs to be 'front', 'back', or 'both'");
		line2=strnword(line2,2);
		CHECKS(line2,"color format: face color");
		er=graphicsreadcolor(&line2,fltv1);
		CHECKS(er!=3,"color values need to be between 0 and 1");
		CHECKS(er!=4,"color name not recognized");
		CHECKS(er!=6,"alpha values need to be between 0 and 1");
		CHECKS(er==0,"format is either 3 numbers or color name, and then optional alpha value");
		er=surfsetcolor(srf,face,fltv1);
		CHECKS(!er,"BUG: error in surfsetcolor");
		CHECKS(!line2,"unexpected text following color"); }

	else if(!strcmp(word,"thickness")) {				// thickness
		CHECKS(srf,"need to enter surface name before thickness");
		itct=sscanf(line2,"%lg",&f1);
		CHECKS(itct==1,"thickness value is missing");
		CHECKS(f1>=0,"thickness value needs to be at least 0");
		er=surfsetedgepts(srf,f1);
		CHECKS(!er,"BUG: error in surfsetedgepts");
		CHECKS(!strnword(line2,2),"unexpected text following thickness"); }

	else if(!strcmp(word,"stipple")) {					// stipple
		CHECKS(srf,"need to enter surface name before stipple");
		itct=sscanf(line2,"%i %i",&i1,&i2);
		CHECKS(itct==2,"stipple format: factor pattern");
		CHECKS(i1>=1,"stipple factor needs to be >=1");
		CHECKS(i2>=0 && i2 <=0xFFFF,"stipple pattern needs to be between 0x00 and 0xFFFF");
		er=surfsetstipple(srf,i1,i2);
		CHECKS(!er,"BUG: error in surfsetstipple");
		CHECKS(!strnword(line2,3),"unexpected text following stipple"); }

	else if(!strcmp(word,"polygon")) {					// polygon
		CHECKS(srf,"need to enter surface name before polygon");
		itct=sscanf(line2,"%s %s",facenm,nm1);
		CHECKS(itct==2,"polygon format: face drawmode");
		face=surfstring2face(facenm);
		CHECKS(face!=PFnone,"in polygon, face name needs to be 'front', 'back', or 'both'");
		dm=surfstring2dm(nm1);
		CHECKS(dm!=DMnone,"in polygon, drawing mode is not recognized");
		er=surfsetdrawmode(srf,face,dm);
		CHECKS(!er,"BUG: error in surfsetdrawmode");
		CHECKS(!strnword(line2,3),"unexpected text following polygon"); }

	else if(!strcmp(word,"shininess")) {				// shininess
		CHECKS(srf,"need to enter surface name before shininess");
		itct=sscanf(line2,"%s %lg",facenm,&f1);
		CHECKS(itct==2,"shininess format: face value");
		face=surfstring2face(facenm);
		CHECKS(face!=PFnone,"face name needs to be 'front', 'back', or 'both'");
		CHECKS(f1>=0 && f1<=128,"shininess value needs to be between 0 and 128");
		er=surfsetshiny(srf,face,f1);
		CHECKS(!er,"BUG: error in surfsetshiny");
		CHECKS(!strnword(line2,3),"unexpected text following shininess"); }

	else if(!strcmp(word,"max_panels")) {					// max_panels
		CHECKS(srf,"need to enter surface name before max_panels");
		itct=sscanf(line2,"%s %i",shapenm,&i1);
		CHECKS(itct==2,"max_panels format: shape number");
		ps=surfstring2ps(shapenm);
		CHECKS(ps<PSMAX,"in max_panels, unknown panel shape");
		CHECKS(dim!=1 || ps<=PSsph,"in max_panels, panel shape is not permitted for a 1-D system");
		CHECKS(i1>=0,"max_panels number needs to be at least 0");
		er=surfsetmaxpanel(srf,dim,ps,i1);
		CHECKS(er!=-1,"out of memory allocating panels");
		CHECKS(!er,"BUG: error in surfsetmaxpanel");
		CHECKS(!strnword(line2,3),"unexpected text following max_panels"); }

	else if(!strcmp(word,"panel")) {							// panel
		CHECKS(srf,"need to enter surface name before panel");
		itct=sscanf(line2,"%s",shapenm);
		CHECKS(itct==1,"in panel, panel shape needs to be entered");
		ps=surfstring2ps(shapenm);
		CHECKS(ps<PSMAX,"in panel, unknown panel shape");
		CHECKS(dim!=1 || ps<=PSsph,"in panel, panel shape is not permitted for a 1-D system");
		line2=strnword(line2,2);
		CHECKS(line2,"panel data missing");
		if(ps==PSrect) {														// panel r ...
			itct=sscanf(line2,"%s",nm1);
			CHECKS(itct==1,"error reading rectangle panel direction");
			line2=strnword(line2,2);
			CHECKS(line2,"panel data missing"); }
		i1=surfpanelparams(ps,dim);
		itct=strreadnd(line2,i1,fltv1,NULL);
		CHECKS(itct==i1,"panel either has too few parameters or has unreadable parameters");
		line2=strnword(line2,i1+1);
		if(line2) {
			itct=sscanf(line2,"%s",nm);
			CHECKS(itct==1,"Error reading panel name");
			line2=strnword(line2,2); }
		else
			nm[0]='\0';
		er=surfaddpanel(srf,dim,ps,nm1,fltv1,nm);
		CHECKS(er!=-1,"out of memory adding panel");
		CHECKS(er!=1,"BUG: error in surfaddpanel");
		CHECKS(er!=2,"BUG: error in surfaddpanel");
		CHECKS(er!=3,"unable to read rectangle direction");
		CHECKS(er!=4,"drawing slices and stacks need to be positive numbers");
		CHECKS(er!=5,"cylinder ends need to be at different locations");
		CHECKS(er!=6,"outward pointing vector cannot be 0 length");
		CHECKS(er!=7,"radius needs to be positive");
		CHECKS(er!=8,"normal vector cannot be 0 length");
		CHECKS(er!=9,"panel name was used previously for a different panel shape");

		CHECKS(!line2,"unexpected text following panel"); }

	else if(!strcmp(word,"jump") || !strcmp(word,"rotate_jump")) {								// jump
		CHECKS(srf,"need to enter surface name before jump");
		itct=sscanf(line2,"%s %s",nm,facenm);
		CHECKS(itct==2,"format for jump: panel1 face1 -> panel2 face2");
		ps=PanelShape(0);
		p=0;
		while(ps<PSMAX && (p=stringfind(srf->pname[ps],srf->npanel[ps],nm))==-1) ps=PanelShape(ps + 1);
		CHECKS(p>=0,"first panel name listed in jump is not recognized");
		face=surfstring2face(facenm);
		CHECKS(face<=PFback,"first face listed in jump needs to be 'front' or 'back'");
		CHECKS(line2=strnword(line2,3),"format for jump: panel1 face1 -> panel2 face2");
		itct=sscanf(line2,"%s",nm);
		CHECKS(itct==1,"format for jump: panel1 face1 -> panel2 face2");
		if(!strcmp(nm,"->")) i1=0;
		else if(!strcmp(nm,"<->")) i1=1;
		else CHECKS(0,"jump operator needs to be -> or <->");
		CHECKS(line2=strnword(line2,2),"format for jump: panel1 face1 -> panel2 face2");
		itct=sscanf(line2,"%s %s",nm,facenm);
		CHECKS(itct==2,"format for jump: panel1 face1 -> panel2 face2");
		p2=stringfind(srf->pname[ps],srf->npanel[ps],nm);
		CHECKS(p2>=0,"second panel name listed in jump is not recognized, or not same shape as first panel");
		face2=surfstring2face(facenm);
		CHECKS(face2<=PFback,"second face listed in jump needs to be 'front' or 'back'");
		er=surfsetjumppanel(srf,srf->panels[ps][p],face,i1,srf->panels[ps][p2],face2);			// cplx
		CHECKS(!er,"BUG: error in surfsetjumppanel");
		CHECKS(!strnword(line2,3),"unexpected text following jump"); }

	else if(!strcmp(word,"neighbors") || !strcmp(word,"neighbours")) {					// neighbors
		CHECKS(srf,"need to enter surface name before neighbors");
		itct=sscanf(line2,"%s",nm);
		CHECKS(itct==1,"format for neighbors: panel neigh1 neigh2 ...");
		pnl=readpanelname(sim,srf,nm);
		CHECKS(pnl,"first panel name listed in neighbors, '%s', not recognized",nm);
		for(i1=0;i1<maxpnllist && (line2=strnword(line2,2));i1++) {
			itct=sscanf(line2,"%s",nm);
			CHECKS(itct==1,"format for neighbors: panel neigh1 neigh2 ...");
			pnllist[i1]=readpanelname(sim,srf,nm);
			CHECKS(pnllist[i1],"neighbor panel name '%s' not recognized",nm); }
		CHECKS(i1<maxpnllist,"too many neighbor panels listed in one line");
		er=surfsetneighbors(pnl,pnllist,i1,1);
		CHECKS(!er,"out of memory allocating neighbor list"); }

	else if(!strcmp(word,"neighbor_action") || !strcmp(word,"neighbour_action")) {
		CHECKS(srf,"need to enter surface name before neighbor_action");
		itct=sscanf(line2,"%s",nm);
		CHECKS(itct==1,"format for neighbor_action: 'hop' or 'stay'");
		i=0;
		if(!strcmp(nm,"hop")) i=1;
		else if(!strcmp(nm,"stay")) i=0;
		else CHECKS(0,"format for neighbor_action: 'hop' or 'stay'");
		surfsetneighhop(srf,i);
		CHECKS(!strnword(line2,2),"unexpected text following neighbor_action"); }

	else if(!strcmp(word,"unbounded_emitter")) {	// unbounded_emitter
		CHECKS(srf,"need to enter surface name before unbounded_emitter");
		itct=sscanf(line2,"%s %s %lf",facenm,nm,&f1);
		CHECKS(itct==3,"format for unbounded_emitter: face species amount position");
		face=surfstring2face(facenm);
		CHECKS(face==PFfront || face==PFback,"face must be 'front' or 'back'");
		i=stringfind(sim->mols->spname,sim->mols->nspecies,nm);
		CHECKS(i>0,"unrecognized species name");
		line2=strnword(line2,4);
		CHECKS(line2,"format for unbounded_emitter: face species amount position");
		itct=strreadnd(line2,dim,fltv1,NULL);
		CHECKS(itct==dim,"format for unbounded_emitter: face species amount position");
		er=surfaddemitter(srf,face,i,f1,fltv1,dim);
		CHECKS(er==0,"out of memory adding emitter");
		CHECKS(!strnword(line2,dim+1),"unexpected text following unbounded_emitter"); }

	else if(!strcmp(word,"action_front")) {				// action_front
		CHECKS(0,"the action_front statement has been replaced with action (remove the underscore)"); }

	else if(!strcmp(word,"action_back")) {				// action_back
		CHECKS(0,"the action_back statement has been replaced with action (remove the underscore)"); }

	else if(!strcmp(word,"action_both")) {				// action_both
		CHECKS(0,"the action_both statement has been replaced with action (remove the underscore)"); }

	else if(!strcmp(word,"polygon_front")) {			// polygon_front
		CHECKS(0,"the polygon_front statement has been replaced with polygon (remove the underscore)"); }

	else if(!strcmp(word,"polygon_back")) {			// polygon_back
		CHECKS(0,"the polygon_back statement has been replaced with polygon (remove the underscore)"); }

	else if(!strcmp(word,"polygon_both")) {			// polygon_both
		CHECKS(0,"the polygon_both statement has been replaced with polygon (remove the underscore)"); }

	else if(!strcmp(word,"color_front") || !strcmp(word,"colour_front")) {				// color_front
		CHECKS(0,"the color_front statement has been replaced with color (remove the underscore)"); }

	else if(!strcmp(word,"color_back") || !strcmp(word,"colour_back")) {					// color_back
		CHECKS(0,"the color_back statement has been replaced with color (remove the underscore)"); }

	else if(!strcmp(word,"color_both") || !strcmp(word,"colour_both")) {		// color, color_both
		CHECKS(0,"the color_both statement has been replaced with color (remove the underscore)"); }

	else {																				// unknown word
		CHECKS(0,"syntax error within surface block: statement not recognized"); }

	return srf;

 failure:
	simParseError(sim,pfp);
	return NULL; }


/* loadsurface */
int loadsurface(simptr sim,ParseFilePtr *pfpptr,char *line2) {
	ParseFilePtr pfp;
	char word[STRCHAR],errstring[STRCHAR];
	int done,pfpcode,firstline2;
	surfaceptr srf;

	pfp=*pfpptr;
	done=0;
	srf=NULL;
	firstline2=line2?1:0;

	while(!done) {
		if(pfp->lctr==0)
			simLog(sim,2," Reading file: '%s'\n",pfp->fname);
		if(firstline2) {
			strcpy(word,"name");
			pfpcode=1;
			firstline2=0; }
		else
			pfpcode=Parse_ReadLine(&pfp,word,&line2,errstring);
		*pfpptr=pfp;
		CHECKS(pfpcode!=3,"%s",errstring);

		if(pfpcode==0);															// already taken care of
		else if(pfpcode==2) {													// end reading
			done=1; }
		else if(!strcmp(word,"end_surface")) {				// end_surface
			CHECKS(!line2,"unexpected text following end_surface");
			return 0; }
		else if(!line2) {															// just word
			CHECKS(0,"unknown word or missing parameter"); }
		else {
			srf=surfreadstring(sim,pfp,srf,word,line2);
			CHECK(srf); }}

	CHECKS(0,"end of file encountered before end_surface statement");	// end of file

 failure:																					// failure
	if(ErrorType!=1) simParseError(sim,pfp);
	*pfpptr=pfp=NULL;
	return 1; }


/* surfupdateparams */
int surfupdateparams(simptr sim) {
	surfacessptr srfss;
	surfaceptr srf;
	int nspecies,i,s;
	double dt,**difc,sum,**difstep;
	enum MolecState ms1,ms2,ms3,ms4;
	enum PanelFace face,face2;
	surfactionptr actdetails;
	
	srfss=sim->srfss;
	
	if(sim->mols) {
		if(sim->mols->condition<SCok) return 2;

		nspecies=sim->mols->nspecies;
		dt=sim->dt;
		difc=sim->mols->difc;
		difstep=sim->mols->difstep;
		
		for(s=0;s<srfss->nsrf;s++) {					// set probabilities
			srf=srfss->srflist[s];
			for(i=1;i<nspecies;i++)
				for(ms1=(MolecState)0;ms1<MSMAX;ms1=(MolecState)(ms1+1))
					for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1))
						if(srf->action[i][ms1][face]==SAmult && srf->actdetails[i][ms1][face]) {
							actdetails=srf->actdetails[i][ms1][face];
							for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1))	// record actual probability of each event
								actdetails->srfprob[ms2]=srfcalcprob(sim,srf,i,ms1,face,ms2); }

			for(i=1;i<nspecies;i++)
				for(ms1=(MolecState)0;ms1<MSMAX;ms1=(MolecState)(ms1+1))
					for(face=(PanelFace)0;face<3;face=(PanelFace)(face+1)) {
						if(srf->action[i][ms1][face]==SAmult && srf->actdetails[i][ms1][face]) {
							actdetails=srf->actdetails[i][ms1][face];
							sum=0;
							for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1))			// find ms1->ms1 probability
								if(!srfsamestate(ms1,face,ms2,NULL)) sum+=actdetails->srfprob[ms2];
								if(sum>1.0) {
									actdetails->srfprob[ms1]=0;
									for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1))
										actdetails->srfprob[ms2]/=sum; 
								}
								else {
									srfsamestate(ms1,face,MSsoln,&ms3);
									actdetails->srfprob[ms3]=1.0-sum; 
								}
								for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1)) {				// record reverse probabilities
									srfreverseaction(ms1,face,ms2,&ms3,&face2,&ms4);
									if(face2!=PFboth && srf->actdetails[i][ms3][face2]) {
										actdetails->srfrevprob[ms2]=srf->actdetails[i][ms3][face2]->srfprob[ms4]; }}

								sum=0;																	// find cumulative probability
								for(ms2=(MolecState)0;ms2<MSMAX1;ms2=(MolecState)(ms2+1)) {
									sum+=actdetails->srfprob[ms2];
									actdetails->srfcumprob[ms2]=sum; }}}}
	
		surfsetemitterabsorption(sim); }
	
	return 0; }


/* surfupdatelists */
int surfupdatelists(simptr sim) {
	surfacessptr srfss;
	surfaceptr srf;
	int i,ll,maxmollist,s,totpanel,pindex,p;
	enum MolecState ms;
	enum SMLflag *newsrfmollist;
	double totarea,*areatable,area;
	panelptr *paneltable;
	enum PanelShape ps;

	srfss=sim->srfss;

	if(sim->mols) {
		if(sim->mols->condition<=SClists) return 2;

		if(sim->mols->nlist>srfss->maxmollist) {
			maxmollist=sim->mols->maxlist;								// allocate srfmollist array
			if(maxmollist>0) {
				newsrfmollist=(enum SMLflag*) calloc(maxmollist,sizeof(enum SMLflag));
				if(!newsrfmollist) return 1; }
			else newsrfmollist=NULL;
			for(ll=0;ll<maxmollist;ll++)
				newsrfmollist[ll]=SMLno;
			free(srfss->srfmollist);
			srfss->srfmollist=newsrfmollist;
			srfss->maxmollist=maxmollist; }
		srfss->nmollist=sim->mols->nlist;

		if(srfss->nmollist) {
			for(i=1;i<sim->mols->nspecies;i++)					// set srfmollist flags
				for(ms=(MolecState)0;ms<MSMAX;ms=(MolecState)(ms+1)) {
					ll=sim->mols->listlookup[i][ms];
					if(molismobile(sim,i,ms))
						srfss->srfmollist[ll] = (SMLflag)(srfss->srfmollist[ll] | SMLdiffuse);
					if(rxnisprod(sim,i,ms,1))
						srfss->srfmollist[ll] = (SMLflag)(srfss->srfmollist[ll] | SMLreact);
					if(ms!=MSsoln)
						srfss->srfmollist[ll] = (SMLflag)(srfss->srfmollist[ll] | SMLsrfbound); }}

		for(s=0;s<srfss->nsrf;s++) {									// area lookup tables
			srf=srfss->srflist[s];
			totarea=surfacearea(srf,sim->dim,&totpanel);
			if(totpanel) {
				areatable=(double*)calloc(totpanel,sizeof(double));
				if(!areatable) return 1;
				paneltable=(panelptr*)calloc(totpanel,sizeof(panelptr));
				if(!paneltable) {free(areatable);return 1;}
				pindex=0;
				area=0;
				for(ps=(PanelShape)0;ps<PSMAX;ps=(PanelShape)(ps+1)) {
					for(p=0;p<srf->npanel[ps];p++) {
						area+=panelarea(srf->panels[ps][p],sim->dim);
						areatable[pindex]=area;
						paneltable[pindex]=srf->panels[ps][p];
						pindex++; }}
				srf->totarea=totarea;
				srf->totpanel=totpanel;
				free(srf->areatable);
				srf->areatable=areatable;
				free(srf->paneltable);
				srf->paneltable=paneltable; }}}

	return 0; }


/* surfupdate */
int surfupdate(simptr sim) {
	int er;
	surfacessptr srfss;

	srfss=sim->srfss;
	if(srfss) {
		if(srfss->condition<=SClists) {
			er=surfupdatelists(sim);
			if(er) return er;
			surfsetcondition(srfss,SCparams,1); }
		if(srfss->condition==SCparams) {
			er=surfupdateparams(sim);
			if(er) return er;
			surfsetcondition(srfss,SCok,1); }}
	return 0; }


/******************************************************************************/
/************************* core simulation functions **************************/
/******************************************************************************/


/* panelside */
enum PanelFace panelside(double* pt,panelptr pnl,int dim,double *distptr,int strict) {
	enum PanelFace face;
	double **point,*front,dist,cylnorm[3];
	int d;

	point=pnl->point;
	front=pnl->front;
	dist=0;

	if(pnl->ps==PSrect) {														// rectangle
		d=(int)front[1];
		dist=front[0]*(pt[d]-point[0][d]); }
	else if(pnl->ps==PStri || pnl->ps==PSdisk) {			// triangle, disk
		for(d=0;d<dim;d++) dist+=(pt[d]-point[0][d])*front[d]; }
	else if(pnl->ps==PSsph || pnl->ps==PShemi) {			// sphere, hemisphere
		for(d=0;d<dim;d++) dist+=(pt[d]-point[0][d])*(pt[d]-point[0][d]);
		dist=front[0]*(sqrt(dist)-point[1][0]); }
	else if(pnl->ps==PScyl) {												// cylinder
		if(dim==2) {
			dist+=(pt[0]-point[0][0])*front[0];
			dist+=(pt[1]-point[0][1])*front[1];
			dist=front[2]*(fabs(dist)-point[2][0]); }
		else {
			dist=Geo_LineNormal3D(point[0],point[1],pt,cylnorm);
			dist=front[2]*(dist-point[2][0]); }}
	else
		dist=0;

	if(dist>0) face=PFfront;
	else if(!strict || dist<0) face=PFback;
	else face=PFnone;

	if(distptr) *distptr=dist;
	return face; }


/* panelnormal */
void panelnormal(panelptr pnl,double *pos,enum PanelFace face,int dim,double *norm) {
	int d;
	double **point,*front;

	point=pnl->point;
	front=pnl->front;

	if(face!=PFback) face=PFfront;

	if(pnl->ps==PSrect) {
		for(d=0;d<dim;d++) norm[d]=0;
		norm[(int)front[1]]=((face==PFfront && front[0]==1) || (face==PFback && front[0]==-1))?1.0:-1.0; }
	else if(pnl->ps==PStri || pnl->ps==PSdisk) {
		if(face==PFfront)
			for(d=0;d<dim;d++) norm[d]=front[d];
		else
			for(d=0;d<dim;d++) norm[d]=-front[d]; }
	else if(pnl->ps==PSsph || pnl->ps==PShemi) {
		Geo_SphereNormal(point[0],pos,((face==PFfront && front[0]==1) || (face==PFback && front[0]==-1))?1:-1,dim,norm); }
	else if(pnl->ps==PScyl) {
		if(dim==2) {
			Geo_LineNormal2D(point[0],point[1],pos,norm);
			if((face==PFback && front[2]==1) || (face==PFfront && front[2]==-1))
				for(d=0;d<dim;d++) norm[d]*=-1; }
		else if(dim==3) {
			Geo_LineNormal3D(point[0],point[1],pos,norm);
			if((face==PFback && front[2]==1) || (face==PFfront && front[2]==-1))
				for(d=0;d<dim;d++) norm[d]*=-1; }}

	return; }


/* lineXpanel */
int lineXpanel(double *pt1,double *pt2,panelptr pnl,int dim,double *crsspt,enum PanelFace *face1ptr,enum PanelFace *face2ptr,double *crossptr,double *cross2ptr,int *veryclose) {
	surfaceptr srf;
	double **point,*front,dist1,dist2;
	double dot,cross,cross2,nrdist,nrpos;
	int intsct,d;
	enum PanelFace face1,face2,facein;

	srf=pnl->srf;
	point=pnl->point;
	front=pnl->front;
	face1=panelside(pt1,pnl,dim,&dist1,0);
	face2=panelside(pt2,pnl,dim,&dist2,0);
	cross=cross2=-1;

	if(pnl->ps==PSrect) {														// rectangle
		
		if(face1==face2) return 0;
		cross=dist1/(dist1-dist2);
		// printf("lineXpanel line 3308, pt1: [0]=%f,[1]=%f,[2]=%f, pt2: [0]=%f,[1]=%f,[2]=%f\n", pt1[0],pt1[1],pt1[2],pt2[0],pt2[1],pt2[2]);
		// printf("lineXpanel line 3309, pname=%s, face1=%d, face2=%d, dist1=%f, dist2=%f, cross=%f\n", pnl->pname, face1,face2, dist1, dist2, cross);
		for(d=0;d<dim;d++) crsspt[d]=pt1[d]+cross*(pt2[d]-pt1[d]);
		if(dim==1) intsct=1;
		else if(dim==2) {
			d=(int)front[2];
			intsct=((point[0][d]<=crsspt[d] && crsspt[d]<=point[1][d]) || (point[1][d]<=crsspt[d] && crsspt[d]<=point[0][d])); }
		else {
			d=(int)front[2];
			// printf("lineXpanel, line 3317, d=%d, point[0][d]=%f, point[1][d]=%f\n", d, point[0][d], point[1][d]);
			intsct=((point[0][d]<=crsspt[d] && crsspt[d]<=point[1][d]) || (point[1][d]<=crsspt[d] && crsspt[d]<=point[0][d]));
			d=(d+1)%3;
			if(d==(int)front[1]) d=(d+1)%3;
			// printf("lineXpanel, line 3321, d=%d, point[1][d]=%f, point[2][d]=%f\n", d, point[1][d], point[2][d]);
			intsct=intsct && ((point[1][d]<=crsspt[d] && crsspt[d]<=point[2][d]) || (point[2][d]<=crsspt[d] && crsspt[d]<=point[1][d])); }}

	else if(pnl->ps==PStri) {												// triangle
		if(face1==face2) return 0;
		cross=dist1/(dist1-dist2);
		for(d=0;d<dim;d++) crsspt[d]=pt1[d]+cross*(pt2[d]-pt1[d]);
		if(dim==1) intsct=1;
		else if(dim==2) {
			//intsct=((point[0][0]<=crsspt[0] && crsspt[0]<=point[1][0]) || (point[1][0]<=crsspt[0] && crsspt[0]<=point[0][0]));
			//intsct=intsct && ((point[0][1]<=crsspt[1] && crsspt[1]<=point[1][1]) || (point[1][1]<=crsspt[1] && crsspt[1]<=point[0][1])); }
			//VCELL
			double v1[2] = {point[1][0] - point[0][0], point[1][1] - point[0][1]};
			double v2[2] = {crsspt[0] - point[0][0], crsspt[1] - point[0][1]};
			intsct = v1[0] * v2[0] + v1[1] * v2[1] >= 0;
			double v3[2] = {crsspt[0] - point[1][0], crsspt[1] - point[1][1]};
			intsct = intsct && (v1[0] * v3[0] + v1[1] * v3[1] <= 0); 
		}
		else {
			intsct=Geo_PtInTriangle2(point,crsspt); }}

	else if(pnl->ps==PSsph || pnl->ps==PShemi) {		// sphere, hemisphere
		facein=front[0]>0?PFback:PFfront;
		if(face1==facein && face1==face2) return 0;
		cross=Geo_LineXSphs(pt1,pt2,point[0],point[1][0],dim,&cross2,&nrdist,&nrpos);
		if(face1==face2 && (nrdist>point[1][0] || nrpos<0 || nrpos>1)) return 0;
		if(face1==facein)
			cross=cross2;
		for(d=0;d<dim;d++) crsspt[d]=pt1[d]+cross*(pt2[d]-pt1[d]);
		if(pnl->ps==PSsph) intsct=1;
		else {
			dot=0;
			for(d=0;d<dim;d++) dot+=(crsspt[d]-point[0][d])*point[2][d];
			intsct=(dot<=0);
			if(!intsct && face1==face2) {
				cross=cross2;
				face1=(face1==PFfront)?PFback:PFfront;
				for(d=0;d<dim;d++) crsspt[d]=pt1[d]+cross*(pt2[d]-pt1[d]);
				dot=0;
				for(d=0;d<dim;d++) dot+=(crsspt[d]-point[0][d])*point[2][d];
				intsct=(dot<=0); }}}

	else if(pnl->ps==PScyl) {									// cylinder
		facein=(int)front[2]>0?PFback:PFfront;
		if(face1==facein && face1==face2) return 0;
		if(dim==2) cross=Geo_LineXCyl2s(pt1,pt2,point[0],point[1],front,point[2][0],&cross2,&nrdist,&nrpos);
		else cross=Geo_LineXCyls(pt1,pt2,point[0],point[1],point[2][0],&cross2,&nrdist,&nrpos);
		if(face1==face2 && (nrdist>point[2][0] || nrpos<0 || nrpos>1)) return 0;
		if(face1==facein)
			cross=cross2;
		for(d=0;d<dim;d++) crsspt[d]=pt1[d]+cross*(pt2[d]-pt1[d]);
		intsct=Geo_PtInSlab(point[0],point[1],crsspt,dim);
		if(!intsct && face1==face2) {
			cross=cross2;
			face1=(face1==PFfront)?PFback:PFfront;
			for(d=0;d<dim;d++) crsspt[d]=pt1[d]+cross*(pt2[d]-pt1[d]);
			intsct=Geo_PtInSlab(point[0],point[1],crsspt,dim); }}

	else if(pnl->ps==PSdisk) {												// disk
		if(face1==face2) return 0;
		cross=dist1/(dist1-dist2);
		for(d=0;d<dim;d++) crsspt[d]=pt1[d]+cross*(pt2[d]-pt1[d]);
		dot=0;
		for(d=0;d<dim;d++) dot+=(crsspt[d]-point[0][d])*(crsspt[d]-point[0][d]);
		intsct=(dot<=point[1][0]*point[1][0]); }

	else
		intsct=0;

	if(face1ptr) *face1ptr=face1;
	if(face2ptr) *face2ptr=face2;
	if(crossptr) *crossptr=cross;
	if(cross2ptr) *cross2ptr=cross2;
	if(veryclose) {
		*veryclose=0;
		if(dist1<VERYCLOSE) *veryclose+=1;
		if(dist2<VERYCLOSE) *veryclose+=2; }
	return intsct; }


/* lineexitpanel */
int lineexitpanel(double *pt1,double *pt2,panelptr pnl,int dim,double *pnledgept,int *exitside) {
	int d;
	double **point,*front,dist,end1[2],end2[2];

	for(d=0;pt1[d]==pt2[d] && d<dim;d++);
	if(d==dim) return 1;						// return 1 if pt1==pt2

	point=pnl->point;
	front=pnl->front;
	*exitside=1;

	if(pnl->ps==PSrect) {
		if(dim==1)
			pnledgept[0]=point[0][0];
		else if(dim==2)
			Geo_LineExitLine2(pt1,pt2,point[0],point[1],pnledgept,exitside);
		else if(dim==3)
			Geo_LineExitRect(pt1,pt2,front,point[0],point[2],pnledgept,exitside); }

	else if(pnl->ps==PStri) {
		if(dim==1)
			pnledgept[0]=point[0][0];
		else if(dim==2)
			Geo_LineExitLine2(pt1,pt2,point[0],point[1],pnledgept,exitside);
		else if(dim==3)
			Geo_LineExitTriangle2(pt1,pt2,point,pnledgept,exitside); }

	else if(pnl->ps==PShemi) {
		if(dim==2)
			Geo_LineExitArc2(pt1,pt2,point[0],point[1][0],point[2],pnledgept,exitside);
		else if(dim==3)
			Geo_LineExitHemisphere(pt1,pt2,point[0],point[1][0],point[2],pnledgept); }

	else if(pnl->ps==PScyl) {
		if(dim==2) {
			dist=0;
			dist+=(pt1[0]-point[0][0])*front[0];
			dist+=(pt1[1]-point[0][1])*front[1];
			dist=(dist>0)?point[2][0]:-point[2][0];
			end1[0]=point[0][0]+dist*front[0];
			end1[1]=point[0][1]+dist*front[1];
			end2[0]=point[1][0]+dist*front[0];
			end2[1]=point[1][1]+dist*front[1];
			Geo_LineExitLine2(pt1,pt2,end1,end2,pnledgept,exitside); }
		else if(dim==3)
			Geo_LineExitCylinder(pt1,pt2,point[0],point[1],point[2][0],pnledgept,exitside); }

	else if(pnl->ps==PSdisk) {
		if(dim==2) {
			end1[0]=point[0][0]+point[1][0]*front[1];
			end1[1]=point[0][1]-point[1][0]*front[0];
			end2[0]=point[0][0]-point[1][0]*front[1];
			end2[1]=point[0][1]+point[1][0]*front[0];
			Geo_LineExitLine2(pt1,pt2,end1,end2,pnledgept,exitside); }
		else if(dim==3)
			Geo_LineExitSphere(pt1,pt2,pnl->point[0],pnl->point[1][0],pnledgept); }

	return 0; }


/* paneledgenormal */
void paneledgenormal(panelptr pnl,double *pnledgept,int dim,int edgenum,double *normal) {
	double **point,*front,normx,normy,normz,dot;

	point=pnl->point;
	front=pnl->front;

	if(pnl->ps==PSrect) {
		if(edgenum==0) edgenum=1;
		if(dim==2) {
			normal[0]=point[edgenum+1][0];
			normal[1]=point[edgenum+1][1]; }
		else {
			normal[0]=point[edgenum+3][0];
			normal[1]=point[edgenum+3][1];
			normal[2]=point[edgenum+3][2]; }}

	else if(pnl->ps==PStri) {
		if(edgenum==0) edgenum=1;
		if(dim==2) {
			normal[0]=point[edgenum+1][0];
			normal[1]=point[edgenum+1][1]; }
		else {
			normal[0]=point[edgenum+2][0];
			normal[1]=point[edgenum+2][1];
			normal[2]=point[edgenum+2][2]; }}

	else if(pnl->ps==PSsph) {
		if(dim==2) {
			normx=pnledgept[0]-point[0][0];
			normy=pnledgept[1]-point[0][1];
			dot=1.0/sqrt(normx*normx+normy*normy);
			normal[0]=-dot*normy;
			normal[1]=dot*normx; }
		else {
			normx=pnledgept[0]-point[0][0];
			normy=pnledgept[1]-point[0][1];
			normz=pnledgept[2]-point[0][2];
			if(normz<normx) {
				dot=1.0/sqrt(normx*normx+normy*normy);
				normal[0]=-dot*normy;
				normal[1]=dot*normx;
				normal[2]=0; }
			else {
				dot=1.0/sqrt(normx*normx+normz*normz);
				normal[0]=-dot*normz;
				normal[1]=0;
				normal[2]=dot*normx; }}}

	else if(pnl->ps==PShemi) {
		if(dim==2) {
			if(edgenum==0) {
				normx=pnledgept[0]-point[0][0];
				normy=pnledgept[1]-point[0][1];
				dot=1.0/sqrt(normx*normx+normy*normy);
				normal[0]=-dot*normy;
				normal[1]=dot*normx; }
			else {
				normal[0]=point[2][0];
				normal[1]=point[2][1]; }}
		else {
			if(edgenum==0) {
				normx=pnledgept[0]-point[0][0];
				normy=pnledgept[1]-point[0][1];
				normz=pnledgept[2]-point[0][2];
				if(normz<normx) {
					dot=1.0/sqrt(normx*normx+normy*normy);
					normal[0]=-dot*normy;
					normal[1]=dot*normx;
					normal[2]=0; }
				else {
					dot=1.0/sqrt(normx*normx+normz*normz);
					normal[0]=-dot*normz;
					normal[1]=0;
					normal[2]=dot*normx; }}
			else {
				normal[0]=point[2][0];
				normal[1]=point[2][1];
				normal[2]=point[2][2]; }}}

	else if(pnl->ps==PScyl) {
		if(edgenum==0) edgenum=1;
		if(dim==2) {
			normal[0]=point[edgenum+2][0];
			normal[1]=point[edgenum+2][1]; }
		else {
			normal[0]=point[edgenum+2][0];
			normal[1]=point[edgenum+2][1];
			normal[2]=point[edgenum+2][2]; }}

	else if(pnl->ps==PSdisk) {
		if(dim==2) {
			if(edgenum==1 || edgenum==0) {
				normal[0]=front[1];
				normal[1]=-front[0]; }
			else {
				normal[0]=-front[1];
				normal[1]=front[0]; }}
		else {
			normx=pnledgept[0]-point[0][0];
			normy=pnledgept[1]-point[0][1];
			normz=pnledgept[2]-point[0][2];
			dot=1.0/sqrt(normx*normx+normy*normy+normz*normz);
			normal[0]=dot*normx;
			normal[1]=dot*normy;
			normal[2]=dot*normz; }}

	return; }


/* ptinpanel */
int ptinpanel(double *pt,panelptr pnl,int dim) {
	surfaceptr srf;
	double **point,*front;
	double len2,dot;
	int intsct,d;

	srf=pnl->srf;
	point=pnl->point;
	front=pnl->front;

	if(pnl->ps==PSrect) {														// rectangle
		if(dim==1) intsct=1;
		else if(dim==2) {
			d=(int)front[2];
			intsct=((point[0][d]<=pt[d] && pt[d]<=point[1][d]) || (point[1][d]<=pt[d] && pt[d]<=point[0][d])); }
		else {
			d=(int)front[2];
			intsct=((point[0][d]<=pt[d] && pt[d]<=point[1][d]) || (point[1][d]<=pt[d] && pt[d]<=point[0][d]));
			d=(d+1)%3;
			if(d==(int)front[1]) d=(d+1)%3;
			intsct=intsct && ((point[1][d]<=pt[d] && pt[d]<=point[2][d]) || (point[2][d]<=pt[d] && pt[d]<=point[1][d])); }}

	else if(pnl->ps==PStri) {												// triangle
		if(dim==1) intsct=1;
		else if(dim==2)
			intsct=Geo_PtInSlab(point[0],point[1],pt,2);
		else
			intsct=Geo_PtInTriangle2(point,pt); }

	else if(pnl->ps==PSsph || pnl->ps==PShemi) {		// sphere, hemisphere
		if(pnl->ps==PSsph) intsct=1;
		else {
			dot=0;
			for(d=0;d<dim;d++) dot+=(pt[d]-point[0][d])*point[2][d];
			intsct=(dot<=0); }}

	else if(pnl->ps==PScyl) {									// cylinder
		intsct=Geo_PtInSlab(point[0],point[1],pt,dim); }

	else if(pnl->ps==PSdisk) {												// disk
		len2=0;
		for(d=0;d<dim;d++) len2+=(pt[d]-point[0][d])*(pt[d]-point[0][d]);
		if(len2<=point[1][0]*point[1][0]) intsct=1;
		else {
			dot=0;
			for(d=0;d<dim;d++) dot+=(pt[d]-point[0][d])*front[d];
			len2-=dot*dot;
			intsct=(len2<=point[1][0]*point[1][0]); }}

	else
		intsct=0;

	return intsct; }


/* surfaction */
enum SrfAction surfaction(surfaceptr srf,enum PanelFace face,int ident,enum MolecState ms,int *i2ptr,enum MolecState *ms2ptr,int site) {
	enum SrfAction act;
	enum MolecState ms2,ms3;
	double r,*srfcumprob;
	int i2;
	surfactionptr actdetails;

	i2=ident;
	ms2=ms;
	act=srf->action[ident][ms][face];
	
	if(act==SAmult) {
		actdetails=srf->actdetails[ident][ms][face];
		srfcumprob=actdetails->srfcumprob;
		r=randCOD();
		ms2=MSnone;
		for(ms3=(MolecState)0;ms3<MSMAX1 && ms2==MSnone;ms3=(MolecState)(ms3+1))
			if(r<srfcumprob[ms3]) ms2=ms3;
		i2=actdetails->srfnewspec[ms2];

		if(i2==0) act=SAabsorb;
		else if(ms==MSsoln) {
			if(face==PFfront) {														// solution, front
				if(ms2==MSsoln) act=SAreflect;
				else if(ms2==MSbsoln) act=SAtrans;
				else act=SAadsorb; }
			else {																	// solution, back
				if(ms2==MSsoln) act=SAtrans;
				else if(ms2==MSbsoln) act=SAreflect;
				else act=SAadsorb; }}
		else {
			if(face==PFnone) {														// bound, none
				if(ms2==ms) act=SAno;
				else if(ms2==MSsoln || ms2==MSbsoln) {
					if(actdetails->srfrevprob[ms2]>0) act=SArevdes;
					else act=SAirrevdes; }
				else
					act=SAflip; }
			else if(face==PFfront) {												// bound, front
				if(ms2==MSsoln) act=SAreflect;
				else if(ms2==MSbsoln) act=SAtrans;
				else act=SAadsorb; }
			else {																	// bound, back
				if(ms2==MSsoln) act=SAtrans;
				else if(ms2==MSbsoln) act=SAreflect;
				else act=SAadsorb; }}}

	if(i2ptr) *i2ptr=i2;
	if(ms2ptr) *ms2ptr=ms2;
	return act; }


/* rxnXsurface */
int rxnXsurface(simptr sim,moleculeptr mptr1,moleculeptr mptr2,int rxn_site_indx1,int rxn_site_indx2) {
	int dim,p,i1,i2,result;
	double *pos1,*pos2,dc1,dc2,rxnpt,crsspt[3],cross,cross2,dsum;
	boxptr bptr;
	panelptr pnl;
	enum PanelFace face1,face2,facein;
	enum MolecState ms1,ms2;

	if(!sim->srfss) return 0;
	dim=sim->dim;
	i1=mptr1->ident;
	i2=mptr2->ident;
	ms1=mptr1->mstate;
	ms2=mptr2->mstate;
	pos1=mptr1->pos;
	pos2=mptr2->pos;
	// dc1=sim->mols->difc[i1][ms1];
	// dc2=sim->mols->difc[i2][ms2];
	dsum=MolCalcDifcSum(sim,mptr1,mptr2,&dc1,&dc2);
	if(dc1==0 && dc2==0) dc1=dc2=1;
	// rxnpt=dc1/(dc1+dc2);
	rxnpt=dc1/dsum;

	result=0;
	for(bptr=pos2box(sim,mptr1->pos);bptr && result==0;bptr=line2nextbox(sim,pos1,pos2,bptr)) {
		for(p=0;p<bptr->npanel && result==0;p++) {
			pnl=bptr->panel[p];
			if(lineXpanel(pos1,pos2,pnl,dim,crsspt,&face1,&face2,&cross,&cross2,NULL)) {
				if((ms1==MSup || ms1==MSdown) && pnl==mptr1->pnl) result=0;
				else if((ms2==MSup || ms2==MSdown) && pnl==mptr2->pnl) result=0;
				else if(face1!=face2) {					// opposite sides of a surface
					if(rxnpt<cross || (rxnpt==cross && face1==PFback)) result=!(surfaction(pnl->srf,face2,i2,ms2,NULL,NULL,rxn_site_indx2)==SAtrans);
					else result=!(surfaction(pnl->srf,face1,i1,ms1,NULL,NULL,rxn_site_indx1)==SAtrans); }
				else {											// across a sphere, cylinder, etc.
					facein=(face1==PFfront)?PFback:PFfront;
					if(rxnpt<cross || (rxnpt==cross && face1==PFback)) result=!(surfaction(pnl->srf,face2,i2,ms2,NULL,NULL,rxn_site_indx2)==SAtrans && surfaction(pnl->srf,facein,i2,ms2,NULL,NULL,rxn_site_indx2)==SAtrans);
					else if(rxnpt<cross2 || (rxnpt==cross2 && face1==PFfront)) result=!(surfaction(pnl->srf,face1,i1,ms1,NULL,NULL,rxn_site_indx1)==SAtrans && surfaction(pnl->srf,face2,i2,ms2,NULL,NULL,rxn_site_indx2)==SAtrans);
					else result=!(surfaction(pnl->srf,face1,i1,ms1,NULL,NULL,rxn_site_indx1)==SAtrans && surfaction(pnl->srf,facein,i1,ms1,NULL,NULL,rxn_site_indx1)==SAtrans); }}}}
	return result; }



/* fixpt2panel */
void fixpt2panel(double *pt,panelptr pnl,int dim,enum PanelFace face,double epsilon) {
	int d,sign;
	double **point,*front,norm[3],dist,factor,dist2;
	enum PanelFace faceat;

	point=pnl->point;
	front=pnl->front;

	faceat=panelside(pt,pnl,dim,&dist,1);
	if((faceat==face || face==PFnone) && fabs(dist)<=epsilon) return;

	if(pnl->ps==PSrect) {
		for(d=0;d<dim;d++) norm[d]=0;
		norm[(int)front[1]]=front[0]; }
	else if(pnl->ps==PStri || pnl->ps==PSdisk) {
		for(d=0;d<dim;d++) norm[d]=front[d]; }
	else if(pnl->ps==PSsph || pnl->ps==PShemi) {
		Geo_SphereNormal(point[0],pt,(int)front[0],dim,norm); }
	else if(pnl->ps==PScyl) {
		if(dim==2) {
			dist2=0;
			for(d=0;d<dim;d++) dist2+=(pt[d]-point[0][d])*front[d];
			sign=((dist2>0 && front[2]==1) || (dist2<0 && front[2]==-1))?1:-1;
			norm[0]=sign*front[0];
			norm[1]=sign*front[1]; }
		else if(dim==3) {
			Geo_LineNormal3D(point[0],point[1],pt,norm);
			if(front[2]==-1)
				for(d=0;d<dim;d++) norm[d]*=-1; }}
	else {
		for(d=0;d<dim;d++) norm[d]=0;
		norm[0]=1; }

	for(d=0;d<dim;d++) pt[d]-=dist*norm[d];
	if(face==PFnone || face==PFboth) return;
	sign=(face==PFfront)?1:-1;
	for(factor=1.0;face!=panelside(pt,pnl,dim,NULL,1);factor*=2.0) {
		for(d=0;d<dim;d++) pt[d]+=sign*factor*DBL_EPSILON*norm[d]; }
	return; }


/* movept2panel */
void movept2panel(double *pt,panelptr pnl,int dim,double margin) {
	double **point,*front,inpt0[3],inpt1[3],inpt2[3],*inpoint[3];
	double lo,hi,dot;
	int d;

	point=pnl->point;
	front=pnl->front;

	if(pnl->ps==PSrect) {														// rectangle
		if(dim==1);
		else if(dim==2) {
			d=(int)front[2];
			lo=((point[0][d]<point[1][d])?point[0][d]:point[1][d])+margin;
			hi=((point[0][d]<point[1][d])?point[1][d]:point[0][d])-margin;
			if(pt[d]<lo) pt[d]=lo;
			else if(pt[d]>hi) pt[d]=hi; }
		else {
			d=(int)front[2];
			lo=((point[0][d]<point[1][d])?point[0][d]:point[1][d])+margin;
			hi=((point[0][d]<point[1][d])?point[1][d]:point[0][d])-margin;
			if(pt[d]<lo) pt[d]=lo;
			else if(pt[d]>hi) pt[d]=hi;
			d=(d+1)%3;
			if(d==(int)front[1]) d=(d+1)%3;
			lo=((point[0][d]<point[3][d])?point[0][d]:point[3][d])+margin;
			hi=((point[0][d]<point[3][d])?point[3][d]:point[0][d])-margin;
			if(pt[d]<lo) pt[d]=lo;
			else if(pt[d]>hi) pt[d]=hi; }}
	else if(pnl->ps==PStri) {												// triangle
		if(dim==1);
		else if(dim==2) {
			Geo_InsidePoints2(point[0],point[1],margin,inpt0,inpt1,dim);
			Geo_NearestSlabPt(inpt0,inpt1,pt,pt,dim); }
		else {
			inpoint[0]=inpt0;
			inpoint[1]=inpt1;
			inpoint[2]=inpt2;
			Geo_InsidePoints32(point,margin,inpoint);
			Geo_NearestTriPt2(inpoint,point+3,front,pt,pt); }}
	else if(pnl->ps==PSsph || pnl->ps==PShemi) {		// sphere, hemisphere
		if(pnl->ps==PSsph);
		else {
			dot=0;
			for(d=0;d<dim;d++) dot+=(pt[d]-point[0][d])*point[2][d];
			if(dot>0)
				for(d=0;d<dim;d++) pt[d]-=dot*point[2][d]+margin; }}
	else if(pnl->ps==PScyl) {												// cylinder
		Geo_InsidePoints2(point[0],point[1],margin,inpt0,inpt1,dim);
		Geo_NearestSlabPt(inpt0,inpt1,pt,pt,dim); }
	else if(pnl->ps==PSdisk)												// disk
		Geo_NearestCylPt(point[0],front,point[1][0]-margin,dim,pt,pt);

	return; }


/* closestpanelpt */
int closestpanelpt(panelptr pnl,int dim,double *testpt,double *pnlpt) {
	double **point,*front;
	double dot,x,y,end1point[2],end2point[2];
	int d,edgenum;

	point=pnl->point;
	front=pnl->front;
	edgenum=0;

	if(pnl->ps==PSrect) {														// rectangle
		if(dim==1) {
			pnlpt[0]=point[0][0];
			edgenum=1; }
		else if(dim==2) {
			d=(int)front[1];
			pnlpt[d]=point[0][d];
			d=(int)front[2];
			if((point[0][d]<testpt[d] && testpt[d]<point[1][d]) || (point[1][d]<testpt[d] && testpt[d]<point[0][d]))
				pnlpt[d]=testpt[d];
			else if(fabs(testpt[d]-point[0][d])<fabs(testpt[d]-point[1][d])) {
				pnlpt[d]=point[0][d];
				edgenum=1; }
			else {
				pnlpt[d]=point[1][d];
				edgenum=2; }}
		else {
			d=(int)front[1];
			pnlpt[d]=point[0][d];
			d=(int)front[2];
			if((point[0][d]<testpt[d] && testpt[d]<point[1][d]) || (point[1][d]<testpt[d] && testpt[d]<point[0][d]))
				pnlpt[d]=testpt[d];
			else if(fabs(testpt[d]-point[0][d])<fabs(testpt[d]-point[1][d])) {
				pnlpt[d]=point[0][d];
				edgenum=4; }
			else {
				pnlpt[d]=point[1][d];
				edgenum=2; }
			d=(d+1)%3;
			if(d==(int)front[1]) d=(d+1)%3;
			if((point[0][d]<testpt[d] && testpt[d]<point[3][d]) || (point[3][d]<testpt[d] && testpt[d]<point[0][d]))
				pnlpt[d]=testpt[d];
			else if(fabs(testpt[d]-point[0][d])<fabs(testpt[d]-point[3][d])) {
				pnlpt[d]=point[0][d];
				edgenum=1; }
			else {
				pnlpt[d]=point[3][d];
				edgenum=3; }}}

	else if(pnl->ps==PStri) {												// triangle
		if(dim==1) {
			pnlpt[0]=point[0][0];
			edgenum=1; }
		else if(dim==2) {
			edgenum=Geo_NearestLineSegPt(point[0],point[1],testpt,pnlpt,2); }
		else {
			edgenum=Geo_NearestTrianglePt2(point,front,testpt,pnlpt); }}

	else if(pnl->ps==PSsph) {												// sphere
		if(dim==1) {
			if(testpt[0]>point[0][0])
				pnlpt[0]=point[0][0]+point[1][0];
			else
				pnlpt[0]=point[0][0]-point[1][0];
			edgenum=1; }
		else {
			Geo_NearestSpherePt(point[0],point[1][0],(int)front[0],dim,testpt,pnlpt);
			edgenum=0; }}

	else if(pnl->ps==PShemi) {											// hemisphere
		dot=0;
		for(d=0;d<dim;d++) dot+=(testpt[d]-point[0][d])*point[2][d];
		if(dot<0) {
			Geo_NearestSpherePt(point[0],point[1][0],(int)front[0],dim,testpt,pnlpt);
			edgenum=0; }
		else if(dim==2) {
			x=point[2][1];
			y=-point[2][0];
			dot=(testpt[0]-point[0][0])*x+(testpt[1]-point[0][1])*y;
			if(dot>0) {
				edgenum=2;
				dot=1; }
			else {
				edgenum=1;
				dot=-1; }
			pnlpt[0]=point[0][0]+dot*point[1][0]*x;
			pnlpt[1]=point[0][1]+dot*point[1][0]*y; }
		else {
			Geo_NearestRingPt(point[0],point[2],point[1][0],3,testpt,pnlpt);
			edgenum=1; }}

	else if(pnl->ps==PScyl) {												// cylinder
		if(dim==2) {
			edgenum=Geo_NearestLineSegPt(point[0],point[1],testpt,pnlpt,2);
			dot=(testpt[0]-point[0][0])*front[0]+(testpt[1]-point[0][1])*front[1];
			dot=(dot>0)?1:-1;
			pnlpt[0]+=dot*point[2][0]*front[0];
			pnlpt[1]+=dot*point[2][0]*front[1]; }
		else {
			edgenum=Geo_NearestCylinderPt(point[0],point[1],point[2][0],3,testpt,pnlpt); }}

	else if(pnl->ps==PSdisk) {										// disk
		if(dim==2) {
			end1point[0]=point[0][0]+point[1][0]*front[1];
			end1point[1]=point[0][1]-point[1][0]*front[0];
			end2point[0]=point[0][0]-point[1][0]*front[1];
			end2point[1]=point[0][1]+point[1][0]*front[0];
			edgenum=Geo_NearestLineSegPt(end1point,end2point,testpt,pnlpt,2); }
		else {
			edgenum=Geo_NearestDiskPt(point[0],front,point[1][0],dim,testpt,pnlpt); }}

	return edgenum; }


/* closestsurfacept */
double closestsurfacept(surfaceptr srf,int dim,double *testpt,double *pnlpt,panelptr *pnlptr) {
	enum PanelShape ps;
	int p,d;
	panelptr pnl,pnlbest;
	double dist2,dist2best,pnlpt1[DIMMAX],pnlptbest[DIMMAX];

	dist2best=DBL_MAX;
	pnlbest=NULL;
	for(ps=PanelShape(0);ps<PSMAX;ps=PanelShape(ps + 1))
		for(p=0;p<srf->npanel[ps];p++) {
			pnl=srf->panels[ps][p];
			closestpanelpt(pnl,dim,testpt,pnlpt1);
			dist2=0;
			for(d=0;d<dim;d++) dist2+=(testpt[d]-pnlpt1[d])*(testpt[d]-pnlpt1[d]);
			if(dist2<dist2best) {
				dist2best=dist2;
				pnlbest=pnl;
				for(d=0;d<dim;d++) pnlptbest[d]=pnlpt1[d]; }}

	if(dist2best==DBL_MAX) return -1;
	if(pnlpt)
		for(d=0;d<dim;d++) pnlpt[d]=pnlptbest[d];
	if(pnlptr)
		*pnlptr=pnlbest;
	return sqrt(dist2best); }


/* movemol2closepanel */
void movemol2closepanel(simptr sim,moleculeptr mptr,int dim,double epsilon,double neighdist,double margin) {
	int nn,p,d,er,exitside,iter,edgenum,newedge;
	double pt[DIMMAX],newedgept[DIMMAX],pnledgept[DIMMAX],oldnormal[DIMMAX],newnormal[DIMMAX],newedgenormal[DIMMAX],newpos[DIMMAX];
	double dot,dist2,len;
	int thetasign;
	panelptr pnl,newpnl;
	enum PanelFace face;
	enum MolecState ms;

	ms=mptr->mstate;
	if(ms==MSfront) face=PFfront;
	else if(ms==MSback) face=PFback;
	else face=PFnone;
	pnl=mptr->pnl;

	fixpt2panel(mptr->pos,pnl,dim,face,epsilon);			// move mptr->pos to panel plane
	if(ptinpanel(mptr->pos,pnl,dim)) return;					// nothing to worry about so return

	if(ptinpanel(mptr->posx,pnl,dim))									// set mptr->via to starting point that was in the panel
		for(d=0;d<dim;d++) mptr->via[d]=mptr->posx[d];
	else
		closestpanelpt(pnl,dim,mptr->pos,mptr->via);

	iter=0;
	while(iter==0 || !ptinpanel(mptr->pos,pnl,dim)) {
		iter++;
		er=lineexitpanel(mptr->via,mptr->pos,pnl,dim,pnledgept,&exitside);
		if(er || iter>20) {
			//printf("time = %g.  serno = %li.  ",sim->time,mptr->serno);	// this code is here for debugging purposes but shouldn't be needed.
			//simLog(sim,7,"BUG in movemol2closepanel %s\n",er?"mptr->via is same as mptr->pos":"iter >1");
			movept2panel(mptr->pos,pnl,dim,margin);
			fixpt2panel(mptr->pos,pnl,dim,face,epsilon);
			break; }
		nn=0;																						// pick a random neighbor for this edge point
		newpnl=NULL;
		newedge=0;
		for(p=0;p<pnl->nneigh;p++) {										// figure out values for newpnl and newedgept
			edgenum=closestpanelpt(pnl->neigh[p],dim,pnledgept,pt);
			dist2=0;
			for(d=0;d<dim;d++) dist2+=(pt[d]-pnledgept[d])*(pt[d]-pnledgept[d]);
			if(dist2<neighdist && coinrandD(1.0/++nn)) {
				newpnl=pnl->neigh[p];
				for(d=0;d<dim;d++) newedgept[d]=pt[d];
				newedge=edgenum; }}

		if(nn) {
			// if(mptr->s_index==0){															// rotate line to account for new normal vector
				for(d=0;d<dim;d++) mptr->pos[d]+=(newedgept[d]-pnledgept[d]);				// translate position to account for any edge mismatch
			//	complex_pos(sim->dim,mptr);				
			// }

			if(dim==2) {
				if(newedge) thetasign=1;
				else thetasign=coinrandD(0.5)?1:-1;
				paneledgenormal(newpnl,newedgept,dim,newedge,newedgenormal);
				len=sqrt((mptr->pos[0]-newedgept[0])*(mptr->pos[0]-newedgept[0])+(mptr->pos[1]-newedgept[1])*(mptr->pos[1]-newedgept[1]));
				mptr->pos[0]=newedgept[0]-thetasign*len*newedgenormal[0];
				mptr->pos[1]=newedgept[1]-thetasign*len*newedgenormal[1]; }
			else {
				panelnormal(newpnl,newedgept,face,dim,newnormal);
				if(newedge) {
					panelnormal(pnl,pnledgept,face,dim,oldnormal);
					paneledgenormal(newpnl,newedgept,dim,newedge,newedgenormal);
					dot=(mptr->pos[0]-newedgept[0])*newedgenormal[0]+(mptr->pos[1]-newedgept[1])*newedgenormal[1]+(mptr->pos[2]-newedgept[2])*newedgenormal[2];
					thetasign=(dot>0)?-1:1;
					Sph_RotateVectWithNormals3D(newedgept,mptr->pos,newpos,oldnormal,newnormal,thetasign);
					if(fabs(dot)<0.01) {												// check result if there's any doubt
						dot=(newpos[0]-newedgept[0])*newedgenormal[0]+(newpos[1]-newedgept[1])*newedgenormal[1]+(newpos[2]-newedgept[2])*newedgenormal[2];
						if(dot>0)																	// rotated wrong way, so do it again
							Sph_RotateVectWithNormals3D(newedgept,mptr->pos,newpos,oldnormal,newnormal,-thetasign); }
					// if(mptr->s_index==0){
						for(d=0;d<dim;d++) mptr->pos[d]=newpos[d];
					//	complex_pos(sim->dim,mptr); 
					// }
				}
				else
					Sph_RotateVectWithNormals3D(newedgept,mptr->pos,mptr->pos,NULL,newnormal,0); }
			pnl=mptr->pnl=newpnl; }
		else {																						// no new panel so bounce off of the edge of current panel
			paneledgenormal(pnl,pnledgept,dim,exitside,newnormal);
			dot=0;
			for(d=0;d<dim;d++) dot+=(mptr->pos[d]-pnledgept[d])*newnormal[d];
			// if(mptr->s_index==0){
				for(d=0;d<dim;d++) mptr->pos[d]-=2.0*newnormal[d]*dot;
			//	complex_pos(sim->dim,mptr);
			// }
			for(d=0;d<dim;d++) newedgept[d]=pnledgept[d]; }

		fixpt2panel(mptr->pos,pnl,dim,face,epsilon);
		for(d=0;d<dim;d++) mptr->via[d]=newedgept[d]; }

	return; }


/* surfacereflect */
void surfacereflect(moleculeptr mptr,panelptr pnl,double *crsspt,int dim,enum PanelFace face) {
	int d,axis;
	double *pos,*front,norm[3],norm2[3],dot;

	pos=mptr->pos;
	front=pnl->front;

	if(mptr->mstate==MSsoln) {
		if(pnl->ps==PSrect) {
			axis=(int)front[1];
			pos[axis]-=2.0*(pos[axis]-crsspt[axis]); 
			// printf("pname=%s, mptr->serno=%d, axis=%d, pos=%f\n", pnl->pname, mptr->serno, axis, pos[axis]); 
		}
		else if(pnl->ps==PStri || pnl->ps==PSdisk) {
			dot=0;
			for(d=0;d<dim;d++) dot+=(pos[d]-crsspt[d])*front[d];
			// if(mptr->s_index==0){
				for(d=0;d<dim;d++) pos[d]-=2.0*front[d]*dot; 
			//	complex_pos(dim,mptr);
			// }
		}
		else if(pnl->ps==PSsph || pnl->ps==PShemi) {
			Geo_SphereNormal(pnl->point[0],crsspt,1,dim,norm);
			dot=0;
			for(d=0;d<dim;d++) dot+=(pos[d]-crsspt[d])*norm[d];
			// if(mptr->s_index==0){
				for(d=0;d<dim;d++) pos[d]-=2.0*front[d]*dot; 
			//	complex_pos(dim,mptr);
			// }
		}		
		else if(pnl->ps==PScyl) {
			if(dim==2) {
				dot=0;
				for(d=0;d<dim;d++) dot+=(pos[d]-crsspt[d])*front[d];
				// if(mptr->s_index==0){
					for(d=0;d<dim;d++) pos[d]-=2.0*front[d]*dot; 
				//	complex_pos(dim,mptr);
				// }
			}
			else {
				Geo_LineNormal3D(pnl->point[0],pnl->point[1],crsspt,norm);
				dot=0;
				for(d=0;d<dim;d++) dot+=(pos[d]-crsspt[d])*norm[d];
				// if(mptr->s_index==0){
					for(d=0;d<dim;d++) pos[d]-=2.0*norm[d]*dot; 
				//	complex_pos(dim,mptr);
				// }
		}}
		else {																			// reflection for surface-bound molecules
			panelnormal(pnl,crsspt,face,dim,norm);	// norm is normal for collision panel
			panelnormal(mptr->pnl,crsspt,PFfront,dim,norm2);			// norm2 is normal for surface-bound panel
			dot=0;
			for(d=0;d<dim;d++) dot+=norm[d]*norm2[d];
			// if(mptr->s_index==0){
				for(d=0;d<dim;d++) norm[d]-=dot*norm2[d];		// make norm in plane of surface-bound panel
			//	complex_pos(dim,mptr);
			// }
			dot=0;
			for(d=0;d<dim;d++) dot+=norm[d]*norm[d];
			dot=sqrt(dot);
			if(dot==0) dot=1;														// should never happen
			for(d=0;d<dim;d++) norm[d]/=dot;						// normalize norm
			dot=0;
			for(d=0;d<dim;d++) dot+=(pos[d]-crsspt[d])*norm[d];
			// if(mptr->s_index==0){	
				for(d=0;d<dim;d++) pos[d]-=2.0*norm[d]*dot; 
			//	complex_pos(dim,mptr);
			// }
		}}

	return; }


/* surfacejump */
int surfacejump(moleculeptr mptr,panelptr pnl,double *crsspt,enum PanelFace face,int dim) {
	double **point,**point2,*front,*front2,dot,cent[3],delta[3];
	panelptr pnl2;
	enum PanelFace face2;
	int d,dir,ps;

	pnl2=pnl->jumpp[face];
	face2=pnl->jumpf[face];
	if(!pnl2) return 0;
	if(!(face2==PFfront || face2==PFback)) return 0;

	point=pnl->point;
	front=pnl->front;
	ps=pnl->ps;
	point2=pnl2->point;
	front2=pnl2->front;

	if(ps==PSrect) {
		Geo_RectCenter(point,cent,dim);
		Geo_RectCenter(point2,delta,dim);
		for(d=0;d<dim;d++) delta[d]-=cent[d];
		dir=(front[0]==front2[0])?1:-1; }
	else if(ps==PStri) {
		Geo_TriCenter(point,cent,dim);
		Geo_TriCenter(point2,delta,dim);
		for(d=0;d<dim;d++) delta[d]-=cent[d];
		dot=0;
		for(d=0;d<dim;d++) dot+=front[d]*front2[d];
		dir=dot>0?1:-1; }
	else if(ps==PSsph || ps==PShemi) {
		for(d=0;d<dim;d++)
			delta[d]=(crsspt[d]-point[0][d])*point2[1][0]/point[1][0]+point2[0][d]-crsspt[d];
		dir=(front[0]==front2[0])?1:-1; }
	else if(ps==PScyl) {
		for(d=0;d<dim;d++)
			delta[d]=(crsspt[d]-point[0][d])*point2[2][0]/point[2][0]+point2[0][d]-crsspt[d];
		dir=(front[2]==front2[2])?1:-1; }
	else if(ps==PSdisk) {
		for(d=0;d<dim;d++) delta[d]=point2[0][d]-point[0][d];
		dot=0;
		for(d=0;d<dim;d++) dot+=front[d]*front2[d];
		dir=dot>0?1:-1; }
	else {
		for(d=0;d<dim;d++) delta[d]=0;
		dir=1; }

	// printf("jump, line 4179, serno=%d, posoffset[0]=%f, posoffset[1]=%f, posoffset[2]=%f\n", mptr->serno, mptr->posoffset[0], mptr->posoffset[1], mptr->posoffset[2]);
	for(d=0;d<dim;d++) {
		crsspt[d]+=delta[d];
		mptr->pos[d]+=delta[d];
		mptr->posoffset[d]-=delta[d]; }
	// printf("jump, line 4184, serno=%d, posoffset[0]=%f, posoffset[1]=%f, posoffset[2]=%f\n", mptr->serno, mptr->posoffset[0], mptr->posoffset[1], mptr->posoffset[2]);
	// printf("jump, line 4185, serno=%d, pos[0]=%f, pos[1]=%f, pos[2]=%f\n", mptr->serno, mptr->pos[0], mptr->pos[1], mptr->pos[2]);
	fixpt2panel(crsspt,pnl2,dim,face2,0);
	dir*=(face!=face2)?1:-1;
	if(dir==-1) surfacereflect(mptr,pnl2,crsspt,dim,face);
	return 1; }

int rotating_jump(moleculeptr mptr, panelptr pnl, double *crsspt, enum PanelFace face, int dim){
	double **point1, **point2, *front1,*front2, dot; 
	double Sc1[3], Sc2[3];
	double pnlpoint1[3]={0,0,0};
	double pnlpoint2[3]={0,0,0};
	panelptr pnl2;
	enum PanelFace face2;
	int d,dir,ps,d1,d2;
	double r_angle;
	double pos[dim];

	pnl2=pnl->jumpp[face];
	face2=pnl->jumpf[face];
	if(!pnl2) return 0;
	if(!(face2==PFfront || face2==PFback)) return 0;
	
	point1=pnl->point;
	front1=pnl->front;
	point2=pnl2->point;
	front2=pnl2->front;
	ps=pnl->ps;	
	
	if(ps==PStri){
		for(d1=0;d1<dim;d1++){
			for(d2=0;d2<3;d2++){
				pnlpoint1[d1]+=point1[d2][d1];
				pnlpoint2[d1]+=point2[d2][d1];
			}
		}
		for(d=0;d<dim;d++){
			pnlpoint1[d]=(double)pnlpoint1[d]/3;
			pnlpoint2[d]=(double)pnlpoint2[d]/3;
		}

		Sph_Cart2Sc(pnlpoint1,Sc1);  // Sc1[0]: r, Sc1[1]: altitude (or inclination), Sc1[2]: azimuth
		Sph_Cart2Sc(pnlpoint2,Sc2);
		}
		r_angle=Sc1[1]+Sc2[1];
		for(d=0;d<dim;d++) pos[d]=mptr->pos[d];

		if(pnlpoint1[0]==0 && pnlpoint2[0]==0 && pnlpoint1[1]!=pnlpoint2[1])  {			// rotating along x axis: x'=x, y'=y*cos(theta)-z*sin(theta), z'=y*sin(theta)+z*cos(theta)
			printf("rotate about x\n");
			// if the to-panel < from-panel, then it's counter-clockwise, then r_angle>0
			r_angle=(pnlpoint2[1]<pnlpoint1[1])?r_angle:-r_angle;
			mptr->pos[1]= pos[1]*cos(r_angle) - pos[2]*sin(r_angle);				
			mptr->pos[2]= pos[1]*sin(r_angle) + pos[2]*cos(r_angle);
			// if(mptr->tot_sunit>1) rotation_update(mptr->rotate_mtrx, 0, r_angle);						// 0: rotate about x		
		}
		if(pnlpoint1[1]==0 && pnlpoint2[1]==0 && pnlpoint1[0]!=pnlpoint2[0]){			// rotating along y axis: x'=x*cos(theta)+z*sin(theta), y'=y, z'=-x*sin(theta)+z*cos(theta)
			printf("rotate about y\n");
			// don't know why for rotating about y, the rule is opposite
			r_angle=(pnlpoint2[0]<pnlpoint1[0])?-r_angle:r_angle;
			mptr->pos[0]= pos[0]*cos(r_angle) + pos[2]*sin(r_angle);
			mptr->pos[2]=-pos[0]*sin(r_angle) + pos[2]*cos(r_angle);
			// if(mptr->tot_sunit>1) rotation_update(mptr->rotate_mtrx, 1, r_angle);						// 1: rotate about y
		}
		
		// printf("rotating_jump, pnlpoint1[0]=%f, pnlpoint1[1]=%f, pnlpoint1[2]=%f, Sc1[0]=%f, Sc1[1]=%f, Sc1[2]=%f\n", pnlpoint1[0], pnlpoint1[1], pnlpoint1[2], Sc1[0], Sc1[1], Sc1[2]);
		// printf("rotating_jump, pnlpoint2[0]=%f, pnlpoint2[1]=%f, pnlpoint2[2]=%f, Sc2[0]=%f, Sc2[1]=%f, Sc2[2]=%f\n", pnlpoint2[0], pnlpoint2[1], pnlpoint2[2], Sc2[0], Sc2[1], Sc2[2]);
		// printf("rotating_jump, serno=%d,r_angle=%f, pos[0]=%f, pos[1]=%f, pos[2]=%f, new_pos[0]=%f, new_pos[1]=%f, new_pos[2]=%f\n", mptr->serno, r_angle, pos[0], pos[1],pos[2],mptr->pos[0],mptr->pos[1],mptr->pos[2]);

		for(d=0;d<dim;d++){
			crsspt[d]+=mptr->pos[d]-pos[d];
		}
		fixpt2panel(crsspt,pnl2,dim,face2,0);
		dir*=(face!=face2)?1:-1;
		if(dir==-1) surfacereflect(mptr,pnl2,crsspt,dim,face);
		return 1;	
}


/* dosurfinteract */
int dosurfinteract(simptr sim,moleculeptr mptr,int ll,int m,panelptr pnl,enum PanelFace face,double *crsspt, int* actptr) {
	int done,dim,d,i,i2,isneigh,p,jump,j;
	enum PanelFace newface;
	enum MolecState ms,ms2;
	enum SrfAction act;																// smoldyn.h line 181	
	double x,norm[DIMMAX],epsilon,neighdist,margin;
	moleculeptr mptr2;
	int cmptcheck;
	int site,s;

	dim=sim->dim;
	done=0;
	i=mptr->ident;
	ms=mptr->mstate;
	epsilon=sim->srfss->epsilon;
	neighdist=sim->srfss->neighdist;
	margin=sim->srfss->margin;

	isneigh=0;																		// check for neighbor
	if(face!=PFnone && mptr->pnl && pnl->nneigh>0)
		for(p=0;p<mptr->pnl->nneigh;p++)
			if(mptr->pnl->neigh[p]==pnl) isneigh=1;

	if(isneigh) {																	// find action, new identity, and new state
		if(mptr->pnl->srf->neighhop==0) act=SAtrans;
		else act=coinrandD(0.5)?SAadsorb:SAtrans;
		i2=i;
		ms2=ms; }
	else if((face==PFfront || face==PFback) && pnl->emitterabsorb[face] && pnl->emitterabsorb[face][i]>0) {
		// only unbound molecules such as camN0C0, free ca2+ are allowed to be absorbed
		if(randCCD()<pnl->emitterabsorb[face][i]) {
			for(s=0;s<sim->mols->spsites_num[mptr->ident];s++){
				if(mptr->sites[s]->bind){
					act=SAreflect;
					break;
			}}
			act=SAabsorb;
		}
		else act=SAreflect;
		i2=i;
		ms2=MSsoln; 
		//if(sim->events)
		//`	fprintf(sim->events,"dosurf time=%f act=%d mptr->serno=%d\n", sim->time,act,mptr->serno);
	}
	else{
		site=0;
		for(j=0;j<sim->mols->spsites_num[i];j++) {		
			if(mptr->sites[j]->bind) {
				site=j;
				break;
			}
		}
		act=surfaction(pnl->srf,face,i,ms,&i2,&ms2,site);
	}
	if(actptr) *actptr=act;
	
	if(act==SAno);														// no action

	else if(act==SAtrans) {												// transmit
		newface=(face==PFfront)?PFback:PFfront;
		fixpt2panel(crsspt,pnl,dim,newface,epsilon);
		if(i2!=i) molchangeident(sim,mptr,ll,m,i2,ms,mptr->pnl); }

	else if(act==SAreflect) {											// reflect, act=0
		surfacereflect(mptr,pnl,crsspt,dim,face);
		fixpt2panel(crsspt,pnl,dim,face,epsilon);
		if(ms!=MSsoln) movemol2closepanel(sim,mptr,dim,epsilon,neighdist,margin);
		if(panelside(mptr->pos,pnl,dim,NULL,1)!=face) fixpt2panel(mptr->pos,pnl,dim,face,0);
		if(i2!=i) molchangeident(sim,mptr,ll,m,i2,ms,mptr->pnl); }

	else if(act==SAabsorb) {											// absorb
		molkill(sim,mptr,ll,m);
		done=1; }

	// periodic jump: flip along one axis then rotate by a fixed angle
	else if(act==SArotate_jump){
		jump=rotating_jump(mptr,pnl,crsspt,face,dim);	
	}

	else if(act==SAjump) {												// jump (just solution-phase molecules), act=3	
		jump=surfacejump(mptr,pnl,crsspt,face,dim);
		// printf("surfact=jump, jump=%d\n", jump);
		if(!jump) {				// transmit
			newface=(face==PFfront)?PFback:PFfront;
			fixpt2panel(crsspt,pnl,dim,newface,epsilon); }
	}

	else if(act==SAport) {												// port
		newface=(face==PFfront)?PFback:PFfront;
    fixpt2panel(crsspt,pnl,dim,newface,epsilon);
		mptr->list=pnl->srf->port[face]->llport;
		if(ll>=0 && m<sim->mols->sortl[ll]) sim->mols->sortl[ll]=m;
		done=1; }

	else if(act==SAadsorb) {											// adsorb
		molchangeident(sim,mptr,ll,m,i2,ms2,pnl);
		for(d=0;d<dim;d++) mptr->pos[d]=crsspt[d];
		if(!ptinpanel(mptr->pos,mptr->pnl,dim))
			movept2panel(mptr->pos,mptr->pnl,dim,margin);
		if(ms2==MSfront) fixpt2panel(mptr->pos,mptr->pnl,dim,PFfront,epsilon);
		else if(ms2==MSback) fixpt2panel(mptr->pos,mptr->pnl,dim,PFback,epsilon);
		done=1; }

	else if(act==SArevdes || act==SAirrevdes) {		// desorb
		for(d=0;d<dim;d++) crsspt[d]=mptr->pos[d];
		fixpt2panel(crsspt,pnl,dim,ms2==MSsoln?PFfront:PFback,epsilon);
		//mptr2=getnextmol(sim->mols);
		mptr2=getnextmol_cplx(sim->mols,1,mptr->ident);
		if(!mptr2) return -1;
		mptr2->ident=i2;
		mptr2->mstate=MSsoln;
		mptr2->serno=mptr->serno;
		x=desorbdist(sim->mols->difstep[i2][MSsoln],act==SArevdes?SPArevAds:SPAirrDes);
		panelnormal(pnl,crsspt,ms2==MSsoln?PFfront:PFback,dim,norm);
		for(d=0;d<dim;d++) {
			mptr2->posx[d]=mptr2->via[d]=crsspt[d];
			//mptr2->pos[d]=crsspt[d]+x*norm[d];
			mptr2->posoffset[d]=mptr->posoffset[d]; }
		// if(mptr2->s_index==0){
			for(d=0;d<dim;d++) mptr2->pos[d]=crsspt[d]+x*norm[d];
			// complex_pos(sim->dim,mptr2);
		// }

		mptr2->box=pos2box(sim,mptr2->pos);
		mptr2->list=sim->mols->listlookup[i2][MSsoln];
		sim->eventcount[ETdesorb]++;
		molkill(sim,mptr,ll,m);
		checksurfaces1mol(sim,mptr2);
		done=1; }

	else if(act==SAflip) {												// on-surface flipping
		molchangeident(sim,mptr,ll,m,i2,ms2,pnl);
		if(ms2==MSfront) fixpt2panel(mptr->pos,mptr->pnl,dim,PFfront,epsilon);
		else if(ms2==MSback) fixpt2panel(mptr->pos,mptr->pnl,dim,PFback,epsilon);
		done=1; }

	return done; }


/* checksurfaces1mol */
int checksurfaces1mol(simptr sim,moleculeptr mptr) {
  int dim,d,done,p,lxp,it,flag;
  boxptr bptr1;
  double crossmin,crossmin2,crssptmin[3],crsspt[3],cross,*via,*pos;
  enum PanelFace face,facemin; 
  panelptr pnl,pnlmin;

  dim=sim->dim;
  via=mptr->via;
  // cplx, this line turns out to be very important
  for(d=0;d<dim;d++) via[d]=mptr->posx[d];
  pos=mptr->pos;
  done=0;
  it=0;
  while(!done) {
    if(++it>50) {
      for(d=0;d<dim;d++) pos[d]=mptr->posx[d];
      simLog(sim,7,"checksurfaces1mol(), SURFACE CALCULATION ERROR: molecule could not be placed after 50 iterations\n");
      break; }
    crossmin=crossmin2=2;
    facemin=PFfront;
    pnlmin=NULL;
    for(bptr1=pos2box(sim,via);bptr1;bptr1=line2nextbox(sim,via,pos,bptr1)) {
		for(p=0;p<bptr1->npanel;p++) {
        	pnl=bptr1->panel[p];
        	if(pnl!=mptr->pnl) {
		  		// printf("checksurfaces1mol, p=%d, pnl->pname=%s\n", p, pnl->pname);
		  		// printf("checksurfaces1mol, time=%f, tot_sunit=%d, via[0]=%f, via[1]=%f, via[2]=%f, pos[0]=%f, pos[1]=%f, pos[2]=%f\n", sim->time, mptr->tot_sunit, via[0], via[1], via[2], pos[0], pos[1], pos[2]);
          		lxp=lineXpanel(via,pos,pnl,dim,crsspt,&face,NULL,&cross,NULL,NULL);
          		if(lxp && cross<=crossmin2) {
            		if(cross<=crossmin) {			// what is cros: dist1/(dist1-dist2)
              			crossmin2=crossmin;
              			crossmin=cross;
              			pnlmin=pnl;
              			for(d=0;d<dim;d++) crssptmin[d]=crsspt[d];
              			facemin=face; }
            		else
              			crossmin2=cross; }}}}
	// printf("checksurfaces1mol line 4341, serno=%d, crossmin=%f, crossmin2=%f, posx[0]=%f, posx[1]=%f, posx[2]=%f, pos[0]=%f, pos[1]=%f, pos[2]=%f\n", mptr->serno, crossmin, crossmin2, mptr->posx[0],mptr->posx[1],mptr->posx[2],mptr->pos[0],mptr->pos[1], mptr->pos[2]);

    if(crossmin<2) {											// a panel was crossed, so deal with it
      flag=(crossmin2!=crossmin && crossmin2-crossmin<VERYCLOSE)?1:0;
      if(flag) {
        for(d=0;d<dim;d++) pos[d]=via[d];
		printf("checksurface1mol line 4349, serno=%d, via[0]=%f, via[1]=%f, via[2]=%f\n", mptr->serno, mptr->via[0], mptr->via[1], mptr->via[2]);
        done=1; }
      else {
#ifdef VCELL
        surfUpdateRate(sim, mptr, facemin, pnlmin);	// this is a very inefficient function.  If we defined suface rate as function then we have to evaluate the rate every time step when surface activity happens
#endif
        	done=dosurfinteract(sim,mptr,sim->mols->listtype[mptr->list]==MLTsystem?mptr->list:-1,-1,pnlmin,facemin,crssptmin, NULL);
        	for(d=0;d<dim;d++) via[d]=crssptmin[d];
        	sim->eventcount[ETsurf]++;
			// printf("checksurfaces1mol, line 4353, time=%f, pnlmin->pname=%s, crssptmin[0]=%f, [1]=%f, [2]=%f\n", sim->time, pnlmin->pname, crssptmin[0], crssptmin[1], crssptmin[2]);
	 	}}
    else{																	// nothing was crossed
      done=1;
	  // printf("checksurfaces1mol, line 4358, crossmin=%f\n", crossmin);	 
	}	
  }
  return 0; }

int checksurfaces(simptr sim,int ll, int reborn){
	moleculeptr *mlist,mptr,mptr_tmp,dif_bind, dif_molec;	
	int s,m,nmol,result,it,m_next,s0;	

	nmol=sim->mols->nl[ll];
	mlist=sim->mols->live[ll];
	if(!sim->srfss) return 0;
	if(!sim->mols) return 0;

	if(!reborn) m=0;
	else m=sim->mols->topl[ll];
	mptr=mlist[0];

	while(m<nmol) {
		mptr=mlist[m];
		dif_bind=dif_molec=NULL;
		result=0;
		m_next=m+mptr->tot_sunit;
		//m_next=m+;
		
		if(sim->mols->difc[mptr->ident][mptr->mstate]==0){
			m=m_next;
			continue;
		}

		if(mptr->pos!=mptr->pos_tmp && mptr->complex_id==-1){
			m=m_next;
			continue;
		}
		mptr_tmp=mptr;
		s0=0;
		if(mptr->complex_id!=-1){
			if(sim->mols->complexlist[mptr->complex_id]->dif_molec){
				dif_molec=sim->mols->complexlist[mptr->complex_id]->dif_molec;
				dif_bind=sim->mols->complexlist[mptr->complex_id]->dif_bind;
				mptr_tmp=dif_bind->to;
				m=m_next;
				continue;
				s0=1;	
			}
		}
		if(reborn){
			if(mptr->complex_id>-1){
			// if(sim->mols->complexlist[mptr->complex_id]->diffuse_updated<=1)
				if(mptr->prev_pos[0]==mptr->posx[0] && mptr->prev_pos[1]==mptr->posx[1] && mptr->prev_pos[2]==mptr->posx[2]) {
					// m+=mptr->tot_sunit;
					m++;
					continue;
		}}}

		for(s=s0,mptr_tmp=mptr;s<mptr->tot_sunit;s++,mptr_tmp=mptr_tmp->to){
			result=checksurfaces_cplx(sim,mptr_tmp,m+mptr_tmp->s_index,ll,reborn);
			if(mptr_tmp->ident!=0 && !boundarytest(sim,mptr_tmp->pos)) 
				printf("time=%f s=%d mptr_tmp->serno=%d result=%d pos0=%f pos1=%f pos2=%f\n", sim->time,s,mptr_tmp->serno,result,mptr_tmp->pos[0],mptr_tmp->pos[1],mptr_tmp->pos[2]);
		}
		s=0;
		mptr_tmp=mptr;
		it=0;
		while(s<mptr->tot_sunit && mptr->tot_sunit>1){
			if(dif_molec){
				s++;
				mptr_tmp=mlist[m+s];
				continue;
			}
			else{
				if(!boundarytest(sim,mptr_tmp->pos)){
					result=checksurfaces_cplx(sim,mptr_tmp,m+s,ll,reborn);
					while(!boundarytest(sim,mptr_tmp->pos)){
						printf("chksurf line4486, not in cmpt, time=%f ident=%d complex_id=%d serno=%d s_index=%d pos0=%f pos1=%f pos2=%f prev_pos0=%f prev_pos1=%f prev_pos2=%f\n", sim->time, mptr_tmp->ident, mptr_tmp->complex_id, mptr_tmp->serno, mptr_tmp->s_index, mptr_tmp->pos[0], mptr_tmp->pos[1], mptr_tmp->pos[2], mptr_tmp->prev_pos[0], mptr_tmp->prev_pos[1], mptr_tmp->prev_pos[2]);
						result=checksurfaces_cplx(sim,mptr_tmp,m+s,ll,reborn);
						it++;
						if(it>50){
							fprintf(sim->events, "surf line4496, it over 50, time=%f, mptr_tmp->serno=%d\n", sim->time, mptr_tmp->serno); 
							result=-2;
							s=mptr->tot_sunit;
							break;
							// return -2;
						}
					}
					s=0;
					mptr_tmp=mlist[m+s];
				}
				else{
					s++;
					mptr_tmp=mlist[m+s];
				}
			}
		}	
  		if(result!=0){
			printf("surf.c, time=%f, m=%d, result=%d\n", sim->time, m, result);
			return result;
		}
		m=m_next;
	}
	

  return result;
}
/* checksurfaces. */
int checksurfaces_cplx(simptr sim, moleculeptr mptr, int m, int ll,int reborn) {
	int act,dim,d,done,p,lxp,it,flag;
	boxptr bptr1;
	double crossmin,crossmin2,crssptmin[3],crsspt[3],cross,*via,*pos;
	double pos_offset[3];
	enum PanelFace face,facemin;
	panelptr pnl,pnlmin;
	// difadj_ptr adj_list;

	dim=sim->dim;
	act=-1;
	cross=0;									// cplx
	crsspt[0]=0; crsspt[1]=0; crsspt[2]=0;		// cplx

	via=mptr->via;
	for(d=0;d<dim;d++) {
		via[d]=mptr->posx[d];
		mptr->prev_pos[d]=mptr->pos[d];
	}
	pos=mptr->pos;
	done=0;
	it=0;
	while(!done) {
		if(++it>50) {
			for(d=0;d<dim;d++) pos[d]=mptr->posx[d];
			simLog(sim,7,"checksurfaces(), SURFACE CALCULATION ERROR: molecule could not be placed after 50 iterations\n");
			break;
		}
		crossmin=crossmin2=2;
		facemin=PFfront;
		pnlmin=NULL;
		for(bptr1=pos2box(sim,via);bptr1;bptr1=line2nextbox(sim,via,pos,bptr1)) {
			for(p=0;p<bptr1->npanel;p++) {
				pnl=bptr1->panel[p];
				if(pnl!=mptr->pnl) {
					lxp=lineXpanel(via,pos,pnl,dim,crsspt,&face,NULL,&cross,NULL,NULL);
					if(lxp && cross<=crossmin2) {
						if(cross<=crossmin) {
							crossmin2=crossmin;
							crossmin=cross;
							pnlmin=pnl;
							for(d=0;d<dim;d++) crssptmin[d]=crsspt[d];
							facemin=face; }
						else
							crossmin2=cross; }}}}
			
		if(crossmin<2) {														// a panel was crossed, so deal with it
			flag=(crossmin2!=crossmin && crossmin2-crossmin<VERYCLOSE)?1:0;
			if(flag) {
				for(d=0;d<dim;d++) pos[d]=via[d];
				done=1; }
			else {
				done=dosurfinteract(sim,mptr,ll,m,pnlmin,facemin,crssptmin,&act);
				for(d=0;d<dim;d++) via[d]=crssptmin[d];
				sim->eventcount[ETsurf]++; }
		}else																	// nothing was crossed
		done=1; 
	}
	if(!boundarytest(sim,mptr->pos) && mptr->ident>0)
		printf("time=%f ident=%d serno=%d act=%d prevpos[0]=%f prevpos[1]=%f prevpos[2]=%f posx[0]=%f posx[1]=%f posx[2]=%f\n",sim->time,mptr->ident,mptr->serno,act,mptr->prev_pos[0],mptr->prev_pos[1],mptr->prev_pos[2],mptr->posx[0],mptr->posx[1],mptr->posx[2]);
	if(mptr->tot_sunit>1  && act!=-1){	
		for(d=0;d<dim;d++) pos_offset[d]=mptr->pos[d]-mptr->prev_pos[d];
		if(pos_offset[0]==0 && pos_offset[1]==0 && pos_offset[2]==0) return 0;

		if(complex_pos(sim,mptr,"pos_surfline4570",pos_offset,1)==-1){
			fprintf(sim->events,"chksurf, time=%f, act=%d, reborn=%d, mptr->serno=%d, pos_offset[0]=%f, pos_offset[1]=%f, pos_offset[2]=%f\n", sim->time, act, reborn, mptr->serno, pos_offset[0], pos_offset[1], pos_offset[2]); 
			return -3;
		}
	}

	return 0; }


/* checksurfacebound. */
int checksurfacebound(simptr sim,int ll) {
	int nmol,m,done;
	moleculeptr mptr,*mlist;

	if(!sim->srfss) return 0;
	if(!sim->mols) return 0;

	nmol=sim->mols->nl[ll];
	mlist=sim->mols->live[ll];
	for(m=0;m<nmol;m++) {
		mptr=mlist[m];
		// printf("before dosurfinteract, time=%f, serno=%d, mptr->pos[0]=%f, pos[1]=%f, pos[2]=%f\n", sim->time, mptr->serno, mptr->pos[0], mptr->pos[1], mptr->pos[2]);
		if(mptr->mstate!=MSsoln) {
			done=dosurfinteract(sim,mptr,ll,m,mptr->pnl,PFnone,mptr->posx, NULL);
			if(done==-1) simLog(sim,10,"Unable to allocate memory in dosurfinteract\n"); }
		// printf("after dosurfinteract, time=%f, serno=%d, mptr->pos[0]=%f, pos[1]=%f, pos[2]=%f\n", sim->time, mptr->serno, mptr->pos[0], mptr->pos[1], mptr->pos[2]);
	}
	return 0; }

