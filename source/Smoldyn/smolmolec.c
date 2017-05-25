/* Steven Andrews, started 10/22/2001.
 This is a library of functions for the Smoldyn program.
 See documentation called Smoldyn_doc1.pdf and Smoldyn_doc2.pdf, and the Smoldyn
 website, which is at www.smoldyn.org.
 Copyright 2003-2013 by Steven Andrews.  This work is distributed under the terms
 of the Gnu Lesser General Public License (LGPL). */

#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include <stdarg.h>
#include "Geometry.h"
#include "math2.h"
#include "random2.h"
#include "Rn.h"
#include "RnSort.h"
#include "smoldyn.h"
#include "smoldynfuncs.h"
#include "smoldynconfigure.h"
#include <iostream>
#include <algorithm>

/******************************************************************************/
/*********************************** Molecules ********************************/
/******************************************************************************/


/******************************************************************************/
/****************************** Local declarations ****************************/
/******************************************************************************/

// enumerated type functions
enum MolListType molstring2mlt(char *string);
char *molmlt2string(enum MolListType mlt,char *string);

// low level utilities
char *molpos2string(simptr sim,moleculeptr mptr,char *string);

// memory management
moleculeptr molalloc(simptr sim, int dim);
void molfree(simptr sim, moleculeptr mptr);
void molfreesurfdrift(double *****surfdrift,int maxspec,int maxsrf);
molssptr molssalloc(molssptr mols,int maxspecies);
int mollistalloc(molssptr mols,int maxlist,enum MolListType mlt);
int molexpandlist(molssptr mols,int dim,int ll,int nspaces,int nmolecs);
int molpatternindexalloc(int **indexptr,int n);
int molpatternalloc(simptr sim,int maxpattern);
complexptr complexalloc(simptr sim, moleculeptr mptr);
void complexfree(complexptr cplxptr);

// data structure output

// structure setup
int molsetmaxspecies(simptr sim,int max);
int molsupdateparams(molssptr mols,double dt);
int molsupdatelists(simptr sim);

// adding and removing molecules
moleculeptr newestmol(molssptr mols);
int molgetexport(simptr sim,int ident,enum MolecState ms);
int molputimport(simptr sim,int nmol,int ident,enum MolecState ms,panelptr pnl,enum PanelFace face);
int moldummyporter(simptr sim);

// core simulation functions
double power(double a, int b);

/******************************************************************************/
/********************************* enumerated types ***************************/
/******************************************************************************/


/* molstring2ms */
enum MolecState molstring2ms(char *string) {
	enum MolecState ans;

	if(!strcmp(string,"solution")) ans=MSsoln;
	else if(!strcmp(string,"fsoln")) ans=MSsoln;
	else if(!strcmp(string,"soln")) ans=MSsoln;
	else if(!strcmp(string,"aq")) ans=MSsoln;
	else if(!strcmp(string,"front")) ans=MSfront;
	else if(!strcmp(string,"back")) ans=MSback;
	else if(!strcmp(string,"up")) ans=MSup;
	else if(!strcmp(string,"down")) ans=MSdown;
	else if(!strcmp(string,"bsoln")) ans=MSbsoln;
	else if(!strcmp(string,"all")) ans=MSall;
	else ans=MSnone;
	return ans;
}


/* molms2string */
char *molms2string(enum MolecState ms,char *string) {
	if(ms==MSsoln) strcpy(string,"solution");
	else if(ms==MSfront) strcpy(string,"front");
	else if(ms==MSback) strcpy(string,"back");
	else if(ms==MSup) strcpy(string,"up");
	else if(ms==MSdown) strcpy(string,"down");
	else if(ms==MSbsoln) strcpy(string,"bsoln");
	else if(ms==MSall) strcpy(string,"all");
	else if(ms==MSsome) strcpy(string,"some");
	else strcpy(string,"none");
	return string; }


/* molstring2mlt */
enum MolListType molstring2mlt(char *string) {
	enum MolListType ans;

	if(!strcmp(string,"system")) ans=MLTsystem;
	else if(!strcmp(string,"port")) ans=MLTport;
	else ans=MLTnone;
	return ans; }


/* molmlt2string */
char *molmlt2string(enum MolListType mlt,char *string) {
	if(mlt==MLTsystem) strcpy(string,"system");
	else if(mlt==MLTport) strcpy(string,"port");
	else strcpy(string,"none");
	return string; }


/******************************************************************************/
/****************************** low level utilities ***************************/
/******************************************************************************/


/* molismobile */
int molismobile(simptr sim,int species,enum MolecState ms) {
	molssptr mols;
	int d,dim,s;
	enum PanelShape ps;

	mols=sim->mols;
	dim=sim->dim;
	if(mols->difc[species][ms]>0) return 1;
	if(mols->difm && mols->difm[species] && mols->difm[species][ms])
		for(d=0;d<dim*dim;d++)
			if(mols->difm[species][ms][d]!=0) return 1;
	if(mols->drift && mols->drift[species] && mols->drift[species][ms])
		for(d=0;d<dim;d++)
			if(mols->drift[species][ms][d]!=0) return 1;
	if(mols->surfdrift && mols->surfdrift[species] && mols->surfdrift[species][ms])
		for(s=0;s<sim->srfss->nsrf;s++)
			if(mols->surfdrift[species][ms][s])
				for(ps=PanelShape(0);ps<PSMAX;ps=PanelShape(ps+1))
					if(mols->surfdrift[species][ms][s][ps])
						for(d=0;d<dim-1;d++)
							if(mols->surfdrift[species][ms][s][ps][d]!=0) return 1;

	return 0; }


/* molwildcardname */
int molwildcardname(molssptr mols,char *name,int channel,int itest) {
	static int i[]={-1,-1};									// i is most recent returned index
	static char nm[][STRCHAR]={"\0","\0"};			// nm is name with wildcards
	static unsigned int *flags=NULL;				// flags lists match results for all species
	static int nspecies=0;
	int i2,ch,match;
	unsigned int mask;

	if(!mols) {										// free memory and reset static variables (mode 1)
		for(ch=0;ch<2;ch++) {
			i[ch]=-1;
			nm[ch][0]='\0';
			free(flags);
			flags=NULL;
			nspecies=0; }
		strEnhWildcardMatch(NULL,NULL);
		return 0; }

	if(name) {							// initialize using a new name (mode 2)
		strncpy(nm[channel],name,STRCHAR);
		i[channel]=1;
		if(nspecies==0);
		else if(nspecies!=mols->nspecies) {
			nspecies=0;
			free(flags); }
		else {
			mask=~(1<<channel);
			for(i2=0;i2<nspecies;i2++)
				flags[i2]&=mask; }
		match=strEnhWildcardMatch(name,"test");
		if(match<0) return match-1; }

	if(itest>0) {						// test specific index (mode 4)
		if(!nspecies) {
			nspecies=mols->nspecies;
			flags=(unsigned int*) calloc(nspecies,sizeof(unsigned int));
			if(!flags) {molwildcardname(NULL,NULL,0,0);return -2;}
			for(i2=0;i2<nspecies;i2++) flags[i2]=0; }
		mask=1<<channel;
		if(!(flags[0]&mask)) {
			flags[0]|=mask;
			for(i2=1;i2<nspecies;i2++)
				if(strEnhWildcardMatch(nm[channel],mols->spname[i2]))
					flags[i2]|=mask; }
		return flags[itest]&mask; }

													// return next match (mode 3)
	if(i[channel]<0) return -1;						// all matches already done
	for(i2=i[channel];i2<mols->nspecies && !strEnhWildcardMatch(nm[channel],mols->spname[i2]);i2++);
	if(i2==mols->nspecies) i[channel]=i2=-1;
	else {											// where to start next search
		i[channel]=i2;
		if(!name) i[channel]++; }
	return i2; }


/* readmolname */
int readmolname(simptr sim,char *str,enum MolecState *msptr,int channel) {
	char nm[STRCHAR],*pareno,*parenc;
	int itct,i;
	enum MolecState ms;

	if(!str) return -1;
	itct=sscanf(str,"%s",nm);

	if(itct!=1) return -1;	// cannot read name
	pareno=strchr(nm,'(');
	if(pareno) {
		*pareno='\0';
		pareno++;
		parenc=strrchr(pareno,')');
		if(parenc && *(parenc+1)=='\0') *parenc='\0';
		else return -2;				// improper close parenthesis
		ms=molstring2ms(pareno);
		if(ms==MSnone) return -3; }		// cannot read state
	else ms=MSsoln;

	if(!strcmp(nm,"all")) i=-5;		// all
	else if(strchr(nm,'*') || strchr(nm,'?')) {	// wildcard character
		if(channel<0) i=-4;
		else i=molwildcardname(sim->mols,nm,channel,0);
		if(i>0) i=-6;						// at least one match
		else i=-4; }						// no match
	else {
		i=stringfind(sim->mols->spname,sim->mols->nspecies,nm);
		if(i<0) i=-4; }		// unknown molecule name
	
	if(msptr) *msptr=ms;
	return i; }

int readmolname_cplx(simptr sim,char *str, int **identptr, int **sitesptr){
	if(!str) return -3;
	char *cptr, molec[STRCHAR], *s,str0[STRCHAR];
	int i,ident_indx,site_indx,count;
	
	count=1;
	sscanf(str,"%s",str0);
	if(!strcmp(str0,"all")){
		identptr[0]=NULL;
		sitesptr[0]=NULL;
		return -5;
	} 

	while(cptr=strchr(str0,'~')){
		count++; 
		*cptr=' ';
		cptr++;
	}
	identptr[0]=(int*) calloc(count,sizeof(int));
	
	if(count==1){
		sitesptr[0]=NULL;
		// ms1=strextract(str,"()");,
		// if(ms1) msptr[0][0]=molstring2ms(ms1);
		// else msptr[0][0]=MSsoln;
			
		if(strchr(str0,'@')!=0){
			//if(strchr(str0,'-')!=0){
			//	sitesptr[0]=(int*) calloc(count+1,sizeof(int));
			//}
			sitesptr[0]=(int*) calloc(count+1,sizeof(int));
			sscanf(str0,"%s",molec);
			cptr=strsplit(molec,"@");
			ident_indx=stringfind(sim->mols->spname,sim->mols->nspecies,molec);
			if(ident_indx<0) return -1;
			else identptr[0][0]=ident_indx;
			if(strchr(cptr,'-')){
				s=strsplit(cptr,"-");
				site_indx=stringfind(sim->mols->spsites_name[ident_indx],sim->mols->spsites_num[ident_indx],cptr);
				if(site_indx<0) return -2;
				else sitesptr[0][0]=site_indx;
				sitesptr[0][1]=stringfind(sim->mols->spname,sim->mols->nspecies,s);
			}
			else{
				site_indx=stringfind(sim->mols->spsites_name[ident_indx],sim->mols->spsites_num[ident_indx],cptr);
				if(site_indx<0) return -2;
				else sitesptr[0][0]=site_indx;
				sitesptr[0][1]=-1;
			}
		}
		else{
			ident_indx=stringfind(sim->mols->spname,sim->mols->nspecies,str0);
			if(ident_indx<0) return -1;
			else identptr[0][0]=ident_indx;
		}
		
	}
	if(count==2){
		s=str0;
		sitesptr[0]=(int*) calloc(count,sizeof(int));
		for(i=0;i<count;i++){
			sscanf(s,"%s",molec);
			s=s+strlen(molec)+1; 
			
			cptr=strsplit(molec,"@");
			ident_indx=stringfind(sim->mols->spname,sim->mols->nspecies,molec);			
			if(ident_indx<0) return -1;
			else identptr[0][i]=ident_indx;
			printf("cptr=%s, molec=%s\n",cptr,molec);
	
			site_indx=stringfind(sim->mols->spsites_name[ident_indx],sim->mols->spsites_num[ident_indx],cptr);
			if(site_indx<0) return -2;
			else sitesptr[0][i]=site_indx;
			printf("identptr[0][i]=%d, sitesptr[0][i]=%d\n", identptr[0][i], sitesptr[0][i]);
		}
	}
	return count;
}

/* molstring2pattern */
int molstring2pattern(simptr sim,char *str,enum MolecState *msptr,char *pat,int mode) {
	char nm[STRCHAR],*pareno,*parenc,*sitestr;
	int itct;
	enum MolecState ms;

	if(!str || !pat) return -1;
	itct=sscanf(str,"%s",nm);
	if(itct!=1) return -1;							// cannot read name

	if(nm[strlen(nm)-1]==')') {
		parenc=nm+strlen(nm)-1;
		*parenc='\0';
		pareno=strrchr(nm,'(');
		if(!pareno) return -2;						// improper open parenthesis
		*pareno='\0';
		pareno++;
		ms=molstring2ms(pareno);
		if(ms==MSnone) return -3; }					// cannot read state
	else ms=MSsoln;

	if(mode==0) pat[0]='\0';
	else if(mode==1 || (mode==2 && strchr(pat,'\n'))) strcat(pat," ");
	else if(mode==2) strcat(pat,"\n");
	strcat(pat,nm);

	if(msptr) *msptr=ms;
	return 0; }


/* molpatternindex */
int molpatternindex(simptr sim,const char *pattern,int isrule,int **indexptr) {
	int npattern,pat,er,i,j,istart,jstart,matchwords,subwords,degenwords,totalwords,nspecies,haswildcard;
	int iword,iword2,somethingnew;
	int **patindex,*index;
	int maxmatch=4,ispecies[4];
	char **patlist,*patstring,matchstr[STRCHAR],*substr,*newline,teststring[STRCHAR],deststring[STRCHAR];

	patlist=sim->mols->patlist;
	patindex=sim->mols->patindex;
	npattern=sim->mols->npattern;

	if(npattern)
		pat=locateVstr(patlist,pattern,npattern,1);			// look for pattern in patlist list
	else pat=-1;

	if(pat<0 || strcmp(patlist[pat],pattern) || isrule!=patindex[pat][PDrule]) {	// this is a new pattern
		pat++;
		if(npattern==sim->mols->maxpattern) {						// expand pattern list if its already full
			er=molpatternalloc(sim,2*npattern+2);
			if(er) return -1;
			patlist=sim->mols->patlist;
			patindex=sim->mols->patindex; }
		patstring=patlist[npattern];										// swap last pattern in list with correct spot
		index=patindex[npattern];
		for(j=npattern-1;j>=pat;j--) {
			patlist[j+1]=patlist[j];
			patindex[j+1]=patindex[j]; }
		sim->mols->npattern++;
		npattern=sim->mols->npattern;

		patlist[pat]=patstring;
		patindex[pat]=index;
		strcpy(patlist[pat],pattern);
		patindex[pat][PDnresults]=0;
		patindex[pat][PDnspecies]=1;
		patindex[pat][PDrule]=isrule; }

	if(patindex[pat][PDnspecies]<sim->mols->nspecies) {			// update pattern list for current species list
		strcpy(matchstr,pattern);												// assign to matchstr and substr and determine newline
		newline=strchr(matchstr,'\n');
		if(newline) {
			*newline='\0';
			substr=newline+1; }
		else substr=NULL;
		istart=patindex[pat][PDnspecies];								// istart is first species to start checking
		jstart=patindex[pat][PDnresults];								// jstart is first result that will be added
		if(istart==1) {																	// setup header
			matchwords=wordcount(matchstr);
			subwords=newline?wordcount(substr):0;
			degenwords=newline?1:0;
			patindex[pat][PDmatch]=matchwords;
			patindex[pat][PDsubst]=subwords;
			patindex[pat][PDdegen]=degenwords; }
		else {																					// read header
			matchwords=patindex[pat][PDmatch];
			subwords=patindex[pat][PDsubst];
			degenwords=patindex[pat][PDdegen]; }
		haswildcard=strpbrk(pattern,"*?&|{}[]")?1:0;
		nspecies=sim->mols->nspecies;
		totalwords=matchwords+subwords+degenwords;

		if(!strcmp(pattern,"all")) {														// pattern == "all"
			if(patindex[pat][PDalloc]<PDMAX+totalwords*nspecies) {
				er=molpatternindexalloc(&patindex[pat],PDMAX+totalwords*nspecies);
				if(er) return -1; }
			for(i=istart,j=jstart;i<nspecies;i++,j++)
				patindex[pat][PDMAX+j*totalwords]=i;
			patindex[pat][PDnresults]=j;
			patindex[pat][PDnspecies]=nspecies;
			sortVii(patindex[pat]+PDMAX,NULL,patindex[pat][PDnresults]); }	// list is sorted

		else if(!haswildcard && jstart==1) {										// no wildcards, so 1 entry, and it's already done
			patindex[pat][PDnspecies]=nspecies; }

		else if(!haswildcard) {																	// pattern has no wildcards, so just one entry
			if(patindex[pat][PDalloc]<PDMAX+totalwords) {
				er=molpatternindexalloc(&patindex[pat],PDMAX+totalwords);
				if(er) return -1; }
			for(iword=0;iword<matchwords;iword++) {								// matchwords
				i=stringfind(sim->mols->spname,nspecies,strnword(matchstr,iword+1));
				if(i<0) return -2;																	// unknown molecule name
				patindex[pat][PDMAX+iword]=i; }
			for(iword=0;iword<subwords;iword++) {									// subwords
				i=stringfind(sim->mols->spname,nspecies,strnword(substr,iword+1));
				if(i<0 && !isrule) return -3;												// ?? need to expand to enable rules
				patindex[pat][PDMAX+matchwords+iword]=i; }
			if(patindex[pat][PDdegen])
				 patindex[pat][PDMAX+matchwords+subwords]=1;				// degeneracy is 1
			patindex[pat][PDnresults]=1;
			patindex[pat][PDnspecies]=nspecies; }

		else if(!newline && matchwords==1) {										// just single match word
			j=jstart;
			for(i=istart;i<nspecies;i++) {
				if(strEnhWildcardMatch(matchstr,sim->mols->spname[i])) {
					if(patindex[pat][PDalloc]<PDMAX+(j+1)) {
						er=molpatternindexalloc(&patindex[pat],PDMAX+2*(j+1));
						if(er) return -1; }
					patindex[pat][PDMAX+j]=i;
					j++;
					patindex[pat][PDnresults]=j; }}
			patindex[pat][PDnspecies]=nspecies;
			sortVii(patindex[pat]+PDMAX,NULL,patindex[pat][PDnresults]); }	// this list is sorted

		else if(!newline) {																				// several match words, no substitute words
			j=jstart;
			if(matchwords>maxmatch) return -4;
			for(iword=0;iword<matchwords;iword++)										// loop to get all species permutations
				for(ispecies[iword]=1;ispecies[iword]<sim->mols->nspecies;ispecies[iword]++) {
					somethingnew=0;
					for(iword2=0;iword2<matchwords && !somethingnew;iword2++)
						if(ispecies[iword2]>=istart) somethingnew=1;			// already done if not somethingnew
					if(somethingnew) {
						strcpy(teststring,sim->mols->spname[ispecies[0]]);
						for(iword2=1;iword2<matchwords;iword2++) {
							strcat(teststring," ");
							strcat(teststring,sim->mols->spname[ispecies[iword2]]); }
						if(strEnhWildcardMatch(matchstr,teststring)) {
							if(patindex[pat][PDalloc]<PDMAX+totalwords*(j+1)) {
								er=molpatternindexalloc(&patindex[pat],PDMAX+2*totalwords*(j+1));
								if(er) return -1; }
							for(iword2=0;iword2<matchwords;iword2++)
								patindex[pat][PDMAX+totalwords*j+iword2]=ispecies[iword2];
							j++;
							patindex[pat][PDnresults]=j; }}}
			patindex[pat][PDnspecies]=nspecies; }

		else {																										// both match and substitute words
			j=jstart;
			if(matchwords>maxmatch) return -4;
			for(iword=0;iword<matchwords;iword++)										// loop to get all species permutations on match side
				for(ispecies[iword]=1;ispecies[iword]<sim->mols->nspecies;ispecies[iword]++) {
					somethingnew=0;
					for(iword2=0;iword2<matchwords && !somethingnew;iword2++)
						if(ispecies[iword2]>=istart) somethingnew=1;			// already done if not somethingnew
					if(somethingnew) {
						strcpy(teststring,sim->mols->spname[ispecies[0]]);
						for(iword2=1;iword2<matchwords;iword2++) {
							strcat(teststring," ");
							strcat(teststring,sim->mols->spname[ispecies[iword2]]); }
						if(strEnhWildcardMatchAndSub(matchstr,teststring,substr,deststring)) {
							if(patindex[pat][PDalloc]<PDMAX+totalwords*(j+1)) {
								er=molpatternindexalloc(&patindex[pat],PDMAX+2*totalwords*(j+1));
								if(er) return -1; }
							for(iword2=0;iword2<matchwords;iword2++)
								patindex[pat][PDMAX+totalwords*j+iword2]=ispecies[iword2];
							for(iword2=0;iword2<subwords;iword2++) {
								i=stringfind(sim->mols->spname,nspecies,strnword(deststring,iword+1));
								if(i<0 && !isrule) return -3;											// ?? need to expand to enable rules, also I'm not sure this should be an error if it's not a rule.
								patindex[pat][PDMAX+totalwords*j+matchwords+iword2]=i; }
							if(patindex[pat][PDdegen])
								patindex[pat][PDMAX+matchwords+subwords]=1;				// ?? degeneracy is given as 1 but it might not really
							j++;
							patindex[pat][PDnresults]=j; }}}
			patindex[pat][PDnspecies]=nspecies; }}

	if(indexptr) *indexptr=patindex[pat];
	return 0; }


/* molstring2index1 */
int molstring2index1(simptr sim,char *str,enum MolecState *msptr,int **indexptr) {
	int er,isall,*index, *site;
	char pattern[STRCHAR];

	er=molstring2pattern(sim,str,msptr,pattern,0);
	if(er) return er;
	isall=strcmp(pattern,"all")?0:1;
	er=molpatternindex(sim,pattern,0,&index);
	if(indexptr) *indexptr=index;
	if(isall && !er) return -5;
	if(er==-1) return -7;
	if(er==-2) return -4;
	if(index[PDnresults]==1 && index[PDmatch]==1)
		return index[PDMAX];
	return 0; }


/* molpos2string. */
char *molpos2string(simptr sim,moleculeptr mptr,char *string) {
	int d,dim,done,p,tryagain,count;
	char *line2;
	double newpos[DIMMAX],crosspt[DIMMAX],dist;
	boxptr bptr;
	panelptr pnl;

	dim=sim->dim;
	done=0;
	dist=0;
	count=0;

	line2=string;												// write position to string
	for(d=0;d<dim;d++) {
		sprintf(line2," %g",mptr->pos[d]);
		line2+=strlen(line2); }

	if(!sim->srfss) done=1;
	while(!done) {
		line2=string;											// read in written position
		for(d=0;d<dim;d++) {
			sscanf(line2,"%lg",&newpos[d]);
			line2=strnword(line2,2); }

		tryagain=0;
		bptr=pos2box(sim,newpos);
		if(bptr!=pos2box(sim,mptr->pos)) tryagain=1;		// check for same box
		for(p=0;p<bptr->npanel && tryagain==0;p++) {		// check for no panels crossed
			pnl=bptr->panel[p];
			if(mptr->pnl!=pnl && lineXpanel(mptr->pos,newpos,pnl,dim,crosspt,NULL,NULL,NULL,NULL,NULL)) tryagain=1; }
		if(!tryagain) done=1;

		if(!done) {
			if(++count>50) {
				simLog(sim,8,"WARNING: unable to write %s molecule position (%s) on the correct side of all surfaces\n",sim->mols->spname[mptr->ident],string);
				return string; }

			if(dist==0) {
				for(d=0;d<dim;d++) dist+=(newpos[d]-mptr->pos[d])*(newpos[d]-mptr->pos[d]);
				dist=50*sqrt(dist); }

			line2=string;												// write position to string
			for(d=0;d<dim;d++) {
				sprintf(line2," %g",mptr->pos[d]+unirandCCD(-dist,dist));
				line2+=strlen(line2); }}}
		
		return string; }


/* molchangeident */
void molchangeident(simptr sim,moleculeptr mptr,int ll,int m,int i,enum MolecState ms,panelptr pnl) {
	int dim,ll2;
	double epsilon;

	if(i==0) {
		molkill(sim,mptr,ll,m);
		return; }

	dim=sim->dim;
	epsilon=sim->srfss?sim->srfss->epsilon:0;

	mptr->ident=i;
	mptr->mstate=ms;
	if(ms==MSsoln || ms==MSbsoln) mptr->pnl=NULL;
	else mptr->pnl=pnl;

	if(ms==MSsoln && !pnl);												// soln -> soln
	else if(ms==MSsoln) {													// surf -> front soln
		fixpt2panel(mptr->posx,pnl,dim,PFfront,epsilon); }
	else if(ms==MSbsoln) {												// surf -> back soln
		mptr->mstate=MSsoln;
		fixpt2panel(mptr->posx,pnl,dim,PFback,epsilon); }
	else if(ms==MSfront)													// any -> front surf
		fixpt2panel(mptr->pos,pnl,dim,PFfront,epsilon);
	else if(ms==MSback)														// any -> back surf
		fixpt2panel(mptr->pos,pnl,dim,PFback,epsilon);
	else																					// any -> up or down
		fixpt2panel(mptr->pos,pnl,dim,PFnone,epsilon);

	ll2=sim->mols->listlookup[i][ms];
	if(ll>=0 && ll2!=ll) {
		mptr->list=ll2;
		if(m<0) sim->mols->sortl[ll]=0;
		else if(m<sim->mols->sortl[ll]) sim->mols->sortl[ll]=m; }
	return; }


/* molssetgausstable */
int molssetgausstable(simptr sim,int size) {
	int er;
	molssptr mols;
	double *newtable;
	if(er) return er;
	mols=sim->mols;

	if(mols->ngausstbl>0 && (mols->ngausstbl==size || size==-1)) return 0;
	if(size<1) size=4096;
	else if(!is2ton(size)) return 3;

	newtable=(double*) calloc(size,sizeof(double));
	CHECKMEM(newtable);
	randtableD(newtable,size,1);
	randshuffletableD(newtable,size);

	if(mols->gausstbl) free(mols->gausstbl);
	mols->ngausstbl=size;
	mols->gausstbl=newtable;
	return 0;
failure:
	simLog(sim,10,"Unable to allocate memory in molssetgausstable");
	return 1; }


/* molsetdifc */
void molsetdifc(simptr sim,int ident,int *index,enum MolecState ms,double difc) {
	int j;
	enum MolecState mslo,mshi;

	if(index) {
		for(j=0;j<index[PDnresults];j++)
			molsetdifc(sim,index[PDMAX+j],NULL,ms,difc); }

	if(ms==MSbsoln) ms=MSsoln;
	else if(ms==MSnone) return;
	if(ms!=MSall) mshi=(enum MolecState)((mslo=ms)+1);
	else {mslo=(enum MolecState)(0);mshi=(enum MolecState)(MSMAX);}

	for(ms=mslo;ms<mshi;ms=(enum MolecState)(ms+1)){
		if(ms!=MSimmobl){
			sim->mols->difc[ident][ms]=difc;
		}else sim->mols->difc[ident][ms]=0;
	}
	molsetcondition(sim->mols,SCparams,0);
	rxnsetcondition(sim,-1,SCparams,0);
	surfsetcondition(sim->srfss,SCparams,0);
	return; }


/* molsetdifm */
int molsetdifm(simptr sim,int ident,int *index,enum MolecState ms,double *difm) {
	int j,d,dim;
	enum MolecState mslo,mshi;
	double *difmat,dm2[DIMMAX*DIMMAX];

	if(index) {
		for(j=0;j<index[PDnresults];j++)
			molsetdifm(sim,index[PDMAX+j],NULL,ms,difm); }

	dim=sim->dim;
	if(!difm) return 0;

	if(ms==MSbsoln) ms=MSsoln;
	else if(ms==MSnone) return 0;
	if(ms!=MSall) mshi=(enum MolecState)((mslo=ms)+1);
	else {mslo=(enum MolecState)(0);mshi=(enum MolecState)(MSMAX);}

	for(ms=mslo;ms<mshi;ms=(enum MolecState)(ms+1)) {
		difmat=sim->mols->difm[ident][ms];
		if(!difmat) {
			difmat=(double*) calloc(sim->dim*sim->dim,sizeof(double));
			CHECKMEM(difmat);
			sim->mols->difm[ident][ms]=difmat; }
		for(d=0;d<sim->dim*sim->dim;d++)
			difmat[d]=difm[d];
		dotMMD(difmat,difmat,dm2,dim,dim,dim);
		sim->mols->difc[ident][ms]=traceMD(dm2,dim)/dim; }

	molsetcondition(sim->mols,SCparams,0);
	return 0;
failure:
	simLog(sim,10,"Unable to allocate memory in molsetdifm");
	return 1; }


/* molsetdrift */
int molsetdrift(simptr sim,int ident,int *index,enum MolecState ms,double *drift) {
	int j,d,dim;
	enum MolecState mslo,mshi;
	double *driftvect;

	if(index) {
		for(j=0;j<index[PDnresults];j++)
			molsetdrift(sim,index[PDMAX+j],NULL,ms,drift); }

	dim=sim->dim;
	if(!drift) return 0;

	if(ms==MSbsoln) ms=MSsoln;
	else if(ms==MSnone) return 0;
	if(ms!=MSall) mshi=(enum MolecState)((mslo=ms)+1);
	else {mslo=(enum MolecState)(0);mshi=(enum MolecState)(MSMAX);}

	for(ms=mslo;ms<mshi;ms=(enum MolecState)(ms+1)) {
		driftvect=sim->mols->drift[ident][ms];
		if(!driftvect) {
			driftvect=(double*) calloc(sim->dim,sizeof(double));
			CHECKMEM(driftvect);
			sim->mols->drift[ident][ms]=driftvect; }
		for(d=0;d<sim->dim;d++)
			driftvect[d]=drift[d]; }

	molsetcondition(sim->mols,SCparams,0);
	return 0;
failure:
	simLog(sim,10,"Unable to allocate memory in molsetdrift");
	return 1; }


/* molsetsurfdrift */
int molsetsurfdrift(simptr sim,int ident,int *index,enum MolecState ms,int surface,enum PanelShape ps,double *drift) {
	int d,i1,j,dim,s1,slo,shi;
	enum MolecState mslo,mshi,ms1;
	enum PanelShape ps1,pslo,pshi;
	molssptr mols;

	if(index) {
		for(j=0;j<index[PDnresults];j++)
			molsetsurfdrift(sim,index[PDMAX+j],NULL,ms,surface,ps,drift); }

	dim=sim->dim;
	mols=sim->mols;

	if(ms!=MSall) mshi=(enum MolecState)((mslo=ms)+1);
	else {mslo=(enum MolecState)(1);mshi=(enum MolecState)(MSMAX);}

	if(surface!=-1) shi=(slo=surface)+1;
	else {slo=0;shi=sim->srfss->nsrf;}

	if(ps!=PSall) pshi=PanelShape((pslo=ps)+1);
	else {pslo=PanelShape(0);pshi=PanelShape(PSMAX);}

	if(!mols->surfdrift) {
		CHECKMEM(mols->surfdrift=(double*****) calloc(mols->maxspecies,sizeof(double****)));
		for(i1=0;i1<mols->maxspecies;i1++) mols->surfdrift[i1]=NULL; }

	if(!mols->surfdrift[ident]) {
		CHECKMEM(mols->surfdrift[ident]=(double****) calloc(MSMAX,sizeof(double***)));
		for(ms1=(enum MolecState)0;ms1<MSMAX;ms1=(enum MolecState)(ms1+1)) mols->surfdrift[ident][ms1]=NULL; }
	for(ms1=mslo;ms1<mshi;ms1=(enum MolecState)(ms1+1)) {
		if(!mols->surfdrift[ident][ms1]) {
			CHECKMEM(mols->surfdrift[ident][ms1]=(double***) calloc(sim->srfss->maxsrf,sizeof(double**)));
			for(s1=0;s1<sim->srfss->maxsrf;s1++) mols->surfdrift[ident][ms1][s1]=NULL; }
		for(s1=slo;s1<shi;s1++) {
			if(!mols->surfdrift[ident][ms1][s1]) {
				CHECKMEM(mols->surfdrift[ident][ms1][s1]=(double**) calloc(PSMAX,sizeof(double*)));
				for(ps1=PanelShape(0);ps1<PSMAX;ps1=PanelShape(ps1+1)) mols->surfdrift[ident][ms1][s1][ps1]=NULL; }
			for(ps1=pslo;ps1<pshi;ps1=PanelShape(ps1+1)) {
				if(!mols->surfdrift[ident][ms1][s1][ps1]) {
					CHECKMEM(mols->surfdrift[ident][ms1][s1][ps1]=(double*) calloc(dim-1,sizeof(double)));
					for(d=0;d<dim-1;d++) mols->surfdrift[ident][ms1][s1][ps1][d]=0; }

				for(d=0;d<dim-1;d++)
					mols->surfdrift[ident][ms1][s1][ps1][d]=drift[d]; }}}
	
	molsetcondition(sim->mols,SCparams,0);
	return 0;
failure:
	simLog(sim,10,"Unable to allocate memory in molsetsurfdrift");
	return 1; }


/* molsetdisplaysize */
void molsetdisplaysize(simptr sim,int ident,int *index,enum MolecState ms,double dsize) {
	int j;
	enum MolecState mslo,mshi;

	if(index) {
		for(j=0;j<index[PDnresults];j++)
			molsetdisplaysize(sim,index[PDMAX+j],NULL,ms,dsize); }

	if(ms==MSbsoln) ms=MSsoln;
	else if(ms==MSnone) return;
	if(ms!=MSall) mshi=(enum MolecState)((mslo=ms)+1);
	else {mslo=(enum MolecState)(0);mshi=(enum MolecState)(MSMAX);}

	for(ms=mslo;ms<mshi;ms=(enum MolecState)(ms+1))
		sim->mols->display[ident][ms]=dsize;

	return; }


/* molsetcolor */
void molsetcolor(simptr sim,int ident,int *index,enum MolecState ms,double *color) {
	int col,j;
	enum MolecState mslo,mshi;

	if(index) {
		for(j=0;j<index[PDnresults];j++)
			molsetcolor(sim,index[PDMAX+j],NULL,ms,color); }

	if(ms==MSbsoln) ms=MSsoln;
	else if(ms==MSnone) return;
	if(ms!=MSall) mshi=(enum MolecState)((mslo=ms)+1);
	else {mslo=(enum MolecState)(0);mshi=(enum MolecState)(MSMAX);}

	for(ms=mslo;ms<mshi;ms=(enum MolecState)(ms+1))
		for(col=0;col<3;col++)
			sim->mols->color[ident][ms][col]=color[col];

	return; }


/* molsetlistlookup */
void molsetlistlookup(simptr sim,int ident,int *index,enum MolecState ms,int ll) {
	int i,j,skip;
	enum MolecState mslo,mshi;

	if(index) {
		for(j=0;j<index[PDnresults];j++)
			molsetlistlookup(sim,index[PDMAX+j],NULL,ms,ll); }

	if(ms==MSbsoln) ms=MSsoln;
	else if(ms==MSnone) return;
	if(ms!=MSall) mshi=(enum MolecState)((mslo=ms)+1);
	else {mslo=(enum MolecState)(0);mshi=(enum MolecState)(MSMAX);}

	if(ident>=0) {
		for(ms=mslo;ms<mshi;ms=(enum MolecState)(ms+1))
			sim->mols->listlookup[ident][ms]=ll; }
	else {
		for(i=0;i<sim->mols->nspecies;i++) {
			for(ms=mslo;ms<mshi;ms=(enum MolecState)(ms+1)) {
				skip=0;
				if(ident==-7 && molismobile(sim,i,ms)==0) skip=1;
				else if(ident==-8 && molismobile(sim,i,ms)==1) skip=1;
				if(!skip) sim->mols->listlookup[i][ms]=ll; }}}
	return; }


/* molsetexist */
void molsetexist(simptr sim,int ident,enum MolecState ms,int exist) {
	if(ident<=0) return;
	if(ms==MSnone) return;
	else if(ms==MSbsoln) ms=MSsoln;
	else if(ms==MSall) {
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) sim->mols->exist[ident][ms]=exist;
		return; }
	sim->mols->exist[ident][ms]=exist;
	return; }


/* molcount */
int molcount(simptr sim,int i,int *index,enum MolecState ms,boxptr bptr,int max) {
	int count,ll,nmol,top,m,j,nresults,uselist;
	moleculeptr *mlist;
	molssptr mols;
	enum MolecState msval;

	mols=sim->mols;
	if(!mols) return 0;
	if(max<0) max=INT_MAX;
	count=0;

	if(i<0 && ms==MSall) {																	// all species, all states
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				if(bptr) count+=bptr->nmol[ll];
				else count+=mols->sortl[ll]; }
		if(!bptr) {
			count+=sim->mols->nd-sim->mols->topd;								// count resurrected molecules
			for(ll=0;ll<mols->nlist;ll++)												// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(mlist[m]->ident>0)
							count++; }}}

	else if(i<0) {																					// all species, one state
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				if(bptr) {
					mlist=bptr->mol[ll];
					top=bptr->nmol[ll]; }
				else {
					mlist=mols->live[ll];
					top=mols->sortl[ll]; }
				for(m=0;m<top && count<max;m++)										// count properly sorted molecules
					if(mlist[m]->ident>0 && mlist[m]->mstate==ms)
						count++; }
		if(!bptr) {
			mlist=mols->dead;																		// count resurrected molecules
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(mlist[m]->ident>0 && mlist[m]->mstate==ms)
					count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(mlist[m]->ident>0 && mlist[m]->mstate==ms)
							count++; }}}

	else if(i>0 && ms==MSall) {															// one species, all states
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				uselist=0;
				for(msval=(enum MolecState)0;msval<MSMAX && !uselist;msval=(enum MolecState)(msval+1))
					if(mols->listlookup[i][msval]==ll) uselist=1;
				if(uselist) {
					if(bptr) {
						mlist=bptr->mol[ll];
						top=bptr->nmol[ll]; }
					else {
						mlist=mols->live[ll];
						top=mols->sortl[ll]; }
					for(m=0;m<top && count<max;m++)										// count properly sorted molecules
						if(mlist[m]->ident==i)
							count++; }}
		if(!bptr) {
			mlist=mols->dead;																			// count resurrected molecules
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(mlist[m]->ident==i)
					count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(mlist[m]->ident==i)
							count++; }}}

	else if(i>0) {																					// single species, single state
		ll=mols->listlookup[i][ms];
		if(bptr) {
			mlist=bptr->mol[ll];
			top=bptr->nmol[ll]; }
		else {
			mlist=mols->live[ll];
			top=mols->sortl[ll]; }
		for(m=0;m<top && count<max;m++)												// count properly sorted molecules
			if(mlist[m]->ident==i && mlist[m]->mstate==ms) count++;
		if(!bptr) {
			mlist=mols->dead;																		// count resurrected molecules
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(mlist[m]->ident==i && mlist[m]->mstate==ms) count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(mlist[m]->ident==i && mlist[m]->mstate==ms) count++; }}}

	else if(index) {																				// list of molecules in index
		nresults=index[PDnresults];
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				uselist=0;
				if(ms==MSall) {
					for(j=0;j<nresults;j++)
						for(msval=(enum MolecState)0;msval<MSMAX && !uselist;msval=(enum MolecState)(msval+1))
							if(mols->listlookup[index[PDMAX+j]][msval]==ll) uselist=1; }
				else {
					for(j=0;j<nresults;j++)
						if(mols->listlookup[index[PDMAX+j]][ms]==ll) uselist=1; }
				if(uselist) {
					if(bptr) {
						mlist=bptr->mol[ll];
						top=bptr->nmol[ll]; }
					else {
						mlist=mols->live[ll];
						top=mols->sortl[ll]; }
					for(m=0;m<top && count<max;m++)									// count properly sorted molecules
						if(locateVi(index+PDMAX,mlist[m]->ident,nresults,0)>=0 && (ms==MSall || mlist[m]->mstate==ms))
							count++; }}
		if(!bptr) {																			// count resurrected molecules
			mlist=mols->dead;
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(locateVi(index+PDMAX,mlist[m]->ident,nresults,0)>=0 && (ms==MSall || mlist[m]->mstate==ms))
					count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(locateVi(index+PDMAX,mlist[m]->ident,nresults,0)>=0 && (ms==MSall || mlist[m]->mstate==ms))
							count++; }}}

	return count; }


/* molcount_cplx */
int molcount_cplx(simptr sim,int i,int *index,enum MolecState ms,boxptr bptr,int max) {
	int count,ll,nmol,top,m,j,nresults,uselist;
	moleculeptr *mlist;
	molssptr mols;
	enum MolecState msval;

	mols=sim->mols;
	if(!mols) return 0;
	if(max<0) max=INT_MAX;
	count=0;

	if(i<0 && ms==MSall) {														// all species, all states
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				if(bptr) count+=bptr->nmol[ll];
				else count+=mols->sortl[ll]; }
		if(!bptr) {
			count+=sim->mols->nd-sim->mols->topd;								// count resurrected molecules
			for(ll=0;ll<mols->nlist;ll++)										// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						//if(mlist[m]->ident>0)
						if(mlist[m]->ident>0 && mlist[m]->s_index==0)
							count++; }}}

	else if(i<0) {																					// all species, one state
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				if(bptr) {
					mlist=bptr->mol[ll];
					top=bptr->nmol[ll]; }
				else {
					mlist=mols->live[ll];
					top=mols->sortl[ll]; }
				for(m=0;m<top && count<max;m++)										// count properly sorted molecules
					if(mlist[m]->ident>0 && mlist[m]->mstate==ms && mlist[m]->s_index==0)
						count++; }
		if(!bptr) {
			mlist=mols->dead;																		// count resurrected molecules
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(mlist[m]->ident>0 && mlist[m]->mstate==ms && mlist[m]->s_index==0)
					count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(mlist[m]->ident>0 && mlist[m]->mstate==ms && mlist[m]->s_index==0)
							count++; }}}

	else if(i>0 && ms==MSall) {															// one species, all states
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				uselist=0;
				for(msval=(enum MolecState)0;msval<MSMAX && !uselist;msval=(enum MolecState)(msval+1))
					if(mols->listlookup[i][msval]==ll) uselist=1;
				if(uselist) {
					if(bptr) {
						mlist=bptr->mol[ll];
						top=bptr->nmol[ll]; }
					else {
						mlist=mols->live[ll];
						top=mols->sortl[ll]; }
					for(m=0;m<top && count<max;m++)										// count properly sorted molecules
						if(mlist[m]->ident==i && mlist[m]->s_index==0)
							count++; }}
		if(!bptr) {
			mlist=mols->dead;																			// count resurrected molecules
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(mlist[m]->ident==i && mlist[m]->s_index==0)
					count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(mlist[m]->ident==i && mlist[m]->s_index==0)
							count++; }}}

	else if(i>0) {																					// single species, single state
		ll=mols->listlookup[i][ms];
		if(bptr) {
			mlist=bptr->mol[ll];
			top=bptr->nmol[ll]; }
		else {
			mlist=mols->live[ll];
			top=mols->sortl[ll]; }
		for(m=0;m<top && count<max;m++)												// count properly sorted molecules
			if(mlist[m]->ident==i && mlist[m]->mstate==ms && mlist[m]->s_index==0) count++;
		if(!bptr) {
			mlist=mols->dead;																		// count resurrected molecules
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(mlist[m]->ident==i && mlist[m]->mstate==ms && mlist[m]->s_index==0) count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(mlist[m]->ident==i && mlist[m]->mstate==ms && mlist[m]->s_index==0) count++; }}}

	else if(index) {																				// list of molecules in index
		nresults=index[PDnresults];
		for(ll=0;ll<mols->nlist;ll++)
			if(mols->listtype[ll]==MLTsystem) {
				uselist=0;
				if(ms==MSall) {
					for(j=0;j<nresults;j++)
						for(msval=(enum MolecState)0;msval<MSMAX && !uselist;msval=(enum MolecState)(msval+1))
							if(mols->listlookup[index[PDMAX+j]][msval]==ll) uselist=1; }
				else {
					for(j=0;j<nresults;j++)
						if(mols->listlookup[index[PDMAX+j]][ms]==ll) uselist=1; }
				if(uselist) {
					if(bptr) {
						mlist=bptr->mol[ll];
						top=bptr->nmol[ll]; }
					else {
						mlist=mols->live[ll];
						top=mols->sortl[ll]; }
					for(m=0;m<top && count<max;m++)									// count properly sorted molecules
						if(locateVi(index+PDMAX,mlist[m]->ident,nresults,0)>=0 && (ms==MSall || mlist[m]->mstate==ms) && mlist[m]->s_index==0)
							count++; }}
		if(!bptr) {																						// count resurrected molecules
			mlist=mols->dead;
			nmol=mols->nd;
			top=mols->topd;
			for(m=top;m<nmol && count<max;m++)
				if(locateVi(index+PDMAX,mlist[m]->ident,nresults,0)>=0 && (ms==MSall || mlist[m]->mstate==ms) && mlist[m]->s_index==0)
					count++;
			for(ll=0;ll<mols->nlist;ll++)													// count molecules that need sorting
				if(mols->listtype[ll]==MLTsystem) {
					mlist=mols->live[ll];
					nmol=mols->nl[ll];
					top=mols->sortl[ll];
					for(m=top;m<nmol && count<max;m++)
						if(locateVi(index+PDMAX,mlist[m]->ident,nresults,0)>=0 && (ms==MSall || mlist[m]->mstate==ms) && mlist[m]->s_index==0)
							count++; }}}

	return count; }


/* MolCalcDifcSum */
//double MolCalcDifcSum(simptr sim, int i1, enum MolecState ms1, boxptr box_tmp1, int i2, enum MolecState ms2, boxptr box_tmp2){
double MolCalcDifcSum(simptr sim, moleculeptr mptr1, moleculeptr mptr2, double *dc1, double *dc2){
	double sum1,sum2, adj_factor1,adj_factor2;
	int i1, i2, s, dif_ident;
	enum MolecState ms1, ms2;
	moleculeptr mptr_bind1, dif_molec1, mptr_bind2, dif_molec2,cplx_dif;
	GSList *r;
	intptr_t r_indx;

	sum1=sum2=0;
	if(mptr1){
		dif_ident=difmolec(sim,mptr1,&dif_molec1);
		if(dif_molec1){
			i1=dif_ident;
			ms1=dif_molec1->mstate;
		}
		else{
			i1=mptr1->ident;
			ms1=mptr1->mstate;
		}
		sum1=sim->mols->difc[i1][ms1]; //*adj_factor1;
		if(dc1) *dc1=sum1;
	}
	if(mptr2){
		dif_ident=difmolec(sim,mptr2,&dif_molec2);
		if(dif_molec2){
			i2=dif_ident;
			ms2=dif_molec1->mstate;
		}
		else{
			i2=mptr1->ident;
			ms2=mptr1->mstate;
		}
		sum2=sim->mols->difc[i2][ms2];//*adj_factor2;
		if(dc2) *dc2=sum2;
	}

	// printf("sum1=%f, sum2=%f, adj_factor1=%f, adj_factor2=%f\n", sum1, sum2, adj_factor1, adj_factor2);	

	return sum1+sum2;
}


/******************************************************************************/
/****************************** memory management *****************************/
/******************************************************************************/

/* molalloc */
moleculeptr molalloc(simptr sim, int dim) {
	moleculeptr mptr;
	int d, d1;

	mptr=NULL;
	CHECKMEM(mptr=(moleculeptr) malloc(sizeof(struct moleculestruct)));
	mptr->serno=0;
	mptr->list=-1;
	mptr->m=-1;

	mptr->pos=NULL;				// potentially double free
	mptr->posx=NULL;
	mptr->via=NULL;
	mptr->posoffset=NULL;
	mptr->ident=0;
	mptr->mstate=MSsoln;
	mptr->box=NULL;
	mptr->pnl=NULL;
	//cplx
	mptr->s_index=-1;
	mptr->to=NULL;				// potentially double free
	mptr->from=NULL;			// potentially double free
	mptr->tot_sunit=0;	
	mptr->pos_tmp=NULL;			// potentially double free
	mptr->theta_init=0;
	mptr->phi_init=0;
	mptr->sdist_tmp=0;
	mptr->sdist_init=0;
	mptr->prev_pos=NULL;
	mptr->complex_id=-1;				// complex_id
	mptr->sites=NULL;
	mptr->sites_val=-1;
	mptr->sites_valx=-1;
	mptr->dif_molec=NULL;
	mptr->dif_site=-1;
	mptr->sim_time=-1;
	mptr->vchannel=NULL;
	mptr->arrival_time=-1;
	mptr->bind_id=-1;

	CHECKMEM(mptr->pos=(double*) calloc(dim,sizeof(double)));
	CHECKMEM(mptr->posx=(double*) calloc(dim,sizeof(double)));
	CHECKMEM(mptr->via=(double*) calloc(dim,sizeof(double)));
	CHECKMEM(mptr->posoffset=(double*) calloc(dim,sizeof(double)));
	CHECKMEM(mptr->prev_pos=(double*) calloc(dim,sizeof(double)));
	
	/*
	for(d=0;d<3;d++){
		for(d1=0;d1<3;d1++) {
			if(d1==d) mptr->rotate_mtrx[d][d1]=1;
			else mptr->rotate_mtrx[d][d1]=0;		
	}}
	*/
	
	for(d=0;d<dim;d++){
		mptr->pos[d]=mptr->posx[d]=mptr->via[d]=mptr->posoffset[d]=mptr->prev_pos[d]=0;	}
	mptr->pos_tmp=mptr->pos;
	return mptr;
	
 failure:
	molfree(sim,mptr);
	simLog(NULL,10,"Unable to allocate memory in molalloc");
	return NULL; }


/* molfree */
void molfree(simptr sim, moleculeptr mptr) {
	int d,k;
	siteptr site_tmp;
	if(!mptr) return;
	
	// printf("mptr->serno=%d\n", mptr->serno);
	if(mptr->to)	mptr->to=NULL;
	if(mptr->from)	mptr->from=NULL;
	if(mptr->dif_molec) mptr->dif_molec=NULL;

	if(mptr->pos) {	
		mptr->pos=mptr->pos_tmp;	
		free(mptr->pos);			
		mptr->pos=NULL;
		mptr->pos_tmp=NULL;
 	}

	if(mptr->via) {
		free(mptr->via);				
		mptr->via=NULL;		
	}

	if(mptr->posx) {
		free(mptr->posx);			
		mptr->posx=NULL;
	}

	if(mptr->posoffset) {
		free(mptr->posoffset);	
		mptr->posoffset=NULL; 
	}
	
	if(mptr->prev_pos) {
		free(mptr->prev_pos);
		mptr->prev_pos=NULL;
	}
	
	if(mptr->sites) {
		for(k=0;k<sim->mols->spsites_num[mptr->ident];k++) {
			site_tmp=mptr->sites[k];
			if(site_tmp->bind) site_tmp->bind=NULL;
			if(site_tmp->value) free(site_tmp->value); site_tmp->value=NULL;
			free(site_tmp);
		}
		free(mptr->sites);
	}			

	if(mptr->vchannel) {
		fclose(mptr->vchannel->voltage_file);
		mptr->vchannel->voltage_file=NULL;
		free(mptr->vchannel);
		mptr->vchannel=NULL;
	}

	if(mptr) free(mptr); mptr=NULL;
	return;
}

/* complexptr */
complexptr complexalloc(simptr sim, moleculeptr mptr){
	complexptr cplxptr;

	cplxptr=NULL;
	CHECKMEM(cplxptr=(complexptr) malloc(sizeof(struct complexstruct)));	
	cplxptr->zeroindx_molec=mptr;	
	cplxptr->dif_molec=NULL;
	cplxptr->dif_bind=NULL;
	cplxptr->dif_bind_site=-1;
	cplxptr->serno=0;
	cplxptr->diffuse_updated=0;
	return cplxptr;
failure:
	simLog(sim,10,"Unable to allocate memory in complexalloc()");
	return NULL;
} 

void complexfree(complexptr cplxptr){
	if(!cplxptr) return;
	if(cplxptr->zeroindx_molec)
		cplxptr->zeroindx_molec=NULL;
	if(cplxptr->dif_molec)
		cplxptr->dif_molec=NULL;
	if(cplxptr->dif_bind)
		cplxptr->dif_bind=NULL;
	free(cplxptr);
	return;
}

/* molexpandsurfdrift */
int molexpandsurfdrift(simptr sim,int oldmaxspec,int oldmaxsrf) {	//?? needs to be called whenever maxspecies or maxsrf increase
	double *****oldsurfdrift;
	int i,s;
	enum MolecState ms;
	enum PanelShape ps;
	
	if(!sim->mols->surfdrift) return 0;
	oldsurfdrift=sim->mols->surfdrift;
	sim->mols->surfdrift=NULL;
	
	for(i=0;i<oldmaxspec;i++)
		if(oldsurfdrift[i])
			for(ms=(enum MolecState)0;ms<MSMAX;ms=(enum MolecState)(ms+1))
				if(oldsurfdrift[i][ms])
					for(s=0;s<oldmaxsrf;s++)
						if(oldsurfdrift[i][ms][s])
							for(ps=PanelShape(0);ps<PSMAX;ps=PanelShape(ps+1))
								if(oldsurfdrift[i][ms][s][ps]) {
									CHECK(molsetsurfdrift(sim,i,NULL,ms,s,ps,oldsurfdrift[i][ms][s][ps])==0); }
	
	molfreesurfdrift(oldsurfdrift,oldmaxspec,oldmaxsrf);
	return 0;
failure:
	return 1; }


/* molfreesurfdrift */
void molfreesurfdrift(double *****surfdrift,int maxspec,int maxsrf) {
	int i,s;
	enum MolecState ms;
	enum PanelShape ps;
	
	if(surfdrift) {
		for(i=0;i<maxspec;i++)
			if(surfdrift[i]) {
				for(ms=(enum MolecState)0;ms<MSMAX;ms=(enum MolecState)(ms+1))
					if(surfdrift[i][ms]) {
						for(s=0;s<maxsrf;s++)
							if(surfdrift[i][ms][s]) {
								for(ps=PanelShape(0);ps<PSMAX;ps=PanelShape(ps+1))
									free(surfdrift[i][ms][s][ps]);
								free(surfdrift[i][ms][s]); }
						free(surfdrift[i][ms]); }
				free(surfdrift[i]); }
		free(surfdrift); }
	return; }


/* molpatternindexalloc */
int molpatternindexalloc(int **indexptr,int n) {
	int j;
	int *index,*newindex;

	index=*indexptr;
	if(n<PDMAX) n=(index?index[PDalloc]*2:PDMAX+1);
	newindex=(int*) calloc(n,sizeof(int));
	if(!newindex) return 1;
	j=0;
	if(index)
		for(;j<index[PDalloc] && j<n;j++)
			newindex[j]=index[j];
	for(;j<n;j++)
		newindex[j]=0;
	newindex[PDalloc]=n;
	free(index);
	*indexptr=newindex;
	return 0; }


/* molpatternalloc */
int molpatternalloc(simptr sim,int maxpattern) {
	int i,er;
	char **newpatlist;
	int **newpatindex;

	newpatlist=(char **) calloc(maxpattern,sizeof(char*));
	if(!newpatlist) return 1;
	newpatindex=(int **) calloc(maxpattern,sizeof(int*));
	if(!newpatindex) return 1;

	for(i=0;i<sim->mols->maxpattern;i++) {
		newpatlist[i]=sim->mols->patlist[i];
		newpatindex[i]=sim->mols->patindex[i]; }
	for(;i<maxpattern;i++) {
		newpatlist[i]=EmptyString();
		if(!newpatlist[i]) return 1;
		newpatlist[i][0]='\0';
		newpatindex[i]=NULL;
		er=molpatternindexalloc(&newpatindex[i],8);
		if(er) return 1; }

	free(sim->mols->patlist);
	free(sim->mols->patindex);
	sim->mols->maxpattern=maxpattern;
	sim->mols->patlist=newpatlist;
	sim->mols->patindex=newpatindex;
	return 0; }


/* molssalloc */
molssptr molssalloc(molssptr mols,int maxspecies) {
	int newmols,i,j,**newexist,**newlistlookup,*newexpand,oldmaxspecies, *newspsites_num, *newvolt_dependent;// *newspdifsites;
	enum MolecState ms;
	char **newspname;
	double **newdifc,**newdifstep,***newdifm,***newdrift,**newdisplay,***newcolor;

	if(maxspecies<1) return NULL;
	maxspecies++;

	newmols=0;

	if(!mols) {
		mols=(molssptr) malloc(sizeof(struct molsuperstruct));
		CHECKMEM(mols);
		newmols=1;

		mols->condition=SCinit;
		mols->sim=NULL;
		mols->maxspecies=0;
		mols->nspecies=1;
		mols->spname=NULL;
		mols->maxpattern=0;
		mols->npattern=0;
		mols->patlist=NULL;
		mols->patindex=NULL;
		mols->difc=NULL;
		mols->difstep=NULL;
		mols->difm=NULL;
		mols->drift=NULL;
		mols->surfdrift=NULL;
		mols->display=NULL;
		mols->color=NULL;
		mols->exist=NULL;
		mols->dead=NULL;
		mols->maxdlimit=-1;
		mols->maxd=0;
		mols->nd=0;
		mols->topd=0;
		mols->maxlist=0;
		mols->nlist=0;
		mols->listlookup=NULL;
		mols->listname=NULL;
		mols->listtype=NULL;
		mols->live=NULL;
		mols->maxl=NULL;
		mols->nl=NULL;
		mols->topl=NULL;
		mols->sortl=NULL;
		mols->diffuselist=NULL;
		mols->serno=1;
		mols->ngausstbl=0;
		mols->gausstbl=NULL;
		mols->expand=NULL; 

		mols->complexlist=NULL;
		mols->ncomplex=0; 		//-1;
		mols->max_complex=0;
		mols->spsites_num=NULL;
		mols->spsites_name=NULL;
		mols->spsites_binding=NULL;
		mols->maxsitecode=0;
		// mols->max_sites=5;
		mols->complex_connect=g_hash_table_new(g_direct_hash,g_direct_equal);
		mols->spdifsites=g_hash_table_new(g_direct_hash,g_direct_equal);
		mols->Mlist=NULL;
		mols->volt_dependent=NULL;

	}
	

	if(maxspecies>mols->maxspecies) {
		oldmaxspecies=mols->maxspecies;

		CHECKMEM(newspname=(char**) calloc(maxspecies,sizeof(char*)));
		for(i=0;i<maxspecies;i++) newspname[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newspname[i]=mols->spname[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newspname[i]=EmptyString()); }
		strncpy(newspname[0],"empty",STRCHAR-1);

		CHECKMEM(newdifc=(double**) calloc(maxspecies,sizeof(double*)));
		for(i=0;i<maxspecies;i++) newdifc[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newdifc[i]=mols->difc[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newdifc[i]=(double*) calloc(MSMAX,sizeof(double)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newdifc[i][ms]=0; }

		CHECKMEM(newdifstep=(double**) calloc(maxspecies,sizeof(double*)));
		for(i=0;i<maxspecies;i++) newdifstep[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newdifstep[i]=mols->difstep[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newdifstep[i]=(double*) calloc(MSMAX,sizeof(double)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newdifstep[i][ms]=0; }

		CHECKMEM(newdifm=(double***) calloc(maxspecies,sizeof(double**)));
		for(i=0;i<maxspecies;i++) newdifm[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newdifm[i]=mols->difm[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newdifm[i]=(double**) calloc(MSMAX,sizeof(double*)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newdifm[i][ms]=NULL; }

		CHECKMEM(newdrift=(double***) calloc(maxspecies,sizeof(double**)));
		for(i=0;i<maxspecies;i++) newdrift[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newdrift[i]=mols->drift[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newdrift[i]=(double**) calloc(MSMAX,sizeof(double*)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newdrift[i][ms]=NULL; }

		CHECKMEM(newdisplay=(double**) calloc(maxspecies,sizeof(double*)));
		for(i=0;i<maxspecies;i++) newdisplay[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newdisplay[i]=mols->display[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newdisplay[i]=(double*) calloc(MSMAX,sizeof(double)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newdisplay[i][ms]=3; }

		CHECKMEM(newcolor=(double ***) calloc(maxspecies,sizeof(double **)));
		for(i=0;i<maxspecies;i++) newcolor[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newcolor[i]=mols->color[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newcolor[i]=(double**) calloc(MSMAX,sizeof(double*)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newcolor[i][ms]=NULL;
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) {
				CHECKMEM(newcolor[i][ms]=(double*) calloc(3,sizeof(double)));
				newcolor[i][ms][0]=newcolor[i][ms][1]=newcolor[i][ms][2]=0; }}

		CHECKMEM(newexist=(int**) calloc(maxspecies,sizeof(int*)));
		for(i=0;i<maxspecies;i++) newexist[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newexist[i]=mols->exist[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newexist[i]=(int*) calloc(MSMAX,sizeof(int)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newexist[i][ms]=0; }

		CHECKMEM(newlistlookup=(int**) calloc(maxspecies,sizeof(int*)));
		for(i=0;i<maxspecies;i++) newlistlookup[i]=NULL;
		for(i=0;i<oldmaxspecies;i++) newlistlookup[i]=mols->listlookup[i];
		for(;i<maxspecies;i++) {
			CHECKMEM(newlistlookup[i]=(int*) calloc(MSMAX,sizeof(int)));
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) newlistlookup[i][ms]=-1; }

		CHECKMEM(newexpand=(int*) calloc(maxspecies,sizeof(int)));
		for(i=0;i<oldmaxspecies;i++) newexpand[i]=mols->expand[i];
		for(;i<maxspecies;i++) newexpand[i]=0;
		
		CHECKMEM(newspsites_num=(int*) calloc(maxspecies,sizeof(int)));
		for(i=0;i<oldmaxspecies;i++) newspsites_num[i]=mols->spsites_num[i];
		for(;i<maxspecies;i++) newspsites_num[i]=0;

		CHECKMEM(newvolt_dependent=(int*) calloc(maxspecies,sizeof(int)));
		for(i=0;i<oldmaxspecies;i++) newvolt_dependent[i]=mols->volt_dependent[i];
		for(;i<maxspecies;i++) newvolt_dependent[i]=0;


		mols->maxspecies=maxspecies;
		free(mols->spname);
		mols->spname=newspname;
		free(mols->difc);
		mols->difc=newdifc;
		free(mols->difstep);
		mols->difstep=newdifstep;
		free(mols->difm);
		mols->difm=newdifm;
		free(mols->drift);
		mols->drift=newdrift;
		free(mols->display);
		mols->display=newdisplay;
		free(mols->color);
		mols->color=newcolor;
		free(mols->exist);
		mols->exist=newexist;
		free(mols->listlookup);
		mols->listlookup=newlistlookup;
		free(mols->expand);
		mols->expand=newexpand;
		free(mols->spsites_num);
		mols->spsites_num=newspsites_num;
		free(mols->volt_dependent);
		mols->volt_dependent=newvolt_dependent;
		//g_hash_table_destroy(mols->spdifsites);
		if(mols->surfdrift && mols->sim->srfss) { CHECK(molexpandsurfdrift(mols->sim,oldmaxspecies,mols->sim->srfss->maxsrf)==0); }}

	return mols;

 failure:
	simLog(NULL,10,"Unable to allocate memory in molssalloc");
	return NULL; }


/* mollistalloc */
int mollistalloc(molssptr mols,int maxlist,enum MolListType mlt) {
	int *maxl,*nl,*topl,*sortl,*diffuselist,ll,m;
	moleculeptr **live,mptr;
	char **listname;
	enum MolListType *listtype;
	// complexptr *complexlist;

	if(maxlist<=0) return -2;
	if(!mols) return -3;
	maxlist+=mols->maxlist;

	listname=NULL;							// allocate new arrays
	listtype=NULL;
	live=NULL;
	maxl=NULL;
	nl=NULL;
	topl=NULL;
	sortl=NULL;
	diffuselist=NULL;

	CHECKMEM(listname=(char**) calloc(maxlist,sizeof(char*)));
	for(ll=0;ll<maxlist;ll++) listname[ll]=NULL;

	CHECKMEM(listtype=(enum MolListType*) calloc(maxlist,sizeof(enum MolListType)));
	for(ll=0;ll<maxlist;ll++) listtype[ll]=MLTnone;

	CHECKMEM(live=(moleculeptr**) calloc(maxlist,sizeof(moleculeptr*)));
	for(ll=0;ll<maxlist;ll++) live[ll]=NULL;

	CHECKMEM(maxl=(int*) calloc(maxlist,sizeof(int)));
	for(ll=0;ll<maxlist;ll++) maxl[ll]=0;

	CHECKMEM(nl=(int*) calloc(maxlist,sizeof(int)));
	for(ll=0;ll<maxlist;ll++) nl[ll]=0;

	CHECKMEM(topl=(int*) calloc(maxlist,sizeof(int)));
	for(ll=0;ll<maxlist;ll++) topl[ll]=0;

	CHECKMEM(sortl=(int*) calloc(maxlist,sizeof(int)));
	for(ll=0;ll<maxlist;ll++) sortl[ll]=0;

	CHECKMEM(diffuselist=(int*) calloc(maxlist,sizeof(int)));
	for(ll=0;ll<maxlist;ll++) diffuselist[ll]=0;

	for(ll=0;ll<mols->maxlist;ll++) {			// copy over existing portions
		listname[ll]=mols->listname[ll];
		listtype[ll]=mols->listtype[ll];
		live[ll]=mols->live[ll];
		maxl[ll]=mols->maxl[ll];
		nl[ll]=mols->nl[ll];
		topl[ll]=mols->topl[ll];
		sortl[ll]=mols->sortl[ll];
		diffuselist[ll]=mols->diffuselist[ll];
	 }

	for(ll=mols->maxlist;ll<maxlist;ll++) {					// listnames and listtypes
		CHECKMEM(listname[ll]=EmptyString());
		listtype[ll]=mlt; }

	for(ll=mols->maxlist;ll<maxlist;ll++) maxl[ll]=1;		// calculate maxl
	for(m=mols->topd;m<mols->nd;m++) {
		mptr=mols->dead[m];
		if(mptr && mptr->list>=mols->maxlist && mptr->list<maxlist) maxl[mptr->list]++; }
	for(ll=mols->maxlist;ll<maxlist;ll++) {
		maxl[ll]*=2;
		if(maxl[ll]>mols->maxd) maxl[ll]=mols->maxd; }

	for(ll=mols->maxlist;ll<maxlist;ll++) {					// allocate live lists
		CHECKMEM(live[ll]=(moleculeptr*) calloc(maxl[ll],sizeof(moleculeptr)));
		for(m=0;m<maxl[ll];m++) live[ll][m]=NULL; 
	}

	if(mols->maxlist) {										// free any old lists
		free(mols->listname);
		free(mols->listtype);
		free(mols->live);
		free(mols->maxl);
		free(mols->nl);
		free(mols->topl);
		free(mols->sortl);
		free(mols->diffuselist); 
	}

	ll=mols->maxlist;
	mols->maxlist=maxlist;									// store new lists
	mols->listname=listname;
	mols->listtype=listtype;
	mols->live=live;
	mols->maxl=maxl;
	mols->nl=nl;
	mols->topl=topl;
	mols->sortl=sortl;
	mols->diffuselist=diffuselist;
	return ll;

 failure:
	if(listname)
		for(ll=0;ll<maxlist;ll++) free(listname[ll]);
	free(listname);
	free(listtype);
	if(live)
		for(ll=mols->maxlist;ll<maxlist;ll++) free(live[ll]);
	free(live);
	free(maxl);
	free(nl);
	free(topl);
	free(sortl);
	free(diffuselist);
	simLog(NULL,10,"Unable to allocate memory in mollistalloc");
	return -1; }


/* molexpandlist */
int molexpandlist(molssptr mols,int dim,int ll,int nspaces,int nmolecs) {
	moleculeptr *newlist,*oldlist;
	int m,nold,maxold,maxnew;

	if(!mols || ll>=mols->nlist) return 2;
	if(ll>=0 && nmolecs>0) return 2;							// can't add molecules to live list

	maxold=ll<0?mols->maxd:mols->maxl[ll];						// maxold is previous allocated size
	nold=ll<0?mols->nd:mols->nl[ll];							// nold is previous number of molecules
	oldlist=ll<0?mols->dead:mols->live[ll];						// oldlist is previous list

	maxnew=nspaces>0?maxold+nspaces:2*maxold+1;					// maxnew is new allocated size
	if(nold+nmolecs>maxnew) return 3;
	// printf("molexpandlist, line 1480, nold=%d, nmolecs=%d, maxnew=%d\n", nold, nmolecs, maxnew);
	
	newlist=(moleculeptr*) calloc(maxnew,sizeof(moleculeptr));
	CHECKMEM(newlist);
	for(m=0;m<maxold;m++) newlist[m]=oldlist[m];
	for(;m<maxnew;m++) newlist[m]=NULL;
	if(ll<0) {
		free(mols->dead);
		mols->dead=newlist;
		mols->maxd=maxnew; }
	else {
		free(mols->live[ll]);
		mols->live[ll]=newlist;
		mols->maxl[ll]=maxnew; }
	
	if(nmolecs) {
		for(m=mols->nd-1;m>=mols->topd;m--) {					// copy resurrected molecules higher on list
			newlist[m+nmolecs]=newlist[m];
			newlist[m]=NULL; }
		for(m=mols->topd;m<mols->topd+nmolecs;m++) {		// create new empty molecules
			newlist[m]=molalloc(mols->sim,dim);
			if(!newlist[m]) return 4; }
		mols->topd+=nmolecs;
		mols->nd+=nmolecs; }
	return 0;
failure:
	simLog(NULL,10,"Unable to allocate memory in molexpandlist");
	return 1; }


/* molssfree */
void molssfree(molssptr mols,int maxsrf) {
	int m,ll,i,maxspecies;
	enum MolecState ms;

	if(!mols) return;
	maxspecies=mols->maxspecies;

	free(mols->expand);
	free(mols->gausstbl);

	// free(mols->spdifsites);
	if(mols->spdifsites) g_hash_table_destroy(mols->spdifsites);

	if(mols->Mlist){
		for(ll=0;ll<mols->nl[ll];ll++) free(mols->Mlist[ll]);
	}

	for(ll=0;ll<mols->maxlist;ll++) {
		if(mols->listname) free(mols->listname[ll]);
		if(mols->live && mols->live[ll]) {
			for(m=0;m<mols->nl[ll];m++)
				molfree(mols->sim,mols->live[ll][m]);
			printf("m=%d, ll=%d \n", m, ll);
			free(mols->live[ll]); 
	}}
	free(mols->diffuselist);
	free(mols->sortl);
	free(mols->topl);
	free(mols->nl);
	free(mols->maxl);
	free(mols->live);
	free(mols->listtype);
	free(mols->listname);
	
	if(mols->complexlist) {
		for(i=0;i<mols->max_complex;i++) 
			complexfree(mols->complexlist[i]);
		free(mols->complexlist);
	}	

	if(mols->listlookup) {
		for(i=0;i<maxspecies;i++) free(mols->listlookup[i]);
		free(mols->listlookup); }

	if(mols->exist) {
		for(i=0;i<maxspecies;i++) free(mols->exist[i]);
		free(mols->exist); }

	if(mols->dead) {
		for(m=0;m<mols->nd;m++) molfree(mols->sim,mols->dead[m]);
		free(mols->dead); }

	if(mols->color) {
		for(i=0;i<maxspecies;i++)
			if(mols->color[i]) {
				for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) free(mols->color[i][ms]);
				free(mols->color[i]); }
		free(mols->color); }

	if(mols->display) {
		for(i=0;i<maxspecies;i++) free(mols->display[i]);
		free(mols->display); }

	molfreesurfdrift(mols->surfdrift,mols->maxspecies,maxsrf);
	
	if(mols->drift) {
		for(i=0;i<maxspecies;i++)
			if(mols->drift[i]) {
				for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) free(mols->drift[i][ms]);
				free(mols->drift[i]); }
		free(mols->drift); }

	if(mols->difm) {
		for(i=0;i<maxspecies;i++)
			if(mols->difm[i]) {
				for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) free(mols->difm[i][ms]);
				free(mols->difm[i]); }
		free(mols->difm); }

	if(mols->difstep) {
		for(i=0;i<maxspecies;i++) free(mols->difstep[i]);
		free(mols->difstep); }

	if(mols->difc) {
		for(i=0;i<maxspecies;i++) free(mols->difc[i]);
		free(mols->difc); }
	/*
	if(mols->dif_cmpt){
		for(i=0;i<maxspecies;i++){
			if(mols->dif_cmpt[i]){	
				difadj_tmp=mols->dif_cmpt[i];
				while(difadj_tmp->next){
					difadj_next=difadj_tmp->next;
					free(difadj_tmp);
					difadj_tmp=difadj_next;
				}			
				free(difadj_tmp);
			}}
		free(mols->dif_cmpt);
	}	
	*/

	if(mols->patindex) {
		for(i=0;i<mols->maxpattern;i++) free(mols->patindex[i]);
		free(mols->patindex); }

	if(mols->patlist) {
		for(i=0;i<mols->maxpattern;i++) free(mols->patlist[i]);
		free(mols->patlist); }

	if(mols->spname) {
		for(i=0;i<maxspecies;i++) free(mols->spname[i]);
		free(mols->spname); }
	
	if(mols->spsites_name)
		for(i=0;i<mols->maxspecies;i++) free(mols->spsites_name[i]);
	
	if(mols->spsites_binding)
		for(i=0;i<mols->maxspecies;i++) free(mols->spsites_binding[i]);
	
	if(mols->spsites_num)	free(mols->spsites_num);
	if(mols->volt_dependent) free(mols->volt_dependent);
	if(mols->complex_connect) {
		g_hash_table_destroy(mols->complex_connect);
	}

	free(mols);
	return; }



/******************************************************************************/
/*************************** data structure output ****************************/
/******************************************************************************/

/* molssoutput */
void molssoutput(simptr sim) {
	int nspecies,i,ll,same,sum,s;
	molssptr mols;
	char string[STRCHAR];
	double maxstep;
	enum MolecState ms;
	enum PanelShape ps;

	simLog(sim,2,"MOLECULE PARAMETERS\n");
	if(!sim || !sim->mols) {
		simLog(sim,2," No molecule superstructure defined\n\n");
		return; }
	mols=sim->mols;
	nspecies=mols->nspecies;

	if(mols->condition!=SCok)
		simLog(sim,7," Molecule superstructure condition: %s\n",simsc2string(mols->condition,string));
	simLog(sim,1," Next molecule serial number: %li\n",mols->serno);
	if(mols->gausstbl) simLog(sim,1," Table for Gaussian distributed random numbers has %i values\n",mols->ngausstbl);
	else simLog(sim,1," Table for Gaussian distributed random numbers has not been set up\n");

	simLog(sim,1," %i species allocated\n",mols->maxspecies-1);
	simLog(sim,2," %i species defined:\n",mols->nspecies-1);
	maxstep=-1;
	for(i=1;i<nspecies;i++) {
		simLog(sim,2," %s:\n",mols->spname[i]);
		simLog(sim,1,"  states used:");
		sum=0;
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
			if(mols->exist[i][ms]) {
				sum++;
				simLog(sim,1," %s",molms2string(ms,string)); }
		if(!sum) simLog(sim,1," none");
		simLog(sim,1,"\n");

		same=1;
		for(ms=(enum MolecState)(0);ms<MSMAX && same;ms=(enum MolecState)(ms+1)) {
			if(mols->difc[i][ms]!=mols->difc[i][MSsoln]) same=0;
			if(mols->difm[i][ms] && !mols->difm[i][MSsoln]) same=0;
			if(!mols->difm[i][ms] && mols->difm[i][MSsoln]) same=0;
			if(mols->drift[i][ms] && !mols->drift[i][MSsoln]) same=0;
			if(!mols->drift[i][ms] && mols->drift[i][MSsoln]) same=0;
			if(mols->listlookup[i][ms]!=mols->listlookup[i][MSsoln]) same=0; }
		if(same) {
			if(mols->difstep[i][MSsoln]>maxstep) maxstep=mols->difstep[i][MSsoln];
			simLog(sim,2,"  all states: difc=%g, rms step=%g",mols->difc[i][MSsoln],mols->difstep[i][MSsoln]);
			if(mols->difm[i][MSsoln]) simLog(sim,2," (anisotropic)");
			if(mols->drift[i][MSsoln]) simLog(sim,2," (drift)");
			if(mols->listname) simLog(sim,2,", list=%s",mols->listname[mols->listlookup[i][MSsoln]]);
			simLog(sim,2,", number=%i\n",molcount(sim,i,NULL,MSall,NULL,-1)); }
		else {
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
				if(mols->exist[i][ms]) {
					if(mols->difstep[i][ms]>maxstep) maxstep=mols->difstep[i][ms];
					simLog(sim,2,"  %s: difc=%g, rms step=%g",molms2string(ms,string),mols->difc[i][ms],mols->difstep[i][ms]);
					if(mols->difm[i][ms]) simLog(sim,2," (anisotropic)");
					if(mols->drift[i][ms]) simLog(sim,2," (drift)");
					if(mols->listname) simLog(sim,2,", list=%s",mols->listname[mols->listlookup[i][ms]]);
					simLog(sim,2,", number=%i\n",molcount(sim,i,NULL,ms,NULL,-1)); }}

		if(mols->surfdrift && mols->surfdrift[i]) {
			simLog(sim,2,"  surface drift:\n");
			for(ms=(enum MolecState)1;ms<MSMAX;ms=(enum MolecState)(ms+1))
				if(mols->surfdrift[i][ms] && sim->srfss) {
					simLog(sim,2,"   %s:",molms2string(ms,string));
					for(s=0;s<sim->srfss->nsrf;s++)
						if(mols->surfdrift[i][ms][s])
							for(ps=PanelShape(0);ps<PSMAX;ps=PanelShape(ps+1))
								if(mols->surfdrift[i][ms][s][ps]) {
									simLog(sim,2," %s,%s",sim->srfss->snames[s],surfps2string(ps,string));
									if(sim->dim==2) simLog(sim,2,"=%g",mols->surfdrift[i][ms][s][ps][0]);
									else simLog(sim,2,"=(%g,%g)",mols->surfdrift[i][ms][s][ps][0],mols->surfdrift[i][ms][s][ps][1]); }
					simLog(sim,2,"\n"); }}

		if(sim->graphss) {
			same=1;
			for(ms=(enum MolecState)(0);ms<MSMAX && same;ms=(enum MolecState)(ms+1)) {
				if(mols->display[i][ms]!=mols->display[i][MSsoln]) same=0;
				if(mols->color[i][ms][0]!=mols->color[i][MSsoln][0]) same=0;
				if(mols->color[i][ms][1]!=mols->color[i][MSsoln][1]) same=0;
				if(mols->color[i][ms][2]!=mols->color[i][MSsoln][2]) same=0; }
			if(same) {
				simLog(sim,2,"  all states:");
				if(mols->display[i][MSsoln])
					simLog(sim,2," color= %g,%g,%g, size=%g\n",mols->color[i][MSsoln][0],mols->color[i][MSsoln][1],mols->color[i][MSsoln][2],mols->display[i][MSsoln]);
				else simLog(sim,2," not displayed to graphics\n"); }
			else {
				for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
					if(mols->exist[i][ms]) {
						simLog(sim,2,"  %s:",molms2string(ms,string));
						if(mols->display[i][ms])
							simLog(sim,2," color= %g,%g,%g, display size= %g\n",mols->color[i][ms][0],mols->color[i][ms][1],mols->color[i][ms][2],mols->display[i][ms]);
						else simLog(sim,2," not displayed to graphics\n"); }}}}

	if(mols->dead==NULL) simLog(sim,1," No dead list allocated\n");
	simLog(sim,1," Dead list: allocated size=%i, number of molecules=%i",mols->maxd,mols->nd);
	if(mols->topd!=mols->nd) simLog(sim,1,", top value=%i",mols->topd);
	simLog(sim,1,"\n");
	if(mols->maxdlimit>=0) simLog(sim,1,"  limited to %i molecules\n",mols->maxdlimit);

	simLog(sim,2," %i molecule lists:\n",mols->nlist);
	for(ll=0;ll<mols->nlist;ll++) {
		if(mols->live[ll]==NULL) simLog(sim,1,"  list %i is not allocated\n",ll);
		simLog(sim,1,"  %s: type=%s, allocated size=%i, number of molecules=%i",mols->listname[ll],molmlt2string(mols->listtype[ll],string),mols->maxl[ll],mols->nl[ll]);
		if(mols->topl[ll]!=mols->nl[ll] && mols->topl!=0) simLog(sim,1,", top value=%i",mols->topl[ll]);
		if(mols->sortl[ll]!=mols->nl[ll]) simLog(sim,1,", sort value=%i",mols->sortl[ll]);
		simLog(sim,1,"\n");
		simLog(sim,2,"%s%s%s",ll==0?"  ":" ",mols->listname[ll],ll==mols->nlist-1?"\n":","); }

	simLog(sim,1," Diffusion molecule lists:");
	for(ll=0;ll<mols->nlist;ll++)
		if(mols->diffuselist[ll]) simLog(sim,1," %s",mols->listname[ll]);
	simLog(sim,1,"\n");

	simLog(sim,2," Overall spatial resolution:");
	if(maxstep==-1 || mols->condition<SCok) simLog(sim,2," not computed\n");
	else simLog(sim,2," %g\n",maxstep);
	simLog(sim,2,"\n");
	return; }


/* writemols */
void writemols(simptr sim,FILE *fptr) {
	int i,d,ll,dim;
	char **spname,string[STRCHAR];
	enum MolecState ms;
	molssptr mols;
	double val0,val1,val2;

	mols=sim->mols;
	if(!mols) return;
	dim=sim->dim;
	spname=mols->spname;
	fprintf(fptr,"# Molecule parameters\n");

	fprintf(fptr,"max_species %i\n",mols->maxspecies-1);
	for(i=1;i<mols->nspecies;i++) fprintf(fptr,"species %s\n",spname[i]);
	fprintf(fptr,"\n");
	if(sim->mols->maxdlimit>=0)
		fprintf(fptr,"max_mol %i\n",sim->mols->maxdlimit);
	fprintf(fptr,"gauss_table_size %i\n\n",mols->ngausstbl);

	for(ll=0;ll<mols->nlist;ll++)
		if(mols->listtype[ll]==MLTsystem)
			fprintf(fptr,"molecule_lists %s\n",mols->listname[ll]);
	fprintf(fptr,"\n");
	
	for(i=1;i<mols->nspecies;i++) {
		val0=mols->difc[i][0];
		for(ms=(enum MolecState)(1);ms<MSMAX && mols->difc[i][ms]==val0;ms=(enum MolecState)(ms+1));
		if(ms==MSMAX) fprintf(fptr,"difc %s(all) %g\n",spname[i],val0);
		else {
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
				if(mols->difc[i][ms]>0)
					fprintf(fptr,"difc %s(%s) %g\n",spname[i],molms2string(ms,string),mols->difc[i][ms]); }
		
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) {
			if(mols->difm[i][ms]) {
				fprintf(fptr,"difm %s(%s)",spname[i],molms2string(ms,string));
				for(d=0;d<dim*dim;d++)
					fprintf(fptr," %g",mols->difm[i][ms][d]);
				fprintf(fptr,"\n"); }}
		
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) {
			if(mols->drift[i][ms]) {
				fprintf(fptr,"drift %s(%s)",spname[i],molms2string(ms,string));
				for(d=0;d<dim;d++)
					fprintf(fptr," %g",mols->drift[i][ms][d]);
				fprintf(fptr,"\n"); }}
		
		if(mols->nlist) {
			ll=mols->listlookup[i][0];
			for(ms=(enum MolecState)(1);ms<MSMAX && mols->listlookup[i][ms]==ll;ms=(enum MolecState)(ms+1));
			if(ms==MSMAX) fprintf(fptr,"mol_list %s(all) %s\n",spname[i],mols->listname[ll]);
			else {
				for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
					fprintf(fptr,"mol_list %s(%s) %s\n",spname[i],molms2string(ms,string),mols->listname[mols->listlookup[i][ms]]); }}
		
		val0=mols->display[i][0];
		for(ms=(enum MolecState)(1);ms<MSMAX && mols->display[i][ms]==val0;ms=(enum MolecState)(ms+1));
		if(ms==MSMAX) fprintf(fptr,"display_size %s(all) %g\n",spname[i],val0);
		else {
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
				fprintf(fptr,"display_size %s(%s) %g\n",spname[i],molms2string(ms,string),mols->display[i][ms]); }
		
		val0=mols->color[i][0][0];
		val1=mols->color[i][0][1];
		val2=mols->color[i][0][2];
		for(ms=(enum MolecState)(1);ms<MSMAX && mols->color[i][ms][0]==val0 && mols->color[i][ms][1]==val1 && mols->color[i][ms][2]==val2;ms=(enum MolecState)(ms+1));
		if(ms==MSMAX) fprintf(fptr,"color %s(all) %g %g %g\n",spname[i],val0,val1,val2);
		else {
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
				fprintf(fptr,"color %s(%s) %g %g %g\n",spname[i],molms2string(ms,string),mols->color[i][ms][0],mols->color[i][ms][1],mols->color[i][ms][2]); }
		fprintf(fptr,"\n"); }
	return; }


/* writemolecules */
void writemolecules(simptr sim,FILE *fptr) {
	int m,ll;
	char **spname,string[STRCHAR];
	molssptr mols;
	moleculeptr mptr;

	mols=sim->mols;
	if(!mols) return;
	spname=mols->spname;
	fprintf(fptr,"# Individual molecules\n");
	
	for(ll=0;ll<mols->nlist;ll++)
		if(mols->listtype[ll]==MLTsystem)
			for(m=0;m<mols->nl[ll];m++) {
				mptr=mols->live[ll][m];
				if(mptr->ident>0) {
					if(mptr->mstate==MSsoln)
						fprintf(fptr,"mol 1 %s",spname[mptr->ident]);
					else {
						fprintf(fptr,"surface_mol 1 %s(%s) %s",spname[mptr->ident],molms2string(mptr->mstate,string),mptr->pnl->srf->sname);
						fprintf(fptr," %s %s",surfps2string(mptr->pnl->ps,string),mptr->pnl->pname); }
					fprintf(fptr,"%s\n",molpos2string(sim,mptr,string)); }}

	return; }


/* checkmolparams */
int checkmolparams(simptr sim,int *warnptr) {
	int dim,i,nspecies,m,ll,warn,error,sum,same;
	molssptr mols;
	moleculeptr mptr;
	wallptr *wlist;
	char **spname,string[STRCHAR];
	double m2[DIMMAX*DIMMAX],diag;
	enum MolecState ms;

	error=warn=0;
	mols=sim->mols;
	if(!mols) {
		if(warnptr) *warnptr=warn;
		return 0; }
	dim=sim->dim;
	nspecies=mols->nspecies;
	wlist=sim->wlist;
	spname=mols->spname;

	if(mols->condition!=SCok) {
		warn++;
		simLog(sim,7," WARNING: molecule structure %s\n",simsc2string(mols->condition,string)); }

	for(ll=0;ll<mols->nlist;ll++) {				// check molecule list sorting
		for(m=0;m<mols->nl[ll];m++) {
			mptr=mols->live[ll][m];
			if(!mptr) {error++;simLog(sim,10," SMOLDYN BUG: NULL molecule in live list %i at %i\n",ll,m);}
			else if(mptr->list!=mols->listlookup[mptr->ident][mptr->mstate]) {error++;simLog(sim,10," SMOLDYN BUG: molecule list value for species %i (%s) is %i but should be %i\n",mptr->ident,molms2string(mptr->mstate,string),mptr->list,mols->listlookup[mptr->ident][mptr->mstate]);}
			else if(mptr->list!=ll) {warn++;simLog(sim,9," WARNING: mis-sorted molecule in live list %i at %i\n",ll,m);}
			else if(!mptr->ident) {warn++;simLog(sim,5," WARNING: empty molecule in live list %i at %i\n",ll,m);} }
		for(;m<mols->maxl[ll];m++) {
			mptr=mols->live[ll][m];
			if(mptr) {error++;simLog(sim,10," SMOLDYN BUG: misplaced molecule in live list %i at %i\n",ll,m);} }}

	for(m=0;m<mols->topd;m++) {
		mptr=mols->dead[m];
		if(!mptr) {error++;simLog(sim,10," SMOLDYN BUG: NULL molecule in dead list at %i\n",m);}
		else if(mptr->list!=-1) {error++;simLog(sim,10," SMOLDYN BUG: mis-sorted molecule in dead list at %i (species %i, serno %li)\n",m,mptr->ident,mptr->serno);}
		else if(mptr->ident) {error++;simLog(sim,10," SMOLDYN BUG: live molecule in dead list at %i\n",m);} }
	for(;m<mols->nd;m++) {
		mptr=mols->dead[m];
		if(!mptr) {error++;simLog(sim,10," SMOLDYN BUG: NULL molecule in resurrected list at %i\n",m);}
		else if(mptr->list==-1) {error++;simLog(sim,10," SMOLDYN BUG: mis-sorted molecule in resurrected list at %i\n",m);}
		else if(!mptr->ident) {error++;simLog(sim,10," BUG: dead molecule in resurrected list at %i\n",m);} }
	for(;m<mols->maxd;m++) {
		mptr=mols->dead[m];
		if(mptr) {error++;simLog(sim,10," SMOLDYN BUG: misplaced molecule in dead list at %i\n",m);} }

	for(ll=0;ll<mols->nlist;ll++)
		for(m=0;m<mols->nl[ll];m++)	{									// check for molecules outside system
			mptr=mols->live[ll][m];
			if(!posinsystem(sim,mptr->pos)) {
				simLog(sim,5," WARNING: molecule #%li, of type '%s', is outside system volume\n",mptr->serno,spname[mptr->ident]);
				warn++; }}

	for(i=1;i<nspecies;i++)														// check for asymmetric diffusion matrices
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
			if(mols->difm[i][ms]) {
				dotMMD(mols->difm[i][ms],mols->difm[i][ms],m2,dim,dim,dim);
				if(!issymmetricMD(m2,dim)) {
					simLog(sim,5," WARNING: diffusion matrix for molecule %s (%s) is asymmetric\n",spname[i],molms2string(ms,string));
					warn++; }}

	for(i=1;i<nspecies;i++) {													// check for unused molecules
		sum=0;
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) sum+=mols->exist[i][ms];
		if(!sum) {
			simLog(sim,5," WARNING: molecule %s is never used\n",spname[i]);
			warn++; }}

	if(sim->graphss && sim->graphss->graphics>1) {		// check for molecules that may not display
		diag=systemdiagonal(sim);
		for(i=1;i<nspecies;i++) {
			same=1;
			for(ms=(enum MolecState)(0);ms<MSMAX && same;ms=(enum MolecState)(ms+1))
				if(mols->display[i][ms]!=mols->display[i][MSsoln]) same=0;
			for(ms=(enum MolecState)(0);ms<MSMAX && (same==0 || ms==MSsoln);ms=(enum MolecState)(ms+1)) {
				if(mols->display[i][ms]>0.1*diag) {
					simLog(sim,5," WARNING: very large display size for molecule %s (%s)\n",spname[i],same?"all":molms2string(ms,string));
					warn++; }
				if(mols->display[i][ms]<0.001*diag) {
					simLog(sim,5," WARNING: very small display size for molecule %s (%s)\n",spname[i],same?"all":molms2string(ms,string));
					warn++; }}}}

	if(warnptr) *warnptr=warn;
	return error; }


/******************************************************************************/
/********************************* structure set up ***************************/
/******************************************************************************/


/* molenablemols */
int molenablemols(simptr sim,int maxspecies) {
	molssptr mols;
	int er;

	if(sim->mols) {									// check for redundant function call
		if(maxspecies==-1) {
			if(sim->mols->nspecies<sim->mols->maxspecies) return 0; }
		else {
			if(maxspecies==sim->mols->maxspecies) return 0;
			if(maxspecies<sim->mols->maxspecies) return 2; }}

	if(maxspecies<0) maxspecies=sim->mols?sim->mols->maxspecies*2+1:5;	// need to initialize or increase maxspecies
	mols=molssalloc(sim->mols,maxspecies);
	if(!mols) return 1;
	sim->mols=mols;
	mols->sim=sim;
	molsetcondition(sim->mols,SClists,0);
	boxsetcondition(sim->boxs,SClists,0);
	er=rxnexpandmaxspecies(sim,maxspecies+1,mols->maxsitecode+1);
	if(er) return 1;
	er=surfexpandmaxspecies(sim,maxspecies+1);
	if(er) return 1;
	rxnsetcondition(sim,-1,SClists,0);
	surfsetcondition(sim->srfss,SClists,0);
	portsetcondition(sim->portss,SClists,0);
	return 0; }


/* molsetcondition */
void molsetcondition(molssptr mols,enum StructCond cond,int upgrade) {
	if(!mols) return;
	if(upgrade==0 && mols->condition>cond) mols->condition=cond;
	else if(upgrade==1 && mols->condition<cond) mols->condition=cond;
	else if(upgrade==2) mols->condition=cond;
	if(mols->sim && mols->condition<mols->sim->condition) {
		cond=mols->condition;
		simsetcondition(mols->sim,cond==SCinit?SClists:cond,0); }
	return; }

/* addmollist */
int addmollist(simptr sim,const char *nm,enum MolListType mlt) {
	int ll,er;
	molssptr mols;

	if(!sim->mols) {
		er=molenablemols(sim,-1);
		if(er) return -1; }
	mols=sim->mols;
	if(!mols || !nm) return -3;
	if(stringfind(mols->listname,mols->nlist,nm)!=-1) return -2;
	if(mols->nlist==mols->maxlist) {
		er=mollistalloc(mols,mols->maxlist+1,mlt);
		if(er<0) return -1; }
	ll=mols->nlist++;
	mols->listtype[ll] = mlt;
	strcpy(mols->listname[ll],nm);
	boxsetcondition(sim->boxs,SClists,0);
	rxnsetcondition(sim,-1,SClists,0);
	surfsetcondition(sim->srfss,SClists,0);
	portsetcondition(sim->portss,SClists,0);
	return ll; }


/* molsetmaxspecies */
int molsetmaxspecies(simptr sim,int max) {
	return molenablemols(sim,max); }

/* molsetmaxmol */
int molsetmaxmol(simptr sim,int max) {
	int er;

	if(!sim->mols) {
		er=molenablemols(sim,-1);
		if(er) return er; }
	if(max>=0 && max<sim->mols->maxd) return 5;
	sim->mols->maxdlimit=max;
	return 0; }

/* moldifsites */
int moldifsites(simptr sim, char *species, char *site_name){
	int er, sp_indx, site_indx;
	molssptr mols;
	char *cptr, str0[STRCHAR];
	
	mols=sim->mols;
	if(!mols->spdifsites) mols->spdifsites=g_hash_table_new(g_direct_hash,g_direct_equal);
	sp_indx=stringfind(mols->spname,mols->nspecies,species);
	printf("sp_indx=%d\n",sp_indx);
	if(sp_indx<0) return -1;
	site_indx=stringfind(mols->spsites_name[sp_indx],mols->spsites_num[sp_indx],site_name);
	printf("site_indx=%d\n",site_indx);
	if(site_indx<0) return -2;
	g_hash_table_insert(mols->spdifsites,GINT_TO_POINTER(sp_indx),g_slist_append((GSList*)g_hash_table_lookup(mols->spdifsites,GINT_TO_POINTER(sp_indx)),GINT_TO_POINTER(site_indx)));

	return 0;
}

int moladdsites(simptr sim, char *species, char *site_name, int sitecode){
	int er, sp_indx, i, site_indx;
	molssptr mols=sim->mols;
	char **newsp_sites_name;
	int *newsp_sites_binding;
	
	er=molenablemols(sim,-1);
	if(er) return -1; 

	sp_indx=stringfind(mols->spname,mols->nspecies,species);
	if(sp_indx<0) return -2;

	if(mols->spsites_name==NULL){ 
		mols->spsites_name=(char***) calloc(mols->maxspecies, sizeof(char**));					// max_sites starts at 1;	
	}
	if(mols->spsites_binding==NULL){
		mols->spsites_binding=(int**) calloc(mols->maxspecies, sizeof(int*));
	}

	if(sitecode+1>mols->spsites_num[sp_indx]){
		newsp_sites_name=(char**) calloc(sitecode+1,sizeof(char*));
		newsp_sites_binding=(int*) calloc(sitecode+1,sizeof(int));
		for(i=0;i<mols->spsites_num[sp_indx];i++){
			newsp_sites_name[i]=mols->spsites_name[sp_indx][i];
			newsp_sites_binding[i]=mols->spsites_binding[sp_indx][i];
		}
		for(;i<sitecode+1;i++){
			newsp_sites_name[i]=EmptyString();
			if(strchr(site_name,'@')) {
				newsp_sites_binding[i]=1;
				site_name++;
				strcpy(newsp_sites_name[i],site_name);
			}
			else{	
				newsp_sites_binding[i]=0;
				strcpy(newsp_sites_name[i],site_name);
			}
		}
		
		free(mols->spsites_name[sp_indx]);
		mols->spsites_name[sp_indx]=newsp_sites_name;
		free(mols->spsites_binding[sp_indx]);
		mols->spsites_binding[sp_indx]=newsp_sites_binding;

		mols->spsites_num[sp_indx]=sitecode+1;
	}

	if(sitecode> mols->maxsitecode) mols->maxsitecode=sitecode;
	return 0;		
}

/* moladdspecies */
int moladdspecies(simptr sim,const char *nm) {
	molssptr mols;
	int found,er;
	char nm_tmp[strlen(nm)];

	er=molenablemols(sim,-1);
	if(er) return -1;
	mols=sim->mols;
	if(!strcmp(nm,"empty")) return -4;
	if(strchr(nm,'?') || strchr(nm,'*')) return -6;

	if(strstr(nm,"<v>")){
		strcpy(nm_tmp,nm);
	 	strsplit(nm_tmp,"<");
		found=stringfind(mols->spname,mols->nspecies,nm_tmp);
		if(found>=0) return -5;
		mols->volt_dependent[mols->nspecies]=1;
		strncpy(mols->spname[mols->nspecies++],nm_tmp,STRCHAR);
	}
	else{
		found=stringfind(mols->spname,mols->nspecies,nm);
		if(found>=0) return -5;
		strncpy(mols->spname[mols->nspecies++],nm,STRCHAR);
	}

	molsetcondition(mols,SClists,0);
	rxnsetcondition(sim,-1,SClists,0);
	surfsetcondition(sim->srfss,SClists,0);
	return mols->nspecies-1; 
}


/* molsetexpansionflag */
int molsetexpansionflag(simptr sim,int i,int flag) {
	int i2;

	if(!sim->mols) return 2;
	if(i==-1) {
		for(i2=1;i2<sim->mols->nspecies;i2++)
			sim->mols->expand[i2]=flag; }
	else if(i<0 || i>=sim->mols->nspecies) return 3;
	else sim->mols->expand[i]=flag;
	return 0; }


/* molsupdateparams */
int molsupdateparams(molssptr mols,double dt) {
	int i,ll;
	enum MolecState ms;

	for(ll=0;ll<mols->nlist;ll++) mols->diffuselist[ll]=0;		// set diffuselist
	for(i=0;i<mols->nspecies;i++)
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) {
			if(molismobile(mols->sim,i,ms))
				mols->diffuselist[mols->listlookup[i][ms]]=1; }

	for(i=0;i<mols->nspecies;i++)					// calculate difstep
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
			mols->difstep[i][ms]=sqrt(2.0*mols->difc[i][ms]*dt);

	return 0; }


/* molsupdatelists */
int molsupdatelists(simptr sim) {
	int i,ll,m,ndif,nfix,ok,er;
	enum MolecState ms;
	molssptr mols;
	moleculeptr mptr;
	
	mols=sim->mols;

	er=molssetgausstable(sim,-1);				// gaussian lookup table
	if(er) return 1;

	for(i=1;i<mols->nspecies;i++)					// set exist values
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
			mols->exist[i][ms]=0;
	for(m=mols->topd;m<mols->nd;m++) {
		mptr=mols->dead[m];
		mols->exist[mptr->ident][mptr->mstate]=1; }
	for(ll=0;ll<mols->nlist;ll++)
		for(m=0;m<mols->nl[ll];m++) {
			mptr=mols->live[ll][m];
			mols->exist[mptr->ident][mptr->mstate]=1; }
	for(i=1;i<mols->nspecies;i++) {
		for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) {
			if(mols->exist[i][ms]==0 && rxnisprod(sim,i,ms,0)) mols->exist[i][ms]=1;
			if(mols->exist[i][ms]==0 && issurfprod(sim,i,ms)) mols->exist[i][ms]=1; }
		if(mols->exist[i][MSsoln]==0 && rxnisprod(sim,i,MSbsoln,0)) mols->exist[i][MSsoln]=1;
		if(mols->exist[i][MSsoln]==0 && issurfprod(sim,i,MSbsoln)) mols->exist[i][MSsoln]=1; }

	for(ll=0;ll<mols->nlist;ll++)					// create system molecule lists if none yet
		if(mols->listtype[ll]==MLTsystem) ll=mols->nlist+1;
	if(ll==mols->nlist && mols->maxd>0 && mols->nspecies>1) {
		ndif=nfix=0;								// fixed and diffuse lists
		for(i=1;i<mols->nspecies;i++)
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1)) {
				if(molismobile(sim,i,ms)) ndif=1;
				else{ 
					// printf("i=%d, ms=%d\n",i,ms);
					nfix=1;
				} 
		}
		if(ndif) {
			ll=addmollist(sim,"diffuselist",MLTsystem);
			if(ll<0) return 1;
			molsetlistlookup(sim,-7,NULL,MSall,ll); }
		if(nfix) {
			ll=addmollist(sim,"fixedlist",MLTsystem);
			if(ll<0) return 1;
			molsetlistlookup(sim,-8,NULL,MSall,ll); }}

	ok=1;															// set any list lookup values that weren't done yet
	for(i=0;i<mols->nspecies && ok;i++)
		for(ms=(enum MolecState)(0);ms<MSMAX && ok;ms=(enum MolecState)(ms+1))
			if(mols->listlookup[i][ms]<0)
				ok=0;
	if(!ok) {
		ll=stringfind(mols->listname,mols->nlist,"unassignedlist");
		if(ll<0) {
			ll=addmollist(sim,"unassignedlist",MLTsystem);
			if(ll<0) return 1; }
		for(i=0;i<mols->nspecies;i++)
			for(ms=(enum MolecState)(0);ms<MSMAX;ms=(enum MolecState)(ms+1))
				if(mols->listlookup[i][ms]<0)
					molsetlistlookup(sim,i,NULL,ms,ll); }

	for(m=mols->topd;m<mols->nd;m++) {		// set molecule list values for molecules in dead list
		mptr=mols->dead[m];
		mptr->list=mols->listlookup[mptr->ident][mptr->mstate]; }
	


	return 0; }


/* molsupdate */
int molsupdate(simptr sim) {
	int er;
	molssptr mols;

	mols=sim->mols;
	if(mols) {
		if(mols->condition<=SClists) {
			er=molsupdatelists(sim);
			if(er) return er;
			molsetcondition(mols,SCparams,1); }
		if(mols->condition==SCparams) {
			er=molsupdateparams(mols,sim->dt);
			if(er) return er;
			molsetcondition(mols,SCok,1); }}
	return 0; }


/******************************************************************************/
/*********************** adding and removing molecules ************************/
/******************************************************************************/

/* molkill */
void molkill(simptr sim,moleculeptr mptr,int ll,int m) {
	int s,dim,d,*sortl,s1;
	moleculeptr mptr_bind;

	if(mptr->dif_molec) return;

	dim=sim->dim;
	sortl=sim->mols->sortl;	

	mptr->ident=0;
	mptr->mstate=MSsoln;
	mptr->list=-1;
	mptr->pos=mptr->pos_tmp;
	for(d=0;d<sim->dim;d++) {
		mptr->posoffset[d]=0;
		mptr->posx[d]=0;
		mptr->prev_pos[d]=0;
		mptr->pos[d]=0;
	}
	mptr->pnl=NULL;
	mptr->s_index=NULL;
	mptr->to=NULL;
	mptr->from=NULL;
	mptr->tot_sunit=0;

	// mptr->react_permit=0;
	mptr->complex_id=-1;
		
	for(s=0;s<sim->mols->spsites_num[mptr->ident];s++){
		if(mptr->sites[s]->bind){
			mptr_bind=mptr->sites[s]->bind;
			if(mptr==mptr_bind->dif_molec)
				mptr_bind->dif_molec=NULL;
			/*
			if(mptr_bind->list!=-1)
				molkill(sim,mptr_bind,mptr_bind->list,mptr_bind->m);
			*/
			for(s1=0;s1<sim->mols->spsites_num[mptr_bind->ident];s1++){
				if(mptr_bind->sites[s1]->bind==mptr){
					mptr_bind->sites[s1]->bind=NULL;
					mptr_bind->sites[s]->value[0]=0;
				}
			}	
			mptr->sites[s]->bind=NULL;
		}
	}

  if(ll<0);
	else if(m<0) sim->mols->sortl[ll]=0;
	else if(m<sim->mols->sortl[ll]) sim->mols->sortl[ll]=m;
	return; }


/* getnextmol */
moleculeptr getnextmol(molssptr mols) {
	//moleculeptr mptr;
	moleculeptr *mptr_cplx, mptr;
	int er,nmol;
	int sunit=1;
	int s;

	if(mols->topd==0) {
		if(mols->maxdlimit>=0 && mols->maxd>=mols->maxdlimit) return NULL;
		//nmol=mols->maxd+1;
		nmol=mols->maxd+sunit;
		if(mols->maxdlimit>=0 && mols->maxd+nmol>mols->maxdlimit)
			nmol=mols->maxdlimit-mols->maxd;
		er=molexpandlist(mols,mols->sim->dim,-1,nmol,nmol);
		if(er) return NULL; }

	mptr_cplx=(moleculeptr*) calloc(sunit, sizeof(moleculeptr));
	for(s=0;s<sunit;s++){
		mptr_cplx[s]=mols->dead[--mols->topd];
		mptr_cplx[s]->serno=mols->serno++;	
		mptr_cplx[s]->s_index=s;
	}
	//printf("not cplx nmol: %d\n", nmol);
	mptr=mptr_cplx[0];
	return mptr; }


/* getnextmol_cplx */
moleculeptr getnextmol_cplx(molssptr mols, int sunit, int ident) {
	moleculeptr mptr_tmp, mptr, *mptr_cplx;
	int er,nmol;
	int s, s_from, s_to;
	double theta_init, phi_init;
	int d,i,k,spsites_num,site_indx;
	char **spsites_name, *site_bind, sitename_tmp[STRCHAR];
	complexptr *complexlist_old, *complexlist_new, complex_tmp=NULL;
	
	if(sunit>1){
		complexlist_old=mols->complexlist;
		complexlist_new=(complexptr*) calloc(mols->max_complex, sizeof(complexptr));
		if(complexlist_new!=NULL){
			if(complexlist_old){
				for(i=0;i<mols->ncomplex;i++) complexlist_new[i]=complexlist_old[i];
				free(complexlist_old);	}
			mols->complexlist=complexlist_new;
		}
		mols->ncomplex++;
	}

	if(mols->topd<sunit) {
		if(mols->maxdlimit>=0 && mols->maxd>=mols->maxdlimit) return NULL;
		//nmol=mols->maxd+1;
		nmol=mols->maxd+sunit;
		if(mols->maxdlimit>=0 && mols->maxd+nmol>mols->maxdlimit)
			nmol=mols->maxdlimit-mols->maxd;
		er=molexpandlist(mols,mols->sim->dim,-1,nmol,nmol);
		//er=molexpandlist(mols,mols->sim->dim,-1,sunit,sunit);		// the second last input is for nspaces, maxnew=nspaces>0?maxold+nspaces:2*maxold+1;
		if(er) return NULL; 
	}

	// rand() between 0 and RAND_MAX

	mptr_cplx=&(mols->dead[mols->topd-1]);
	theta_init=(double)rand()/RAND_MAX*PI/2;
	phi_init=(double)rand()/RAND_MAX*PI/2;

	spsites_num=mols->spsites_num[ident];
	for(s=0;s<sunit;s++){
		mptr_tmp=mols->dead[--mols->topd];
		mptr_tmp->serno=mols->serno++;	
		mptr_tmp->ident=ident;
		mptr_tmp->s_index=sunit-1-s;
		mptr_tmp->tot_sunit=sunit;
		mptr_tmp->theta_init=theta_init;
		mptr_tmp->phi_init=phi_init;
		mptr_tmp->ident=ident;
		mptr_tmp->bind_id=ident;
		if(spsites_num>0){
			mptr_tmp->sites=(siteptr*) calloc(spsites_num,sizeof(siteptr));
			for(k=0;k<spsites_num;k++){
				mptr_tmp->sites[k]=(siteptr) calloc(1,sizeof(sitestruct));
				mptr_tmp->sites[k]->value=(int*) malloc(sizeof(int));
				mptr_tmp->sites[k]->value[0]=0;
				mptr_tmp->sites[k]->value_tmp=mptr_tmp->sites[k]->value;
				mptr_tmp->sites[k]->time=-1;
		}}
		if(sunit>1) {
			if(mptr_tmp->s_index==0){
				complex_tmp=complexalloc(mols->sim, mptr_cplx[s]);
				if(complex_tmp){
					complex_tmp->serno=mols->ncomplex;
					mols->complexlist[mols->ncomplex-1]=complex_tmp;			// the first element will start from 1 if offsetting it
					complex_tmp->layer=1; // 2;
			}}
			mptr_tmp->complex_id=mols->ncomplex-1;
		}	
		if(mols->volt_dependent[ident]==1){
			mptr_tmp->vchannel=(vchnlptr) calloc(1,sizeof(vchnlstruct));	
			mptr_tmp->vchannel->voltage_file=fopen(mols->sim->vfile,"r");
			mptr_tmp->vchannel->vtime=-1;
			mptr_tmp->vchannel->vtime_n=-1;
			mptr_tmp->vchannel->molec_gen=-1;
			//mptr_tmp->vchannel->voltage_file=fopen("/home/neuro/Documents/from_axon/dat_files/soma_v_10Hz.txt","r");
			//mptr_tmp->vchannel->voltage=0;
			mptr_tmp->vchannel->mptr=mptr_tmp;
		}
	}

	for(s=0;s<sunit;s++){	
		mptr_tmp=mptr_cplx[-s];
		if(sunit>1){
			// adjusted for index and s differences
			s_from=(s+1)>sunit?(s+1)%sunit:(s+1);			
			s_to=(s-1)<0?(s-1+sunit):(s-1);
			mptr_tmp->from=mptr_cplx[-s_from];
			mptr_tmp->to=mptr_cplx[-s_to];
		}
		if(mols->spsites_name)
			spsites_name=mols->spsites_name[ident];

		for(k=0;k<spsites_num;k++){			// default unbound, sites_state=1;
			if(mols->spsites_binding[ident][k]==1)
				mptr_tmp->sites[k]->site_type=1;
			else mptr_tmp->sites[k]->site_type=0;
			mptr_tmp->sites[k]->bind=NULL;
			mptr_tmp->sites[k]->indx=k;
			// mptr_cplx[s]->sites[k]->value=(int*) malloc(sizeof(int));
			// mptr_cplx[s]->sites[k]->value[0]=0;

			if(strstr(spsites_name[k],"n.")){
				// printf("molec.c line2569, ident=%d, k=%d, spsites_name[k]:%s\n", ident, k, spsites_name[k]);
				strcpy(sitename_tmp,spsites_name[k]);
				site_bind=strsplit(sitename_tmp,".");							
				site_indx=stringfind(mols->spsites_name[mptr_tmp->to->ident],mols->spsites_num[mptr_tmp->to->ident],site_bind);
				mptr_tmp->sites[k]->value=mptr_tmp->to->sites[site_indx]->value;
		}
	}}
	mptr=mptr_cplx[-sunit+1];
	return mptr; 
}

/* newestmol */
moleculeptr newestmol(molssptr mols) {
	return mols->dead[mols->topd-1]; }

/* addmol */
int addmol(simptr sim,int nmol,int ident,double *poslo,double *poshi,int sort) {
	int m,d;
	moleculeptr mptr;

	for(m=0;m<nmol;m++) {
		// printf("addmol!\n");
		mptr=getnextmol(sim->mols);
		if(!mptr) return 3;
		mptr->ident=ident;
		mptr->mstate=MSsoln;
		mptr->list=sim->mols->listlookup[ident][MSsoln];
		if(poslo==poshi)
			for(d=0;d<sim->dim;d++)
				mptr->posx[d]=mptr->pos[d]=poslo[d];
		else
			for(d=0;d<sim->dim;d++)
				mptr->posx[d]=mptr->pos[d]=unirandOOD(poslo[d],poshi[d]);
		if(sim->boxs && sim->boxs->nbox)
			mptr->box=pos2box(sim,mptr->pos);
		else mptr->box=NULL; }
	molsetexist(sim,ident,MSsoln,1);
	if(sort)
		if(molsort(sim,1)) return 1;
	return 0; }

/* addmol_cplx */
int addmol_cplx(simptr sim,int num_mol,int sunit,int bind_num,int *ident,int *sites,double *poslo,double *poshi,int sort, int ident_free, int sites_free) {
	int m,d;
	moleculeptr mptr, mptr_tmp, mptr_bound;

	if(sunit>1) sim->mols->max_complex+=num_mol;
	for(m=0;m<num_mol;m++) {
		mptr=getnextmol_cplx(sim->mols,sunit,ident[0]);
		if(!mptr) return 3;
		mptr->mstate=MSsoln;
		mptr->list=sim->mols->listlookup[ident[0]][MSsoln];
		if(bind_num==1 && sites!=NULL) mptr->sites[0]->value[0]=1;
		if(poslo==poshi) {
			for(d=0;d<sim->dim;d++)
				mptr->posx[d]=mptr->prev_pos[d]=mptr->pos[d]=poslo[d];
			complex_pos_init(sim,mptr,NULL); }
		else {
			for(d=0;d<sim->dim;d++)
				mptr->posx[d]=mptr->prev_pos[d]=mptr->pos[d]=unirandOOD(poslo[d],poshi[d]);
			complex_pos_init(sim,mptr,NULL);	 }

		if(sim->boxs && sim->boxs->nbox) mptr->box=pos2box(sim,mptr->pos);
		else mptr->box=NULL;
		mptr_tmp=mptr;
		while(mptr_tmp->to!=NULL && mptr_tmp->s_index < mptr_tmp->to->s_index) {
			mptr_tmp->to->mstate=MSsoln;
			mptr_tmp->to->list=sim->mols->listlookup[ident[0]][MSsoln];
			if(sim->boxs && sim->boxs->nbox) mptr_tmp->box=pos2box(sim,mptr_tmp->pos);
			else mptr_tmp->box=NULL;
			// complex_pos_init(sim,mptr_tmp,NULL);	
			mptr_tmp=mptr_tmp->to;
		}
		if(bind_num>1) {
			if(!sites) return 4;
			mptr_bound=getnextmol_cplx(sim->mols,1,ident[1]);
			if(!mptr_bound) return 3;
			mptr_bound->mstate=MSsoln;
			mptr_bound->list=sim->mols->listlookup[ident[1]][MSsoln];

			// mptr_bound->pos_tmp=mptr_bound->pos;
			mptr_bound->pos=mptr->pos;
			for(d=0;d<sim->dim;d++)	mptr_bound->posx[d]=mptr_bound->prev_pos[d]=mptr_bound->pos[d];
			mptr->sites[0]->value[0]=1;
			mptr_bound->sites[1]->value[0]=1;
			mptr->sites[0]->bind=mptr_bound;
			mptr_bound->sites[1]->bind=mptr;
			mptr_bound->dif_molec=mptr;
			mptr_bound->sites_val=molecsites_state(sim->mols,mptr_bound);
			if(sim->events){
				fprintf(sim->events,"rxn_time=%f start bound_ident=%s bound_serno=%d bound_state=%d  pos[0]=%f pos[1]=%f pos[2]=%f complex_id=%d\n", sim->time,sim->mols->spname[mptr_bound->ident], mptr_bound->serno, mptr_bound->sites_val, mptr_bound->pos[0], mptr_bound->pos[1], mptr_bound->pos[2], mptr_bound->complex_id);}
		}		
		mptr->sites_val=molecsites_state(sim->mols,mptr);
		if(sim->events)	{
			fprintf(sim->events,"rxn_time=%f start ident=%s serno=%d sites_val=%d pos[0]=%f pos[1]=%f pos[2]=%f complex_id=%d\n", sim->time,sim->mols->spname[mptr->ident],mptr->serno,mptr->sites_val,mptr->pos[0],mptr->pos[1],mptr->pos[2], mptr->complex_id);}
	
	}	
	molsetexist(sim,ident[0],MSsoln,1);
	if(bind_num>1) molsetexist(sim,ident[1],MSsoln,1);
	if(sites_free) if(sites) free(sites);	
	if(ident_free) if(ident) free(ident);
	if(sort)
		if(molsort(sim,1)) return 1;
	return 0; }


/* addsurfmol */
int addsurfmol(simptr sim,int nmol,int ident,enum MolecState ms,double *pos,panelptr pnl,int surface,enum PanelShape ps,char *pname){
	int dim,m,d,totpanel,panel;
	moleculeptr mptr;
	int s,slo,shi,pslo,pshi,p,plo,phi,pindex;
	double *areatable,area,mpos[DIMMAX],totarea;
	panelptr *paneltable;
	surfaceptr srf;

	dim=sim->dim;

	if(pnl || (surface>=0 && ps!=PSall && pname && strcmp(pname,"all"))) {			// add to a specific panel
		if(!pnl) {
			srf=sim->srfss->srflist[surface];
			panel=stringfind(srf->pname[ps],srf->npanel[ps],pname);
			if(panel<0) return 2;
			pnl=srf->panels[ps][panel]; }
		for(m=0;m<nmol;m++) {
			mptr=getnextmol(sim->mols);
			if(!mptr) return 3;
			mptr->ident=ident;
			mptr->mstate=ms;
			mptr->list=sim->mols->listlookup[ident][ms];
			mptr->pnl=pnl;
			if(pos)
				for(d=0;d<dim;d++) mpos[d]=pos[d];
			else
				panelrandpos(pnl,mpos,dim);
			if(ms==MSfront) fixpt2panel(mpos,pnl,dim,PFfront,0);
			else if(ms==MSback) fixpt2panel(mpos,pnl,dim,PFback,0);
			if(mptr->s_index==0){
				for(d=0;d<dim;d++) mptr->pos[d]=mptr->prev_pos[d]=mptr->posx[d]=mpos[d];
				// complex_pos(sim->dim,mptr);
			}	
			if(sim->boxs && sim->boxs->nbox) mptr->box=pos2box(sim,mpos);
			else mptr->box=NULL; }}

	else {
		totarea=surfacearea2(sim,surface,ps,pname,&totpanel);		// create area lookup tables
		if(totpanel<1) return 2;
		areatable=(double*) calloc(totpanel,sizeof(double));
		if(!areatable) return 1;
		paneltable=(panelptr*) calloc(totpanel,sizeof(panelptr));
		if(!paneltable) {free(areatable);return 1; }

		slo=(surface>=0)?surface:0;
		shi=(surface>=0)?surface+1:sim->srfss->nsrf;
		pslo=(ps!=PSall)?ps:0;
		pshi=(ps!=PSall)?ps+1:PSMAX;

		pindex=0;																						// fill in area lookup tables
		area=0;
		for(s=slo;s<shi;s++)
			for(ps=PanelShape(pslo);ps<pshi;ps=PanelShape(ps + 1)) {
				srf=sim->srfss->srflist[s];
				if(!pname || !strcmp(pname,"all")) {plo=0;phi=srf->npanel[ps];}
				else if((panel=stringfind(srf->pname[ps],srf->npanel[ps],pname))<0) plo=phi=0;
				else {plo=panel;phi=panel+1;}
				for(p=plo;p<phi;p++) {
					area+=surfacearea2(sim,s,ps,srf->pname[ps][p],NULL);
					areatable[pindex]=area;
					paneltable[pindex]=srf->panels[ps][p];
					pindex++; }}

		for(m=0;m<nmol;m++) {															// place molecules
			mptr=getnextmol(sim->mols);
			if(!mptr) {free(paneltable);free(areatable);return 3;}
			mptr->ident=ident;
			mptr->mstate=ms;
			mptr->list=sim->mols->listlookup[ident][ms];
			pindex=intrandpD(totpanel,areatable);
			pnl=paneltable[pindex];
			mptr->pnl=pnl;
			panelrandpos(pnl,mpos,dim);
			if(ms==MSfront) fixpt2panel(mpos,pnl,dim,PFfront,0);
			else if(ms==MSback) fixpt2panel(mpos,pnl,dim,PFback,0);
			if(mptr->s_index==0){
				for(d=0;d<dim;d++) mptr->pos[d]=mptr->prev_pos[d]=mptr->posx[d]=mpos[d];
				// complex_pos(sim->dim,mptr);
			}				
			if(sim->boxs && sim->boxs->nbox) mptr->box=pos2box(sim,mpos);
			else mptr->box=NULL; }

		free(paneltable);
		free(areatable); }

	molsetexist(sim,ident,ms,1);
	return 0; 
}

/* addsurfmol_cplx*/
int addsurfmol_cplx(simptr sim,int nmol,int sunit,int bind_num,int *ident,int *sites,enum MolecState ms,double *pos,panelptr pnl,int surface,enum PanelShape ps,char *pname) {
	int dim,m,d,totpanel,panel;
	moleculeptr mptr;
	int s,slo,shi,pslo,pshi,p,plo,phi,pindex;
	double *areatable,area,mpos[DIMMAX],totarea;
	panelptr *paneltable;
	surfaceptr srf;

	dim=sim->dim;

	if(pnl || (surface>=0 && ps!=PSall && pname && strcmp(pname,"all"))) {			// add to a specific panel
		if(!pnl) {
			srf=sim->srfss->srflist[surface];
			panel=stringfind(srf->pname[ps],srf->npanel[ps],pname);
			if(panel<0) return 2;
			pnl=srf->panels[ps][panel]; }
		for(m=0;m<nmol;m++) {
			mptr=getnextmol_cplx(sim->mols,sunit,ident[0]);
			if(!mptr) return 3;
			mptr->mstate=ms;
			mptr->list=sim->mols->listlookup[ident[0]][ms];
			if(bind_num==1 && sites!=NULL) mptr->sites[0]->value[0]=1;
			mptr->pnl=pnl;
			if(pos)
				for(d=0;d<dim;d++) mpos[d]=pos[d];
			else
				panelrandpos(pnl,mpos,dim);
			if(ms==MSfront) fixpt2panel(mpos,pnl,dim,PFfront,0);
			else if(ms==MSback) fixpt2panel(mpos,pnl,dim,PFback,0);
			if(mptr->s_index==0){
				for(d=0;d<dim;d++) mptr->pos[d]=mptr->prev_pos[d]=mptr->posx[d]=mpos[d];
				complex_pos_init(sim,mptr,NULL);
			}
			if(sim->boxs && sim->boxs->nbox) mptr->box=pos2box(sim,mpos);
			else mptr->box=NULL; }}

	else {
		totarea=surfacearea2(sim,surface,ps,pname,&totpanel);		// create area lookup tables
		if(totpanel<1) return 2;
		areatable=(double*) calloc(totpanel,sizeof(double));
		if(!areatable) return 1;
		paneltable=(panelptr*) calloc(totpanel,sizeof(panelptr));
		if(!paneltable) {free(areatable);return 1; }

		slo=(surface>=0)?surface:0;
		shi=(surface>=0)?surface+1:sim->srfss->nsrf;
		pslo=(ps!=PSall)?ps:0;
		pshi=(ps!=PSall)?ps+1:PSMAX;

		pindex=0;																						// fill in area lookup tables
		area=0;
		for(s=slo;s<shi;s++)
			for(ps=PanelShape(pslo);ps<pshi;ps=PanelShape(ps + 1)) {
				srf=sim->srfss->srflist[s];
				if(!pname || !strcmp(pname,"all")) {plo=0;phi=srf->npanel[ps];}
				else if((panel=stringfind(srf->pname[ps],srf->npanel[ps],pname))<0) plo=phi=0;
				else {plo=panel;phi=panel+1;}
				for(p=plo;p<phi;p++) {
					area+=surfacearea2(sim,s,ps,srf->pname[ps][p],NULL);
					areatable[pindex]=area;
					paneltable[pindex]=srf->panels[ps][p];
					pindex++; }}

		for(m=0;m<nmol;m++) {															// place molecules
			mptr=getnextmol_cplx(sim->mols,sunit,ident[0]);
			if(!mptr) {free(paneltable);free(areatable);return 3;}
			mptr->mstate=ms;
			mptr->list=sim->mols->listlookup[ident[0]][ms];
			pindex=intrandpD(totpanel,areatable);
			pnl=paneltable[pindex];
			mptr->pnl=pnl;
			panelrandpos(pnl,mpos,dim);
			if(ms==MSfront) fixpt2panel(mpos,pnl,dim,PFfront,0);
			else if(ms==MSback) fixpt2panel(mpos,pnl,dim,PFback,0);
			if(mptr->s_index==0){
				for(d=0;d<dim;d++) mptr->pos[d]=mptr->prev_pos[d]=mptr->posx[d]=mpos[d];
				complex_pos_init(sim,mptr,NULL);
			}
			if(sim->boxs && sim->boxs->nbox) mptr->box=pos2box(sim,mpos);
			else mptr->box=NULL;
			if(sim->events){
				fprintf(sim->events,"rxn_time=%f start ident=%s serno=%d sites_val=%d pos[0]=%f pos[1]=%f pos[2]=%f complex_id=%d\n", sim->time,sim->mols->spname[mptr->ident],mptr->serno,mptr->sites_val,mptr->pos[0],mptr->pos[1],mptr->pos[2],mptr->complex_id);		
			}
		 }

		free(paneltable);
		free(areatable); }

	molsetexist(sim,ident[0],ms,1);
	if(bind_num>1);

	//if(ms) free(ms);
	if(sites) free(sites);
	if(ident) free(ident);

	return 0; }

/* addcompartmol */
int addcompartmol(simptr sim,int nmol,int ident,compartptr cmpt) {
	int d,dim,m,er;
	moleculeptr mptr;

	if(cmpt->npts==0 && cmpt->ncmptl==0) return 2;
	dim=sim->dim;
	for(m=0;m<nmol;m++) {
		mptr=getnextmol(sim->mols);
		if(!mptr) return 3;
		mptr->ident=ident;
		mptr->mstate=MSsoln;
		mptr->list=sim->mols->listlookup[ident][MSsoln];
		er=compartrandpos(sim,mptr->pos,cmpt);
		if(er) return 2;
		for(d=0;d<dim;d++) mptr->posx[d]=mptr->prev_pos[d]=mptr->pos[d];
		if(sim->boxs && sim->boxs->nbox) mptr->box=pos2box(sim,mptr->pos);
		else mptr->box=NULL; }
	molsetexist(sim,ident,MSsoln,1);
	return 0; }

/* addcompartmol_cplx*/
int addcompartmol_cplx(simptr sim,int num_mol,int sunit,int bind_num,int *ident,int *sites_val,compartptr cmpt){
	int d,dim,m,k,er;
	moleculeptr mptr, mptr_tmp, mptr_bound;
	int sunit_i;

	if(cmpt->npts==0 && cmpt->ncmptl==0) return 2;
	dim=sim->dim;
	if(sunit>1) sim->mols->max_complex+=num_mol;
	for(m=0;m<num_mol;m++) {
		mptr=getnextmol_cplx(sim->mols,sunit,ident[0]);
		if(!mptr) return 3;
		mptr->mstate=MSsoln;
		mptr->list=sim->mols->listlookup[ident[0]][MSsoln];
		if(bind_num==1 && sites_val!=NULL) {
			mptr->sites[sites_val[0]]->value[0]=1;		// phosphorylation or actin binding set to 1	
			if(sites_val[1]>0)
				mptr->bind_id=sites_val[1];
			if(mptr->sites[sites_val[0]]->site_type==0) {
				mptr->dif_site=sites_val[0];
				if(sites_val[1]>0){ 
					mptr->sites[sites_val[0]]->bind=mptr;	// bind to itself;
					sim->mols->complexlist[mptr->complex_id]->dif_molec=mptr;
					sim->mols->complexlist[mptr->complex_id]->dif_bind_site=sites_val[0];
					sim->mols->complexlist[mptr->complex_id]->dif_bind=mptr;
				}
			}
		}
		mptr_tmp=mptr;
		sunit_i=1;
		while(mptr_tmp->to!=NULL && mptr_tmp->s_index!=mptr_tmp->to->s_index && sunit_i<mptr->tot_sunit){
			mptr_tmp->to->mstate=MSsoln;
			mptr_tmp->to->list=sim->mols->listlookup[ident[0]][MSsoln];
			sunit_i++;
			if(bind_num==1 && sites_val!=NULL) {
				if(sites_val[1]==-1){
					mptr_tmp->to->sites[sites_val[0]]->value[0]=1; 
				}
			}
			mptr_tmp=mptr_tmp->to;
		}
		er=compartrandpos(sim,mptr->pos,cmpt);
		if(er) return 2;
		for(d=0;d<dim;d++) mptr->posx[d]=mptr->prev_pos[d]=mptr->pos[d];
		if(sim->boxs && sim->boxs->nbox) mptr->box=pos2box(sim,mptr->pos);
		else mptr->box=NULL; 
		
		complex_pos_init(sim,mptr,cmpt);	
		mptr_tmp=mptr;
		sunit_i=1;
		while(mptr_tmp->to!=NULL && mptr_tmp->s_index!=mptr_tmp->to->s_index && sunit_i<=mptr->tot_sunit) {		
			for(d=0;d<dim;d++) mptr_tmp->to->posx[d]=mptr_tmp->to->prev_pos[d]=mptr_tmp->to->pos[d];
			if(sim->boxs && sim->boxs->nbox) mptr_tmp->to->box=pos2box(sim,mptr_tmp->to->pos);
			else mptr_tmp->to->box=NULL;
			mptr_tmp=mptr_tmp->to;
			sunit_i++;
		}	
		if(bind_num>1){
			if(!sites_val) return 4;											// no binding sites information
			mptr_bound=getnextmol_cplx(sim->mols,1,ident[1]);
			if(!mptr_bound) return 3;
			mptr_bound->mstate=MSsoln;
			mptr_bound->list=sim->mols->listlookup[ident[1]][MSsoln];
			mptr_bound->pos=mptr->pos;
			for(d=0;d<dim;d++)  mptr_bound->pos_tmp[d]=mptr_bound->pos[d];
			for(d=0;d<dim;d++) mptr_bound->posx[d]=mptr_bound->prev_pos[d]=mptr_bound->pos[d];
			if(sim->boxs && sim->boxs->nbox) mptr_bound->box=mptr->box;	
			else mptr_bound->box=NULL;
		
			mptr->sites[sites_val[0]]->value[0]=1;
			mptr_bound->sites[sites_val[1]]->value[0]=1;
			mptr->sites[sites_val[0]]->bind=mptr_bound;
			mptr_bound->sites[sites_val[1]]->bind=mptr;
			mptr_bound->dif_molec=mptr;
			mptr_bound->sites_val=molecsites_state(sim->mols,mptr_bound);	
			if(sim->events){
				fprintf(sim->events,"rxn_time=%f start bound_ident=%s bound_serno=%d sites_val=%d pos[0]=%f pos[1]=%f pos[2]=%f complex_id=%d\n", sim->time,sim->mols->spname[mptr_bound->ident],mptr_bound->serno,mptr_bound->sites_val,mptr_bound->pos[0],mptr_bound->pos[1],mptr_bound->pos[2], mptr_bound->complex_id);
			}
		}
		mptr->sites_val=molecsites_state(sim->mols,mptr);
		if(sim->events){
			fprintf(sim->events,"rxn_time=%f start ident=%s serno=%d sites_val=%d pos[0]=%f pos[1]=%f pos[2]=%f complex_id=%d\n", sim->time,sim->mols->spname[mptr->ident],mptr->serno,mptr->sites_val,mptr->pos[0],mptr->pos[1],mptr->pos[2],mptr->complex_id);		
		}
	}	
	molsetexist(sim,ident[0],MSsoln,1);
	if(bind_num>1)	molsetexist(sim,ident[1],MSsoln,1);
	if(ident) free(ident);
	if(sites_val) free(sites_val); 	// alllocated in readmolname_cplx

	return 0; }

/******************************************************************************/
/*************************** core simulation functions ************************/
/******************************************************************************/


/* molsort */
int molsort(simptr sim,int onlydead2live) {
	molssptr mols;
	int nlist,*maxl,*nl,*topl,*sortl,m,ll,ll2;
	moleculeptr *dead,**live,*mlist,mptr;
	enum MolListType *listtype;
	boxptr bptr;
	int *list_new;
	int list_maxlen=0;

	if(!sim->mols) return 0;
	mols=sim->mols;
	dead=mols->dead;
	nlist=mols->nlist;
	listtype=mols->listtype;
	live=mols->live;
	maxl=mols->maxl;
	nl=mols->nl;
	topl=mols->topl;
	sortl=mols->sortl;

	if(!onlydead2live) {
   		for(ll=0;ll<nlist;ll++)								// reset topl indicies  // onlydead2live=0
    	 	topl[ll]=nl[ll];

    	for(ll=0;ll<nlist;ll++) {							// sort live lists
      		mlist=live[ll];
	  		// printf("molsort, time=%f, sortl[%d]=%d, topl[%d]=%d\n", sim->time, ll, sortl[ll], ll, topl[ll]);
      		for(m=sortl[ll];m<topl[ll];m++) {
        		if(mlist[m]->list!=ll) {
          			mptr=mlist[m];
					mptr->m=m;
          			if(mptr->list==-1) {							// move to dead list
						// printf("molsort() line 2751 m=%d, serno=%d\n", m, mptr->serno);
						if(mptr->box) boxremovemol(mptr,ll);
            			dead[mols->nd++]=dead[mols->topd];
            			dead[mols->topd++]=mptr;
						mptr->m=mols->topd;
            			mlist[m]=NULL;
		 			}
          			else {													// move to another live list
            			ll2=mptr->list;
            			bptr=mptr->box;
						// printf("molsort() line 2769 m=%d, ll2=%d, serno=%d,", m, ll2, mptr->serno);
			
         	  	 		if(mptr->box) boxremovemol(mptr,ll);
            			if(nl[ll2]==maxl[ll2])
              			if(molexpandlist(mols,sim->dim,ll2,-1,0)) {
                			simLog(sim,10,"out of memory in molsort\n");return 1;}
            			live[ll2][nl[ll2]++]=mptr;
						mptr->m=nl[ll2];
            			mlist[m]=NULL;
            			if(listtype[ll2]==MLTsystem) {
              			if(bptr) mptr->box=bptr;
              			else mptr->box=pos2box(sim,mptr->pos);
              			if(boxaddmol(mptr,ll2)) {
                			simLog(sim,10,"out of memory in molsort\n");return 1;}}}

         		 mlist[m]=mlist[--topl[ll]];				// compact original live list
          		 mlist[topl[ll]]=mlist[--nl[ll]];
          		 mlist[nl[ll]]=NULL;
          		 m--; }}}
	}

	for(m=mols->topd;m<mols->nd;m++) {		// move molecules from resurrected to reborn
		// printf("molsort, time=%f, topd=%d, nd=%d\n", sim->time, mols->topd, mols->nd);
		mptr=dead[m];
		ll2=mptr->list;

		if(nl[ll2]==maxl[ll2])
			if(molexpandlist(mols,sim->dim,ll2,-1,0)) {
				simLog(sim,10,"out of memory in molsort\n");return 1;}
		live[ll2][nl[ll2]++]=mptr;
		mptr->m=nl[ll2];
		dead[m]=NULL;
		if(listtype[ll2]==MLTsystem) {
				if(boxaddmol(mptr,ll2)) {
				simLog(sim,10,"out of memory in molsort\n");return 1;}}}
	mols->nd=mols->topd;

  	if(!onlydead2live) {
    	for(ll=0;ll<nlist;ll++)								// reset sortl indicies
      	sortl[ll]=nl[ll]; }
	
	/*	
	// create Mlist, indices for shuffling	
	if(!mols->Mlist){
		mols->Mlist=(int**) calloc(mols->nlist+1,sizeof(int*));			// the last row records the corresponding list length
		mols->Mlist[mols->nlist]=(int*) calloc(mols->nlist,sizeof(int));	
		for(ll=0;ll<mols->nlist;ll++) 
			mols->Mlist[mols->nlist][ll]=0;
	}

	for(ll=0;ll<mols->nlist;ll++){
		if(mols->nl[ll]!=mols->Mlist[mols->nlist][ll]){	
			list_new=(int*) calloc(mols->nl[ll],sizeof(int));
			for(m=0;m<mols->nl[ll];m++){
				list_new[m]=m;
			}
			std::random_shuffle(&list_new[0],&list_new[mols->nl[ll]-1]);
			free(mols->Mlist[ll]);
			mols->Mlist[ll]=list_new;
			mols->Mlist[mols->nlist][ll]=mols->nl[ll];
		}
		else if(mols->nl[ll]>0){
			std::random_shuffle(&(sim->mols->Mlist[ll][0]),&(sim->mols->Mlist[ll][mols->nl[ll]-1]));
		}
	}
	*/
	return 0; }


/* moldosurfdrift */
void moldosurfdrift(simptr sim,moleculeptr mptr,double dt) {
	int i,s,axis;
	enum MolecState ms;
	enum PanelShape ps;
	double *****surfdrift,vect[3],drift1,drift2,*pt1,*pt2,dist,unit0[3],unit1[3],unit2[3],top[3];
	panelptr pnl;

	i=mptr->ident;
	ms=mptr->mstate;
	pnl=mptr->pnl;
	s=pnl->srf->selfindex;
	ps=pnl->ps;
	surfdrift=sim->mols->surfdrift;
	vect[0]=vect[1]=vect[2]=0;

	if(surfdrift[i][ms][s] && surfdrift[i][ms][s][ps]) {
		if(sim->dim==2) {
			drift1=surfdrift[i][ms][s][ps][0]*dt;
			if(ps==PSrect)
				vect[(int)(pnl->front[2])]=drift1;
			else if(ps==PStri || ps==PScyl) {
				vect[0]=-drift1*pnl->front[1];
				vect[1]=drift1*pnl->front[0]; }
			else if(ps==PSsph || ps==PShemi) {
				vect[0]=-drift1*(mptr->pos[1]-pnl->point[0][1])/pnl->point[1][0];
				vect[1]=drift1*(mptr->pos[0]-pnl->point[0][0])/pnl->point[1][0]; }
			else if(ps==PSdisk) {
				pt1=mptr->pos;
				pt2=pnl->point[0];
				dist=sqrt((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+(pt2[1]-pt1[1])*(pt2[1]-pt1[1]));
				if(dist>VERYCLOSE) {
					vect[0]=drift1*(pt2[0]-pt1[0])/dist;
					vect[1]=drift1*(pt2[1]-pt1[1])/dist; }
				else {
					vect[0]=-drift1*pnl->front[1];
					vect[1]=drift1*pnl->front[0]; }}
			mptr->pos[0]+=vect[0];
			mptr->pos[1]+=vect[1]; }

		else {
			drift1=surfdrift[i][ms][s][ps][0]*dt;
			drift2=surfdrift[i][ms][s][ps][1]*dt;
			if(ps==PSrect) {
				vect[(int)(pnl->front[2])]=drift1;
				axis=0;
				if(axis==(int)(pnl->front[1]) || axis==(int)(pnl->front[2])) axis++;
				if(axis==(int)(pnl->front[1]) || axis==(int)(pnl->front[2])) axis++;
				vect[axis]=drift2; }
			else if(ps==PStri) {
				Geo_TriUnitVects(pnl->point[0],pnl->point[1],pnl->point[2],unit0,unit1,unit2);
				vect[0]=drift1*unit1[0]+drift2*unit2[0];
				vect[1]=drift1*unit1[1]+drift2*unit2[1];
				vect[2]=drift1*unit1[2]+drift2*unit2[2]; }
			else if(ps==PSsph) {
				top[0]=pnl->point[0][0];
				top[1]=pnl->point[0][1];
				top[2]=pnl->point[0][2]+pnl->point[1][0];
				Geo_SphereUnitVects(pnl->point[0],top,mptr->pos,(int)(pnl->front[0]),unit0,unit1,unit2);
				vect[0]=drift1*unit1[0]+drift2*unit2[0];
				vect[1]=drift1*unit1[1]+drift2*unit2[1];
				vect[2]=drift1*unit1[2]+drift2*unit2[2]; }
			else if(ps==PScyl) {
				Geo_CylUnitVects(pnl->point[0],pnl->point[1],mptr->pos,(int)(pnl->front[2]),unit0,unit1,unit2);
				vect[0]=drift1*unit1[0]+drift2*unit2[0];
				vect[1]=drift1*unit1[1]+drift2*unit2[1];
				vect[2]=drift1*unit1[2]+drift2*unit2[2]; }
			else if(ps==PShemi) {
				top[0]=pnl->point[0][0]-pnl->point[2][0];
				top[1]=pnl->point[0][1]-pnl->point[2][1];
				top[2]=pnl->point[0][2]-pnl->point[2][2];
				Geo_SphereUnitVects(pnl->point[0],top,mptr->pos,(int)(pnl->front[0]),unit0,unit1,unit2);
				vect[0]=drift1*unit1[0]+drift2*unit2[0];
				vect[1]=drift1*unit1[1]+drift2*unit2[1];
				vect[2]=drift1*unit1[2]+drift2*unit2[2]; }
			else if(ps==PSdisk) {
				Geo_DiskUnitVects(pnl->point[0],pnl->front,mptr->pos,unit0,unit1,unit2);
				vect[0]=drift1*unit1[0]+drift2*unit2[0];
				vect[1]=drift1*unit1[1]+drift2*unit2[1];
				vect[2]=drift1*unit1[2]+drift2*unit2[2]; }
			mptr->pos[0]+=vect[0];
			mptr->pos[1]+=vect[1];
			mptr->pos[2]+=vect[2]; }}

	return; }


/* diffuse */
int diffuse(simptr sim) {
	molssptr mols;
	int ll,m,d,nmol,dim,i,ngtablem1;
	enum MolecState ms;
	double flt1;
	double v1[DIMMAX],v2[DIMMAX],**difstep,***difm,***drift,epsilon,margin,neighdist,*gtable,dt;
	moleculeptr *mlist;
	moleculeptr mptr, mptr_tmp;
	int incmpt_flag=0;	
	int incmpt_posx_flag=0;
	int m_next;
	complexptr cplx;

	if(!sim->mols) return 0;
	dim=sim->dim;
	mols=sim->mols;
	ngtablem1=mols->ngausstbl-1;
	gtable=mols->gausstbl;
	difstep=mols->difstep;
	difm=mols->difm;
	drift=mols->drift;
	dt=sim->dt;
	flt1=sqrt(2.0*dt);
	epsilon=(sim->srfss)?sim->srfss->epsilon:0;
	margin=(sim->srfss)?sim->srfss->margin:0;
	neighdist=(sim->srfss)?sim->srfss->neighdist:0;
	double offset[dim];
	int updated_flag;									// to check whether a bound molec has been updated or not, to prevent double update
	double adj_factor;
	char vstr[STRCHAR], *vstr1;
	double volt,vtime;

	for(ll=0;ll<mols->nlist;ll++){
		mlist=mols->live[ll];
		nmol=mols->nl[ll];
		for(m=0;m<nmol;m++){
			mptr=mlist[m];
			mptr->sites_valx=mptr->sites_val=molecsites_state(sim->mols,mptr);
			if(mols->diffuselist[ll]){
				for(d=0;d<dim;d++) mptr->posx[d]=mptr->pos[d];	
				incmpt_posx_flag=boundarytest(sim,mptr->posx);
				if(incmpt_posx_flag==0){
					printf("diffuse, incmpt_flag=%d incmpt_posx_flag=%d complex_id=%d m=%d time=%f sites_val=%d ident=%d serno=%d s_index=%d posx[0]=%f posx[1]=%f posx[2]=%f pos[0]=%f pos[1]=%f pos[2]=%f prev_pos[0]=%f prev_pos[1]=%f prev_pos[2]=%f\n", incmpt_flag, incmpt_posx_flag, mptr->complex_id, m, sim->time, mptr->sites_val, mptr->ident, mptr->serno, mptr->s_index, mptr->posx[0], mptr->posx[1], mptr->posx[2],mptr->pos[0], mptr->pos[1], mptr->pos[2], mptr->prev_pos[0], mptr->prev_pos[1], mptr->prev_pos[2]);
					return -1;
				}
				if(mptr->tot_sunit>1 && mptr->s_index==0) sim->mols->complexlist[mptr->complex_id]->diffuse_updated=0;
			}
	
			if(sim->mols->volt_dependent[mptr->ident]==1){
				if(mptr->vchannel->vtime_n< sim->time){
					if(fgets(vstr,STRCHAR,mptr->vchannel->voltage_file)){
						vstr1=strsplit(vstr,"\t");
						sscanf(vstr1,"%lf",&volt);
						sscanf(vstr,"%lf",&vtime);	
						printf("volt=%f vtime=%f sim->time=%f\n",volt,vtime,sim->time);
						mptr->vchannel->voltage=volt;
						mptr->vchannel->vtime=sim->time;
						mptr->vchannel->vtime_n=vtime;
					}
			}}
			else sim->mols->volt_dependent[mptr->ident]=0;		
		}
	}
	
	for(ll=0;ll<mols->nlist;ll++)
		if(mols->diffuselist[ll]){
			mlist=mols->live[ll];
			nmol=mols->nl[ll];		
			m=0;
			mptr=mlist[0];
			for(m=0;m<nmol;m+=mptr->tot_sunit){
				updated_flag=0;
				mptr=mlist[m];
				i=mptr->ident;
				ms=mptr->mstate;
			
				if(mptr->complex_id!=-1){ 
					cplx=sim->mols->complexlist[mptr->complex_id];
					if(cplx->dif_molec)
						mptr=cplx->dif_bind;
				}
				else cplx=NULL;
				
				if(mptr->pnl && mols->surfdrift && mols->surfdrift[i] && mols->surfdrift[i][ms]){
					if(mptr->s_index==0){
						moldosurfdrift(sim,mptr,dt);						// surface drift
						// complex_pos(sim,mptr,"posx",);
				}}
				
				/*
				if(drift[i][ms]){  											// drift 
					if(mptr->s_index==0){															
						for(d=0;d<dim;d++) {
							mptr->prev_pos[d]=mptr->pos[d];
							mptr->pos[d]+=drift[i][ms][d]*dt; 
							offset[d]=mptr->pos[d]-mptr->posx[d];
						}
						complex_pos(sim,mptr,"posx",&offset[0],1);
				}}
				*/
				if(!difm[i][ms]){															// isotropic diffusion
					if(sim->interfaces) {
						if(sim->interfaces[mptr->ident]){	
							if(posincompart(sim,mptr->pos,sim->interfaces[mptr->ident]->cmpt)) {
								theta,difc=f(sim->interfaces[mptr->ident]->cmpt, mptr->pos)
								exacthittingtime
								mptr_new_pos=sbm()



					}
					}else { difc_type=1;}	
				
					
					if(difc_type==1) {
						if(mptr->tot_sunit==1 && mptr->pos==mptr->pos_tmp) {
							for(d=0;d<dim;d++) {
									mptr->prev_pos[d]=mptr->pos[d];	
									mptr->pos[d]+=difstep[i][ms]*gtable[randULI()&ngtablem1];	
							}
						}
						else if(mptr->complex_id!=-1){
							if(cplx->dif_molec){
								//i=cplx->dif_molec->ident;
								i=cplx->dif_molec->bind_id;
								ms=cplx->dif_molec->mstate;
							}
							updated_flag=mols->complexlist[mptr->complex_id]->diffuse_updated;
							if(updated_flag==1);
							else{	
								for(d=0;d<dim;d++) {
									mptr->prev_pos[d]=mptr->pos[d];			
									mptr->pos[d]+=difstep[i][ms]*gtable[randULI()&ngtablem1]; 
									offset[d]=mptr->pos[d]-mptr->posx[d];
								}
								if(complex_pos(sim,mptr,"posx", &offset[0],1)==-1){
									printf("time=%f, m=%d, mptr->ident=%d, mptr->s_index=%d, mptr->serno=%d, pos0=%f, pos1=%f, pos2=%f, offset[0]=%f, offset[1]=%f, offset[2]=%f\n", sim->time, m, mptr->ident, mptr->s_index,mptr->serno, mptr->pos[0], mptr->pos[1], mptr->pos[2], offset[0], offset[1], offset[2]);
									return -1;
					}}}}
				}else{																	    // anisotropic diffusion
					for(d=0;d<dim;d++)
						v1[d]=flt1*gtable[randULI()&ngtablem1];
					dotMVD(difm[i][ms],v1,v2,dim,dim);
					for(d=0;d<dim;d++) mptr->pos[d]+=v2[d]; }

				if(mptr->mstate!=MSsoln){													// surface-bound molecules
					if(dim>1)
						movemol2closepanel(sim,mptr,dim,epsilon,neighdist,margin);
					else
						mptr->pos[0]=mptr->posx[0]; }
				}}									// 1D surface-bound molecules aren't allowed to move

	return 0; }

/*complex_pos_init for camkii ring strucutre, 2 layers*/
void complex_pos_init(simptr sim, moleculeptr mptr,compartptr cmpt){
	// if(sim->time>0) return;
	if(mptr->s_index!=0) return;
	if(mptr->complex_id==-1) return;
	if(!sim->mols->complexlist) return;

	moleculeptr mptr_tmp;
	double r=8; // 80 angstrom = 8 nm
	// double r=0;
	int sunit_i;
	int total_sunit=mptr->tot_sunit;
	double theta_tmp, k, g, h;
	double x,y,z, x0, y0, z0, x1,y1,z1;
	int d;
	double tmp_pos[3];
	int layer=sim->mols->complexlist[mptr->complex_id]->layer;
	int group_sunit=total_sunit/layer;
	
	/* solve for (x0,y0,z0) at the center of the ring based on one cornor point
		(x0,y0,z0) is not the origin, but the center of ring, a variable
	*/
	k=tan(mptr->phi_init);
	g=cos(mptr->theta_init);
	x=mptr->pos[0];
	y=mptr->pos[1];
	z=mptr->pos[2];
	x0=x + sqrt(r*r*(1-g*g)/(1+k*k));	
	y0=y + k*(x0-x);
	z0=z + r*g;

	// layer2, with a shifted ring center, which is perpendicular to the plain phi=phi_init
	if(layer>1){
		h=10;
		x1=x0+h*sin(mptr->phi_init);
		y1=y0-h*cos(mptr->phi_init);
		z1=z0;
	}

	//#define PI 3.1415926 
	// PI is actually defined in lib/math2.h
	mptr_tmp=mptr->to;
	while(mptr_tmp!=NULL && mptr_tmp->s_index<total_sunit && mptr_tmp!=mptr){	
		theta_tmp=2*PI/(total_sunit/layer)*(mptr_tmp->s_index % (total_sunit/layer))+mptr_tmp->theta_init;

		// all subunits have the same phi
		if(mptr_tmp->s_index<group_sunit){
			mptr_tmp->pos[2]= z0-r*cos(theta_tmp);
			mptr_tmp->pos[1]= y0-r*sin(theta_tmp)*sin(mptr_tmp->phi_init);
			mptr_tmp->pos[0]= x0-r*sin(theta_tmp)*cos(mptr_tmp->phi_init);
			//printf("theta_tmp=%f, mptr_tmp->s_index=%d\n",theta_tmp,mptr_tmp->s_index);
		}
		/*		
		else if(mptr_tmp->s_index==group_sunit){
			mptr_tmp->pos[0]=x1-sqrt(r*r*(1-g*g)/(1+k*k));
			mptr_tmp->pos[1]=y1-k*(x1-mptr_tmp->pos[0]);
			mptr_tmp->pos[2]=z1-r*g;
			printf("theta_tmp=%f, mptr_tmp->s_index=%d, s_index=%d, dist^2=%f\n",theta_tmp,mptr_tmp->s_index,mptr_tmp->from->from->from->from->from->from->s_index, molec_distance(sim,mptr_tmp->pos,mptr_tmp->from->from->from->from->from->from->pos));

		}
		*/
		else{
			mptr_tmp->pos[2]= z1-r*cos(theta_tmp);
			mptr_tmp->pos[1]= y1-r*sin(theta_tmp)*sin(mptr_tmp->phi_init);
			mptr_tmp->pos[0]= x1-r*sin(theta_tmp)*cos(mptr_tmp->phi_init);
			printf("theta_tmp=%f, mptr_tmp->s_index=%d, s_index=%d, dist^2=%f\n",theta_tmp,mptr_tmp->s_index,mptr_tmp->from->from->from->from->from->from->s_index, molec_distance(sim,mptr_tmp->pos,mptr_tmp->from->from->from->from->from->from->pos));

		}


		if(cmpt==NULL){
			mptr_tmp=mptr_tmp->to;
			// sunit_i++;
		}
		else{
			for(d=0;d<sim->dim;d++) 
				tmp_pos[d]=mptr_tmp->pos[d];
			while(!posincompart(sim,tmp_pos,cmpt)){
				compartrandpos(sim,tmp_pos,cmpt);
				z=tmp_pos[2];
				y=tmp_pos[1];
				x=tmp_pos[0];
				x0=x + sqrt(r*r*(1-g*g)/(1+k*k));	
				y0=y + k*(x0-x);
				z0=z + r*g;	
				theta_tmp=mptr->theta_init;
				mptr_tmp=mptr;
				// sunit_i=0;		// then the following code makes sunit_i=1
				if(layer>1){
					x1=x0+h*sin(mptr->phi_init);
					y1=y0+h*cos(mptr->phi_init);
					z1=z0;
				}
			}
			mptr_tmp->pos[2]=tmp_pos[2];
			mptr_tmp->pos[1]=tmp_pos[1];
			mptr_tmp->pos[0]=tmp_pos[0];
			mptr_tmp=mptr_tmp->to;
			//sunit_i++;
		}
	}
	// printf("assign position done, mptr_tmp->s_index=%d, mptr_tmp->to->s_index=%d\n", mptr_tmp->s_index, mptr_tmp->to->s_index);
	
	if(total_sunit>1){	
		for(sunit_i=1;sunit_i<=total_sunit;sunit_i++){
			mptr_tmp->sdist_init=molec_distance(sim,mptr_tmp->pos,mptr_tmp->to->pos);
			mptr_tmp=mptr_tmp->to;
	}}
		
	return;
}

// Q: do I need to convert everything like diffusion to the spherical coordinate?
// when loop around the subunits, use "from" (serno increase) is better than "to" (serno decreases), because when loop through the molecule list, serno number decreases? 
int complex_pos(simptr sim, moleculeptr mptr, char* pos_to_update, double* offset, int dist_calculation){			
	moleculeptr mptr_tmp, mptr_bind;
	int sunit_i, site,sunit;
	int d, d1;
	int total_sunit=mptr->tot_sunit;
    double sunit_dist, offset_tmp[3]={0,0,0};
	int dim=sim->dim;
	int s,k;
	GHashTable* complex_connect=sim->mols->complex_connect;
	complexptr complex_bind, complex_tmp;
	char pos_tmp[STRCHAR];
	int sindex_tmp, sindex_bind;
	
	if(!sim->mols->complexlist) return 0;
	complex_tmp=sim->mols->complexlist[mptr->complex_id];
	if(offset[0]!=0 || offset[1]!=0 || offset[2]!=0){
		for(sunit_i=1,mptr_tmp=mptr;sunit_i<total_sunit && mptr_tmp->to && mptr_tmp!=mptr_tmp->to; sunit_i++, mptr_tmp=mptr_tmp->to){
				for(d=0;d<dim;d++){
					// if allow sunit_i == total_sunit, mptr_tmp positions will be updated twice
					if(strstr(pos_to_update,"posx")){
						mptr_tmp->to->prev_pos[d]=mptr_tmp->to->posx[d];
						mptr_tmp->to->pos[d]=mptr_tmp->to->posx[d]+offset[d];
					}
					else{
						mptr_tmp->to->prev_pos[d]=mptr_tmp->to->pos[d];
						mptr_tmp->to->pos[d]=mptr_tmp->to->pos[d]+ offset[d];	// r_offset[d];
					}
		}}
		complex_tmp->diffuse_updated++;
		/*
		if(pos_to_update){ 
			if(complex_connect){
				mptr_tmp=mptr;
				for(s=0;s<total_sunit;s++){
					for(k=0;k<sim->mols->spsites_num[mptr->ident];k++){
						if(!mptr_tmp->sites[k]->bind) continue;
						if(mptr_tmp->sites[k]->bind->complex_id==-1) continue;			// not a complex bound at the site
						mptr_bind=mptr_tmp->sites[k]->bind;
						complex_bind=sim->mols->complexlist[mptr_bind->complex_id];
						// mptr_bind may equal mptr;
						printf("molec.c: 3363, time=%f, pos_to_update: %s, mptr_tmp->complex_id=%d, mptr_bind->complex_id=%d, tmp_sindex=%d, bind_sindx=%d\n", sim->time, pos_to_update, mptr_tmp->complex_id, mptr_bind->complex_id, mptr_tmp->s_index, mptr_bind->s_index);
						mptr_bind->sites[k]->bind=NULL;
						printf("molec.c: 3365, unbind, mptr_bind->serno=%d, s_index=%d\n", mptr_bind->serno, mptr_bind->s_index);
						printf("molec.c: 3366, mptr_bind, complex_bind->diffuse_updated=%d, pos0=%f, pos1=%f, pos2=%f, prevpos0=%f, prevpos1=%f, prevpos2=%f, posx0=%f, posx1=%f, posx2=%f\n", complex_bind->diffuse_updated, mptr_bind->pos[0], mptr_bind->pos[1], mptr_bind->pos[2], mptr_bind->prev_pos[0], mptr_bind->prev_pos[1], mptr_bind->prev_pos[2], mptr_bind->posx[0], mptr_bind->posx[1], mptr_bind->posx[2]);
						printf("molec.c: 3367, mptr_tmp, complex_tmp->diffuse_updated=%d, pos0=%f, pos1=%f, pos2=%f, prevpos0=%f, prevpos1=%f, prevpos2=%f, posx0=%f, posx1=%f, posx2=%f\n", complex_tmp->diffuse_updated, mptr_tmp->pos[0], mptr_tmp->pos[1],mptr_tmp->pos[2], mptr_tmp->prev_pos[0], mptr_tmp->prev_pos[1], mptr_tmp->prev_pos[2], mptr_tmp->posx[0], mptr_tmp->posx[1], mptr_tmp->posx[2]);
						strcpy(pos_tmp,pos_to_update);
						strcat(pos_tmp,"_line3369");
						printf("mptr_bind->pos==mptr_tmp->pos: %d, mptr_tmp->pos==mptr_bind->pos:%d\n", mptr_bind->pos==mptr_tmp->pos, mptr_tmp->pos==mptr_bind->pos);
						if(mptr_bind->prev_pos[0]!=mptr_tmp->prev_pos[0] || mptr_bind->prev_pos[1]!=mptr_tmp->prev_pos[1] || mptr_bind->pos[2]!=mptr_tmp->pos[2]){
							if(complex_tmp->diffuse_updated > complex_bind->diffuse_updated)	
								for(d=0;d<sim->dim;d++) mptr_bind->prev_pos[d]=mptr_bind->pos[d]-offset[d];
							else printf("molec.c line3374, complex_tmp->diffuse_updated=%d, complex_bind->diffuse_updated=%d\n", complex_tmp->diffuse_updated, complex_bind->diffuse_updated);
						}
						complex_pos(sim,mptr_bind,pos_tmp,offset,dist_calculation);	
						mptr_bind->sites[k]->bind=mptr_tmp;
						printf("bind\n");
				}
				mptr_tmp=mptr_tmp->to;
			}}}
		*/
		}
	
		if(dist_calculation){					// wait for all sunits pos updated to update sunits dist
		for(sunit_i=1,mptr_tmp=mptr;sunit_i<=total_sunit && mptr_tmp->to && mptr_tmp!=mptr_tmp->to; sunit_i++, mptr_tmp=mptr_tmp->to){
			if(mptr_tmp->tot_sunit==1) mptr_tmp->sdist_tmp=0; 
			else{
				mptr_tmp->sdist_tmp=molec_distance(sim,mptr_tmp->pos,mptr_tmp->to->pos);
				if(mptr_tmp->sdist_tmp - mptr_tmp->sdist_init>1 || mptr_tmp->sdist_tmp-mptr_tmp->sdist_init<-1){
					printf("molec.c line3388, %s, complex_tmp->diffuse_updated=%d, mptr_tmp->serno=%d, sites_val=%d, s_index=%d, to->serno=%d, to->sites_val=%d, to->s_index=%d\n", pos_to_update, complex_tmp->diffuse_updated, mptr_tmp->serno, mptr_tmp->sites_val, mptr_tmp->s_index, mptr_tmp->to->serno, mptr_tmp->to->sites_val, mptr_tmp->to->s_index); 
					printf("molec.c mptr_tmp pos0=%f, pos1=%f, pos2=%f, prevpos0=%f, prevpos1=%f, prevpos2=%f, posx0=%f, posx1=%f, posx2=%f\n", mptr_tmp->pos[0], mptr_tmp->pos[1], mptr_tmp->pos[2], mptr_tmp->prev_pos[0], mptr_tmp->prev_pos[1], mptr_tmp->prev_pos[2], mptr_tmp->posx[0], mptr_tmp->posx[1], mptr_tmp->posx[2]); 
					printf("molec.c to->mptr_tmp pos0=%f, pos1=%f, pos2=%f, prevpos0=%f, prevpos1=%f, prevpos2=%f, posx0=%f, posx1=%f, posx2=%f\n", mptr_tmp->to->pos[0], mptr_tmp->to->pos[1], mptr_tmp->to->pos[2], mptr_tmp->to->prev_pos[0], mptr_tmp->to->prev_pos[1], mptr_tmp->to->prev_pos[2], mptr_tmp->to->posx[0], mptr_tmp->to->posx[1], mptr_tmp->to->posx[2]); 
					return -1;
		}}}}
	return 0;
}

double power(double a, int b){
	int k=0;
	double tmp=1;

	while(k<b){
		tmp=tmp*a;		
		k++;
	}
	return tmp;
}
 

/* actually, distance squared, only for connected molecs */
double molec_distance(simptr sim, double *pos1, double *pos2){
	int d;
	int dim=sim->dim;
	double dist=0;
	
	for(d=0;d<dim;d++) {
			dist+=power(pos1[d]-pos2[d],2);
	}

	return dist;	
}

int cplxpos_updated_test(simptr sim, moleculeptr mptr, double *offset, int *sunit){
	int d,s;
	int updated_flag=1;
	moleculeptr mptr_tmp;

	mptr_tmp=mptr;
	for(s=0;s<mptr->tot_sunit;s++){
		for(d=0;d<sim->dim;d++){
			if(fabs(mptr->pos[d]-mptr->prev_pos[d]-offset[d])>1e-2){ 
				*sunit=mptr_tmp->s_index;
				return 0;
			}
		}
		mptr_tmp=mptr_tmp->to;
	}
	return updated_flag;
}


/* test if a position inside or outside */
int boundarytest(simptr sim, double *pos){
	compartptr cmptptr;
	int poscmpt=0;		
	int c;
	
	if(!sim->cmptss) return 1;
	if(pos[0]==0 && pos[1]==0 && pos[2]==0) return 1;

	// for(c=0;c<sim->cmptss->ncmpt;c++){
		cmptptr=sim->cmptss->cmptlist[0];			// cmptlist[0] is the compelete compartment
		poscmpt=posincompart(sim,pos,cmptptr);
		if(poscmpt) return poscmpt;
	// }
	return poscmpt;	
}


double ExactHittingTime_sbm(simptr sim, moleculeptr mptr,double difc){
	//z=x/sqrt(2.0*difc);
	//y=z+sqrt(sim->dt*eta);
	x=mptr->pos[2];
	x_tmp=x+sqrt(2.0*difc*sim->dt*fabs(rnd));

	if(sgn(z-interface_pos)!=sgn(y-interface_pos)){
		tau=sim->dt*rnd/(1.0+rnd);
		return tau;
	}	
	else{
		eta=invGaussian(fabs(z)/fabs(y),z**2/sim->dt);
		u=randCOD();
		cross_prob=exp(-2.0*(interface_pos-x)*(interface_pos-x_tmp)/(2.0*difc*sim->dt));
		if(u<cross_prob){
			eta=invGaussian(fabs(z)/fabs(y),z**2/sim->dt);
			tau=sim->dt*eta/(1.0+eta);
			return tau;
		}
		else{
			return sim->dt+time
		}
	}
}

double sgn(double x){
	if(x>0)
		return 1;
	else if(x<0)
		return -1;
	else 
		return 0;
}


double SBM(simptr sim, moleculeptr mptr, cmpt){
	double s,rnd_u,rnd_g,theta
	s=ExactHittingTime_sbm(sim,mptr,&y);

	if(s<sim->time+sim->dt){
		rnd_u=randCOD()
		// gsl_ran_gaussian_ziggurat(const gsl_rng *r, double sigma)
		gsl_ran_gaussian_ziggurat(&rnd_g,1.0);

		// interface naturall seperates two compartment, each with distinct D; maybe a linkedlist structure
		difc1=interface->cmpt->difc 		// the current difc
		difc2=interface->cmpt->next->difc 	// the next difc	
		theta=(sqrt(difc2)-sqrt(difc1))/(sqrt(difc2)+sqrt(difc1));

		if(u<(1.0+theta)/2.0){		
			for(d=0;d<sim->dim-1;d++)
				mptr->pos[d]=mptr->pos[d]+sqrt(2*difc1*s)*fabs(G2)+sqrt(2*difc2*(sim->dt-s))*fabs(G3);			
			mptr->pos[2]=z_interface+sqrt(2*difc1*(sim->dt-s))*fabs(G2);
		}
		else{ 
			for(d=0;d<sim->dim-1;d++)	
				mptr->pos[d]=mptr->pos[d]+sqrt(2*difc
			mptr->pos[2]=z_interface+sqrt(2*difc2*(sim->dt-s))*fabs(G2);
		}
	}
	else
		return x_interface+y; 
}

double invGaussian(double mu, double lambda){
	double rnd,test,y,x;
	
	gsl_ran_gaussian_ziggurat(&rnd,1.0);
	y=rnd*rnd;
	x=mu+(mu*mu*y)/(2.0*lambda)-(mu/(2.0*lambda))*sqrt(4.0*mu*lambda*y + mu*mu*y*y);	
	test=randCOD();
	if(test<=(mu)/(mu+x))
		return x;
	else 
		return (mu*mu)/x;
}


}
