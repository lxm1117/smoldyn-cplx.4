/* Steven Andrews, started 10/22/2001.
 This is the header for data structures for the Smoldyn program.
 See documentation called Smoldyn_doc1.pdf and Smoldyn_doc2.pdf, and the Smoldyn
 website, which is at www.smoldyn.org.
 Copyright 2003-2013 by Steven Andrews.  This work is distributed under the terms
 of the Gnu Lesser General Public License (LGPL). */

#ifndef __smoldyn_h__
#define __smoldyn_h__

#include <time.h>
#include <stdio.h>
#include "List.h"
#include "smoldynconfigure.h"			// generated by CMake from smoldynconfigure.h.in
#include <glib.h>
#include <map>
#include <utility>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#ifdef OPTION_NSV
  #include "nsvc.h"
#endif

#ifdef OPTION_VTK
  #include "vtkwrapper.h"
#endif

#ifdef OPTION_VCELL
	#include <string>
	using namespace std;
#endif

/********************************** General *********************************/

#define DIMMAX 3					// maximum system dimensionality
#define VERYCLOSE 1.0e-12			// distance that's safe from round-off error

enum StructCond {SCinit,SClists,SCparams,SCok};

/********************************* Molecules ********************************/

// #define MSMAX 5
// #define MSMAX1 6
#define MSMAX1 9
#define MSMAX 9
enum MolecState {MSsoln,MSfront,MSback,MSup,MSdown,MSbsoln,MSall,MSnone,MSsome,MSimmobl};
enum MolListType {MLTsystem,MLTport,MLTnone};
#define PDMAX 7
enum PatternData {PDalloc,PDnresults,PDnspecies,PDmatch,PDsubst,PDdegen,PDrule};

typedef struct sitestruct{
	int site_type;							// 0 for states, 1 for binding sites
	struct moleculestruct *bind;	
	int* value;
	int* value_tmp;							// value pointer to its original place
	int indx;
	double time;							// record the latest checked timestamp for a reaction event;
} *siteptr;

typedef struct vchnlstruct{
	struct moleculestruct *mptr; 		
	FILE *voltage_file;
	double voltage;
	double vtime;
	double vtime_n;
	int molec_gen;							// for calcium channel to generate ca2+ in particular 
	// double half_n;						// v1/2 of activation gate			
	// double slope_n;
	// double half_h;
	// double slope_h;							
	// double n_rate;
	// double h_rate;
	// double block_f;
	// double block_r;

} *vchnlptr;

typedef struct moleculestruct {
	// long int serno;							// serial number
	int serno;
	int list;									// destination list number (ll)
	int m;
	double *pos;								// dim dimensional vector for position [d]
	double *posx;								// dim dimensional vector for old position [d]
	double *via;								// location of last surface interaction [d]
	double *posoffset;							// position offset arising from jumps [d]
	int ident;									// species of molecule; 0 is empty (i)
	enum MolecState mstate;			// physical state of molecule (ms)
	struct boxstruct *box;			// pointer to box which molecule is in
	struct panelstruct *pnl;		// panel that molecule is bound to if any
	int s_index;					// subunit index
	struct moleculestruct *from; 	// pointer in
	struct moleculestruct *to;	  	// pointer out
	int tot_sunit;
	double *pos_tmp;
	double theta_init;
	double phi_init;
	double sdist_init;				// distance between current subunit and its 'to' neighbor
	double sdist_tmp;
	double* prev_pos;				// record positions before the latest updated pos, usd for calculating pos_offset
	int complex_id;					// >0 if belongs to a complex; -1 otherwise
	siteptr *sites;
	int sites_val;
	int sites_valx;
	struct moleculestruct *dif_molec;
	int dif_site;
	// double adj_prob;				// accumulated probability of an upcoming rxn in a time step accounted for preccedingly occurred reactions on a molecule
	double sim_time;
	struct vchnlstruct *vchannel;	
	double arrival_time;
	int bind_id;
} *moleculeptr;

typedef struct complexstruct{
	int serno;
	moleculeptr zeroindx_molec;				// the zero index molecule to start with
	moleculeptr dif_molec;
	moleculeptr dif_bind;
	int dif_bind_site;
	int diffuse_updated;
	int layer;
} *complexptr;

/*
typedef struct difadjstruct{
	int molec_ident;
	struct compartstruct *cmpt; 		
	double adj_factor;
	struct difadjstruct *next;
} *difadj_ptr;
*/

typedef struct molsuperstruct {
	enum StructCond condition;				// structure condition
	struct simstruct *sim;					// simulation structure
	int maxspecies;							// maximum number of species
	int nspecies;							// number of species, including empty mols.
	char **spname;							// names of molecular species [i]
	int *volt_dependent;					// whether a species has voltage-dependent sensing or not
	int maxpattern;							// maximum number of patterns
	int npattern;							// actual number of patterns
	char **patlist;							// list of patterns [pat]
	int **patindex;							// species indicies for patterns [pat][j]
	double **difc;							// diffusion constants [i][ms]
	double **difstep;						// rms diffusion step [i][ms]
	double ***difm;							// diffusion matrix [i][ms][d]
	double ***drift;						// drift vector [i][ms][d]
	double *****surfdrift;					// surface drift [i][ms][s][ps][d]
	double **display;						// display size of molecule [i][ms] 
	double ***color;						// RGB color vector [i][ms]
	int **exist;							// flag for if molecule could exist [i][ms]
	moleculeptr *dead;						// list of dead molecules [m]
	int maxdlimit;							// maximum allowed size of dead list
	int maxd;								// size of dead molecule list
	int nd;									// total number of molecules in dead list
	int topd;								// index for dead list; above are resurrected
	int maxlist;							// allocated number of live lists
	int nlist;								// number of live lists
	int **listlookup;						// lookup table for live lists [i][ms]
	char **listname;						// names of molecule lists [ll]
	enum MolListType *listtype;				// types of molecule lists [ll]
	moleculeptr **live;						// live molecule lists [ll][m]
	int *maxl;								// size of molecule lists [ll]
	int *nl;								// number of molecules in live lists [ll]
	int *topl;								// live list index; above are reborn [ll]
	int *sortl;								// live list index; above need sorting [ll]
	int *diffuselist;						// 1 if any listed molecs diffuse [ll]
	long int serno;							// serial number for next resurrected molec.
	int ngausstbl;							// number of elements in gausstbl
	double *gausstbl;						// random numbers for diffusion
	int *expand;							// whether species expand with libmzr [i]

	complexptr *complexlist;				
	int ncomplex;
	int max_complex;
	char ***spsites_name;					// all sites regardless which species	
	int *spsites_num;						// max sites for each species
	int **spsites_binding;					// whether a pariticular site allows bindings or not
	int maxsitecode;
	// int max_sites;
	GHashTable* spdifsites;			
	GHashTable* complex_connect;		
	int **Mlist;							// indices for shuffling molecular list

} *molssptr;

/*********************************** Walls **********************************/

typedef struct wallstruct {
	int wdim;										// dimension number of perpendicular to wall
	int side;										// low side of space (0) or high side (1)
	double pos;										// position of wall along dim axis
	char type;										// properties of wall
	struct wallstruct *opp; 						// pointer to opposite wall
	} *wallptr;

/********************************* Reactions ********************************/
#define HashNum 97
#define MAXORDER 3
#define MAXPRODUCT 16
enum RevParam {RPnone,RPirrev,RPconfspread,RPbounce,RPpgem,RPpgemmax,RPpgemmaxw,RPratio,RPunbindrad,RPpgem2,RPpgemmax2,RPratio2,RPoffset,RPfixed, RPnone1};

#ifdef OPTION_VCELL
class ValueProvider;
typedef ValueProvider* valueproviderptr;
#endif
typedef struct std::pair<std::string,std::string> str_pair;
//typedef struct std::pair<str_pair,std::string> str2_pair;
typedef struct std::map<str_pair,double> radius_map;

typedef struct rctmolecule{						// unit molecule involved in a reaction
	int ident;						
	int *states;								// could be reacting sites or product sites		
	int *sites_indx;
	int sites_num;
	int states_num;	
	enum MolecState rctstate;
	void *difc;
} *rct_mptr;	

typedef struct prdmolecule{
	int ident;
	int *sites_indx;
	int *sites_val;
	int sites_num;
	int site_bind;
	enum MolecState prdstate;
} *prd_mptr;

typedef struct rxnstruct {
	struct rxnsuperstruct *rxnss;				// pointer to superstructure
	char *rname;								// pointer to name of reaction
	int *permit;								// permissions for reactant states [ms]
	int order;
	int nprod;									// number of products
	int molec_num;
	rct_mptr *rct;
	prd_mptr *prd;

	long int *prdserno;							// list of product serno rules [prd]
	listptrli logserno;							// list of serial nums for logging reaction
	char *logfile;								// filename for logging reaction
	double rate;								// requested reaction rate
	double bindrad;								// squared binding radius, if appropriate
	double k_on;
	double k_off;
	double prob;								// reaction probability
	double tau;									// characteristic reaction time
	double phi;									// for 2nd order reactions;
	double miu;
	enum RevParam rparamt;						// type of parameter in rpar
	double rparam;								// parameter for reaction of products
	double unbindrad;						// unbinding radius, if appropriate
	double **prdpos;							// product position vectors [prd][d]
	int disable;								// 1 if reaction is disabled
	struct compartstruct *cmpt;					// compartment reaction occurs in, or NULL
	struct surfacestruct *srf;					// surface reaction on, or NULL
} *rxnptr;

typedef struct rxnsuperstruct {
	enum StructCond condition;		// structure condition
	struct simstruct *sim;			// simulation structure
	// int order;					// still needs order variable, for rxns with 0 order; no need, bc rxnss are grouped by molec_num
	int molec_num;
	int maxspecies;					// maximum number of species 
	int maxsitecode;
	int maxlist;					// copy of maximum number of molecule lists
	int *nrxn;						// number of rxns for each reactant set [i]
	GHashTable *table;				// lookup table for reaction numbers [i][j]
	GHashTable *entrylist; 			// contains all entry keys to the table
	GHashTable *bindrad_eff;		// using rxn entries as keys
	GHashTable *rxnaff;
	GHashTable *probrng_l;			// using rxn entries as keys
	GHashTable *probrng_h;
	//GHashTable *rmaps;				// rxn->rmap now moves to rxnss
	//GHashTable *rmaps_adj;
	GHashTable *radius;
	GHashTable *rxnr_ptr;
	GHashTable *rxn_ord1st;
	int **binding;					// for 2 moleculed reactions only
	
	int maxrxn;						// allocated number of reactions
	int totrxn;						// total number of reactions listed
	char **rname;					// names of reactions [r]
	rxnptr *rxn;					// list of reactions [r]
	// int *rxnmollist;				// live lists that have reactions [ll]
} *rxnssptr;

/********************************* Surfaces *********************************/

#define PSMAX 6															// maximum number of panel shapes
enum PanelFace {PFfront,PFback,PFnone,PFboth};
enum PanelShape {PSrect,PStri,PSsph,PScyl,PShemi,PSdisk,PSall,PSnone};
enum SrfAction {SAreflect,SAtrans,SAabsorb,SAjump,SAport,SAmult,SAno,SAnone,SAadsorb,SArevdes,SAirrevdes,SAflip,SArotate_jump};  // cplx
// act=0: reflect
enum DrawMode {DMno=0,DMvert=1,DMedge=2,DMve=3,DMface=4,DMvf=5,DMef=6,DMvef=7,DMnone};
enum SMLflag {SMLno=0,SMLdiffuse=1,SMLreact=2,SMLsrfbound=4};

typedef struct surfactionstruct {
	int *srfnewspec;						// surface convert mol. species [ms]
	double *srfrate;						// surface action rate [ms]
#ifdef OPTION_VCELL
	valueproviderptr* srfRateValueProvider;	//rate for surface actions: asorption, desorption, transmission...etc.
#endif
	double *srfprob;						// surface action probability [ms]
	double *srfcumprob;						// surface cumulative probability [ms]
	int *srfdatasrc;						// surface data source [ms]
	double *srfrevprob;						// probability of reverse action [ms]
} *surfactionptr;

typedef struct panelstruct {
	char *pname;							// panel name (reference, not owned)
	enum PanelShape ps;						// panel shape
	struct surfacestruct *srf;				// surface that owns this panel
	int npts;								// number of defining points
	double **point;							// defining points, [number][d]
	double front[DIMMAX];					// front parameters, which depend on the shape
	struct panelstruct *jumpp[2];			// panel to jump to, if appropriate [face]
	enum PanelFace jumpf[2];				// face to jump to, if appropriate [face]
	int maxneigh;							// maximum number of neighbor panels
	int nneigh;								// number of neighbor panels
	struct panelstruct **neigh;		// list of neighbor panels [p]
	double *emitterabsorb[2];		// absorption for emitters [face][i]
} *panelptr;

typedef struct surfacestruct {
	char *sname;							// surface name (reference, not owned)
	struct surfacesuperstruct *srfss;		// owning surface superstructure
	int selfindex;							// index of self
	enum SrfAction ***action;			// action for molecules [i][ms][face]
	// enum SrfAction ****action;				// srf->action[i][sitecode][ms][face]	
	surfactionptr ***actdetails;			// action details [i][ms][face]
	int neighhop;							// whether molecules hop between neighbors
	double fcolor[4];						// RGBA color vector for front
	double bcolor[4];						// RGBA color vector for back
	double edgepts;							// thickness of edge for drawing
	unsigned int edgestipple[2];			// edge stippling [factor,pattern]
	enum DrawMode fdrawmode;				// polygon drawing mode for front
	enum DrawMode bdrawmode;				// polygon drawing mode for back
	double fshiny;							// front shininess
	double bshiny;							// back shininess
	int maxpanel[PSMAX];				// allocated number of panels [ps]
	int npanel[PSMAX];					// actual number of panels [ps]
	char **pname[PSMAX];				// names of panels [ps][p]
	panelptr *panels[PSMAX];			// list of panels [ps][p]
	struct portstruct *port[2];			// port, if any, for each face [face]
	double totarea;						// total surface area
	int totpanel;						// total number of panels
	double *areatable;					// cumulative panel areas [pindex]
	panelptr *paneltable;				// sequential list of panels [pindex]
	int *maxemitter[2];					// maximum number of emitters [face][i]
	int *nemitter[2];					// number of emitters [face][i]
	double **emitteramount[2];			// emitter amounts [face][i][emit]
	double ***emitterpos[2];			// emitter positions [face][i][emit][d]
	} *surfaceptr;

typedef struct surfacesuperstruct {
	enum StructCond condition;		// structure condition
	struct simstruct *sim;			// simulation structure
	int maxspecies;					// maximum number of molecular species
	int maxsrf;						// maximum number of surfaces
	int nsrf;						// number of surfaces
	double epsilon;					// max deviation of surface-point from surface
	double margin;					// panel margin away from edge
	double neighdist;				// neighbor distance value
	char **snames;					// surface names [s]
	surfaceptr *srflist;			// list of surfaces [s]
	int maxmollist;					// number of molecule lists allocated
	int nmollist;					// number of molecule lists used
	enum SMLflag *srfmollist;		// flags for molecule lists to check [ll]
	} *surfacessptr;


typedef struct interface_struct{
	struct simstruct *sim;
	struct compartstruct *cmpt;
	double pos;
	double difc1;
	double difc2;  
	double side1;
	double side2;
	int species;
} *interfaceptr;
/*********************************** Boxes **********************************/

typedef struct boxstruct {
	int *indx;									// dim dimensional index of the box [d]
	int nneigh;									// number of neighbors in list
	int midneigh;								// logical middle of neighbor list
	struct boxstruct **neigh;					// all box neighbors, using sim. accuracy
	int *wpneigh;								// wrapping code of neighbors in list
	int nwall;									// number of walls in box
	wallptr *wlist;							// list of walls that cross the box
	int maxpanel;								// allocated number of panels in box
	int npanel;									// number of surface panels in box
	panelptr *panel;						// list of panels in box
	int *maxmol;								// allocated size of live lists [ll]
	int *nmol;									// number of molecules in live lists [ll]
	moleculeptr **mol;					// lists of live molecules in the box [ll][m]
	// molec_list *mol;	
	double* difadj;
	} *boxptr;

typedef struct boxsuperstruct {
	enum StructCond condition;	// structure condition
	struct simstruct *sim;			// simulation structure
	int nlist;									// copy of number of molecule lists
	double mpbox;								// requested number of molecules per box
	double boxsize;							// requested box width
	double boxvol;							// actual box volumes
	int nbox;										// total number of boxes
	int *side;									// number of boxes on each side of space
	double *min;								// position vector for low corner of space
	double *size;								// length of each side of a box
	boxptr *blist; 							// actual array of boxes
	} *boxssptr;

/******************************* Compartments *******************************/

enum CmptLogic {CLequal,CLequalnot,CLand,CLor,CLxor,CLandnot,CLornot,CLnone};

typedef struct compartstruct {
	struct compartsuperstruct *cmptss;	// compartment superstructure
	char *cname;						// compart. name (reference, not owned)
	int nsrf;							// number of bounding surfaces
	surfaceptr *surflist;				// list of bounding surfaces [s]
	int npts;								// number of inside-defining points
	double **points;						// list of inside-defining points [k][d]
	int ncmptl;								// number of logic compartments
	struct compartstruct **cmptl;	// list of logic compartments [cl]
	enum CmptLogic *clsym;			// compartment logic symbol [cl]
	double volume;							// volume of compartment
	int maxbox;								// maximum number of boxes in compartment
	int nbox;								// number of boxes inside compartment
	boxptr *boxlist;						// list of boxes inside compartment [b]
	double *boxfrac;						// fraction of box volume that's inside [b]
	double *cumboxvol;						// cumulative cmpt. volume of boxes [b]
	double *difadj;
	} *compartptr;

typedef struct compartsuperstruct {
	enum StructCond condition;	// structure condition
	struct simstruct *sim;			// simulation structure
	int maxcmpt;								// maximum number of compartments
	int ncmpt;									// actual number of compartments
	char **cnames;							// compartment names [c]
	compartptr *cmptlist;				// list of compartments [c]
	} *compartssptr;

/*********************************** Ports **********************************/

typedef struct portstruct {
	struct portsuperstruct *portss;	// port superstructure
	char *portname;							// port name (reference, not owned)
	surfaceptr srf;							// porting surface (ref.)
	enum PanelFace face;				// active face of porting surface
	int llport;									// live list number for buffer
	} *portptr;

typedef struct portsuperstruct {
	enum StructCond condition;	// structure condition
	struct simstruct *sim;			// simulation structure
	int maxport;								// maximum number of ports
	int nport;									// actual number of ports
	char **portnames;						// port names
	portptr *portlist;					// list of ports
	} *portssptr;


/*********************************** lattice **********************************/

enum LatticeType {LATTICEnone,LATTICEnsv,LATTICEpde};

typedef struct latticestruct {
  struct latticesuperstruct *latticess;	// lattice superstructure
  char *latticename;          // lattice name (reference, not owned)
  enum LatticeType type;      // type of lattice
  double min[DIMMAX];         // lower spatial boundaries
  double max[DIMMAX];         // upper spatial boundaries
  double dx[DIMMAX];          // lattice lengthscale (subvolume width)
  char btype[DIMMAX];         // boundary type (r)eflective or (p)eriodic
  portptr port;               // interface port (ref.)
  int maxreactions;           // maximum number of reactions
  int nreactions;             // number of reactions
  rxnptr *reactionlist;       // list of reactions
  int *reactionmove;          // 0 or 1 for moving reactions
  int maxsurfaces;           // maximum number of surface
  int nsurfaces;             // number of surface
  surfaceptr *surfacelist;       // list of surface
  int maxspecies;             // maximum number of species
  int nspecies;               // number of species
  int *species_index;					// species indecies
  int *maxmols;               // allocated size of molecule list [lat.species]
  int *nmols;                 // number of individual molecules [lat.species]
  double*** mol_positions;    // molecule positions [lat.species][nmols][dim]
#ifdef OPTION_NSV
  NextSubvolumeMethod* nsv;		// nsv class
  NextSubvolumeMethod* pde;		// pde class
#else
  void *nsv;
  void *pde;
#endif
	} *latticeptr;

typedef struct latticesuperstruct {
  enum StructCond condition;	// structure condition
  struct simstruct *sim;			// simulation structure
  int maxlattice;							// maximum number of lattices
  int nlattice;								// actual number of lattices
  char **latticenames;				// lattice names
  latticeptr *latticelist;		// list of lattices
	} *latticessptr;

/********************************* Filaments ********************************/

typedef struct filamentstruct {
	struct filamentsuperstruct *filss;	// filament superstructure
	char *fname;								// filament name
	double color[4];						// filament color
	double edgepts;							// thickness of edge for drawing
	unsigned int edgestipple[2];	// edge stippling [factor, pattern]
	enum DrawMode drawmode;			// polygon drawing mode
	double shiny;								// shininess
	int nmax;										// number of monomers allocated [1,inf)
	int n;											// number of monomers
	int front;									// front index
	int back;										// back index
	double **px;								// Coords. for monomer ends [nmax+1][3]
	double *pl;									// monomer length [nmax]
	double **pa;								// relative ypr angles [nmax][3]
	double **pd;								// relative dcm [nmax][9]
	double **po;								// absolute monomer orientation [9]
	double *pthk;								// thickness of monomer [nmax], [0,inf)
	double lstd;								// minimum energy monomer length
	double astd[3];							// minimum energy bend angle
	double lk;									// force constant for length
	double ak[3];								// force constant for angle
	double kT;									// thermodynamic temperature, [0,inf)
	double treadrate;						// treadmilling rate constant
	char surf;									// character for surface shape
	double spar[2];							// parameters of surface
	} *filamentptr;

typedef struct filamentsuperstruct {
	enum StructCond condition;	// structure condition
	struct simstruct *sim;			// simulation structure
	int maxfil;									// maximum number of filaments
	int nfil;										// actual number of filaments
	char **fnames;							// filament names
	filamentptr *fillist;				// list of filaments
	} *filamentssptr;


/******************************** BioNetGen *********************************/

typedef struct bngstruct {
	struct bngsuperstruct *bngss; // bng superstructure
	char *bngname;              // bng name
  int bngindex;               // index of this bng structure
	double unirate;							// multiplier for unimolecular rxn rates
	double birate;							// multiplier for bimolecular rxn rates
  
	int maxparams;              // maximum number of numeric parameters
	int nparams;                // actual number of numeric parameters
	char **paramnames;          // names of parameters [index]
  char **paramstrings;        // strings for parameter values [index]
	double *paramvalues;        // actual parameter values [index]

  int maxmonomer;             // maximum number of monomers
  int nmonomer;               // actual number of monomers
  char **monomernames;        // names of monomers [index]
  int *monomercount;          // monomer count work space [index]
	double *monomerdifc;				// diffusion coefficient of monomer [index]
	double *monomerdisplaysize;	// display size of monomer [index]
	double **monomercolor;			// color of monomer [index][RGB]
	enum MolecState *monomerstate; // default monomer state [index]

	int maxbspecies;            // maximum number of bng species
	int nbspecies;              // actual number of bng species
	char **bsplongnames;        // complete bng species names [index]
  char **bspshortnames;       // shortened bng species names [index]
	enum MolecState *bspstate;	// default species state [index]
	char **bspcountstr;         // strings for initial bng species counts [index]
  double *bspcount;           // actual initial bng species counts [index]
  int *spindex;               // smoldyn index of this species [index]

	int maxbrxns;               // maximum number of bng reactions
	int nbrxns;                 // acutal number of bng reactions
	char **brxnreactstr;        // strings for reactants [index]
	char **brxnprodstr;         // strings for products [index]
  char **brxnratestr;         // strings for reaction rates [index]
  int **brxnreact;            // reactant bng species indicies [index][rct]
  int **brxnprod;             // product bng species indicies [index][prd]
  int *brxnorder;             // order of bng reaction [index]
  int *brxnnprod;             // number of products of bng reaction [index]
  rxnptr *brxn;               // pointer to this reaction [index]
} *bngptr;

typedef struct bngsuperstruct {
	enum StructCond condition;	// structure condition
	struct simstruct *sim;			// simulation structure
	char *BNG2path;							// path and name of BNG2.pl executable
	int maxbng;                 // maximum number of bng networks
	int nbng;                   // actual number of bng networks
	char **bngnames;            // names of bng networks
	bngptr *bnglist;            // list of bng networks
} *bngssptr;

/********************************* Graphics ********************************/

#define MAXLIGHTS 8						// must be ? GL_MAX_LIGHTS
enum LightParam {LPambient,LPdiffuse,LPspecular,LPposition,LPon,LPoff,LPauto,LPnone};

typedef struct graphicssuperstruct {
	enum StructCond condition;	// structure condition
	struct simstruct *sim;			// simulation structure
	int graphics;								// graphics: 0=none, 1=opengl, 2=good opengl
	int runmode;								// 0=Smoldyn, 1=Libsmoldyn
	int currentit;							// current number of simulation time steps
	int graphicit;							// number of time steps per graphics update
	unsigned int graphicdelay;	// minimum delay (in ms) for graphics updates
	int tiffit;									// number of time steps per tiff save
	double framepts;						// thickness of frame for graphics
	double gridpts;							// thickness of virtual box grid for graphics
	double framecolor[4];				// frame color [c]
	double gridcolor[4];				// grid color [c]
	double backcolor[4];				// background color [c]
	double textcolor[4];				// text color [c]
	int maxtextitems;						// allocated size of item list
	int ntextitems;							// actual size of item list
	char **textitems;						// items to display with text [item]
	enum LightParam roomstate;	// on, off, or auto (on)
	double ambiroom[4];					// global ambient light [c]
	enum LightParam lightstate[MAXLIGHTS];	// on, off, or auto (off) [lt]
	double ambilight[MAXLIGHTS][4];		// ambient light color [lt][c]
	double difflight[MAXLIGHTS][4];		// diffuse light color [lt][c]
	double speclight[MAXLIGHTS][4];		// specular light color [lt][c]
	double lightpos[MAXLIGHTS][3];		// light positions [lt][d]
	} *graphicsssptr;

/******************************** Simulation *******************************/

#define ETMAX 10
enum SmolStruct {SSmolec,SSwall,SSrxn,SSsurf,SSbox,SScmpt,SSport,SSfilament,SScmd,SSsim,SScheck,SSall,SSnone};
enum EventType {ETwall,ETsurf,ETdesorb,ETrxn0,ETrxn1,ETrxn2intra,ETrxn2inter,ETrxn2wrap,ETimport,ETexport};

typedef int (*diffusefnptr)(struct simstruct *);
typedef int (*surfaceboundfnptr)(struct simstruct *,int);	
typedef int (*surfacecollisionsfnptr)(struct simstruct *,int,int);
// typedef int (*surfacecollisionsfnptr)(struct simstruct *, int, int, int); 		// cplx
typedef int (*assignmols2boxesfnptr)(struct simstruct *,int,int);
typedef int (*zeroreactfnptr)(struct simstruct *);
typedef int (*unimolreactfnptr)(struct simstruct *);
typedef int (*bimolreactfnptr)(struct simstruct *,int);
typedef int (*checkwallsfnptr)(struct simstruct *,int,int,boxptr);


typedef struct simstruct {
	enum StructCond condition;	// structure condition
	char *vfile;
	FILE *events;							// record reaction events	
	FILE *logfile;							// file to send output
	char *filepath;							// configuration file path
	char *filename;							// configuration file name
	char *flags;								// command-line options from user
	time_t clockstt;						// clock starting time of simulation
	double elapsedtime;					// elapsed time of simulation
	long int randseed;					// random number generator seed
	int eventcount[ETMAX];			// counter for simulation events
	int dim;										// dimensionality of space.
	double accur;								// accuracy, on scale from 0 to 10
	double time;								// current time in simulation
	double tmin;								// simulation start time
	double tmax;								// simulation end time
	double tbreak;							// simulation break time
	double dt;									// simulation time step
	rxnssptr rxnss[MAXORDER+MAXORDER-1];		// reaction superstructures, used to be rxnss[MAXORDER]
	molssptr mols;							// molecule superstructure
	wallptr *wlist;							// list of walls
	surfacessptr srfss;					// surface superstructure
	boxssptr boxs;							// box superstructure
	compartssptr cmptss;				// compartment superstructure
	portssptr portss;						// port superstructure
	latticessptr latticess;						// port superstructure
	bngssptr bngss;							// bionetget superstructure
	filamentssptr filss;				// filament superstructure
	void* cmds;									// command superstructure, output_root, output_filename, etc.
	graphicsssptr graphss;						// graphics superstructure
	diffusefnptr diffusefn;											// function for molecule diffusion
	surfaceboundfnptr surfaceboundfn;				// function for surface-bound molecules
	surfacecollisionsfnptr surfacecollisionsfn; 	// function for surface collisons
	assignmols2boxesfnptr assignmols2boxesfn;		// function that assigns molecs to boxes
	zeroreactfnptr zeroreactfn;									// function for zero order reactions
	unimolreactfnptr unimolreactfn;								// function for first order reactions
	bimolreactfnptr bimolreactfn;								// function for second order reactions
	checkwallsfnptr checkwallsfn;								// function for molecule collisions with walls
	int multibinding;
	interfaceptr interface;	
	gsl_rng *r;

#ifdef OPTION_VCELL
	VolumeSamplesPtr volumeSamplesPtr;
	ValueProviderFactory* valueProviderFactory;
	AbstractMesh* mesh;
#endif
	} *simptr;


/*********************************** VCell *********************************/


#ifdef OPTION_VCELL

	struct CompartmentIdentifierPair {
		char name[128];
		unsigned char pixel;//the compartmentID
	};

	typedef struct VolumeSamples {
		int num[3];//number of mesh points in X, Y,Z
		double size[3];//actual size in X, Y, Z (e.g. in micron)
		double origin[3];//origin of the X, Y, Z
		unsigned char* volsamples;//compartmentID for each mesh point center
		int nCmptIDPair; // number of compartments
		CompartmentIdentifierPair* compartmentIDPairPtr;//ID vs. comptName pairs.
	}* VolumeSamplesPtr;

	class ValueProviderFactory;
	class AbstractMesh;


	class ValueProvider {
	public:
		virtual ~ValueProvider(){};
		virtual double getConstantValue()=0;
		virtual double getValue(double t, double x, double y, double z, rxnptr rxn)=0;
		virtual double getValue(double t, double x, double y, double z, rxnptr rxn, char* panelName)=0;
		virtual double getValue(double t, double x, double y, double z, surfactionptr actiondetails, char* panelName)=0;
	};

	class ValueProviderFactory {
	public:
		virtual ~ValueProviderFactory(){};
		virtual ValueProvider* createValueProvider(string& rateExp)=0;
		void setSimptr(simptr sim){this->sim = sim;}
		simptr getSimptr(){return this->sim;}
	private:
		simptr sim;
	};

	class AbstractMesh{
	public:
		virtual ~AbstractMesh(){};
		virtual void getCenterCoordinates(int volIndex, double* coords)=0;
		virtual void getDeltaXYZ(double* delta)=0;
		virtual void getNumXYZ(int* num)=0;
	};

#endif


#endif
