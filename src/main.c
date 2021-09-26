#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <fftw3.h>
#include "iniparser.h"

/***************** INIPARSER ***************/
int  parse_ini_file(char * ini_name);

/********* UNIVERSAL CONSTANTS *************/
const double EPS = 8.85418782E-12;    // Vacuum permittivity
const double K = 1.38065E-23;            // Boltzmann Constant
const double EV_TO_K = 11604.52;         // Conversion of EV to K
const double PI = 3.14159265359;



/*************** USER CHOICE ******************/
short int solverType; // Poisson SOlver
short int loadType;   // Particle and velocity sampling
short int TSI;        // Two stream instability
/************ VARIABLE INITIALIZATION **********/
int nParticlesI;      // Number of simulation ions
int nParticlesE; // Number of simulation electrons

int numxCells;             // Total number of cells alonx x
int numyCells;             // Total number of cells alonx y
int nTimeSteps;          // Total time steps (default)
int dumpPeriod;         // Data dump period
double massI;           // Ion mass
double massE;            // Electron mass
double chargeE;         // Electron charge
// int probLoc;  //VDF end location
double driftE;          // TODO
double driftI;
// double EPS;
double offsetE;
double offsetI;
int E_per_cell;
int I_per_cell;


/* Simulation Parameters*/
double density;  // Plasma Density
double ion_spwt;
double electron_spwt;
double stepSize;          // Cell Spacing
double timeStep;        // Time steps

double TE; // electron temperature in eV
double TI;  // ion temperature in eV
double v_te; // electron thermal velocity

double Lambda_D;
double omega_pe;
double Bx,By,Bz;
double B[3];

// double massE 9.10938215e-31		// kg, electron mass
// #define chargeE 1.602176565e-19		// C, electron charge
// #define massI 1.660538921e-27	// kg, atomic mass unit
// #define TE 1.0 // in eV
// #define TI 0.026 // in eV



// double Lambda_D = sqrt((EPS_0*TE*EV_TO_K)/(chargeE*chargeE*density));
// double omega_pe = sqrt((chargeE*chargeE*density)/(massE*EPS_0));


/* Data structure to hold domain information*/
struct Domain
{
	int nix;				/*number of nodes along x*/
  int niy;				/*number of nodes along y*/
  int niz;
	double x0;			/*mesh origin*/
  double y0;
  double z0;
  double dx;			/*cell spacing*/
  double dy;
  double dz;
  double xl;			/*domain length*/
  double yl;
  double zl;
  double xmax;		/*domain max position*/
  double ymax;
  double zmax;

/*data structures*/
	double *phi;		/*potential*/
	double *efx;			/*electric field*/
  double *efy;
  // double *efz;
	double *rho;		/*charge density*/
	double *ndi;		/*ion density*/
	double *nde;		/*electron density*/
  fftw_complex *phik;
  fftw_complex *rhok;
};

/* Data structure for particle storage **/
struct Particle
{
	double x;			/*position*/
	double vx;			/*velocity*/
  double y;
  double vy;
  double z;
  double vz;
};

/* Data structure to hold species information*/
struct Species
{
	double mass;			/*particle mass in kg*/
	double charge;			/*particle charge in Coulomb*/
	double spwt;			/*species specific weight*/
  double temp;

	int np;					/*number of particles*/
	int np_alloc;			/*size of the allocated data array*/
  struct Particle *part;			/*array holding particles*/
};

/** FUNCTION PROTOTYPES **/
void ScatterSpecies(struct Species *species, double *nd);
void ComputeRho(double *rho, struct Species *ions, struct Species *electrons);
bool SolvePotential(double *phi, double *rho);
bool PseudospectralSolvePotential(double *phi, double *rho, fftw_complex *phik, fftw_complex *rhok);
void ComputeEF(double *phi, double *efx, double *efy);
void PushSpecies(struct Species *species, double *efx, double *efy);
void BorisPushSpecies(struct Species *species, double *efx, double *efy, double B[]);
double* CrossProduct(double v1[3], double v2[3]);
void RewindSpecies(struct Species *species, double *efx, double *efy);
void AddParticle(struct Species *species, double x, double vx, double y, double vy, double z, double vz);
double XtoL(double pos);
void scatter(double lx, double ly, double value, double *field);
double gather(double lx, double ly, double *field);
void WriteResults(int ts);
double ComputeKE(struct Species *species);	/*computes kinetic energy*/
double rnd();
double sampleVel(double T, double mass);
void Write_Particle(FILE *file, int ts, struct Species *species);
void Write_density(FILE *file, int ts);

/* GLOBAL VARIABLES */
struct Domain domain;
FILE *file_res;
FILE *file_pot;
FILE *f1;
FILE *f2;
FILE *f3;
FILE *f4;


/*Parsing Input file*/
int parse_ini_file(char * ini_name)
{
    dictionary  *   ini ;

    ini = iniparser_load(ini_name);
    if (ini==NULL) {
        fprintf(stderr, "cannot parse file: %s\n", ini_name);
        return -1 ;
    }
    // iniparser_dump(ini, stderr); // Comment out to fix issues with iniparser

    /*Get Simulation Parameters */
    nTimeSteps  = iniparser_getint(ini,"time:nTimeSteps",-1);
    timeStep    = iniparser_getdouble(ini,"time:timeStep",-1.0);
    stepSize    = iniparser_getdouble(ini,"grid:stepSize",-1.0);
    numxCells   = iniparser_getint(ini,"grid:numxCells",-1);
    numyCells   = iniparser_getint(ini,"grid:numyCells",-1);

    /* SPECIES INFO */
    nParticlesI = iniparser_getint(ini,"population:nParticlesI",-1);
    nParticlesE = iniparser_getint(ini,"population:nParticlesE",-1);
    massI    = iniparser_getdouble(ini,"population:massI",-1.0);
    massE    = iniparser_getdouble(ini,"population:massE",-1.0);
    chargeE  = iniparser_getdouble(ini,"population:chargeE",-1.0);
    density  = iniparser_getdouble(ini,"population:density",-1.0);
    TE     = iniparser_getdouble(ini,"population:TE",-1.0);
    TI     = iniparser_getdouble(ini,"population:TI",-1.0);
    driftE   = iniparser_getdouble(ini,"population:driftE",-1.0);
    driftI   = iniparser_getdouble(ini,"population:driftI",-1.0);
    offsetE  = iniparser_getdouble(ini,"population:offsetE",-1.0);
    offsetI  = iniparser_getdouble(ini,"population:offsetI",-1.0);

    Bx  = iniparser_getdouble(ini,"field:Bx",-1.0);
    By  = iniparser_getdouble(ini,"field:By",-1.0);
    Bz  = iniparser_getdouble(ini,"field:Bz",-1.0);

    B[0]  = Bx;
    B[1]  = By;
    B[2]  = Bz;

    /* DIAGNOSTICS */
    // probLoc = iniparser_getint(ini,"diagnostics:probLoc",-1);
    dumpPeriod            = iniparser_getint(ini,"diagnostics:dumpPeriod",-1);

    /* USER CHOICE */
    solverType            = iniparser_getint(ini,"solver:solverType",-1);
    loadType              = iniparser_getint(ini,"population:loadType",-1);
    TSI                   = iniparser_getint(ini,"population:TSI",-1);

    /* Normalization */ //TO BE ADDED AS A SEPERATE FUNCTION
    // EPS             = EPS_un;//EPS_un;
    v_te     = sqrt(2*K*TE*EV_TO_K/massE); // electron thermal vel
    omega_pe = sqrt((chargeE*chargeE*density)/(massE*EPS));
    Lambda_D = sqrt((EPS*K*TE*EV_TO_K)/(density*chargeE*chargeE));


    printf("********** IMPORTANT PLASMA QUANTITIES ***********\n");
    printf("omega_pe: %f\n",omega_pe);
    printf("Lambda_D: %f\n",Lambda_D);

    printf("*************** Input Sanity Check ***************\n");
    bool SFLAG = true;
    if (stepSize >= 1) {
      printf("ERROR, stepSize is bigger than Debye length.\n");
      SFLAG = false;
    }
    if (timeStep > 0.5) {
      printf("ERROR, timeStep is too big. The recommended value: <%f s\n",(0.5/omega_pe));
      SFLAG = false;
    }
    if (solverType != 1 && solverType != 2) {
      printf("ERROR, Wrong Solver Type. The recommended value: 1 or 2\n");
      printf("solverType: %d",solverType);
      SFLAG = false;
    }
    if (loadType != 1 && loadType != 2 && loadType !=3) {
      printf("ERROR, Wrong Load Type. The recommended value: 1 or 2\n");
      printf("loadType: %d\n",loadType);
      SFLAG = false;
    }
    if (SFLAG==true) {
      printf("STATUS, Input parameters are compatible.\n");
    }
    else {
      printf("ERROR, Input parameters are incompatible.\n");
      exit (EXIT_FAILURE);
    }


    iniparser_freedict(ini);
    return 0;
}






/* --------- main -------------*/

int main(int argc, char *argv[])
{
  /************* INIPARSER ********************/
  if(argc<2) {
    printf("ERROR, at least one argument expected (the input file).\n");
    exit (EXIT_FAILURE);
  }
  parse_ini_file(argv[1]);


	int i,p;
	int ts;						// time step
  srand(time(NULL));
	clock_t start = clock();	// grab starting clock time

	/*---1) initialize domain---*/
	domain.nix = numxCells+1;			// number of nodes
	domain.dx = stepSize;				// cell spacing
	domain.x0 = 0;					// origin
	domain.xl = (domain.nix-1)*domain.dx;//domain length
	domain.xmax = domain.x0+domain.xl;	//max position

  domain.niy = numxCells+1;			// number of nodes
	domain.dy = stepSize;				// cell spacing
	domain.y0 = 0;					// origin
	domain.yl = (domain.niy-1)*domain.dy;//domain length
	domain.ymax = domain.y0+domain.yl;	//max position

  domain.niz = numxCells+1;			// number of nodes
	domain.dz = stepSize;				// cell spacing
	domain.z0 = 0;					// origin
	domain.zl = (domain.niz-1)*domain.dz;//domain length
	domain.zmax = domain.z0+domain.zl;	//max position

	/*allocate data structures, remember these need to be cleared first in C++*/
   domain.phi = (double*) malloc((domain.nix*domain.niy*domain.niz) * sizeof(double));
   domain.rho = (double*) malloc((domain.nix*domain.niy*domain.niz) * sizeof(double));
   domain.efx = (double*) malloc((domain.nix*domain.niy*domain.niz) * sizeof(double));
   domain.efy = (double*) malloc((domain.nix*domain.niy*domain.niz) * sizeof(double));
   // domain.efz = (double*) malloc((domain.nix*domain.niy*domain.niz) * sizeof(double));
   domain.nde = (double*) malloc((domain.nix*domain.niy) * sizeof(double));
   domain.ndi = (double*) malloc((domain.nix*domain.niy) * sizeof(double));
   // domain.rhok = (fftw_complex*) malloc((domain.nix*(domain.niy/2+1)) * sizeof(fftw_complex));
   // domain.phik = (fftw_complex*) malloc((domain.nix*(domain.niy/2+1)) * sizeof(fftw_complex));

	// /*save pointers so we can write phi instead of domain.phi*/
	double *phi = domain.phi;
	double *rho = domain.rho;
	double *efx = domain.efx;
  double *efy = domain.efy;
  // double *efz = domain.efz;
	double *nde = domain.nde;
	double *ndi = domain.ndi;
  // fftw_complex *phik = domain.phik;
  // fftw_complex *rhok = domain.rhok;

  struct Species ions;
	struct Species electrons;

	/*set material data*/
	ions.mass = massI;
	ions.charge = chargeE;
	ions.spwt = density*(domain.xl*domain.yl)/nParticlesI;
  ions.temp = TI;
	ions.np = 0;
	ions.np_alloc = nParticlesI;
	ions.part = (struct Particle*)malloc(nParticlesI * sizeof(struct Particle));

	electrons.mass = massE;	// electrons
	electrons.charge = -chargeE;
	electrons.spwt = density*(domain.xl*domain.yl)/nParticlesE;
  electrons.temp = TE;
	electrons.np = 0;
	electrons.np_alloc = nParticlesE;
	electrons.part = (struct Particle*)malloc(nParticlesE * sizeof(struct Particle));

  double delta_ions = domain.xl/nParticlesI;
	for (p=0;p<nParticlesI;p++)
	{
    double x = domain.x0 + ((rnd()+rnd()+rnd())/3)*(domain.nix-1)*domain.dx;
		double vx = sampleVel(ions.temp*EV_TO_K, ions.mass);	/*initially zero x velocity*/
    double y = domain.y0 + ((rnd()+rnd()+rnd())/3)*(domain.niy-1)*domain.dy;
		double vy = sampleVel(ions.temp*EV_TO_K, ions.mass);	/*initially zero y velocity*/
    double z = 0.0;
    double vz = 0.0;
    if (x>domain.xl) x = x - domain.xl;
    if (x<0) x = x + domain.xl;
    if (y>domain.yl) y = y - domain.yl;
    if (y<0) y = y + domain.yl;
		AddParticle(&ions,x,vx,y,vy,z,vz);
	}

	/*now do the same for electrons*/
	double delta_electrons = domain.xl/nParticlesE;
	for (p=0;p<nParticlesE;p++)
	{
		double x = domain.x0 + ((rnd()+rnd()+rnd())/3)*(domain.nix-1)*domain.dx;
		double vx = sampleVel(electrons.temp*EV_TO_K, electrons.mass);
    double y = domain.x0 + ((rnd()+rnd()+rnd())/3)*(domain.niy-1)*domain.dy;
		double vy = sampleVel(electrons.temp*EV_TO_K, electrons.mass);
    double z = 0.0;
    double vz = 0.0;
    if (x>domain.xl) x = x - domain.xl;
    if (x<0) x = x + domain.xl;
    if (y>domain.yl) y = y - domain.yl;
    if (y<0) y = y + domain.yl;
		AddParticle(&electrons,x,vx,y,vy,z,vz);
	}

  ScatterSpecies(&ions,ndi);
	ScatterSpecies(&electrons,nde);

  ComputeRho(rho, &ions, &electrons);
  // // PseudospectralSolvePotential(phi,rho,phik,rhok);
	SolvePotential(phi,rho);
	ComputeEF(phi,efx,efy);
  //
  RewindSpecies(&ions,efx,efy);
	RewindSpecies(&electrons,efx,efy);

  file_res = fopen("results.dat","w");
	fprintf(file_res,"VARIABLES = x y nde ndi rho phi efx efy\n");
	file_pot = fopen("potential.dat","w");
	WriteResults(0);

  char Name[50];


  for (ts = 1; ts<=nTimeSteps; ts++)
	{
		//compute number density
		ScatterSpecies(&ions,ndi);
		ScatterSpecies(&electrons,nde);

		ComputeRho(rho,&ions, &electrons);
    // PseudospectralSolvePotential(phi,rho,phik,rhok);
		SolvePotential(phi,rho);
		ComputeEF(phi,efx,efy);

		//move particles
    // PushSpecies(&electrons,efx,efy);
		// PushSpecies(&ions,efx,efy);
    BorisPushSpecies(&electrons,efx,efy,B);
		BorisPushSpecies(&ions,efx,efy,B);


		//write diagnostics
		if (ts%10==0)
		{
			double max_phi = phi[(0*domain.niy+0)*domain.niz+0];
      for (int k=0;k<domain.niz;k++)
      for (int j=0;j<domain.niy;j++)
			for (int i=0;i<domain.nix;i++)
      {
        if (phi[(i*domain.niy+j)*domain.niz+k]>max_phi) max_phi=phi[(i*domain.niy+j)*domain.niz+k];
      }
    //
    //
		// 	//compute kinetic energy
			double ke_ions = ComputeKE(&ions);
			double ke_electrons = ComputeKE(&electrons);

			printf("TS: %i \t delta phi: %.3g\n",ts,max_phi-phi[0]);

    // for(int i=0; i<domain.nix; i++)
    // for(int j=0; j<domain.niy; j++)
    // {
    //   fprintf(file_pot, "%g \t %g\n", timeStep*ts, phi[i]);
    // }

      sprintf(Name,"i%d.dat",ts);
      f1 = fopen(Name, "w");

      sprintf(Name,"e%d.dat",ts);
      f2 = fopen(Name, "w");

      sprintf(Name,"data%d.dat",ts);
      f3 = fopen(Name, "w");

      // sprintf(Name,"nde%d.dat",ts);
      // f4 = fopen(Name, "w");

      Write_Particle(f1,ts,&ions);
      Write_Particle(f2,ts,&electrons);
      Write_density(f3,ts);

		}
    //
		// //save animation data
		// if (ts%100==0) WriteResults(ts);
    // WriteResults(ts);
	}

  fclose(file_res);
	fclose(file_pot);

  free(phi);
	free(rho);
  free(efx);
  free(efy);
	free(nde);
  free(ndi);
	// /*free particles*/
	free(ions.part);
	free(electrons.part);

	return 0;
}


/*adds new particle to the species*/
void AddParticle(struct Species *species, double x, double vx, double y, double vy, double z, double vz)
{
	/*abort if we ran out of space to store this particle*/
	if (species->np > species->np_alloc-1)
	{
		printf("Too many particles!\n"); exit(-1);
	}

	/*store position and velocity of this particle*/
	species->part[species->np].x = x;
	species->part[species->np].vx = vx;
  species->part[species->np].y = y;
	species->part[species->np].vy = vy;
  species->part[species->np].z = z;
	species->part[species->np].vz = vz;

  // printf("%g\t%g\t%g\t%g\t%g\t%g\n", x,vx,y,vy,z,vz);
	/*increment particle counter*/
	species->np++;
}

double sampleVel(double T, double mass)
{
  double v_th = sqrt(2*K*T/mass);
  double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd())/3;
  return vt/v_te;
}


double XtoL(double pos)
{
	double li = (pos-domain.x0)/domain.dx;
	return li;
}
//
/* scatters scalar value onto a field at logical coordinate lc*/
void scatter(double lx, double ly, double value, double *field)
{
  int i = (int) lx;
  int j = (int) ly;
  double di = lx - i;
  double dj = ly - j;

  field[i*domain.niy+j] += value*(1-di)*(1-dj);
  field[(i+1)*domain.niy+j] += value*(di)*(1-dj);
  field[i*domain.niy+j+1] += value*(1-di)*(dj);
  field[(i+1)*domain.niy+j+1] += value*(di)*(dj);
  // printf("%g \t%g \t%g\n",lx,ly,field[i*domain.niy+j]);
}

double gather(double lx, double ly, double *field)
{
	int i = (int)lx;
	double di = lx-i;

  int j = (int)ly;
  double dj = ly - j;

  double val = field[i*domain.niy+j]*(1-di)*(1-dj)+field[(i+1)*domain.niy+j]*di*(1-dj)+field[i*domain.niy+j+1]*(1-di)*dj+field[(i+1)*domain.niy+j+1]*di*dj;
  return val;
}


/*scatter particles of species to the mesh*/
void ScatterSpecies(struct Species *species, double *den)
{
	/*initialize densities to zero*/
	memset(den,0,sizeof(double)*domain.nix*domain.niy);

	/*scatter particles to the mesh*/
	for (int p=0;p<species->np;p++)
		{
			double lcx = XtoL(species->part[p].x);
      double lcy = XtoL(species->part[p].y);
			scatter(lcx, lcy, species->spwt, den);
		}

	/*apply periodic boundaries*/
  for(int j=0; j<domain.niy; j++)
    {
    den[0*domain.niy+j] += den[(domain.nix-1)*domain.niy+j];
    den[(domain.nix-1)*domain.niy+j]=den[0*domain.niy+j];
    }

  for(int i=0; i<domain.nix; i++)
    {
    den[i*domain.niy+0] += den[i*domain.niy+domain.niy-1];
    den[i*domain.niy+domain.niy-1] = den[i*domain.niy+0];
    }
  for(int i=0; i<domain.nix; i++)
  for(int j=0; j<domain.niy; j++)
  {
      den[i*domain.niy+j] /=domain.dx*domain.dy;
  }

  for(int i=0; i<domain.nix; i++)
  for(int j=0; j<domain.niy; j++)
  {
      den[i*domain.niy+j] /=density;
      // printf("%g\t%g\t%g\n",i*domain.dx, j*domain.dy, den[i*domain.niy+j]);
  }
}

/*computes charge density by adding ion and electron data*/
void ComputeRho(double *rho, struct Species *ions, struct Species *electrons)
{
	// double *rho = domain.rho;

    // for(int k=1; k<domain.niz-1; k++)
    for(int j=1; j<domain.niy-1; j++)
    for(int i=1; i<domain.nix-1; i++)
        {
        rho[(i*domain.niy+j)*domain.niz+0]=(-domain.ndi[i*domain.niy+j] + domain.nde[i*domain.niy+j]);
        }
    // for(int k=1; k<domain.niz; k++)
    for(int j=0; j<domain.niy; j++)
      {
      rho[(0*domain.niy+j)*domain.niz+0] += rho[((domain.nix-1)*domain.niy+j)*domain.niz+0];
      rho[((domain.nix-1)*domain.niy+j)*domain.niz+0]=rho[(0*domain.niy+j)*domain.niz+0];
      }

    for(int i=0; i<domain.nix; i++)
      {
      rho[(i*domain.niy+0)*domain.niz+0] += rho[(i*domain.niy+domain.niy-1)*domain.niz+0];
      rho[(i*domain.niy+domain.niy-1)*domain.niz+0] = rho[(i*domain.niy+0)*domain.niz+0];
      }
    //
    // for(int j=0; j<domain.niy; j++)
    // for(int i=0; i<domain.nix; i++)
    // {
    //   printf("%g\t%g\t%g\n",i*domain.dx, j*domain.dy, rho[0+domain.niz*(j+domain.niy*i)]);
    // }
}

/* solves potential using the Gauss Seidel Method, returns true if converged
periodic boundaries version*/
bool SolvePotential(double *phi, double *rho)
{
	double L2;
	double dx2 = domain.dx*domain.dx;	/*precompute*/
  double dy2 = domain.dy*domain.dy;

	/*initialize to zero, this is only needed because of periodic b.c. without a dirichlet*/
	memset(phi, 0, sizeof(double)*domain.nix*domain.niy);

	/*solve potential, identical to lesson 2*/
	for (int solver_it=0;solver_it<50000;solver_it++)
	{
		/*Gauss Seidel method, phi[i-1]-2*phi[i]+phi[i+1] = -dx^2*rho[i]/eps_0*/
		for (int k=0;k<domain.niz;k++)
    for (int j=0;j<domain.niy;j++)
    for (int i=0;i<domain.nix;i++)
		{
			int p = i-1; if(p<0) p=domain.nix-2;
			int q = i+1; if(q>domain.nix-1) q=1;
      int r = j-1; if(r<0) r=domain.niy-2;
      int s = j+1; if(s>domain.niy-1) s=1;
			// phi[i] = 0.5*(phi[im] + phi[ip] + dx2*rho[i]/EPS_0);
      double g = 0.5*(1/((1/dx2)+(1/dy2)))*( ((phi[(p*domain.niy+j)*domain.niz+k]+phi[(q*domain.niy+j)*domain.niz+k])/dx2) + ((phi[(i*domain.niy+r)*domain.niz+k]+phi[(i*domain.niy+s)*domain.niz+k])/dy2) + (rho[(i*domain.niy+j)*domain.niz+k]) );
      phi[(i*domain.niy+j)*domain.niz+k] = phi[(i*domain.niy+j)*domain.niz+k] + 1.4*(g - phi[(i*domain.niy+j)*domain.niz+k]);
    }

		/*check for convergence*/
		if (solver_it%25==0)
		{
			double sum = 0;
      for (int k=0;k<domain.niz;k++)
      for (int j=0;j<domain.niy;j++)
      for (int i=0;i<domain.nix;i++)
			{
        int p = i-1; if(p<0) p=domain.nix-2;
  			int q = i+1; if(q>domain.nix-1) q=1;
        int r = j-1; if(r<0) r=domain.niy-2;
        int s = j+1; if(s>domain.niy-1) s=1;
				// double R = -rho[i]/EPS_0 - (phi[im] - 2*phi[i] + phi[ip])/dx2;
        double R = 0.25*(phi[(p*domain.niy+j)*domain.niz+k]+phi[(q*domain.niy+j)*domain.niz+k]+phi[(i*domain.niy+r)*domain.niz+k]+phi[(i*domain.niy+s)*domain.niz+k]+(dx2*rho[(i*domain.niy+j)*domain.niz+k]))-phi[(i*domain.niy+j)*domain.niz+k];
				sum+=R*R;
			}
			L2 = sqrt(sum)/(domain.nix*domain.niy);
			if (L2<1e-5) {return true;}
		}

    // for(int j=0; j<domain.niy; j++)
    // for(int i=0; i<domain.nix; i++)
    // {
    //   printf("%g\t%g\t%g\n",i*domain.dx, j*domain.dy, phi[0+domain.niz*(j+domain.niy*i)]);
    // }
	}
	printf("Gauss-Seidel solver failed to converge, L2=%.3g!\n",L2);
	return false;
}

// bool PseudospectralSolvePotential(double *phi, double *rho, fftw_complex *phik, fftw_complex *rhok)
// {
//   int Nx = domain.nix, Ny = domain.niy, Nh = domain.niy/2+1;
//   int i,j;
//   fftw_plan f,b;
//   double Lx = domain.xl, Ly = domain.yl;
//   double kx, ky;
//   fftw_complex *rhok_dum, *phik_dum;
//
//   rhok_dum = (fftw_complex*) malloc(Nx*Ny * sizeof(fftw_complex));
//   phik_dum = (fftw_complex*) malloc(Nx*Ny * sizeof(fftw_complex));
//
//   f = fftw_plan_dft_r2c_2d(Nx,Ny,&rho[0*Ny+0],&rhok[0*Nh+0],FFTW_ESTIMATE);
//   fftw_execute(f);
//   fftw_destroy_plan(f);
//   fftw_cleanup();
//
//   for(i=0; i<Nx; i++)
//   for(j=0; j<Nh; j++)
//   {
//     rhok_dum[i*Nh+j][0] = rhok[i*Nh+j][0];
//     rhok_dum[i*Nh+j][1] = rhok[i*Nh+j][1];
//   }
//
//   for(j=0; j<Nh; j++)
//   {
//     ky = 2.0*PI*j/Ly;
//     for(i=0; i<Nx/2; i++)
//     {
//       kx = 2.0*PI*i/Lx;
//       phik[i*Nh+j][0] = rhok[i*Nh+j][0]/(kx*kx + ky*ky);
//       phik[i*Nh+j][1] = rhok[i*Nh+j][1]/(kx*kx + ky*ky);
//     }
//
//     for(i=Nx/2+1; i<Nx; i++)
//     {
//       kx = 2.0*PI*(Nx-i)/Lx;
//       phik[i*Nh+j][0] = rhok[i*Nh+j][0]/(kx*kx + ky*ky);
//       phik[i*Nh+j][1] = rhok[i*Nh+j][1]/(kx*kx + ky*ky);
//     }
//   }
//
//   phik[0*Nh+0][0] = 0.0;
//   phik[0*Nh+0][1] = 0.0;
//
//
//   b =fftw_plan_dft_c2r_2d(Nx,Ny,&phik[0*Nh+0],&phi[0*Ny+0],FFTW_ESTIMATE);
//   fftw_execute(b);
//   fftw_destroy_plan(b);
//   fftw_cleanup();
//
//   for(i=0; i<Nx; i++)
//   for(j=0; j<Ny; j++)
//   {
//     phi[i*Ny+j] = phi[i*Ny+j]/(Nx*Ny);
//     // printf("%f\t%f\n",rho[i*Ny+j],phi[i*Ny+j]/(Nx*Ny));
//   }
//
//   fftw_free(phik_dum);
//   fftw_free(rhok_dum);
//
//   return true;
// }

/* computes electric field by differentiating potential*/
void ComputeEF(double *phi, double *efx, double *efy)
{
  memset(efx,0,sizeof(double)*domain.nix*domain.niy*domain.niz);
  memset(efy,0,sizeof(double)*domain.nix*domain.niy*domain.niz);

  for(int k=0; k<domain.niz; k++)
  {
    for(int j=1; j<domain.niy-1; j++)
    for(int i=1; i<domain.nix-1; i++)
    {
      efx[(i*domain.niy+j)*domain.niz+k] = (phi[((i-1)*domain.niy+j)*domain.niz+k]-phi[((i+1)*domain.niy+j)*domain.niz+k])/(2*domain.dx);
      efy[(i*domain.niy+j)*domain.niz+k] = (phi[(i*domain.niy+j-1)*domain.niz+k]-phi[(i*domain.niy+j+1)*domain.niz+k])/(2*domain.dy);
    }

    /*Apply one sided difference at the boundary nodes*/
    for(int j=0; j<domain.niy; j++)                                    // calculate electric field at nodes
   {
      // efx[0*domain.niy+j] = -(phi[1*domain.niy+j]-phi[(domain.nix-2)*domain.niy+j])/(2*domain.dx);
      // efx[(domain.nix-1)*domain.niy+j] = efx[0*domain.niy+j];
      efx[(0*domain.niy+j)*domain.niz+k] = -(phi[(1*domain.niy+j)*domain.niz+k] - phi[(0*domain.niy+j)*domain.niz+k])/(2*domain.dx);
      efx[((domain.nix-1)*domain.niy+j)*domain.niz+k] = -(phi[((domain.nix-1)*domain.niy+j)*domain.niz+k] - phi[((domain.nix-2)*domain.niy+j)*domain.niz+k])/(2*domain.dx);

   }

   for(int i=0; i<domain.nix; i++)
   {
      // efy[i*domain.niy+0] = -(phi[i*domain.niy+1]-phi[i*domain.niy+domain.niy-2])/(2*domain.dy);
      // efy[i*domain.niy+domain.niy-1] = efy[i*domain.niy+0];

      efy[(i*domain.niy+0)*domain.niz+k] = -(phi[(i*domain.niy+1)*domain.niz+k] - phi[(i*domain.niy+0)*domain.niz+k])/(2*domain.dy);
      efy[(i*domain.niy+domain.niy-1)*domain.niz+k] = -(phi[(i*domain.niy+domain.niy-1)*domain.niz+k] - phi[(i*domain.niy+domain.niy-2)*domain.niz+k])/(2*domain.dy);
   }
  }

  // for(int j=0; j<domain.niy; j++)
  // for(int i=0; i<domain.nix; i++)
  // {
  //   printf("%g\t%g\t%g\t%g\n",i*domain.dx, j*domain.dy, efx[0+domain.niz*(j+domain.niy*i)], efy[0+domain.niz*(j+domain.niy*i)]);
  // }

}

void RewindSpecies(struct Species *species, double *efx, double *efy)
{
	/*precompute q/m*/
	double qm = species->charge / species->mass;

	for (int p=0;p<species->np;p++)
	{
		struct Particle *part = &species->part[p];

		/*compute particle node position*/
		double lcx = XtoL(part->x);
    double lcy = XtoL(part->y);

		/*gather electric field onto particle position*/
		double part_efx = gather(lcx,lcy,efx);
    double part_efy = gather(lcy,lcy,efy);

    double wl = Lambda_D*omega_pe;

		/*advance velocity*/
		part->vx -= 0.5*(1/(wl*wl))*(qm*TE)*part_efx*timeStep;
    part->vy -= 0.5*(1/(wl*wl))*(qm*TE)*part_efy*timeStep;
    part->vz = 0.0;
	}
}

/* moves particles of a single species*/
void PushSpecies(struct Species *species, double *efx, double *efy)
{
	/*precompute q/m*/
	double qm = species->charge / species->mass;

	/*loop over particles*/
	for (int p=0;p<species->np;p++)
	{
		/*grab pointer to this particle*/
		struct Particle *part = &species->part[p];

		/*compute particle node position*/
		double lcx = XtoL(part->x);
    double lcy = XtoL(part->y);

		/*gather electric field onto particle position*/
		double part_efx = gather(lcx,lcy,efx);
    double part_efy = gather(lcx,lcy,efy);

    double wl = Lambda_D*omega_pe;

		/*advance velocity*/
		part->vx += 0.5*(1/(wl*wl))*(qm*TE)*part_efx*timeStep;
    part->vy += 0.5*(1/(wl*wl))*(qm*TE)*part_efy*timeStep;
    part->vz = 0.0;

		/*advance position*/
		part->x += timeStep*part->vx;
    part->y += timeStep*part->vy;
    part->z = 0.0;

		/*apply periodic boundaries*/
		if (part->x < domain.x0)	part->x+=domain.xl;
		if (part->x > domain.xmax)	part->x-=domain.xl;
    if (part->y < domain.y0)	part->y+=domain.yl;
		if (part->y > domain.ymax)	part->y-=domain.yl;
	}
}

void BorisPushSpecies(struct Species *species, double *efx, double *efy, double B[])
{
  double qm = species->charge / species->mass;
  double wl = Lambda_D*omega_pe;
  double t_mag2 = 0.0;
  int dim;

  double *v_minus;
  double *v_prime;
  double *v_plus;

  double *t;
  double *s;

  v_minus = (double*) malloc(3 * sizeof(double));
  v_prime = (double*) malloc(3 * sizeof(double));
  v_plus = (double*) malloc(3 * sizeof(double));
  t = (double*) malloc(3 * sizeof(double));
  s = (double*) malloc(3 * sizeof(double));

  for(int p=0; p<species->np; p++)
  {
    struct Particle *part = &species->part[p];

    /*compute particle node position*/
		double lcx = XtoL(part->x);
    double lcy = XtoL(part->y);

		/*gather electric field onto particle position*/
		double part_efx = gather(lcx,lcy,efx);
    double part_efy = gather(lcx,lcy,efy);

    for(dim=0; dim<3; dim++)
    {
      t[dim] = qm*B[dim]*0.5*timeStep;
      t_mag2 += t[dim]*t[dim];
    }

    for(dim=0; dim<3; dim++)
    {
      s[dim] = 2*t[dim]/(1+t_mag2);
    }

    v_minus[0] = part->vx + (TE/(wl*wl))*qm*part_efx*0.5*timeStep;
    v_minus[1] = part->vy + (TE/(wl*wl))*qm*part_efy*0.5*timeStep;
    v_minus[2] = 0.0;

    double *v_minus_cross_t = CrossProduct(v_minus,t);

    v_prime[0] = v_minus[0] + v_minus_cross_t[0];
    v_prime[1] = v_minus[1] + v_minus_cross_t[1];
    v_prime[2] = 0.0;

    double *v_prime_cross_s = CrossProduct(v_prime,t);

    v_plus[0] = v_minus[0] + v_prime_cross_s[0];
    v_plus[1] = v_minus[1] + v_prime_cross_s[1];
    v_plus[2] = 0.0;

    /*advance velocity and position*/
    part->vx = v_plus[0] + (TE/(wl*wl))*qm*part_efx*0.5*timeStep;
    part->x += timeStep*part->vx;
    part->vy = v_plus[1] + (TE/(wl*wl))*qm*part_efy*0.5*timeStep;
    part->y += timeStep*part->vy;
    part->vz = v_plus[2];
    part->z += timeStep*part->vz;

    /*apply periodic boundaries*/
		if (part->x < domain.x0)	part->x+=domain.xl;
		if (part->x > domain.xmax)	part->x-=domain.xl;
    if (part->y < domain.y0)	part->y+=domain.yl;
		if (part->y > domain.ymax)	part->y-=domain.yl;
  }

  free(v_minus);
  free(v_prime);
  free(v_plus);
  free(t);
  free(s);
}

double* CrossProduct(double v1[3], double v2[3])
{
  static double r[3];
  r[0] = v1[1]*v2[2] - v1[2]*v2[1];
  r[1] = -v1[0]*v2[2] + v1[2]*v2[0];
  r[3] = v1[0]*v2[1] - v1[1]*v2[0];

  return r;
}



void WriteResults(int ts)
{
	//fprintf(file_res,"ZONE I=%d T=ZONE_%06d\n",domain.nix,ts);
  for (int k=0;k<domain.niz;k++)
  for (int j=0;j<domain.niy;j++)
  for (int i=0;i<domain.nix;i++)
	{
		fprintf(file_res,"%g %g %g %g %g %g %g %g %g\n",i*domain.dx, j*domain.dy, k*domain.dz,
									domain.nde[i*domain.niy+j],domain.ndi[i*domain.niy+j],
									domain.rho[(i*domain.niy+j)*domain.niz+k],domain.phi[(i*domain.niy+j)*domain.niz+k],domain.efx[(i*domain.niy+j)*domain.niz+k],domain.efy[(i*domain.niy+j)*domain.niz+k]);
	}

	fflush(file_res);
}

/*computes species kinetic energy in electron volt*/
double ComputeKE(struct Species *species)
{
	double ke=0;
	for (int p=0;p<species->np;p++)
		ke += species->part[p].vx*species->part[p].vx + species->part[p].vy*species->part[p].vy;

	/*we now have sum of v^2, multiply by 0.5*mass*/
	ke *= 0.5*(species->spwt*species->mass);

	/*convert to electron volts, 1eV=chargeE joules*/
	ke /=chargeE;

	return ke;
}

double rnd()
{
  double a = rand() % 10;
  a = 0.1*a;
  return a;
}

void Write_Particle(FILE *file, int ts, struct Species *species)
{
    for(int p=0; p<species->np;p++)
    {
        fprintf(file,"%g \t %g \t %g \t %g\n",species->part[p].x, species->part[p].y, species->part[p].vx, species->part[p].vy);
    }
    fflush(file_res);
}

void Write_density(FILE *file, int ts)
{
  for (int i=0;i<domain.nix;i++)
  for (int j=0;j<domain.niy;j++)
    {
      fprintf(file,"%g \t%g \t%g \t%g\n",i*domain.dx, j*domain.dy,
                    domain.nde[i*domain.niy+j],domain.ndi[i*domain.niy+j]);
    }
    fflush(file_res);
}
