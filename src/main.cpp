# include <iostream>
# include <complex>
# include <fftw3.h>
# include <cmath>
# include <cstdlib>
# include <vector>
# include <iterator>
# include <list>
# include <ctime>
# include <random>
# include <cstring>
# include <fstream>
#include <chrono>

extern "C" {
  #include "iniparser.h"
}

using namespace std;

/***************** HDF5 *******************/
#include <string>
#include "H5Cpp.h"
using namespace H5;
// DATA FILE NAME
const H5std_string FILE_NAME( "output/data.h5" );
const int RANK = 2; // 2 for 2D; 3 for 3D
// DATA FILE CREATION
H5File *file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
// DATA GROUP CREATION
Group* groupE = new Group( file->createGroup( "/particle.e" ));
Group* groupI = new Group( file->createGroup( "/particle.i" ));
Group* groupT = new Group( file->createGroup( "/timedata" ));
Group* groupP = new Group( file->createGroup( "/phi" ));
Group* groupR = new Group( file->createGroup( "/rho" ));
Group* groupEFx = new Group( file->createGroup( "/efx" ));
Group* groupEFy = new Group( file->createGroup( "/efy" ));
Group* groupDE = new Group( file->createGroup( "/den.e" ));
Group* groupDI = new Group( file->createGroup( "/den.i" ));
// Group* groupPK = new Group( file->createGroup( "/phik" ));
// Group* groupRK = new Group( file->createGroup( "/rhok" ));

// // Create new dataspace for attribute
//   DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
//   // Create attribute and write to it
//   Attribute varIN = dataset.createAttribute(LenX, datatype, attr_dataspace);
//   varIN.write(datatype, stepSize*numxCells);

/***************** INIPARSER ***************/
int  parse_ini_file(char * ini_name);

/************ RANDOM NUMBER GEN *********/
std::mt19937 mt_gen(0);
std::uniform_real_distribution<double> rnd_dist(0,1.0);
double rnd()
{
    return rnd_dist(mt_gen);
}

/********************************************/
/********* UNIVERSAL CONSTANTS *************/
const double EPS = 8.85418782E-12;    // Vacuum permittivity
const double K = 1.38065E-23;            // Boltzmann Constant
const double EV_TO_K = 11604.52;         // Conversion of EV to K
const double PI = 3.14159265359;
// const double chargeE = 1.602176565E-19; // Charge of an electron
// const double AMU = 1.660538921E-27;

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

double vthE; // electron temperature in eV
double vthI;  // ion temperature in eV
double v_te; // electron thermal velocity

double Lambda_D;
double omega_pe;

double B[3] = {0.0, 0.00001, 0.0001};

/* Class Domain: Hold the domain parameters*/
class Domain
{
public:
    int nix;      // Number of nodes
    double x0;   // initial position
    double dx;   // cell spacing
    double xl;   // domain length
    double xmax; // domain maximum position

    int niy;      // Number of nodes along y
    double y0;   // initial position in y
    double dy;   // cell spacing along y
    double yl;   // domain length along y
    double ymax; // domain maximum position along y

    /* Field Data structures */
    double *phi; // Electric Potential
    double *efx;  // Electric field along x
    double *efy;  // Electric field along y
    double *rho; // Charge Density
    fftw_complex *rhok;
    fftw_complex *phik;

};


/* Class Particle: Hold particle position, velocity and particle identity*/
class Particle
{
public:
    double xpos;  // particle position
    double ypos;  // particle position
    double xvel;  // particle velocity
    double yvel;  // particle velocity

    int id;  // hold particle identity

    // Add a constructor
    Particle(double x, double y, double u, double v):xpos(x), ypos(y), xvel(u), yvel(v){};
};

/* Class Species: Hold species data*/
class Species
{
public:
    // Use linked list for the particles
    list<Particle> part_list;
    double mass;
    double charge;
    double spwt;
    string name;

    int NUM;
    double Temp;
    double *den;
    double *xvel;
    double *yvel;

    void add(Particle part)
    {
        part.id=part_id++;
        part_list.push_back(part);
    }

    // Add a constructor
    Species(string name, double mass, double charge, double spwt, int NUM, double Temp)
    {
        setName(name);
        setMass(mass);
        setCharge(charge);
        setSpwt(spwt);
        setNum(NUM);
        setTemp(Temp);
    }

    // Define the constructor functions
    void setName(string name){this->name = name;}
    void setMass(double mass){this->mass = mass;}
    void setCharge(double charge){this->charge = charge;}
    void setSpwt(double spwt){this->spwt = spwt;}
    void setNum (int NUM){this->NUM = NUM;}
    void setTemp(double Temp){this->Temp = Temp;}

private:
    int part_id = 0;
};

// Define Domain and File as the global variable
Domain domain;


/********** HDF5 *********/
// HDF5 GROUP NAMES
H5std_string gNamePartE = "/particle.e/*";
H5std_string gNamePartI = "/particle.i/*";
H5std_string gNameTE = "/timedata/energy";
H5std_string gNamePhi = "/phi/*";
H5std_string gNameRho = "/rho/*";
H5std_string gNameEFx = "/efx/*";
H5std_string gNameEFy = "/efy/*";
H5std_string gNameDenE = "/den.e/*";
H5std_string gNameDenI = "/den.i/*";
// H5std_string gNamePhik = "/phik/*";
// H5std_string gNameRhok = "/rhok/*";


DataSpace *dataspace;
DataType datatype(H5::PredType::NATIVE_DOUBLE);
DataType dataInttype(H5::PredType::NATIVE_INT);
DataSet* dataset;
Attribute* attr;

// Define Helper functions
void init(Species *species, double xvel,double offset);
void scatterSpecies(Species *species);
void scatterSpeciesVel(Species *species);
void computeRho(double *rho, Species *ions, Species *electrons);
void computeEF(double *phi, double *efx, double *efy);
void pushSpecies(Species *species, double *efx, double *efy);
void rewindSpecies(Species *species, double *efx, double *efy);
void Write_ts(int ts, Species *ions, Species *electrons);
void writeSpecies(int ts, Species *species, H5std_string groupPart, H5std_string gNameDen);
void Write_VDF(FILE *file, int ts, double vdfLocStart, double vdfLocEnd, Species *species);

void writeKE(double energy[][2]); // 2 for KE of two species
void writePot(int ts, double *phi);
void writeRho(int ts, double *rho);
// void writePotk(int ts, fftw_complex *phik);
// void writeRhok(int ts, fftw_complex *rhok);
void writeEFx(int ts, double *efx);
void writeEFy(int ts, double *efy);
void writeAttributes(H5std_string groupPart, double data);
void writeIntAttributes(H5std_string attrName, int data);

void AddSources(Species *species);
void Inlet(Species *species);
// void computePE(double Time);

double computeKE(Species *species);
double XtoL(double xpos);
double YtoL(double ypos);
double gather(double lx, double ly, double *field);
void scatter(double lx, double ly, double value, double *field);
double sampleVel(double T, double mass);

bool solvePotential(double *phi, double *rho);
//bool solvePotentialDirect(double *phi, double *rho);
bool spectralPotentialSolver(double *phi, double *rho, fftw_complex *phik, fftw_complex *rhok);

double* CrossProduct(double v1[], double v2[]);
void BorispushSpecies(Species *species, double *efx, double *efy, double B[]);

/* HDF5 initiate */

// int h5_initiate(){
//
// }
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
    vthE     = iniparser_getdouble(ini,"population:vthE",-1.0);
    vthI     = iniparser_getdouble(ini,"population:vthI",-1.0);
    driftE   = iniparser_getdouble(ini,"population:driftE",-1.0);
    driftI   = iniparser_getdouble(ini,"population:driftI",-1.0);
    offsetE  = iniparser_getdouble(ini,"population:offsetE",-1.0);
    offsetI  = iniparser_getdouble(ini,"population:offsetI",-1.0);


    /* DIAGNOSTICS */
    // probLoc = iniparser_getint(ini,"diagnostics:probLoc",-1);
    dumpPeriod            = iniparser_getint(ini,"diagnostics:dumpPeriod",-1);

    /* USER CHOICE */
    solverType            = iniparser_getint(ini,"solver:solverType",-1);
    loadType              = iniparser_getint(ini,"population:loadType",-1);
    TSI                   = iniparser_getint(ini,"population:TSI",-1);

    /* Normalization */ //TO BE ADDED AS A SEPERATE FUNCTION
    // EPS             = EPS_un;//EPS_un;
    v_te     = sqrt(2*K*vthE*EV_TO_K/massE); // electron thermal vel
    omega_pe = sqrt((chargeE*chargeE*density)/(massE*EPS));
    Lambda_D = sqrt((EPS*K*vthE*EV_TO_K)/(density*chargeE*chargeE));



    cout << "********** IMPORTANT PLASMA QUANTITIES ***********" << '\n';
    cout<< "omega_pe: "<<omega_pe <<endl;
    cout<< "Lambda_D: "<<Lambda_D <<endl;

    cout << "*************** Input Sanity Check ***************" << '\n';
    bool SFLAG = true;
    if (stepSize >= 1) {
      cout<<"ERROR, stepSize is bigger than Debye length."<<endl;
      SFLAG = false;
    }
    if (timeStep > 0.5) {
      cout<<"ERROR, timeStep is too big. The recommended value: <"<<(0.5/omega_pe)<<" s"<<endl;
      SFLAG = false;
    }
    if (solverType != 1 && solverType != 2) {
      cout<<"ERROR, Wrong Solver Type. The recommended value: 1 or 2"<<endl;
      cout << "solverType: " << solverType <<endl;
      SFLAG = false;
    }
    if (loadType != 1 && loadType != 2 && loadType !=3) {
      cout<<"ERROR, Wrong Load Type. The recommended value: 1 or 2"<<endl;
      cout << "loadType: " << loadType <<endl;
      SFLAG = false;
    }
    if (SFLAG==true) {
      cout<<"STATUS, Input parameters are compatible."<<endl;
    }
    else {
      cout<<"ERROR, Input parameters are incompatible."<<endl;
      exit (EXIT_FAILURE);
    }


    iniparser_freedict(ini);
    return 0;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* MAIN FUNCTION ***************************/
int main(int argc, char *argv[])
{
    auto start = chrono::steady_clock::now(); //timer start
    /************* INIPARSER ********************/
    if(argc<2) {
      cout<<"ERROR, at least one argument expected (the input file)."<<endl;
      exit (EXIT_FAILURE);
    }
    parse_ini_file(argv[1]);

    /*********** HDF5 ATTRIBUTES ***************/
    writeAttributes("Lx", numxCells*stepSize);
    writeAttributes("Ly", numyCells*stepSize);
    writeIntAttributes("dp", dumpPeriod);
    writeIntAttributes("Nt", nTimeSteps);
    writeAttributes("dt", timeStep/omega_pe);
    writeIntAttributes("Nx", numxCells+1);
    writeIntAttributes("Ny", numyCells+1);
    writeAttributes("density", density);
    writeAttributes("massI", massI);
    writeAttributes("vthE", vthE);
    writeAttributes("vthI", vthI);

    /*********** **** ***************/
    // double Time = 0;

    // Energy time array
    // Note: (Define all of the time dependent arrays here to write to file in the end)
    double energy[int(nTimeSteps/dumpPeriod)+1][2];


    //omega_pe = sqrt((chargeE*chargeE*density)/(massE*EPS));;
    //Lambda_D = sqrt((EPS*K*vthE*EV_TO_K)/(density*chargeE*chargeE));;


    /*Construct the domain parameters*/
    domain.nix = numxCells+1;
    domain.dx = stepSize;
    domain.x0 = 0;
    domain.xl = (domain.nix-1)*domain.dx;
    domain.xmax = domain.x0 + domain.xl;


    /*Construct the domain parameters for y*/
    domain.niy = numyCells+1;
    domain.dy = stepSize;
    domain.y0 = 0;
    domain.yl = (domain.niy-1)*domain.dy;
    domain.ymax = domain.y0 + domain.yl;

    /*Allocate memory to the domain data structures (Field variables)*/
    domain.phi = new double[domain.nix*domain.niy];
    domain.efx = new double[domain.nix*domain.niy];
    domain.efy = new double[domain.nix*domain.niy];
    domain.rho = new double[domain.nix*domain.niy];
    domain.rhok = (fftw_complex*) fftw_malloc(domain.nix*(domain.niy/2+1) * sizeof(fftw_complex));
    domain.phik = (fftw_complex*) fftw_malloc(domain.nix*(domain.niy/2+1) * sizeof(fftw_complex));



    /*Redifine the field variables */
    double *phi = domain.phi;
    double *efx = domain.efx;
    double *efy = domain.efy;
    double *rho = domain.rho;
    fftw_complex *rhok = domain.rhok;
    fftw_complex *phik = domain.phik;

    /* Clear the domain fields*/

    memset(phi, 0,sizeof(double)*domain.nix*domain.niy);
    memset(efx, 0,sizeof(double)*domain.nix*domain.niy);
    memset(efy, 0,sizeof(double)*domain.nix*domain.niy);
    memset(rho, 0,sizeof(double)*domain.nix*domain.niy);


    /*Species Info: Create vector to hold the data*/
    vector <Species> species_list;

    ion_spwt = (density*numxCells*numyCells*domain.dx*domain.dy)/(nParticlesI);
    electron_spwt = (density*numxCells*numyCells*domain.dx*domain.dy)/(nParticlesE);


    /* Add singly charged ions and electrons */
    /*********************************************/
    /* Create the species lists*/
    species_list.emplace_back("Ion",massI,chargeE,ion_spwt, nParticlesI, vthI);
    species_list.emplace_back("Electrons",massE,-chargeE,electron_spwt, nParticlesE, vthE);

    /*Assign the species list as ions and electrons*/
    Species &ions = species_list[0];
    Species &electrons = species_list[1];

    /*initiate the species density and velocity fields*/


    ions.den = new double[domain.nix*domain.niy];
    ions.xvel = new double[domain.nix*domain.niy];
    ions.yvel = new double[domain.nix*domain.niy];
    electrons.den = new double[domain.nix*domain.niy];
    electrons.xvel = new double[domain.nix*domain.niy];
    electrons.yvel = new double[domain.nix*domain.niy];


    /*clear densities and velocities*/

    memset(ions.den,0,sizeof(double)*domain.nix*domain.niy);
    memset(ions.xvel,0,sizeof(double)*domain.nix*domain.niy);
    memset(ions.yvel,0,sizeof(double)*domain.nix*domain.niy);

    memset(electrons.den,0,sizeof(double)*domain.nix*domain.niy);
    memset(electrons.xvel,0,sizeof(double)*domain.nix*domain.niy);
    memset(electrons.yvel,0,sizeof(double)*domain.nix*domain.niy);


    /*initialize electrons and ions */
    init(&ions,driftI,offsetI);
    init(&electrons,driftE,offsetE);


    cout<< "*********** Parameters ***********"<<endl;
    for(auto &p:species_list)
        cout<< p.name << " mass: " << p.mass<< " charge: " << p.charge << " spwt: " << p.spwt << " Num of particles: " << p.NUM <<endl;
    cout<< "vdriftE: " << sqrt(2*(driftE*chargeE)/massE) <<" vdriftI: " << sqrt(2*(driftI*chargeE)/massI)<<endl;
    cout<< "density: " << density <<endl;
    cout<< "************ Simulation Parameters **********"<<endl;
    cout<< "Nx: " << numxCells << " Ny: " << numyCells <<endl;
    cout<< "Total timesteps: " << nTimeSteps <<endl;
    cout<< "timeStep: " << timeStep << " stepSize: " << domain.dx <<endl;

    cout<< "********** Beginning of Simulation  **********"<<endl;
    /***************************************************************************/

    /*Compute Number Density*/
    scatterSpecies(&ions);
    scatterSpecies(&electrons);

    /*Compute charge density, solve for potential
     and compute the electric field*/

    computeRho(rho, &ions, &electrons);

    /* Poisson Solver */
    if (solverType==1) {
      spectralPotentialSolver(phi, rho, phik, rhok);
    }
    else if (solverType==2) {
      solvePotential(phi, rho);
    }
    computeEF(phi,efx,efy);

    rewindSpecies(&ions,efx,efy);
    rewindSpecies(&electrons,efx,efy);


    /*TIME LOOP*/

    int ti =0; // Time dependent parameter index (Energy[ti])

    for (int ts=0; ts<nTimeSteps+1; ts++)
    {
      //Compute number density
      // scatterSpecies(&ions);
      scatterSpecies(&electrons);

      //Compute velocities
      // scatterSpeciesVel(&ions);  //TODO
      scatterSpeciesVel(&electrons);

      //Compute charge density
      computeRho(rho, &ions, &electrons);

      if (solverType==1) {
        spectralPotentialSolver(phi, rho, phik, rhok);
      }
      else if (solverType==2) {
        solvePotential(phi, rho);
      }

      computeEF(phi, efx, efy);

      //move particles
      // pushSpecies(&ions, efx, efy);  // TODO
      // pushSpecies(&electrons, efx, efy);
      // BorispushSpecies(&ions, efx, efy, B);
      BorispushSpecies(&electrons, efx, efy, B);
      //Write diagnostics
      if(ts%dumpPeriod== 0)
      {
          double max_phi = phi[0];
          for(int i=0; i<domain.nix; i++)
            for(int j=0; j<domain.niy; j++)
              {
                if(phi[i*domain.niy+j]>max_phi) max_phi=phi[i*domain.niy+j];
              }

          printf("TS: %i \t delta_phi: %.3g\n", ts, max_phi-phi[0]);

          writeSpecies(ts, &ions, gNamePartI, gNameDenI);
          writeSpecies(ts, &electrons, gNamePartE, gNameDenE);

          writePot(ts, phi);
          writeRho(ts, rho);
          writeEFx(ts, efx);
          writeEFy(ts, efy);
          // stop writing


          // writePotk(ts, phik);
          // writeRhok(ts, rhok);
          //computePE(Time);

          energy[ti][0] = computeKE(&ions);
          energy[ti][1] = computeKE(&electrons);
          // stop writing
          ti++; // increase time dependent parameter index

      }
      // WritePotOsc(Time,probLoc);
      // Time += timeStep;
    }

    writeKE(energy); // stop writing

    /*free up memory*/
    delete[] phi;
    delete[] rho;
    delete[] efx;
    delete[] efy;
    delete[] ions.den;
    delete[] electrons.den;
    delete[] ions.xvel;
    delete[] ions.yvel;
    delete[] electrons.xvel;
    delete[] electrons.yvel;
	delete[] phik;
	delete[] rhok;


    // /* HDF5*/
    // delete dataset;
    // delete dataspace;
    // cout << "/* OK */" << '\n';
    // /*
    delete groupE;
    delete groupI;
    delete groupT;
    delete groupP;
    delete groupR;
    delete groupDE;
    delete groupDI;
    delete groupEFx;
    delete groupEFy;
    // delete groupPK;
    // delete groupRK;
    // delete file;
    // */
    //****** END OF TIMER ******/

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Total time taken by PICSP: "<< chrono::duration <double> (diff).count() << " s" << endl;

    return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* HELPER FUNCTIONS ***************************/

/*initialize the particle data : initial positions and velocities of each particle*/
void init(Species *species, double xvel, double offset)
{
   //xvel = xvel/v_te;
   //yvel = yvel/v_te;
   // int part_per_cell = species->NUM/(numxCells*numyCells);
   // int xperiod = numxCells*part_per_cell;
   // int yperiod = numyCells*part_per_cell;
   double delta_x = domain.xl/species->NUM;
   double delta_y = domain.yl/species->NUM;
   // double delta_x = domain.xl/double(xperiod);
   // double delta_y = domain.yl/double(yperiod);
   double theta = 2*PI/domain.xl;
   double x=0.0, y=0.0, u=0.0, v=0.0;
   // int xcount=0, ycount=0;
    // sample particle positions and velocities
    for(int p=0; p<species->NUM; p++)
    {
        if (loadType==1) {
          // Maxwellian Loading
          x = domain.x0 + rnd()*(domain.nix-1)*domain.dx;
          u = sampleVel(species->Temp*EV_TO_K, species->mass);
          y = domain.y0 + rnd()*(domain.niy-1)*domain.dy;
          v = sampleVel(species->Temp*EV_TO_K, species->mass);

          // Periodic boundary
          if(x<0) x = x + domain.xl;
          if(x>domain.xl) x = x - domain.xl;
          if(y<0) y = y + domain.yl;
          if(y>domain.yl) y = y - domain.yl;

          // Add to the list
          species->add(Particle(x,y,u,v));

        }
        else if (loadType==2) {
          // Custom Loading
          x = domain.x0 + (p)*delta_x + offset*sin(x*theta);//
          // x = domain.x0 + xcount*delta_x;// + offset*delta_x; //sin(x*theta);
          y = domain.y0+(domain.niy-1)*domain.dy/2;//delta_y;// + rnd()*(domain.niy-1)*domain.dy;
          // if (p%xperiod==0) {
          //   xcount = 0;
          //   ycount = ycount+1;
          //  }
          //  xcount++;
          // TSI == 1 for two stream instability
          // if(TSI==0) {u = 0.0;}
          // else {u = xvel*pow(-1,p) + sampleVel(species->Temp*EV_TO_K, species->mass);}
          u = xvel*pow(-1,p) + sampleVel(species->Temp*EV_TO_K, species->mass);
          v = 0.0;

          // Periodic boundary
          if(x<0) x = x + domain.xl;
          if(x>domain.xl) x = x - domain.xl;

          if(y<0) y = y + domain.yl;
          if(y>domain.yl) y = y - domain.yl;

          // Add to the list
          species->add(Particle(x,y,u,v));

        }
        else if (loadType==3) {
          // Custom Loading
          // x = domain.x0 + p*domain.dx + offset*domain.dx;
          // x = domain.x0 + p*delta_x + offset*delta_x;
          x = domain.x0 + (32*domain.dx) + offset*domain.dx;
          y = domain.y0 + 32*domain.dy;
          cout << "Domain.dx: " <<domain.dx << endl;

          // TSI == 1 for two stream instability
          u = 0.0;
          v = 0.0;

          // Periodic boundary
          if(x<0) x = x + domain.xl;
          if(x>domain.xl) x = x - domain.xl;

          if(y<0) y = y + domain.yl;
          if(y>domain.yl) y = y - domain.yl;

          // Add to the list
          species->add(Particle(x,y,u,v));

        }
    }
}

/*Module to add further species in each time step: p is a macroparticle*/
/*
 void AddSources(Species *species)
 {
 for(int p=0; p<7; p++)
 {
 double x = domain.x0 + rnd()*(domain.ni-1)*domain.dx;
 //double x = (domain.xl - domain.x0)/2;
 //double x = (domain.x0+(domain.ni/100)*domain.dx) + rnd()*domain.xl-((domain.ni/100)*domain.dx);
 double v = sampleVel(species->Temp*EV_TO_K, species->mass);
 // Add to the list
 species->add(Particle(x,v));
 }
 }*/

/*Sample Velocity (According to Birdsall)*/
double sampleVel(double T, double mass)
{
    //double v_te = sqrt(2*K*vthE*EV_TO_K/massE);
    double v_th = sqrt(2*K*T/mass);
    double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
    return 0.5*vt/v_te;
}

/*Covert the physical coordinate to the logical coordinate*/
double XtoL(double xpos)
{
    double lx = (xpos-domain.x0)/domain.dx;
    return lx;
}
double YtoL(double ypos)
{
   double ly = (ypos - domain.y0)/domain.dy;                // converts physical coordinate to logical coordinate along y
   return ly;
}

/*scatter the particle data to the mesh and collect the densities at the mesh */
void scatter(double lx, double ly, double value, double *field)
{
   //double dxdy = domain.dx*domain.dy;

   int i = (int) lx;
   int j = (int) ly;
   double di = lx - i;
   double dj = ly - j;            // scatter particle density at nodes

   field[i*domain.niy+j] += value*(1-di)*(1-dj); //dxdy;
   field[(i+1)*domain.niy+j] += value*(di)*(1-dj); //dxdy;
   field[i*domain.niy+j+1] += value*(1-di)*(dj); //dxdy;
   field[(i+1)*domain.niy+j+1] += value*(di)*(dj); //dxdy;

   // if(di<domain.dx/2 && dj<domain.dy/2) { field[i*domain.niy+j] += value;}
   // else if (di>domain.dx/2 && dj<domain.dy/2) { field[(i+1)*domain.niy+j] += value;}
   // else if (di<domain.dx/2 && dj>domain.dy/2) { field[i*domain.niy+(j+1)] += value;}
   // else if (di>domain.dx/2 && dj>domain.dy/2) { field[(i+1)*domain.niy+(j+1)] += value;}
}

/* Gather field values at logical coordinates*/
double gather(double lx, double ly, double *field)
{
   int i = (int) lx;
   double di = lx - i;
                                                         // collect electric field from nodes at particle position
   int j = (int) ly;
   double dj = ly - j;

   // double val;

  // if(di<domain.dx/2 && dj<domain.dy/2) { val = field[i*domain.niy+j];}
  // else if (di>domain.dx/2 && dj<domain.dy/2) { val = field[(i+1)*domain.niy+j];}
  // else if (di<domain.dx/2 && dj>domain.dy/2) { val = field[i*domain.niy+(j+1)];}
  // else if (di>domain.dx/2 && dj>domain.dy/2) { val = field[(i+1)*domain.niy+(j+1)];}


   double val = field[i*domain.niy+j]*(1-di)*(1-dj)+field[(i+1)*domain.niy+j]*di*(1-dj)+field[i*domain.niy+j+1]*(1-di)*dj+field[(i+1)*domain.niy+j+1]*di*dj;
   return val;
}

/*Scatter the particles to the mesh for evaluating densities*/
void scatterSpecies(Species *species)
{
    /*grab a pointer to the species density data and change
     the density field using the pointer*/

    double *field = species->den;

    /*clear the field*/
    //memset(field,0,sizeof(double)*domain.ni);

    /*scatter particles to the mesh*/
    for(auto &p:species->part_list)
    {
        double lcx = XtoL(p.xpos);
        double lcy = YtoL(p.ypos);
        scatter(lcx,lcy,species->spwt,field);
    }

    for(int j=0; j<domain.niy; j++)
      {
      field[0*domain.niy+j] += field[(domain.nix-1)*domain.niy+j];
      field[(domain.nix-1)*domain.niy+j]=field[0*domain.niy+j];
      }

   for(int i=0; i<domain.nix; i++)
      {
      field[i*domain.niy+0] += field[i*domain.niy+domain.niy-1];
      field[i*domain.niy+domain.niy-1] = field[i*domain.niy+0];
      }

    /*divide by cell volume*/
    for(int i=0; i<domain.nix; i++)
    for(int j=0; j<domain.niy; j++)
        field[i*domain.niy+j] /=domain.dx*domain.dy;

    for(int i=0; i<domain.nix; i++)
    for(int j=0; j<domain.niy; j++)
    {
      field[i*domain.niy+j] /=density;  //Added
      //cout << species->name << " density :" << field[i*domain.niy+j] << endl;
    }

    for(int j=0; j<domain.niy; j++)
    {
      field[0+j] *=2.0;
      field[domain.niy*(domain.nix-1)+j] *= 2.0;
    }
}

/*Scatter the particles to the mesh for evaluating velocities*/
void scatterSpeciesVel(Species *species)
{
    /*grab a pointer to the species velocity field and change
     the velocity field using the pointer*/
    double *field1 = species->xvel;
    double *field2 = species->yvel;

    /*clear the field*/
   // memset(field,0,sizeof(double)*domain.ni);

    /*scatter particles to the mesh*/
    for(auto &p:species->part_list)
    {
        double lcx = XtoL(p.xpos);
        double lcy = YtoL(p.ypos);
        scatter(lcx,lcy,species->spwt*p.xvel,field1);
        scatter(lcx,lcy,species->spwt*p.yvel,field2);
    }


    for(int j=0; j<domain.niy; j++)
      {
      field1[0*domain.niy+j] += field1[(domain.nix-1)*domain.niy+j];
      field1[(domain.nix-1)*domain.niy+j]=field1[0*domain.niy+j];
      field2[0*domain.niy+j] += field2[(domain.nix-1)*domain.niy+j];
      field2[(domain.nix-1)*domain.niy+j]=field2[0*domain.niy+j];
      }

    for(int i=0; i<domain.nix; i++)
      {
      field1[i*domain.niy+0] += field1[i*domain.niy+domain.niy-1];
      field1[i*domain.niy+domain.niy-1] = field1[i*domain.niy+0];
      field2[i*domain.niy+0] += field2[i*domain.niy+domain.niy-1];
      field2[i*domain.niy+domain.niy-1] = field2[i*domain.niy+0];
      }


    /*divide by cell volume*/
    for(int i=0; i<domain.nix; i++)
    for(int j=0; j<domain.niy; j++){
        field1[i*domain.niy+j] /=(species->den[i*domain.niy+j]*domain.dx*domain.dy);
        field2[i*domain.niy+j] /=(species->den[i*domain.niy+j]*domain.dx*domain.dy);
        }

    //field[0][0] *=2.0;
    //field[domain.nix-1] *= 2.0;
    //field[domain.niy-1] *= 2.0;
}

//*******************************************************
void pushSpecies(Species *species, double *efx, double *efy)
{
    // compute charge to mass ratio
    double qm = species->charge/species->mass;
    list<Particle>::iterator it = species->part_list.begin();

    // loop over particles
    while (it!=species->part_list.end())
    {
        // grab a reference to the pointer
        Particle &part = *it;

        // compute particle node position
        double lcx = XtoL(part.xpos);
        double lcy = YtoL(part.ypos);

        // gather electric field onto particle position
        double part_efx = gather(lcx,lcy,efx);
        double part_efy = gather(lcx,lcy,efy);


        double wl = Lambda_D*omega_pe;

        // advance velocity & particle position
        part.xvel += 0.5*(1/(wl*wl))*(qm*vthE)*timeStep*part_efx;
        part.xpos += timeStep*part.xvel;
        if(loadType==3)
        {
        part.yvel = part.yvel;
        part.ypos = part.ypos;
        }
        else
        {
          part.yvel += 0.5*(1/(wl*wl))*(qm*vthE)*timeStep*part_efy;
          part.ypos += timeStep*part.yvel;
        }



      if(part.xpos < domain.x0) //|| part.xpos >= domain.xmax)
      {
        // it = species->part_list.erase(it);
         //continue;
         part.xpos += domain.xl;
      }
      else if(part.xpos >= domain.xmax)
      {
         part.xpos -= domain.xl;
      }
      else if(part.ypos < domain.y0)
      {
         part.ypos += domain.yl;
      }
      else if(part.ypos >= domain.ymax)
      {
         part.ypos -= domain.yl;
      }



        // Remove the particles leaving the domain
        //if(part.xpos < domain.x0 || part.xpos >= domain.xmax)
        //{
          //  it = species->part_list.erase(it);

            /* Encountering Steady state*/
            //part.pos = (domain.xl - domain.x0)/2; // relocate the particle in the middle of the domain
            //part.pos = domain.x0 + rnd()*(domain.ni - 1)*domain.dx;
            //cout << (domain.x0+(domain.ni/100)*domain.dx) << endl;
            /*
             part.pos = (domain.x0+(domain.ni/100)*domain.dx) + rnd()*domain.xl-((domain.ni/100)*domain.dx);
             part.vel = sampleVel(species->Temp*EV_TO_K, species->mass);
             species->add(Particle(part.pos,part.vel));*/

            //continue;
        //}
        else
            it++;
    }
}





void BorispushSpecies(Species *species, double *efx, double *efy, double B[])
{
    // compute charge to mass ratio
    double qm = species->charge/species->mass;
    list<Particle>::iterator it = species->part_list.begin();
    double wl = Lambda_D*omega_pe;

	double t_mag2;
    int dim;

    double *v_minus = new double[3];
    double *v_prime = new double[3];
    double *v_plus = new double[3];

    double *t = new double[3];
    double *s = new double[3];


    // loop over particles
    while (it!=species->part_list.end())
    {
        // grab a reference to the pointer
        Particle &part = *it;

        // compute particle node position
        double lcx = XtoL(part.xpos);
        double lcy = YtoL(part.ypos);

        // gather electric field onto particle position
        double part_efx = gather(lcx,lcy,efx);
        double part_efy = gather(lcx,lcy,efy);



        for(dim=0; dim<3; dim++)
          {
            t[dim] = qm*B[dim]*0.5*timeStep;
          }

        t_mag2 = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];

        for(dim=0; dim<3; dim++)
          {
            s[dim] = 2*t[dim]/(1+t_mag2);
          }

        v_minus[0] = part.xvel + (vthE/(wl*wl))*qm*part_efx*0.5*timeStep;
        v_minus[1] = part.yvel + (vthE/(wl*wl))*qm*part_efy*0.5*timeStep;
        v_minus[2] = 0.0;

        double *v_minus_cross_t = CrossProduct(v_minus, t);

        v_prime[0] = v_minus[0] + v_minus_cross_t[0];
        v_prime[1] = v_minus[1] + v_minus_cross_t[1];
        v_prime[2] = 0.0;

        double *v_prime_cross_s = CrossProduct(v_prime, s);
        v_plus[0] = v_minus[0] + v_prime_cross_s[0];
        v_plus[1] = v_minus[1] + v_prime_cross_s[1];


      // advance velocity & particle position
        part.xvel = v_plus[0] + (vthE/(wl*wl))*qm*part_efx*0.5*timeStep;
        part.xpos += timeStep*part.xvel;
        if(loadType==2)
        {
          part.yvel = part.yvel;
          part.ypos = part.ypos;
        }
        else if(loadType==3)
        {
          part.yvel = part.yvel;
          part.ypos = part.ypos;
        }
        else
        {
          part.yvel = v_plus[1] + (vthE/(wl*wl))*qm*part_efy*0.5*timeStep;
          part.ypos += timeStep*part.yvel;
        }


      if(part.xpos < domain.x0)
      {
         part.xpos += domain.xl;
      }
      else if(part.xpos >= domain.xmax)
      {
         part.xpos -= domain.xl;
      }
      else if(part.ypos < domain.y0)
      {
         part.ypos += domain.yl;
      }
      else if(part.ypos >= domain.ymax)
      {
         part.ypos -= domain.yl;
      }
      else
          it++;
    }

	delete[] v_minus;
	delete[] v_prime;
	delete[] v_plus;
	delete[] t;
	delete[] s;
}


double* CrossProduct(double v1[3], double v2[3])
{
  // double *r = new double[3];
  static double r[3];
  r[0] = v1[1]*v2[2] - v1[2]*v2[1];
  r[1] = -v1[0]*v2[2] + v1[2]*v2[0];
  r[2] = v1[0]*v2[1] - v1[1]*v2[0];

  return r;
}





//*********************************************************
/*Rewind particle velocities by -0.5*timeStep */
void rewindSpecies(Species *species, double *efx, double *efy)
{
    // compute charge to mass ratio
    double qm = species->charge/species->mass;
    for(auto &p:species->part_list)
    {
        // compute particle node position
        double lcx = XtoL(p.xpos);
        double lcy = XtoL(p.ypos);
        // gather electric field onto the particle position
        double part_efx = gather(lcx,lcy,efx);
        double part_efy = gather(lcx,lcy,efy);
        //rewind velocity
        double wl = Lambda_D*omega_pe;
        p.xvel -= 0.5*(1/(wl*wl))*(qm*vthE)*timeStep*part_efx;
        p.yvel -= 0.5*(1/(wl*wl))*(qm*vthE)*timeStep*part_efy;
    }

}

/* Compute the charge densities */
void computeRho(double *rho, Species *ions, Species *electrons)
{
    //double *rho = domain.rho;
    //memset(rho,0,sizeof(double)*domain.ni);

    for(int i=1; i<domain.nix-1; i++)
    for(int j=1; j<domain.niy-1; j++){

        rho[i*domain.niy+j]=(-ions->den[i*domain.niy+j] +
        electrons->den[i*domain.niy+j]);

        //rho[i*domain.niy+j]=ions->charge*ions->den[i*domain.niy+j] +
        //electrons->charge*electrons->den[i*domain.niy+j];
        }

    for(int j=0; j<domain.niy; j++)
      {
      rho[0*domain.niy+j] += rho[(domain.nix-1)*domain.niy+j];
      rho[(domain.nix-1)*domain.niy+j]=rho[0*domain.niy+j];
      }

    for(int i=0; i<domain.nix; i++)
      {
      rho[i*domain.niy+0] += rho[i*domain.niy+domain.niy-1];
      rho[i*domain.niy+domain.niy-1] = rho[i*domain.niy+0];
      }




    /*Reduce numerical noise by setting the densities to zero when less than 1e8/m^3*/
//     if(false){
//         for(int i=0; i<domain.nix; i++)
//         for(int j=0; j<domain.niy; j++)
//             if(fabs(rho[i*domain.niy+j])<1e8*chargeE) rho[i*domain.niy+j]=0;
//     }
}

/* Potential Solver: 1. Gauss-Seidel 2. Direct-Solver*/
bool solvePotential(double *phi, double *rho)
{
   double dx2 = domain.dx*domain.dx;
   // double dy2 = domain.dy*domain.dy;
   double L2;
   double beta = 1.75;

   for(int solver=0; solver<200000; solver++)
   {                                                     // solve electric potential at nodes using
      for(int i=0; i<domain.nix; i++)                  // Gauss - Seidel method
      for(int j=0; j<domain.niy; j++){

           int p=i-1; if(p<0) p=domain.nix-2;
           int q=i+1; if(q>domain.nix-1) q=1;
           int r=j-1; if(r<0) r=domain.niy-2;
           int s=j+1; if(s>domain.niy-1) s=1;

   // double g = 0.5*(1/((1/dx2)+(1/dy2)))*( ((phi[p*domain.niy+j]+phi[q*domain.niy+j])/dx2) + ((phi[i*domain.niy+r]+phi[i*domain.niy+s])/dy2) + (rho[i*domain.niy+j]) );

   // double g = beta*0.25*( ((phi[p*domain.niy+j]+phi[q*domain.niy+j]) + (phi[i*domain.niy+r]+phi[i*domain.niy+s])) - dx2*(rho[i*domain.niy+j]) ) + (1-beta)*phi[i*domain.niy+j];
   //
   // phi[i*domain.niy+j] = g;
   phi[i*domain.niy+j] = beta*0.25*( ((phi[p*domain.niy+j]+phi[q*domain.niy+j]) + (phi[i*domain.niy+r]+phi[i*domain.niy+s])) - dx2*(rho[i*domain.niy+j]) ) + (1-beta)*phi[i*domain.niy+j];


  // phi[i*domain.niy+j] = phi[i*domain.niy+j] + 1.9*(g - phi[i*domain.niy+j]);


   }


   if(solver%100==0)
     {                                                   // check solver convergence
         double sum = 0;
         for(int i=0; i<domain.nix;i++){
         for(int j=0; j<domain.niy;j++){
         int p = i-1; if(p<0) p = domain.nix-2;
         int q = i+1; if(q>domain.nix-1) q = 1;
         int r = j-1; if(r<0) r = domain.niy-2;
         int s = j+1; if(s>domain.niy-1) s = 1;

         //double R = -(rho[i][j]/EPS0) - (phi[p][j]-2*phi[i][j]+phi[q][j])/dx2 - (phi[i][r]-2*phi[i][j]+phi[i][s])/dy2;

         // double R = 0.25*(phi[p*domain.niy]+phi[q*domain.niy]+phi[i*domain.niy]+phi[i*domain.niy]+(dx2*rho[i*domain.niy]))-phi[i*domain.niy];
         // double R = 0.25*(phi[p*domain.niy+j]+phi[q*domain.niy+j]+phi[i*domain.niy+r]+phi[i*domain.niy+s]+(dx2*rho[i*domain.niy+j]))-phi[i*domain.niy+j];
         double R = ((phi[p*domain.niy+j]+phi[q*domain.niy+j]+phi[i*domain.niy+r]+phi[i*domain.niy+s]-4*phi[i*domain.niy+j])/dx2) - (rho[i*domain.niy+j]);

         sum = sum + (R*R);
         }
         }

         L2 = sqrt(sum)/(domain.nix*domain.niy);
         // cout<<L2<<endl;
         if(L2<1e-2){return true;}

      }
   }

   printf("Gauss-Seidel solver failed to converge, L2=%.3g!\n",L2);
	return false;
}


bool spectralPotentialSolver(double *phi, double *rho, fftw_complex *phik, fftw_complex *rhok)
{
   int Nx = domain.nix, Ny = domain.niy, Nh = (Ny/2) + 1;
   int i, j;
   fftw_complex *rhok_dum, *phik_dum;
   // fftw_complex *rhok, *phik, *rhok_dum, *phik_dum;

   fftw_plan p;
   fftw_plan b;

   // rhok = (fftw_complex*) fftw_malloc(Nx*Nh * sizeof(fftw_complex));
   // phik = (fftw_complex*) fftw_malloc(Nx*Nh * sizeof(fftw_complex));
   // Suggested by Rupak
   rhok_dum = (fftw_complex*) fftw_malloc(Nx*Nh * sizeof(fftw_complex));
   phik_dum = (fftw_complex*) fftw_malloc(Nx*Nh * sizeof(fftw_complex));


   // rhok = new fftw_complex[Nx*Nh];
   // phik = new fftw_complex[Nx*Nh];
   // rhok_dum = new fftw_complex[Nx*Nh];


   double Lx = domain.xl, Ly = domain.yl;
   double kx, ky;// Kx, Ky;


   //int num_thrd = omp_get_num_threads();


   // Take FFT of rho
   //fftw_init_threads();

    p = fftw_plan_dft_r2c_2d(Nx, Ny,  &rho[0*Ny+0], &rhok[0*Nh+0], FFTW_ESTIMATE);
   //fftw_plan_with_nthreads(4);
   fftw_execute(p);
   fftw_destroy_plan(p);
   fftw_cleanup();


   // Write dummy file

   for(j=0;j<Nh;j++)
   for(i=0;i<Nx;i++)
   {
      rhok_dum[i*Nh+j][0] = rhok[i*Nh+j][0];
      rhok_dum[i*Nh+j][1] = rhok[i*Nh+j][1];
   }

   // Find phik = rhok / (kx^2 + ky^2)

   for(j=0;j<Nh;j++)
   {
      ky = 2.0*PI*j/Ly;

      //#pragma openmp parallel for
      for(i=0;i<Nx/2;i++)
      {
         kx = 2.0*PI*i/Lx;

         phik[i*Nh+j][0] = rhok_dum[i*Nh+j][0]/(kx*kx+ky*ky);
         phik[i*Nh+j][1] = rhok_dum[i*Nh+j][1]/(kx*kx+ky*ky);

      }
      //#pragma openmp parallel for
      for(i=Nx/2+1;i<Nx;i++)
      {
         kx = 2.0*PI*(Nx-i)/Lx;

         phik[i*Nh+j][0] = rhok_dum[i*Nh+j][0]/(kx*kx+ky*ky);
         phik[i*Nh+j][1] = rhok_dum[i*Nh+j][1]/(kx*kx+ky*ky);

      }
   }

   phik[0][0] = 0;
   phik[0][1] = 0;

   // Write dummy file
   for(j=0;j<Nh;j++)
   for(i=0;i<Nx;i++)
   {
      phik_dum[i*Nh+j][0] = phik[i*Nh+j][0];
      phik_dum[i*Nh+j][1] = phik[i*Nh+j][1];
   }


   // Take iFFT of phik
   //fftw_init_threads();

   b = fftw_plan_dft_c2r_2d(Nx, Ny, &phik[0*Nh+0], &phi[0*Ny+0], FFTW_ESTIMATE);
   //fftw_plan_with_nthreads(4);
   fftw_execute(b);
   fftw_destroy_plan(b);
   fftw_cleanup();
   //fftw_cleanup_threads();


   for(j=0;j<Ny;j++)
   for(i=0;i<Nx;i++)
   {
      phi[i*Ny+j] /= double(Nx*Ny);
   }

   delete[] rhok_dum;
   delete[] phik_dum;
   return true;
}





/* Potential Direct Solver */

// bool solvePotentialDirect(double *x, double *rho)
// {
//     /* Set coefficients, precompute them*/
//     int ni = domain.ni;
//     double dx2 = domain.dx*domain.dx;
//     double *a = new double[ni];
//     double *b = new double[ni];
//     double *c = new double[ni];
//
//     /*Centtral difference on internal nodes*/
//     for(int i=1; i<ni-1; i++)
//     {
//         a[i] = 1; b[i] = -2; c[i] = 1;
//     }
//
//     /*Apply dirichlet boundary conditions on boundaries*/
//     a[0]=0; b[0]=1; c[0]=0;
//     a[ni-1]=0; b[ni-1]=1; c[ni-1]=0;
//
//     /*multiply R.H.S.*/
//     for (int i=1; i<ni-1; i++)
//         x[i]=-rho[i]*dx2/EPS;
//
//     x[0] = 0;
//     x[ni-1] = 0;
//
//     /*Modify the coefficients*/
//     c[0] /=b[0];
//     x[0] /=b[0];
//
//     for(int i=1; i<ni; i++)
//     {
//         double id = (b[i]-c[i-1]*a[i]);
//         c[i] /= id;
//         x[i] = (x[i]-x[i-1]*a[i])/id;
//     }
//
//     /* Now back substitute */
//     for(int i=ni-2; i>=0; i--)
//         x[i] = x[i] - c[i]*x[i+1];
//
//     return true;
// }

/*Compute electric field (differentiating potential)*/
void computeEF(double *phi, double *efx, double *efy)
{
    /*Apply central difference to the inner nodes*/

    for(int j=1; j<domain.niy-1; j++)
    for(int i=1; i<domain.nix-1; i++)
    {
      efx[i*domain.niy+j] = (phi[(i-1)*domain.niy+j]-phi[(i+1)*domain.niy+j])/(2*domain.dx);
      efy[i*domain.niy+j] = (phi[i*domain.niy+(j-1)]-phi[i*domain.niy+(j+1)])/(2*domain.dy);
    }



    /*Apply one sided difference at the boundary nodes*/
    for(int j=0; j<domain.niy; j++)                                    // calculate electric field at boundary nodes
   {
      //efx[0][j] = -(phi[1][j]-phi[domain.nix-2][j])/(2*domain.dx);
      //efx[domain.nix-1][j] = efx[0][j];
      efx[0*domain.niy+j] = -(phi[1*domain.niy+j] - phi[0*domain.niy+j])/(domain.dx);
      efx[(domain.nix-1)*domain.niy+j] = -(phi[(domain.nix-1)*domain.niy+j] - phi[(domain.nix-2)*domain.niy+j])/(domain.dx);

   }

   for(int i=0; i<domain.nix; i++)
   {
      //efy[i][0] = -(phi[i][1]-phi[i][domain.niy-2])/(2*domain.dy);
      //efy[i][domain.niy-1] = efy[i][0];

      efy[i*domain.niy+0] = -(phi[i*domain.niy+1] - phi[i*domain.niy+0])/(domain.dy);
      efy[i*domain.niy+(domain.niy-1)] = -(phi[i*domain.niy+(domain.niy-1)] - phi[i*domain.niy+(domain.niy-2)])/(domain.dy);
   }
}

/* Write the Output results*/
void writeSpecies(int ts, Species *species, H5std_string groupPart, H5std_string gNameDen)
{
    /////////// Particle Data ///////////
    hsize_t nP = species->part_list.size();
    hsize_t  dims[2] = {nP,4};
    groupPart.replace(groupPart.begin()+12,groupPart.end(),to_string(ts));
    dataspace = new DataSpace(RANK, dims); // create new dspace
    dataset = new DataSet(file->createDataSet(groupPart,datatype, *dataspace));
    int i = 0;
    // Array to store particle data as contiguous memory for HDF5
    double varray[species->part_list.size()][4];
    // iterate over particles
    for(auto& p: species->part_list)
    {
        varray[i][0] = p.xpos;
        varray[i][1] = p.ypos;
        varray[i][2] = p.xvel;
        varray[i][3] = p.yvel;
        i++;
    }
    dataset->write(varray, datatype);
    delete dataset;
    delete dataspace;
    ///////////// Density Data //////////////
    hsize_t nx = domain.nix;
    hsize_t ny = domain.niy;
    hsize_t  dimsp[2] = {nx,ny};
    gNameDen.replace(gNameDen.begin()+7,gNameDen.end(),to_string(ts));
    dataspace = new DataSpace(RANK, dimsp); // create new dspace
    dataset = new DataSet(file->createDataSet(gNameDen,datatype, *dataspace));
    dataset->write(species->den, datatype);
    delete dataset;
    delete dataspace;
}

/* ********************************** */

void writeKE(double energy[][2])
{
    hsize_t nT = int(nTimeSteps/dumpPeriod)+1;
    hsize_t  dimstime[2] = {nT,2};
    dataspace = new DataSpace(RANK, dimstime); // create second dspace
    dataset = new DataSet(file->createDataSet(gNameTE,datatype, *dataspace));
    dataset->write(energy, datatype);
    delete dataset;
    delete dataspace;
}

double computeKE(Species *species)
{
    double ke = 0;
    for (auto &p:species->part_list)
    {
        ke += p.xvel*p.xvel + p.yvel*p.yvel;
    }
    /*Multiply 0.5*mass for all particles*/
    ke += 0.5*(species->spwt*species->mass);

    /*Convert the kinetic energy in eV units*/
    ke /= chargeE;
    return ke;
}

void writePot(int ts, double *phi)
{
  hsize_t nx = domain.nix;
  hsize_t ny = domain.niy;
  hsize_t  dimsp[2] = {nx,ny};
  gNamePhi.replace(gNamePhi.begin()+5,gNamePhi.end(),to_string(ts));
  dataspace = new DataSpace(RANK, dimsp); // create new dspace
  dataset = new DataSet(file->createDataSet(gNamePhi,datatype, *dataspace));
  dataset->write(phi, datatype);
  delete dataset;
  delete dataspace;
}

void writeRho(int ts, double *rho)
{
  hsize_t nx = domain.nix;
  hsize_t ny = domain.niy;
  hsize_t  dimsp[2] = {nx,ny};
  gNameRho.replace(gNameRho.begin()+5,gNameRho.end(),to_string(ts));
  dataspace = new DataSpace(RANK, dimsp); // create new dspace
  dataset = new DataSet(file->createDataSet(gNameRho,datatype, *dataspace));
  dataset->write(rho, datatype);
  delete dataset;
  delete dataspace;
}

void writeEFx(int ts, double *efx)
{
  hsize_t nx = domain.nix;
  hsize_t ny = domain.niy;
  hsize_t  dimsp[2] = {nx,ny};
  gNameEFx.replace(gNameEFx.begin()+5,gNameEFx.end(),to_string(ts));
  dataspace = new DataSpace(RANK, dimsp); // create new dspace
  dataset = new DataSet(file->createDataSet(gNameEFx,datatype, *dataspace));
  dataset->write(efx, datatype);
  delete dataset;
  delete dataspace;
}

void writeEFy(int ts, double *efy)
{
  hsize_t nx = domain.nix;
  hsize_t ny = domain.niy;
  hsize_t  dimsp[2] = {nx,ny};
  gNameEFy.replace(gNameEFy.begin()+5,gNameEFy.end(),to_string(ts));
  dataspace = new DataSpace(RANK, dimsp); // create new dspace
  dataset = new DataSet(file->createDataSet(gNameEFy,datatype, *dataspace));
  dataset->write(efy, datatype);
  delete dataset;
  delete dataspace;
}


// void writePotk(int ts, fftw_complex *phik)
// {
  // hsize_t nx = domain.nix;
  // hsize_t ny = domain.niy;
  // hsize_t  dimsp[2] = {nx,ny};
  // gNamePhik.replace(gNamePhik.begin()+6,gNamePhik.end(),to_string(ts));
  // dataspace = new DataSpace(RANK, dimsp); // create new dspace
  // dataset = new DataSet(file->createDataSet(gNamePhik,datatype, *dataspace));
  // dataset->write(phik, datatype);
  // delete dataset;
  // delete dataspace;
// }

// void writeRhok(int ts, fftw_complex *rhok)
// {
  // hsize_t nx = domain.nix;
  // hsize_t ny = domain.niy;
  // hsize_t  dimsp[2] = {nx,ny};
  // gNameRhok.replace(gNameRhok.begin()+6,gNameRhok.end(),to_string(ts));
  // dataspace = new DataSpace(RANK, dimsp); // create new dspace
  // dataset = new DataSet(file->createDataSet(gNameRhok,datatype, *dataspace));
  // dataset->write(rhok, datatype);
  // delete dataset;
  // delete dataspace;
// }

// void computePE(double Time)
// {
//    double pe = 0;
//    for(int j=0; j<domain.niy; j++)
//    for(int i=0; i<domain.nix; i++)
//    {
//      pe += 0.5*domain.phi[i*domain.niy+j]*domain.rho[i*domain.niy+j];
//    }
//
// }

void writeAttributes(H5std_string attrName, double data)
{
  dataspace = new DataSpace(H5S_SCALAR); // create new dspace
  attr = new Attribute(file->createAttribute(attrName, datatype, *dataspace));
  // double data = numxCells*stepSize;
  attr->write(datatype, &data);
  delete attr;
  delete dataspace;
}

void writeIntAttributes(H5std_string attrName, int data)
{
  dataspace = new DataSpace(H5S_SCALAR); // create new dspace
  attr = new Attribute(file->createAttribute(attrName, dataInttype, *dataspace));
  // double data = numxCells*stepSize;
  attr->write(dataInttype, &data);
  delete attr;
  delete dataspace;
}
