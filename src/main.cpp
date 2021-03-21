# include <iostream>
# include <cmath>
# include <cstdlib>
# include <vector>
# include <iterator>
# include <list>
# include <ctime>
# include <random>
# include <cstring>
# include <fstream>

extern "C" {
  #include "iniparser.h"
}

using namespace std;

/*Iniparser function*/
int  parse_ini_file(char * ini_name);

/* Random Number Generator */
std::mt19937 mt_gen(0);
std::uniform_real_distribution<double> rnd_dist(0,1.0);
double rnd()
{
    return rnd_dist(mt_gen);
}


/* Define universal constants */
const double EPS_un = 8.85418782E-12;    // Vacuum permittivity
const double K = 1.38065E-23;        // Boltzmann Constant
// const double chargeE = 1.602176565E-19; // Charge of an electron
const double AMU = 1.660538921E-27;
const double EV_TO_K = 11604.52;
const double pi = 3.14159265359;


/* CHANGED TYPE FROM CONST TO VAR FOR INPUT DATA CONTROL  */
int nParticlesI;      // Number of simulation ions
int nParticlesE; // Number of simulation electrons

int numxCells;             // Total number of cells alonx x
int numyCells;             // Total number of cells alonx y
int nTimeSteps;          // Total time steps (default)
double massI;  // Ion mass
double massE;  // Electron mass
double chargeE; // Electron charge
// double vdfLocStart;  //VDF start location
// double vdfLocEnd;  //VDF end location
// int probLoc;  //VDF end location
double velocity; //5e5           // TODO
double EPS;

/* Simulation Parameters*/
double density; // Plasma Density
/*Calculate the specific weights of the ions and electrons*/
double ion_spwt;
double electron_spwt;
double stepSize;         // Cell Spacing
double timeStep;        // Time steps

double thermalVelocityE; // electron temperature in eV
double thermalVelocityI;  // ion temperature in eV

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

};


/* Class Particle: Hold particle position, velocity and particle identity*/
class Particle
{
public:
    double xpos;  // particle position
    double xvel;  // particle velocity

    double ypos;  // particle position
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
FILE *file_res;
FILE *file_ke;
FILE *file_phi;
FILE *f1;
FILE *f2;
// FILE *f3;
FILE *file_res1;
FILE *file_res2;
FILE *file_sp;

// Define Helper functions
void Init(Species *species, double xvel,double yvel);
void ScatterSpecies(Species *species);
void ScatterSpeciesVel(Species *species);
void ComputeRho(double *rho, Species *ions, Species *electrons);
void ComputeEF(double *phi, double *efx, double *efy);
void PushSpecies(Species *species, double *efx, double *efy);
void RewindSpecies(Species *species, double *efx, double *efy);
void Write_ts(int ts, Species *ions, Species *electrons);
void Write_Particle(FILE *file, int ts, Species *species);
void Write_VDF(FILE *file, int ts, double vdfLocStart, double vdfLocEnd, Species *species);
void WriteKE(double Time, Species *ions, Species *electrons);
void Write_Single_Particle(Species *species);
void AddSources(Species *species);
void Inlet(Species *species);
void Write_pot(double Time);
void ComputePE(double Time);

double ComputeKE(Species *species);
double XtoL(double xpos);
double YtoL(double ypos);
double gather(double lx, double ly, double *field);
void scatter(double lx, double ly, double value, double *field);
double SampleVel(double T, double mass);

bool SolvePotential(double *phi, double *rho);
//bool SolvePotentialDirect(double *phi, double *rho);



/*Parsing Input file*/

int parse_ini_file(char * ini_name)
{
    dictionary  *   ini ;

    ini = iniparser_load(ini_name);
    if (ini==NULL) {
        fprintf(stderr, "cannot parse file: %s\n", ini_name);
        return -1 ;
    }
    iniparser_dump(ini, stderr);

    /*Get Simulation Parameters */
    nTimeSteps  = iniparser_getint(ini,"time:nTimeSteps",-1);
    double timeStep_unorm    = iniparser_getdouble(ini,"time:timeStep",-1.0);
    double stepSize_unorm    = iniparser_getdouble(ini,"grid:stepSize",-1.0);
    numxCells    = iniparser_getint(ini,"grid:numxCells",-1);
    numyCells    = iniparser_getint(ini,"grid:numyCells",-1);

    /* NUM OF COM PARTICLE */
    nParticlesI = iniparser_getint(ini,"population:nParticlesI",-1);
    nParticlesE = iniparser_getint(ini,"population:nParticlesE",-1);
    double massI_unorm  =  iniparser_getdouble(ini,"population:massI",-1.0);
    double massE_unorm   =  iniparser_getdouble(ini,"population:massE",-1.0);
    double chargeE_unorm   =  iniparser_getdouble(ini,"population:chargeE",-1.0);
    double density_unorm  = iniparser_getdouble(ini,"population:density",-1.0);
    double thermalVelocityE_unorm = iniparser_getdouble(ini,"population:thermalVelocityE",-1.0);
    double thermalVelocityI_unorm = iniparser_getdouble(ini,"population:thermalVelocityI",-1.0);
    double velocity_unorm = iniparser_getdouble(ini,"population:driftE",-1.0);
    /* DIAGNOSTICS */
    // vdfLocStart = iniparser_getdouble(ini,"diagnostics:vdfLocStart",-1.0);
    // vdfLocEnd = iniparser_getdouble(ini,"diagnostics:vdfLocEnd",-1.0);
    // probLoc = iniparser_getint(ini,"diagnostics:probLoc",-1);

    /* Normalization */ //TO BE ADDED AS A SEPERATE FUNCTION
    EPS = EPS_un/EPS_un;
    double omega_pe = sqrt((chargeE_unorm*chargeE_unorm*density_unorm)/(massE_unorm*EPS_un));
    double Lambda_D = sqrt((EPS_un*K*thermalVelocityE_unorm*EV_TO_K)/(density_unorm*chargeE_unorm*chargeE_unorm));
    chargeE = chargeE_unorm/chargeE_unorm;
    massI = massI_unorm/massE_unorm;
    massE = massE_unorm/massE_unorm;
    velocity = velocity_unorm; ///velocity_unorm;
    density  = density_unorm/density_unorm;
    timeStep = timeStep_unorm*omega_pe;
    stepSize = stepSize_unorm/Lambda_D;
    thermalVelocityE = thermalVelocityE_unorm/thermalVelocityE_unorm;
    thermalVelocityI = thermalVelocityI_unorm/thermalVelocityE_unorm;
    /*Calculate the specific weights of the ions and electrons*/
    ion_spwt = (density*numxCells*numyCells*stepSize*stepSize)/(nParticlesI);
    electron_spwt = (density*numxCells*numyCells*stepSize*stepSize)/(nParticlesE);
    cout<< "omega_pe: "<<omega_pe <<endl;


    iniparser_freedict(ini);
    return 0 ;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* MAIN FUNCTION ***************************/
int main(int argc, char *argv[])
{
    //iniparser_load(*ini);
    parse_ini_file(argv[1]);


    double Time = 0;

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


    /*Redifine the field variables */
    double *phi = domain.phi;
    double *efx = domain.efx;
    double *efy = domain.efy;
    double *rho = domain.rho;

    /* Clear the domain fields*/

     memset(phi,0,sizeof(double)*domain.nix*domain.niy);
     memset(efx, 0,sizeof(double)*domain.nix*domain.niy);
     memset(efy, 0,sizeof(double)*domain.nix*domain.niy);
     memset(rho,0,sizeof(double)*domain.nix*domain.niy);


    /*Species Info: Create vector to hold the data*/
    vector <Species> species_list;

    ion_spwt = (density*numxCells*numyCells*stepSize*stepSize)/(nParticlesI);
    electron_spwt = (density*numxCells*numyCells*stepSize*stepSize)/(nParticlesE);

    /* Add singly charged ions and electrons */
    /*********************************************/
    /* Create the species lists*/
    species_list.emplace_back("Ion",massI,chargeE,ion_spwt, nParticlesI, thermalVelocityI);
    species_list.emplace_back("Electrons",massE,-chargeE,electron_spwt, nParticlesE, thermalVelocityE);

    /*Assign the species list as ions and electrons*/
    Species &ions = species_list[0];
    Species &electrons = species_list[1];

    /*Initiate the species density and velocity fields*/


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




    /*Initialize electrons and ions */
    Init(&ions,0,0);
    Init(&electrons,velocity,0);

    for(auto &p:species_list)
        cout<< p.name << '\n' << p.mass<< '\n' << p.charge << '\n' << p.spwt << '\n' << p.NUM << endl <<endl;
        cout<< "nParticlesI: " << nParticlesI << " nParticlesE: " << nParticlesE <<endl;
        cout<< "numxCells: " << numxCells << " numyCells: " << numyCells <<endl;
        cout<< "nTimeSteps: " << nTimeSteps <<endl;
        cout<< "massI: " << massI << " massE: " << massE <<endl;
        cout<< "velocity: " << velocity <<endl;
        cout<< "density: " << density <<endl;
        cout<< "timeStep: " << timeStep << " stepSize: " << stepSize <<endl;
    /***************************************************************************/

    /*Compute Number Density*/
    ScatterSpecies(&ions);
    ScatterSpecies(&electrons);

    /*Compute charge density, solve for potential
     and compute the electric field*/

    ComputeRho(rho, &ions, &electrons);
    SolvePotential(phi, rho);
    ComputeEF(phi,efx,efy);

    RewindSpecies(&ions,efx,efy);
    RewindSpecies(&electrons,efx,efy);

    /*------------- Print Output ---------------*/

    /*create a folder named output and*/
    /*delete the previous output folder: print statement is just to show*/
    printf("rm -rf output/*\n");
    system("rm -rf output");

    /*create an output folder*/
    system("mkdir output");

    /*create a seperate directory for phase-space data inside output*/
    system("mkdir output/phase_space");

    /*create a seperate directory for VDF data inside output*/
   system("mkdir output/vdf_output");

    char NAmassE[50];
    // char NAmassI[50];

    file_res = fopen("output/results.dat","w");
    file_ke = fopen("output/ke.dat","w");
    file_sp = fopen("output/part.dat","w");
    // file_phi = fopen("output/phi.dat","w");
    file_res2 = fopen("output/pe.dat","w");


    /*MAIN LOOP*/

    for (int ts=0; ts<nTimeSteps+1; ts++)
    {
        //Compute number density
        ScatterSpecies(&ions);
        ScatterSpecies(&electrons);

        //Compute velocities
        // ScatterSpeciesVel(&ions);  //TODO
        ScatterSpeciesVel(&electrons);

        //Compute charge density
        ComputeRho(rho, &ions, &electrons);

        SolvePotential(phi, rho);
        // SolvePotentialDirect(phi, rho);
        ComputeEF(phi, efx, efy);

        //move particles
        // PushSpecies(&ions, efx, efy);  // TODO
        PushSpecies(&electrons, efx, efy);

        //Write diagnostics
        if(ts%50== 0)
        {
            sprintf(NAmassE,"output/phase_space/i%d.dat",ts);
            f1 = fopen(NAmassE,"w");

            sprintf(NAmassE,"output/phase_space/e%d.dat",ts);
            f2 = fopen(NAmassE,"w");

            // //Added by SAYAN 14/08/2019 for VDF data
            // sprintf(NAmassE,"output/vdf_output/i%d.dat",ts);
            // f3 = fopen(NAmassE,"w");

            sprintf(NAmassE, "output/ef%d.dat",ts);
            file_res1 = fopen(NAmassE,"w");

            ///////////////////////////////////////
            double max_phi = phi[0];
            for(int i=0; i<domain.nix; i++)
              for(int j=0; j<domain.niy; j++)
                {
                  if(phi[i*domain.niy+j]>max_phi) max_phi=phi[i*domain.niy+j];
                }

            //Compute kinetic energy
            //double ke_ions = ComputeKE(&ions)/(ions.NUN*ions.spwt);
            //double ke_electrons = ComputeKE(&electrons)/(electrons.NUN*electrons.spwt);

            printf("TS: %i \t delta_phi: %.3g\n", ts, max_phi-phi[0]);


            WriteKE(Time, &ions, &electrons);
            //Write_ts(ts,&ions,&electrons);

            Write_Particle(f1,ts, &ions);
            Write_Particle(f2,ts, &electrons);
            Write_Single_Particle(&electrons);
            ComputePE(Time);

            // Write_VDF(f3,ts,vdfLocStart,vdfLocEnd, &ions);  //Added by SAYAN 14/08/2019
            Write_pot(Time);

            fclose(f1);
            fclose(f2);
            // fclose(f3);
        }
        // WritePotOsc(Time,probLoc);
        Time += timeStep;
    }

    /*free up memory*/
     delete[] phi;
     delete[] rho;
     delete[] efx;
     delete[] efy;
     delete[] ions.den;
     delete[] ions.xvel;
     delete[] ions.yvel;
     delete[] electrons.den;
     delete[] electrons.xvel;
     delete[] electrons.yvel;


    /*copy other diagnostics to output directory*/
//    system("cp *dat output/");
    /*clean home directory*/
//    system("rm -f *dat");

    return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* HELPER FUNCTIONS ***************************/

/*Initialize the particle data : initial positions and velocities of each particle*/
void Init(Species *species, double xvel, double yvel)
{

   double delta_x = domain.xl/species->NUM;
   double delta_y = domain.yl/species->NUM;
   double theta = 2*pi/domain.xl;
    // sample particle positions and velocities
    for(int p=0; p<species->NUM; p++)
    {
        double x = domain.x0 + (p+0.5)*delta_x + 0.1*sin(theta*x); //domain.x0 + rnd()*(domain.nix-1)*domain.dx;

        double u = xvel*pow(-1,p); // SampleVel(species->Temp*EV_TO_K, species->mass);

        double y = (p+0.5)*delta_y; //domain.y0 + rnd()*(domain.niy-1)*domain.dy; //(p+0.5)*delta_y;
        double v = yvel; // + SampleVel(species->Temp*EV_TO_K, species->mass);                    //yvel; //SampleVel(species->Temp*EV_TO_K, species->mass);

        if(x<0) x = x + domain.xl;
        if(x>domain.xl) x = x - domain.xl;

        if(y<0) y = y + domain.yl;
        if(y>domain.yl) y = y - domain.yl;

        // Add to the list
        species->add(Particle(x,y,u,v));
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
 double v = SampleVel(species->Temp*EV_TO_K, species->mass);
 // Add to the list
 species->add(Particle(x,v));
 }
 }*/

/*Sample Velocity (According to Birdsall)*/
double SampleVel(double T, double mass)
{
    double v_th = sqrt(2*K*T/mass);
    return v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
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
   double dxdy = domain.dx*domain.dy;

   int i = (int) lx;
   int j = (int) ly;
   double di = lx - i;
   double dj = ly - j;            // scatter particle density at nodes

   field[i*domain.niy+j] += value*(1-di)*(1-dj)/dxdy;
   field[(i+1)*domain.niy+j] += value*(di)*(1-dj)/dxdy;
   field[i*domain.niy+j+1] += value*(1-di)*(dj)/dxdy;
   field[(i+1)*domain.niy+j+1] += value*(di)*(dj)/dxdy;
}

/* Gather field values at logical coordinates*/
double gather(double lx, double ly, double *field)
{
   int i = (int) lx;
   double di = lx - i;
                                                         // collect electric field from nodes at particle position
   int j = (int) ly;
   double dj = ly - j;

   double val = field[i*domain.niy+j]*(1-di)*(1-dj)+field[(i+1)*domain.niy+j]*di*(1-dj)+field[i*domain.niy+j+1]*(1-di)*dj+field[(i+1)*domain.niy+j+1]*di*dj;
   return val;
}

/*Scatter the particles to the mesh for evaluating densities*/
void ScatterSpecies(Species *species)
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
    /*for(int i=0; i<domain.nix; i++)
    for(int j=0; j<domain.niy; j++)
        field[i][j] /=domain.dx*domain.dy;*/

    //field[0] *=2.0;
    //field[domain.ni-1] *= 2.0;
}

/*Scatter the particles to the mesh for evaluating velocities*/
void ScatterSpeciesVel(Species *species)
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
    /*for(int i=0; i<domain.nix; i++)
    for(int j=0; j<domain.niy; j++)
        field[i][j] /=(species->den[i][j]*domain.dx*domain.dy);*/

    //field[0][0] *=2.0;
    //field[domain.nix-1] *= 2.0;
    //field[domain.niy-1] *= 2.0;
}

//*******************************************************
void PushSpecies(Species *species, double *efx, double *efy)
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

        // advance velocity
        part.xvel += timeStep*qm*part_efx;
        part.yvel += timeStep*qm*part_efy;

        // Advance particle position
        part.xpos += timeStep*part.xvel;
        part.ypos += timeStep*part.yvel;


     /* if(part.xpos < domain.x0) part.xpos += domain.xl;
      if(part.ypos < domain.y0) part.ypos += domain.yl;

      if(part.xpos > domain.xmax) part.xpos -= domain.xl;
      if(part.ypos > domain.ymax) part.ypos -= domain.yl;*/

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
             part.vel = SampleVel(species->Temp*EV_TO_K, species->mass);
             species->add(Particle(part.pos,part.vel));*/

            //continue;
        //}
        else
            it++;
    }
}
//*********************************************************
/*Rewind particle velocities by -0.5*timeStep */
void RewindSpecies(Species *species, double *efx, double *efy)
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
        //advance velocity
        p.xvel -= 0.5*timeStep*qm*part_efx;
        p.yvel -= 0.5*timeStep*qm*part_efy;
    }
}

/* Compute the charge densities */
void ComputeRho(double *rho, Species *ions, Species *electrons)
{
    //double *rho = domain.rho;
    //memset(rho,0,sizeof(double)*domain.ni);

    for(int i=1; i<domain.nix-1; i++)
    for(int j=1; j<domain.niy-1; j++){
        rho[i*domain.niy+j]=ions->charge*ions->den[i*domain.niy+j] +
        electrons->charge*electrons->den[i*domain.niy+j];
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
    if(false){
        for(int i=0; i<domain.nix; i++)
        for(int j=0; j<domain.niy; j++)
            if(fabs(rho[i*domain.niy+j])<1e8*chargeE) rho[i*domain.niy+j]=0;
    }
}

/* Potential Solver: 1. Gauss-Seidel 2. Direct-Solver*/
bool SolvePotential(double *phi, double *rho)
{
   double dx2 = domain.dx*domain.dx;
   double dy2 = domain.dy*domain.dy;
   double L2;


   for(int solver=0; solver<200000; solver++)
   {                                                     // solve electric potential at nodes using
      for(int i=0; i<domain.nix; i++)                    // Gauss - Seidel method
      for(int j=0; j<domain.niy; j++){

           int p=i-1; if(p<0) p=domain.nix-2;
           int q=i+1; if(q>domain.nix-1) q=1;
           int r=j-1; if(r<0) r=domain.niy-2;
           int s=j+1; if(s>domain.niy-1) s=1;

   double g = 0.5*(1/((1/dx2)+(1/dy2)))*( ((phi[p*domain.niy+j]+phi[q*domain.niy+j])/dx2) + ((phi[i*domain.niy+r]+phi[i*domain.niy+s])/dy2) + (rho[i*domain.niy+j]/EPS) );

  // double R1 = 0.25*(phi[p][j]+phi[q][j]+phi[i][r]+phi[i][s]+(dx2*rho[i][j]/EPS))-phi[i][j];
  phi[i*domain.niy+j] = phi[i*domain.niy+j] + 1.4*(g - phi[i*domain.niy+j]);

   //phi[i][j] += R1;
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

         double R = 0.25*(phi[p*domain.niy+j]+phi[q*domain.niy+j]+phi[i*domain.niy+r]+phi[i*domain.niy+s]+(dx2*rho[i*domain.niy+j]/EPS))-phi[i*domain.niy+j];

         sum = sum + (R*R);
         }
         }

         L2 = sqrt(sum)/(domain.nix*domain.niy);
         if(L2<1e-2){return true;}
      }
   }

   printf("Gauss-Seidel solver failed to converge, L2=%.3g!\n",L2);
	return false;
}

/* Potential Direct Solver */

// bool SolvePotentialDirect(double *x, double *rho)
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
void ComputeEF(double *phi, double *efx, double *efy)
{
    /*Apply central difference to the inner nodes*/
    for(int i=1; i<domain.nix-1; i++)
    for(int j=1; j<domain.niy-1; j++)
    {
      efx[i*domain.niy+j] = (phi[(i-1)*domain.niy+j]-phi[(i+1)*domain.niy+j])/(2*domain.dx);
      efy[i*domain.niy+j] = (phi[i*domain.niy+j-1]-phi[i*domain.niy+j+1])/(2*domain.dy);
    }

    /*Apply one sided difference at the boundary nodes*/
    for(int j=0; j<domain.niy; j++)                                    // calculate electric field at nodes
   {
      //efx[0][j] = -(phi[1][j]-phi[domain.nix-2][j])/(2*domain.dx);
      //efx[domain.nix-1][j] = efx[0][j];
      efx[0*domain.niy+j] = -(phi[1*domain.niy+j] - phi[0*domain.niy+j])/(2*domain.dx);
      efx[(domain.nix-1)*domain.niy+j] = -(phi[(domain.nix-1)*domain.niy+j] - phi[(domain.nix-2)*domain.niy+j])/(2*domain.dx);

   }

   for(int i=0; i<domain.nix; i++)
   {
      //efy[i][0] = -(phi[i][1]-phi[i][domain.niy-2])/(2*domain.dy);
      //efy[i][domain.niy-1] = efy[i][0];

      efy[i*domain.niy+0] = -(phi[i*domain.niy+1] - phi[i*domain.niy+0])/(2*domain.dy);
      efy[i*domain.niy+domain.niy-1] = -(phi[i*domain.niy+domain.niy-1] - phi[i*domain.niy+domain.niy-2])/(2*domain.dy);
   }
}


/*Write the output with time*/
/*void Write_ts(int ts, Species *ions, Species *electrons)
{
    for(int j=0; j<domain.niy; j++){
    for(int i=0; i<domain.nix; i++)
    {
        fprintf(file_res,"%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", i*domain.dx, j*domain.dy, ions->den[i][j], electrons->den[i][j], ions->xvel[i][j], electrons->xvel[i][j], domain.rho[i][j], domain.phi[i][j], domain.efx[i][j],domain.efy[i][j]);

    }
    }
    //fprintf(file_res,"%g \t %g \t %g\n",ts*timeStep, gamma_i[domain.ni-1], gamma_e[domain.ni-1]);
    fflush(file_res);
}*/

// void WritePotOsc(double Time, int probLoc)
// {
//
//     fprintf(file_phi,"%g \t %g\n",Time, domain.phi[probLoc]);
//     fflush(file_phi);
// }

/* Write the Output results*/
void Write_Particle(FILE *file, int ts, Species *species)
{
    for(auto& p: species->part_list)
    {
        fprintf(file,"%g \t %g \t %g \t %g\n",p.xpos, p.ypos, p.xvel, p.yvel);
    }
    fflush(file_res);
}

/******* ADDED BY SAYAN 14/08/2019  *******/
// void Write_VDF(FILE *file, int ts, double vdfLocStart, double vdfLocEnd, Species *species)
// {
//     for(auto& p: species->part_list)
//     {
//         if (p.xpos >= vdfLocStart && p.xpos <= vdfLocEnd)
//         {
//         fprintf(file,"%g \t %g\n",p.xvel, p.yvel);
//         }
//     }
//     fflush(file_res);
// }

/* ********************************** */

void Write_Single_Particle(Species *species)
{
    list<Particle>::iterator it=species->part_list.begin();
    Particle &part = *it;
    for(int i=0; i<10; i++){it++;}
    fprintf(file_sp,"%g \t %g \t %g \t %g\n",part.xpos,part.ypos,part.xvel,part.yvel);
    fflush(file_sp);
}

void WriteKE(double Time, Species *ions, Species *electrons)
{
    double ke_ions = ComputeKE(ions);
    double ke_electrons = ComputeKE(electrons);

    fprintf(file_ke,"%g \t %g \t %g\n",Time, ke_ions, ke_electrons);

    fflush(file_ke);
}

double ComputeKE(Species *species)
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

void Write_pot(double Time)
{

    fprintf(file_res1,"%g\t %g\n",Time, domain.efx[4*domain.niy+4]);

    fflush(file_res1);
}

void ComputePE(double Time)
{
   double pe = 0;
   for(int j=0; j<domain.niy; j++)
   for(int i=0; i<domain.nix; i++)
   {
     pe += 0.5*domain.phi[i*domain.niy+j]*domain.rho[i*domain.niy+j];
   }
   fprintf(file_res2,"%g\t %g\n",Time, pe);

}
