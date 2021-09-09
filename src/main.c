#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

/*constants*/
#define EPS_0 8.85418782e-12  	// F/m, vacuum permittivity
#define K	1.38065e-23			// J/K, Boltzmann constant
#define ME 9.10938215e-31		// kg, electron mass
#define QE 1.602176565e-19		// C, electron charge
#define AMU  1.660538921e-27	// kg, atomic mass unit

/*simulation parameters, these could come from an input file*/
#define PLASMA_DEN	1e14		// plasma density to load
#define NUM_IONS 100000		// number of ions
#define NUM_ELECTRONS 200000	    // number of electrons
#define DX 1e-4					// cell spacing
#define NC 199				// number of cells
#define NUM_TS	5000			// number of time steps
#define DT 1e-12				// time step size

/* Data structure to hold domain information*/
struct Domain
{
	int ni;				/*number of nodes*/
	double x0;			/*mesh origin*/
	double dx;			/*cell spacing*/
	double xl;			/*domain length*/
	double xmax;		/*domain max position*/

/*data structures*/
	double *phi;		/*potential*/
	double *ef;			/*electric field*/
	double *rho;		/*charge density*/
	double *ndi;		/*ion density*/
	double *nde;		/*electron density*/
};

/* Data structure for particle storage **/
struct Particle
{
	double x;			/*position*/
	double v;			/*velocity*/
};

/* Data structure to hold species information*/
struct Species
{
	double mass;			/*particle mass in kg*/
	double charge;			/*particle charge in Coulomb*/
	double spwt;			/*species specific weight*/

	int np;					/*number of particles*/
	int np_alloc;			/*size of the allocated data array*/
  struct Particle *part;			/*array holding particles*/
};

/** FUNCTION PROTOTYPES **/
void ScatterSpecies(struct Species *species, double *nd);
void ComputeRho(struct Species *ions, struct Species *electrons);
bool SolvePotential(double *phi, double *rho);
void ComputeEF(double *phi, double *ef);
void PushSpecies(struct Species *species, double *ef);
void RewindSpecies(struct Species *species, double *ef);
void AddParticle(struct Species *species, double x, double v);
double XtoL(double pos);
void scatter(double lc, double value, double *field);
double gather(double lc, double *field);
void WriteResults(int ts);
double ComputeKE(struct Species *species);	/*computes kinetic energy*/

/* GLOBAL VARIABLES */
struct Domain domain;
FILE *file_res;
FILE *file_pot;
/* --------- main -------------*/

int main()
{
	int i,p;
	int ts;						// time step

	clock_t start = clock();	// grab starting clock time

	/*---1) initialize domain---*/
	domain.ni = NC+1;			// number of nodes
	domain.dx = DX;				// cell spacing
	domain.x0 = 0;					// origin
	domain.xl = (domain.ni-1)*domain.dx;//domain length
	domain.xmax = domain.x0+domain.xl;	//max position

	/*allocate data structures, remember these need to be cleared first in C++*/
   domain.phi = (double*)malloc(domain.ni * sizeof(double));	//potential
   domain.rho = (double*)malloc(domain.ni * sizeof(double));
   domain.ef = (double*)malloc(domain.ni * sizeof(double));
   domain.nde = (double*)malloc(domain.ni * sizeof(double));
   domain.ndi = (double*)malloc(domain.ni * sizeof(double));

	// /*save pointers so we can write phi instead of domain.phi*/
	double *phi = domain.phi;
	double *rho = domain.rho;
	double *ef = domain.ef;
	double *nde = domain.nde;
	double *ndi = domain.ndi;

   struct Species ions;
	struct Species electrons;

	/*set material data*/
	ions.mass = 8*AMU;
	ions.charge = QE;
	ions.spwt = PLASMA_DEN*domain.xl/NUM_IONS;
	ions.np = 0;
	ions.np_alloc = NUM_IONS;
	ions.part = (struct Particle*)malloc(NUM_IONS * sizeof(struct Particle));

	electrons.mass = ME;	// electrons
	electrons.charge = -QE;
	electrons.spwt = PLASMA_DEN*domain.xl/NUM_ELECTRONS;
	electrons.np = 0;
	electrons.np_alloc = NUM_ELECTRONS;
	electrons.part = (struct Particle*)malloc(NUM_ELECTRONS * sizeof(struct Particle));

  double delta_ions = domain.xl/NUM_IONS;
	for (p=0;p<NUM_IONS;p++)
	{
		double x = domain.x0 + p*delta_ions;
		double v = 0;	/*initially zero velocity*/
		AddParticle(&ions,x,v);
	}

	/*now do the same for electrons*/
	double delta_electrons = domain.xl/NUM_ELECTRONS;
	for (p=0;p<NUM_ELECTRONS;p++)
	{
		double x = domain.x0 + p*delta_electrons;
		double v = 0;
		AddParticle(&electrons,x,v);
	}

  ScatterSpecies(&ions,ndi);
	ScatterSpecies(&electrons,nde);

  ComputeRho(&ions, &electrons);
	SolvePotential(phi,rho);
	ComputeEF(phi,ef);

  RewindSpecies(&ions,ef);
	RewindSpecies(&electrons,ef);

  file_res = fopen("results.dat","w");
	fprintf(file_res,"VARIABLES = x nde ndi rho phi ef\n");
	file_pot = fopen("potential.dat","w");
	WriteResults(0);

  for (ts = 1; ts<=NUM_TS; ts++)
	{
		//compute number density
		ScatterSpecies(&ions,ndi);
		ScatterSpecies(&electrons,nde);

		ComputeRho(&ions, &electrons);
		SolvePotential(phi,rho);
		ComputeEF(phi,ef);

		//move particles
		PushSpecies(&electrons,ef);
		PushSpecies(&ions,ef);

		//write diagnostics
		if (ts%500==0)
		{
			//max phi
			double max_phi = phi[0];
			for (int i=0;i<domain.ni;i++)
				if (phi[i]>max_phi) max_phi=phi[i];

			//compute kinetic energy
			double ke_ions = ComputeKE(&ions);
			double ke_electrons = ComputeKE(&electrons);

			printf("TS: %i, ions: (np=%d, KE=%.3g)\t electrons: (np=%d, KE=%.3g)\t delta phi: %.3g\n",
				ts,ions.np,ke_ions,electrons.np,ke_electrons,max_phi-phi[0]);

	      fprintf(file_pot, "%g \t %g\n", DT*ts, phi[5]);
		}

		//save animation data
		if (ts%500==0) WriteResults(ts);
	}

  fclose(file_res);
	fclose(file_pot);

  free(phi);
	free(rho);
  free(ef);
	free(nde);
  free(ndi);
	// /*free particles*/
	free(ions.part);
	free(electrons.part);

	return 0;
}




/*adds new particle to the species*/
void AddParticle(struct Species *species, double x, double v)
{
	/*abort if we ran out of space to store this particle*/
	if (species->np > species->np_alloc-1)
	{
		printf("Too many particles!\n"); exit(-1);
	}

	/*store position and velocity of this particle*/
	species->part[species->np].x = x;
	species->part[species->np].v = v;

	/*increment particle counter*/
	species->np++;
}

double XtoL(double pos)
{
	double li = (pos-domain.x0)/domain.dx;
	return li;
}

/* scatters scalar value onto a field at logical coordinate lc*/
void scatter(double lc, double value, double *field)
{
	int i = (int)lc;
	double di = lc-i;

	field[i] += value*(1-di);
	field[i+1] += value*(di);
}

double gather(double lc, double *field)
{
	int i = (int)lc;
	double di = lc-i;

	/*gather field value onto particle position*/
	double val = field[i]*(1-di) + field[i+1]*(di);
	return val;
}


/*scatter particles of species to the mesh*/
void ScatterSpecies(struct Species *species, double *den)
{
	/*initialize densities to zero*/
	memset(den,0,sizeof(double)*domain.ni);

	/*scatter particles to the mesh*/
	for (int p=0;p<species->np;p++)
		{
			double lc = XtoL(species->part[p].x);
			scatter(lc, species->spwt, den);
		}

	/*apply periodic boundaries*/
	den[0]+=den[domain.ni-1];
	den[domain.ni-1]=den[0];

	/*divide by cell volume*/
	for (int i=0;i<domain.ni;i++)
		den[i]/=domain.dx;
}

/*computes charge density by adding ion and electron data*/
void ComputeRho(struct Species *ions, struct Species *electrons)
{
	double *rho = domain.rho;

	for (int i=0;i<domain.ni;i++){
		rho[i] = ions->charge*domain.ndi[i] + electrons->charge*domain.nde[i];
		//printf("%g\n",rho[i]);
	}
	/*remove numerical noise by setting zero when densities less than 1e8/m^3*/
	if (false)
	{
		for (int i=0;i<domain.ni;i++)
			if (fabs(rho[i])<1e8*QE) rho[i]=0;
	}

}

/* solves potential using the Gauss Seidel Method, returns true if converged
periodic boundaries version*/
bool SolvePotential(double *phi, double *rho)
{
	double L2;
	double dx2 = domain.dx*domain.dx;	/*precompute*/

	/*initialize to zero, this is only needed because of periodic b.c. without a dirichlet*/
	memset(phi, 0, sizeof(double)*domain.ni);

	/*solve potential, identical to lesson 2*/
	for (int solver_it=0;solver_it<5000;solver_it++)
	{
		/*Gauss Seidel method, phi[i-1]-2*phi[i]+phi[i+1] = -dx^2*rho[i]/eps_0*/
		for (int i=0;i<domain.ni;i++)
		{
			int im = i-1;	if (im<0) im=domain.ni-2;
			int ip = i+1;   if (ip>domain.ni-1) ip=1;
			phi[i] = 0.5*(phi[im] + phi[ip] + dx2*rho[i]/EPS_0);
		}

		/*check for convergence*/
		if (solver_it%25==0)
		{
			double sum = 0;
			for (int i=0;i<domain.ni;i++)
			{
				int im = i-1;	if (im<0) im=domain.ni-2;
				int ip = i+1;   if (ip>domain.ni-1) ip=1;
				double R = -rho[i]/EPS_0 - (phi[im] - 2*phi[i] + phi[ip])/dx2;
				sum+=R*R;
			}
			L2 = sqrt(sum)/domain.ni;
			if (L2<1e-6) {return true;}
		}
	}
	printf("Gauss-Seidel solver failed to converge, L2=%.3g!\n",L2);
	return false;
}

/* computes electric field by differentiating potential*/
void ComputeEF(double *phi, double *ef)
{
	for (int i=1;i<domain.ni-1;i++)
		ef[i] = -(phi[i+1]-phi[i-1])/(2*domain.dx);	//central difference

	/*periodic boundaries*/
	ef[0] = -(phi[1]-phi[domain.ni-2])/(2*domain.dx);
	ef[domain.ni-1] = ef[0];
}

void RewindSpecies(struct Species *species, double *ef)
{
	/*precompute q/m*/
	double qm = species->charge / species->mass;

	for (int p=0;p<species->np;p++)
	{
		struct Particle *part = &species->part[p];

		/*compute particle node position*/
		double lc = XtoL(part->x);

		/*gather electric field onto particle position*/
		double part_ef = gather(lc,ef);

		/*advance velocity*/
		part->v -= 0.5*DT*qm*part_ef;
	}
}

/* moves particles of a single species*/
void PushSpecies(struct Species *species, double *ef)
{
	/*precompute q/m*/
	double qm = species->charge / species->mass;

	/*loop over particles*/
	for (int p=0;p<species->np;p++)
	{
		/*grab pointer to this particle*/
		struct Particle *part = &species->part[p];

		/*compute particle node position*/
		double lc = XtoL(part->x);

		/*gather electric field onto particle position*/
		double part_ef = gather(lc,ef);

		/*advance velocity*/
		part->v += DT*qm*part_ef;

		/*advance position*/
		part->x += DT*part->v;

		/*apply periodic boundaries*/
		if (part->x < domain.x0)	part->x+=domain.xl;
		if (part->x > domain.xmax)	part->x-=domain.xl;
	}
}
void WriteResults(int ts)
{
	fprintf(file_res,"ZONE I=%d T=ZONE_%06d\n",domain.ni,ts);
	for (int i=0;i<domain.ni;i++)
	{
		fprintf(file_res,"%g %g %g %g %g %g\n",i*domain.dx,
									domain.nde[i],domain.ndi[i],
									domain.rho[i],domain.phi[i],domain.ef[i]);
	}

	fflush(file_res);
}

/*computes species kinetic energy in electron volt*/
double ComputeKE(struct Species *species)
{
	double ke=0;
	for (int p=0;p<species->np;p++)
		ke += species->part[p].v*species->part[p].v;

	/*we now have sum of v^2, multiply by 0.5*mass*/
	ke *= 0.5*(species->spwt*species->mass);

	/*convert to electron volts, 1eV=QE joules*/
	ke /=QE;

	return ke;
}
