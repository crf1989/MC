#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct
{
  double s[3]; /* sx, sy, sz */
  int NNeighbor[6]; /* Nearest neighbors */
  int NNNeighbor[12]; /* Next nearest neighbors */
} site;

site* lattice; 

/* parameters */
const double J1 = -0.037*1.0;
const double J2 = 0.069*1.0;
const double g = 2;
const double u = 0.6717139; /* this is \mu_B/k_B */
const double S = 3.968626966596886; /* this spin S is in fact sqrt(S(S+1)),
					i.e. sqrt(3.5*4.5) */
const double PI = 3.141592653589793;

double T = 0; /* temperature */
double B = 0; /* magnetic field, Bz */
int L = 16; /* the side length of cubic lattice */
int N; /* N is the number of sites in lattice */

double energy; 
double total_energy;
double energy_square;
double spin[3]; 
double total_spin[3];

/* this function generate a uniform random number between 0 and 1 */
double drand ()
{
  return ((double)random ())/RAND_MAX;
}

/* this function generate a uniform random number between -1 and 1 */
double urand ()
{
  return 2*(drand ()-0.5);
}

/* this function generate a random vector of length S */
void genvec (double* x, double* y, double* z)
{
  double theta = asin (urand ()) + PI/2; /* theta \in [0, PI] */
  double phi = 2*PI * drand (); /* phi \in [0, 2PI] */
  
  *x = S * sin (theta) * cos (phi);
  *y = S * sin (theta) * sin (phi);
  *z = S * cos (theta);
}

void position (int m, int* i, int* j, int* k)
{
  *k = m % L;
  *j = (m % (L*L)) / L;
  *i = m / (L*L);
}

int Index (int i, int j, int k)
{
  if (i < 0) i += L;
  if (i >= L) i -= L;
  if (j < 0) j += L;
  if (j >= L) j -= L;
  if (k < 0) k += L;
  if (k >= L) k -= L;
  
  return i*(L*L) + j*L + k;
}

void init ()
{
  srand (time (0));
  N = L*L*L;
  lattice = malloc (N * sizeof(site));
  
  for (int i = 0; i < N; ++i)
    {
      genvec (&lattice[i].s[0], &lattice[i].s[1], &lattice[i].s[2]);
      /* set the nearest neighbors */
      int a, b, c;
      position (i, &a, &b, &c);
      lattice[i].NNeighbor[0] = Index (a-1, b, c);
      lattice[i].NNeighbor[1] = Index (a+1, b, c);
      lattice[i].NNeighbor[2] = Index (a, b-1, c);
      lattice[i].NNeighbor[3] = Index (a, b+1, c);
      lattice[i].NNeighbor[4] = Index (a, b, c-1);
      lattice[i].NNeighbor[5] = Index (a, b, c+1);
      /* set the next nearest neighbors */
      lattice[i].NNNeighbor[0] = Index (a-1, b-1, c);
      lattice[i].NNNeighbor[1] = Index (a-1, b+1, c);
      lattice[i].NNNeighbor[2] = Index (a+1, b-1, c);
      lattice[i].NNNeighbor[3] = Index (a+1, b+1, c);
      lattice[i].NNNeighbor[4] = Index (a-1, b, c-1);
      lattice[i].NNNeighbor[5] = Index (a-1, b, c+1);
      lattice[i].NNNeighbor[6] = Index (a+1, b, c-1);
      lattice[i].NNNeighbor[7] = Index (a+1, b, c+1);
      lattice[i].NNNeighbor[8] = Index (a, b-1, c-1);
      lattice[i].NNNeighbor[9] = Index (a, b-1, c+1);
      lattice[i].NNNeighbor[10] = Index (a, b+1, c-1);
      lattice[i].NNNeighbor[11] = Index (a, b+1, c+1);
    }
}

void metropolis_step ()
{
  int m = rand() % N;
  double s_new[3];
  genvec (&s_new[0], &s_new[1], &s_new[2]);

  /* energy change due to magnetic field */
  double delta_e = -g*u*B * (s_new[2]-lattice[m].s[2]);
  
  /* energy change due to nearest neighbors */
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 3; ++j)
      delta_e += -J1 * (s_new[j]-lattice[m].s[j]) * lattice[lattice[m].NNeighbor[i]].s[j];

  /* energy change due to next nearest neighbors */
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 3; ++j)
      delta_e += -J2 * (s_new[j]-lattice[m].s[j]) * lattice[lattice[m].NNNeighbor[i]].s[j];

  if (delta_e <= 0 || drand () < exp (-delta_e/T))
    {
      for (int i = 0; i < 3; ++i)
	{
	  spin[i] += s_new[i] - lattice[m].s[i];
	  lattice[m].s[i] = s_new[i];
	}
      energy += delta_e;
    }
  total_energy += energy;
  energy_square += energy * energy;
  for (int i = 0; i < 3; ++i)
    total_spin[i] += spin[i];
}

void monte_carlo_step ()
{
  for (int i = 0; i < N; ++i)
    metropolis_step ();
}

void thermalize (int steps)
{
  for (int i = 0; i < steps; ++i)
    monte_carlo_step ();
  
  /* initial the total energy and spins */
  energy = 0;
  for (int i = 0; i < 3; ++i)
    spin[i] = 0;
  
  for (int i = 0; i < N; ++i)
    {
      energy += -g*u*B * lattice[i].s[2];
      
      for (int j = 0; j < 6; ++j)
	for (int k = 0; k < 3; ++k)
	  energy += -J1 * lattice[i].s[k] * lattice[lattice[i].NNeighbor[j]].s[k];

      for (int j = 0; j < 12; ++j)
	for (int k = 0; k < 3; ++k)
	  energy += -J2 * lattice[i].s[k] * lattice[lattice[i].NNNeighbor[j]].s[k];

      for (int j = 0; j < 3; ++j)
	spin[j] += lattice[i].s[j];
    }

  energy /= 2; /* double counted */
  total_energy = energy;
  energy_square = energy * energy;
  for (int i = 0; i < 3; ++i)
    total_spin[i] = spin[i];
}

int main (int argc, char* argv[])
{
  B = atof (argv[1]);
  init ();

  printf ("#B = %g, N = L^3 = %d\n", B, N);
  double k = 1.3806488;
  double Na = 6.02114129;

  int steps = 1000000;
  for (T = 1; T <= 15; T += 0.5)
    {
      thermalize (1000);
      for (int i = 0; i < steps; ++i)
	monte_carlo_step ();
      
      total_energy /= ((double)steps)*N;
      total_energy *= total_energy;
      energy_square /= ((double)steps)*N;

      printf ("%g\t%g\n", T, k*Na*(energy_square-total_energy)/(T*T*N));
    }

  return 0;
}
