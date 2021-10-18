#include <iostream>
#include <vector>
#include <numeric>
#include <random>
#include <fstream>
#include <string>
#include <omp.h>
#include "constantes.h"
#include "Body.h"
class potencial {
public:
  double value;
  bool ocupation;
  bool electrode;
  bool pivot;
};
typedef std::vector<potencial> data_t;
typedef std::vector<double> data_q;
void print(data_q probability,std::string filename);
void initial_conditions(data_t & data, int nx, int ny);
void boundary_conditions(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff);
void boundary_conditions1(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
void boundary_conditions2(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
void boundary_conditions3(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
void boundary_conditions4(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
void boundary_conditions_wn(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff, double amp, int seed);
void boundary_conditions_pn(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff, double amp, int seed);
bool evolve(data_t & data, int nx, int ny, int nsteps, int ns_est);
void relaxation_step_pivot(data_t & data, int nx, int ny);
void stabilization_step(data_t & data, int nx, int ny);
bool relaxation_step(data_t & data, int nx, int ny);
void print_screen(const data_t & data, int nx, int ny);
void start_gnuplot(void);
void print_gnuplot(const data_t & data, int nx, int ny);
void Get_Q(Body * N, data_q & Q, int nx, int ny, double l, int Nmax);
void deposited_particles(data_t & data, int nx, int ny, Body * N, double l, int Nmax);
void Get_EF(Body * N, int nx, int ny, int Nmax, data_t & data, double Delta, double gamma, bool interacciones);
bool Update_boundary(Body * N, int nx, int ny, int Nmax, data_t & data);
void update_and_check_pos2(Body * N, int nx, int ny, double l, int Nmax, data_t & data, double mu, double sigma, double dt, double coefx, double coefv, int seed, bool interacciones);
void print_fractal (int nx, int ny, data_t & data, std::string filename);
void evolve_system(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, double coefx, double coefv, int seed, double V_diff, bool interacciones);
void PEFRL(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, int seed, double V_diff, bool interacciones);
void evolve_opt(data_t & data, int nx, int ny, double l, int Nmax, Body * N, double V_diff);
void print_potential (int nx, int ny, data_t & data, std::string filename);
void print_potential_size (int nx, int ny, data_t & data, std::string filename, int t);
