#include <iostream>
#include <vector>
#include <numeric>
#include "constantes.h"
#include "Body.h"
class potencial {
public:
  double value;
  bool ocupation;
  bool electrode;
};
typedef std::vector<potencial> data_t;
typedef std::vector<double> data_q;

void initial_conditions(data_t & data, int nx, int ny);
void boundary_conditions(data_t & data, int nx, int ny, Body * N, double l, int Nmax);
void evolve(data_t & data, int nx, int ny, int nsteps, data_q & Q);
void relaxation_step(data_t & data, int nx, int ny, data_q & Q);
void print_screen(const data_t & data, int nx, int ny);
void start_gnuplot(void);
void print_gnuplot(const data_t & data, int nx, int ny);
void Get_Q(Body * N, data_q & Q, int nx, int ny, double l, int Nmax);
void Get_EF(Body * N, int nx, int ny, double l, int Nmax, data_t & data, double gamma);
void Update_boundary(Body * N, int nx, int ny, double l, int Nmax, data_t & data);
void update_and_check_pos(Body * N, int nx, int ny, double l, int Nmax, data_t & data, double mu, double sigma, double dt);
void print_fractal (int nx, int ny, data_t & data);
