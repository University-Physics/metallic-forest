#include <iostream>
#include <vector>


class particula {
public:
  double x, y, q, Vx, Vy;
};
class potencial {
public:
  double value;
  bool ocupation;
};
typedef std::vector<potencial> data_t;
  typedef std::vector<particula> data_p;
 typedef std::vector<double> data_q;

void initial_conditions(data_t & data, int nx, int ny);
void boundary_conditions(data_t & data, int nx, int ny, data_p & N, double l);
void evolve(data_t & data, int nx, int ny, int nsteps, data_q & Q);
void relaxation_step(data_t & data, int nx, int ny, data_q & Q);
void print_screen(const data_t & data, int nx, int ny);
void start_gnuplot(void);
void print_gnuplot(const data_t & data, int nx, int ny);
void Get_Q(data_p & N, data_q & Q, int nx, int ny, double l);
