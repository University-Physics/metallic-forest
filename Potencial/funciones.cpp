#include"constantes.h"
#include"encabezado.h"

void initial_conditions(data_t & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            data[ix*ny + iy].value = 1.0;
        }
    }
}

void boundary_conditions(data_t & data, int nx, int ny, data_p & N, double l)
{
    int ix, iy;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {
        data[ix*ny + iy].value = 100.0;
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {
        data[ix*ny + iy].value = 0.0;
    }
    // first row
    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {
        data[ix*ny + iy].value = 0.0;
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {
        data[ix*ny + iy].value = 0.0;
    }
    for(int ii=0;ii<N.size();ii++)
    {
      if(N[ii].Vx==0 && N[ii].Vy==0 && N[ii].q<0)
	{
      data[int(nx*N[ii].x/l)*ny+int(ny*N[ii].y/l)].value=100;
      data[int(nx*N[ii].x/l)*ny+int(ny*N[ii].y/l)].ocupation=true;
	}
	}
    //new
    //ix = nx/2;
    //for(int iy = ny/3; iy <= 2*ny/3; ++iy) {
    //    data[ix*ny + iy] = -50.0;
    //}
}

void evolve(data_t & data, int nx, int ny, int nsteps, data_q & Q)
{
    start_gnuplot();
    for(int istep = 0; istep < nsteps; ++istep) {
      relaxation_step(data, nx, ny,Q);
        //print_screen(data, nx, ny);
        print_gnuplot(data, nx, ny);
    }
}

void relaxation_step(data_t & data, int nx, int ny, data_q & Q)
{
    // recorrer toda la matriz y aplicar el algoritmo,
    // teniendo cuidado con no modificar las condiciones de
    // frontera
    for(int ix = 1; ix < nx-1; ++ix) {
        for(int iy = 1; iy < ny-1; ++iy) {
            // check that this cell is NOT a boundary condition or a border
            //if ( (ix == nx/2) && (ny/3 <= iy) && (iy <= 2*ny/3) ) continue;
            // update the cell
	  if(data[ix*ny + iy].ocupation==false)
	    {
            data[ix*ny + iy].value = (data[(ix+1)*ny + iy].value + data[(ix-1)*ny + iy].value + data[ix*ny + iy+1].value + data[ix*ny + iy-1].value)/4.0+Q[ix*ny+iy];
	    }
        }
    }

}

void print_screen(const data_t & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            std::cout << data[ix*ny + iy].value << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


void start_gnuplot(void)
{
    std::cout << "set pm3d\n";
    std::cout << "set contour base\n";
    std::cout << "set term gif animate\n";
    std::cout << "set output 'anim.gif'\n";
}

void print_gnuplot(const data_t & data, int nx, int ny)
{
    std::cout << "splot '-' w l lt 3 \n";
    for(int ix = 0; ix < nx; ++ix) {
        double x = XMIN + ix*DELTA;
        for(int iy = 0; iy < ny; ++iy) {
            double y = YMIN + iy*DELTA;
            std::cout << x << "  " << y << "  " << data[ix*ny + iy].value << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "e\n";
}

void Get_Q(data_p & N, data_q & Q, int nx, int ny, double l)
{
  for(int ii=0;ii<N.size();ii++)
    {
      Q[int(nx*N[ii].x/l)*ny+int(ny*N[ii].y/l)]+=N[ii].q;
    }
}
