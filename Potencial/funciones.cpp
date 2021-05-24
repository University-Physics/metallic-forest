#ifndef Vector3D
#include "encabezado.h" 
#endif

 void initial_conditions(data_t & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
	  data[ix*ny + iy].value = 0.0;
	  data[ix*ny+iy].ocupation= false;
	  data[ix*ny+iy].electrode= false;
        }
    }
}

 void boundary_conditions(data_t & data, int nx, int ny, Body * N, double l, int Nmax)
{
    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {
        data[ix*ny + iy].value = -1.0;
	data[ix*ny + iy].electrode = true;
	data[ix*ny + iy].ocupation = true;
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {
      data[ix*ny + iy].value = 1.0;
      data[ix*ny + iy].electrode = true;
      data[ix*ny + iy].ocupation = false;
    }
    // first row
    
    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {
        data[ix*ny + iy].value = 1.0;
	data[ix*ny + iy].electrode = true;
	data[ix*ny + iy].ocupation = false;
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {
      data[ix*ny + iy].value = 1.0;
      data[ix*ny + iy].electrode = true;
      data[ix*ny + iy].ocupation = false;
    }
    
    for(int ii=0;ii<Nmax;ii++)
    {
      aux=N[ii].getR();
      if(N[ii].getoc()==true && data[int(nx*aux[0]/l)*ny + int(ny*aux[1]/l)].electrode == false)
	{
	  data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].value=-1.0;
	  data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].ocupation=true;
	}
    }
}

void evolve(data_t & data, int nx, int ny, int nsteps, int ns_est)
{
  //start_gnuplot();
    for(int istep = 0; istep < nsteps; ++istep) {
      relaxation_step(data, nx, ny);
      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
    for(int istep = 0; istep < nsteps; ++istep) {
      relaxation_step_pivot(data, nx, ny);
      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
    for(int istep = 0; istep < ns_est; ++istep){
      stabilization_step(data, nx, ny);
      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
    for(int istep = 0; istep < nsteps; ++istep) {
      relaxation_step(data, nx, ny);
      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
}

void stabilization_step(data_t & data, int nx, int ny)
{
  for(int ix = 1; ix < nx-1; ++ix) {
        for(int iy = 1; iy < ny-1; ++iy) {
            // check that this cell is NOT a boundary condition or a border
            //if ( (ix == nx/2) && (ny/3 <= iy) && (iy <= 2*ny/3) ) continue;
            // update the cell
	  if(data[ix*ny + iy].ocupation==false && data[ix*ny + iy].pivot==false)
	    {
            data[ix*ny + iy].value = (data[(ix+1)*ny + iy].value + data[(ix-1)*ny + iy].value + data[ix*ny + iy+1].value + data[ix*ny + iy-1].value)/4.0;
	    }
        }
    }
}

 void relaxation_step_pivot(data_t & data, int nx, int ny)
{
    // recorrer toda la matriz y aplicar el algoritmo,
    // teniendo cuidado con no modificar las condiciones de
    // frontera
  for(int ix = 1; ix < nx/10; ++ix) {
    for(int iy = 1; iy < ny/10; ++iy) {
            // check that this cell is NOT a boundary condition or a border
            //if ( (ix == nx/2) && (ny/3 <= iy) && (iy <= 2*ny/3) ) continue;
            // update the cell
      if(data[(ix*10)*ny + (iy*10)].ocupation==false)
	{
	  data[(ix*10)*ny + (iy*10)].value = (data[((ix+1)*10)*ny + (iy*10)].value + data[(ix-1)*10*ny + (iy*10)].value + data[(ix*10)*ny + (iy+1)*10].value + data[(ix*10)*ny + 10*(iy-1)].value)/4.0;
	  data[(ix*10)*ny+(iy*10)].pivot=true;
	}
    }
  }

}

void relaxation_step(data_t & data, int nx, int ny)
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
            data[ix*ny + iy].value = (data[(ix+1)*ny + iy].value + data[(ix-1)*ny + iy].value + data[ix*ny + iy+1].value + data[ix*ny + iy-1].value)/4.0;
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

void Get_Q(Body * N, data_q & Q, int nx, int ny, double l, int Nmax)
{
    Vector3D aux;
    double qaux;
  for(int ii=0;ii<Nmax;ii++)
  {
      aux=N[ii].getR();
      qaux=N[ii].getQ();
      aux[0]-=DELTA;
      aux[1]-=DELTA;
      Q[(int((nx-2)*aux[0]/l)+1)*ny+int((ny-2)*aux[1]/l)+1]+=qaux;
    }
}

void Get_EF(Body * N, int nx, int ny, int Nmax, data_t & data, double Delta, double gamma)
{
  Vector3D aux, aux1, aux2, Faux, dr;
  double q1, q2;
  int auxx, auxy;
  bool auxoc;
  
  for(int ii=0;ii<Nmax;ii++)
    {
      if(N[ii].getoc()==true)
	{
	  N[ii].resetForce();
	}
      else{
	aux=N[ii].getR();
	aux1=N[ii].getV();
	q1=N[ii].getQ();
	N[ii].resetForce();
	auxx=int(aux[0]/DELTA);
	auxy=int(aux[1]/DELTA);
	if(1<=auxx<nx-1 && 1<=auxy<ny-1)
	  {
	    /*
	    for(int jj = 0; jj<ii; jj++)
	      {
		
		if(N[jj].getoc()==false){
		  aux2=N[jj].getR();
		  dr=aux2-aux1;
		  q2=N[jj].getQ();
		  Faux=q1*q2*dr/(norm(dr)*norm(dr)*norm(dr));
		  N[ii].addForce(Faux);
		  N[jj].addForce((-1)*Faux);
		}
	      }
	    */
	    Faux[0]=q1*(data[(auxx-1)*ny+auxy].value-data[(auxx+1)*ny+auxy].value)/(2*Delta)-gamma*aux1[0];
	    Faux[1]=q1*(data[auxx*ny+(auxy-1)].value-data[auxx*ny+(auxy+1)].value)/(2*Delta)-gamma*aux1[1];
	    N[ii].addForce(Faux);
	  }
      }
    }

  return;
}

bool Update_boundary(Body * N, int nx, int ny, int Nmax, data_t & data)
{
  Vector3D aux, aux1, Vaux;
  double qaux;
  int auxx, auxy;
  bool cond=false;
  for(int ii=0;ii<Nmax;ii++)
    {
      aux=N[ii].getR();
      aux1=N[ii].getV();
      qaux=N[ii].getQ();
      auxx=int(aux[0]/DELTA);
      auxy=int(aux[1]/DELTA);
      if( qaux>0 && 1<=auxx<nx-1 && 1<=auxy<ny-1)
	{
	  if(data[(auxx+1)*ny+auxy].ocupation == true ||data[(auxx-1)*ny+auxy].ocupation == true || data[auxx*ny+(auxy+1)].ocupation == true ||  data[auxx*ny+(auxy-1)].ocupation == true)
	    {
	      cond=true;
	      Vaux[0]=0;
	      Vaux[1]=0;
	      N[ii].setV(Vaux);
	      N[ii].setoc(true);
	    }
	}
    }

  return cond;
}
void update_and_check_pos2(Body * N, int nx, int ny, int Nmax, data_t & data, double mu, double sigma, double dt, double coefx, double coefv, int seed)
{
  Vector3D Rnew, Vnew, move;
  double move_x, move_y, d_auxx, d_auxy;
  int auxx, auxy;
  CRandom rand(seed);
  Get_EF(N, nx, ny, Nmax, data, DELTA, 1);
  for(int ii=0;ii<Nmax;ii++)
    {
      if(N[ii].getoc()==false)
	{
	  N[ii].moveV(dt, coefx);
	  N[ii].moveR(dt, coefv);
	  move_x = rand.gauss(mu, sigma);
	  move_y = rand.gauss(mu, sigma);
	  move.load(move_x, move_y, 0);
	  N[ii].setR(N[ii].getR() + (coefx)* move);
	  Rnew=N[ii].getR();
	  Vnew=N[ii].getV();
	  auxx=int((Rnew[0]-DELTA)/DELTA);
          auxy=int((Rnew[1]-DELTA)/DELTA);
	  d_auxx=(Rnew[0]-DELTA)-auxx*DELTA;
	  d_auxy=(Rnew[1]-DELTA)-auxy*DELTA;
	  if(auxx<0)auxx=(nx-2)-((-auxx)%(nx-2));
          if(auxy<0)auxy=(ny-2)-((-auxy)%(ny-2));
          auxx=(auxx)%(nx-2);
	  auxy=(auxy)%(ny-2);
	  Rnew[0]=(auxx+1)*DELTA+d_auxx;
	  Rnew[1]=(auxy+1)*DELTA+d_auxy;
	  N[ii].setV(Vnew);
	  N[ii].setR(Rnew);
	}
    }
}

void evolve_system(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, double coefx, double coefv, int seed)
{
  bool cond;
  //update positions and boundaries
  update_and_check_pos2(N, nx, ny, Nmax, data, mu, sigma, dt, coefx, coefv, seed);
  cond=Update_boundary(N, nx, ny, Nmax, data);
  //Obtain potential using fast algorithm
  if(cond==true)
    {
      evolve_opt(data, nx, ny, l, Nmax, N);
    }
}

void PEFRL(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, int seed)
{
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Zi, 0.0, seed);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Xi, coef1, seed);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, coef2, Lambda, seed);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Xi, Lambda, seed);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Zi, coef1, seed);
}

void print_fractal (int nx, int ny, data_t & data)
{
  std::ofstream myfile;
  myfile.open ("data.txt");
  for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
	  myfile<< ix << "  "<<iy<<"   "<< data[ix*ny + iy].ocupation << "\n";
        }
        myfile << "\n";
    }
    myfile << "\n";
    myfile.close();
}

void print_potential (int nx, int ny, data_t & data)
{
  std::ofstream myfile;
  myfile.open ("pot.txt");
  for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
	  myfile<< ix << "\t"<<iy<<"\t"<< data[ix*ny + iy].value << "\n";
        }
    }
    myfile << "\n";
    myfile.close();
}

void evolve_opt(data_t & data, int nx, int ny, double l, int Nmax, Body * N)
{
  //doesn't reset values in order to obtain faster convergence and increase performance
  boundary_conditions(data, nx, ny, N ,l, Nmax);
  // evolve
  evolve(data, nx, ny, NSTEPS/10, NSTEPS/10);
}
