#ifndef Vector3D
#include "encabezado.h" 
#endif

 void initial_conditions(data_t & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            data[ix*ny + iy].value = 1.0;
        }
    }
}

 void boundary_conditions(data_t & data, int nx, int ny, Body * N, double l, int Nmax)
{
    int ix, iy;
    Vector3D aux, aux1;
    double qaux;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {
        data[ix*ny + iy].value = 100.0;
	data[ix*ny + iy].ocupation = true;
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
    for(int ii=0;ii<Nmax;ii++)
    {
        aux=N[ii].getR();
        aux1=N[ii].getV();
        qaux=N[ii].getQ();
      if(aux1[0]==0 && aux1[1]==0 && qaux>0)
	{
      data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].value=100;
      data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].ocupation=true;
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
  //start_gnuplot();
    for(int istep = 0; istep < nsteps; ++istep) {
      relaxation_step(data, nx, ny,Q);
        //print_screen(data, nx, ny);
        //print_gnuplot(data, nx, ny);
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

void Get_Q(Body * N, data_q & Q, int nx, int ny, double l, int Nmax)
{
    Vector3D aux;
    double qaux;
  for(int ii=0;ii<Nmax;ii++)
  {
      aux=N[ii].getR();
      qaux=N[ii].getQ();
      Q[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)]+=qaux;
    }
}

void Get_EF(Body * N, int nx, int ny, double l, int Nmax, data_t & data, double Delta, double gamma)
{
  Vector3D aux, aux1, Faux;
  double qaux;
  int auxx, auxy;
  for(int ii=0;ii<Nmax;ii++)
    {
      aux=N[ii].getR();
      aux1=N[ii].getV();
      qaux=N[ii].getQ();
      N[ii].resetForce();
      if(aux1[0]==0 && aux1[1]==0 && qaux>0)
	{
	  
	}
      else{
	  auxx=int(nx*aux[0]/l);
	  auxy=int(ny*aux[1]/l);
	  Faux[0]=-(data[auxx*ny+auxy].value-data[(auxx+1)*ny+auxy].value)/Delta-gamma*aux1[0];
	  Faux[1]=-(data[auxx*ny+auxy].value-data[auxx*ny+(auxy+1)].value)/Delta-gamma*aux1[1];
	  N[ii].addForce(Faux);
      }
    }

  return;
}

void Update_boundary(Body * N, int nx, int ny, double l, int Nmax, data_t & data)
{
  Vector3D aux, aux1, Vaux;
  double qaux;
  int auxx, auxy;
  for(int ii=0;ii<Nmax;ii++)
    {
      aux=N[ii].getR();
      aux1=N[ii].getV();
      qaux=N[ii].getQ();
      auxx=int(nx*aux[0]/l);
      auxy=int(ny*aux[1]/l);
      if( qaux>0)
	{
	  if(data[(auxx+1)*ny+auxy].ocupation == true ||data[(auxx-1)*ny+auxy].ocupation == true || data[auxx*ny+(auxy+1)].ocupation == true ||  data[auxx*ny+(auxy-1)].ocupation == true)
	    {
	      Vaux[0]=0;
	      Vaux[1]=0;
	      N[ii].setV(Vaux);
	      data[(auxx)*ny+auxy].ocupation = true;
	    }
	}
    }

  return;
}

void update_and_check_pos(Body * N, int nx, int ny, double l, int Nmax, data_t & data, double mu, double sigma, double dt)
{
  Vector3D Rold, Rnew, Vold, Vnew, Dr, move;
  double qaux, move_x, move_y;
  int auxx, auxy;
  std::vector<int> cond (4,0);
  CRandom rand(42);
  for(int ii=0;ii<Nmax;ii++)
    {
      Rold=N[ii].getR();
      Vold=N[ii].getV();
      qaux=N[ii].getQ();
      Get_EF(N, nx, ny, l, Nmax, data, DELTA, 0.1);
      if(Vold[0]!=0 && Vold[1]!=0)
	{
	  N[ii].moveV(dt, 1);
	  N[ii].moveR(dt, 1);
	  move_x = rand.gauss(mu, sigma);
	  move_y = rand.gauss(mu, sigma);
	  move.load(move_x, move_y, 0);
	  N[ii].setR(N[ii].getR() + move);
	  Rnew=N[ii].getR();
	  Vnew=N[ii].getV();
	  if(Rnew[0]>=l) cond[0]=1;
	  else if(Rnew[1]>=l) cond[1]=1;
	  else if(Rnew[0]<=0) cond[2]=1;
	  else if(Rnew[1]<=0) cond[3]=1;
	  int b =std::accumulate(cond.begin(), cond.end(), 0);
	  if(b>0)
	    {
	      if(b>1){
		N[ii].setV(-1*Vnew);
		N[ii].setR(Rold);
	      }
	      else{
		for(int jj=0; jj<4; jj++)
		  {
		    if(cond[jj]==1 && jj%2!=0)
		      {
			Vnew[0]=-1*Vnew[0];
			N[ii].setV(Vnew);
			N[ii].setR(Rold);
		      }
		    else
		      {
			Vnew[1]=-1*Vnew[1];
			N[ii].setV(Vnew);
			N[ii].setR(Rold);
		      }
		  }
	      }
	    }
	  if (qaux>0 && cond[2]==1)
	    {
	      double daux=0.01;
	      Dr=Rnew-Rold;
	      for (int kk=0; kk<100; kk++)
		{
		  auxx=int(nx*(Rold[0]+kk*daux*Dr[0])/l);
		  auxy=int(ny*(Rold[1]+kk*daux*Dr[0])/l);
		  if(auxx <=nx && auxy <= ny && data[(auxx)*ny+auxy].ocupation == true)
		    {
		      Rold[0]=Rold[0]+(kk-1)*daux*Dr[0];
		      Rold[1]=Rold[1]+(kk-1)*daux*Dr[1];
		      Vold[0]=0;
		      Vold[1]=0;
		      auxx=int(nx*(Rold[0])/l);
		      auxy=int(ny*(Rold[1])/l);
		      N[ii].setR(Rold);
		      N[ii].setV(Vold);
		      data[(auxx)*ny+auxy].ocupation == true;
		      break;
		    }
		}
	    }
	}
    }
}

void print_fractal (int nx, int ny, data_t & data)
{
  for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
	  std::cout << ix << "  "<<iy<<"   "<< data[ix*ny + iy].ocupation << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}
