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
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {
        data[ix*ny + iy].value = 1.0;
	data[ix*ny + iy].ocupation = true;
	data[ix*ny + iy].electrode = true;
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {
      data[ix*ny + iy].value = 0.0;
      data[ix*ny + iy].electrode = true;
    }
    // first row
    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {
        data[ix*ny + iy].value = 0.0;
	data[ix*ny + iy].electrode = true;
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {
      data[ix*ny + iy].value = 0.0;
      data[ix*ny + iy].electrode = true;
    }
    for(int ii=0;ii<Nmax;ii++)
    {
      aux=N[ii].getR();
      if(N[ii].getoc()==true && data[int(nx*aux[0]/l)*ny + int(ny*aux[1]/l)].electrode == false)
	{
	  data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].value=1.0;
	  data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].ocupation=true;
	}
    }
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
      Q[(int((nx-2)*aux[0]/l)+1)*ny+int((ny-2)*aux[1]/l)+1]+=qaux;
    }
}

void Get_EF(Body * N, int nx, int ny, double l, int Nmax, data_t & data, double Delta, double gamma)
{
  Vector3D aux, aux1, Faux;
  double qaux;
  int auxx, auxy;
  bool auxoc;
  for(int ii=0;ii<Nmax;ii++)
    {
      if(auxoc==true)
	{
	  N[ii].resetForce();
	}
      else{
	aux=N[ii].getR();
	aux1=N[ii].getV();
	auxoc=N[ii].getoc();
	qaux=N[ii].getQ();
	N[ii].resetForce();
	auxx=int((nx-2)*aux[0]/l)+1;
	auxy=int((ny-2)*aux[1]/l)+1;
	if(1<=auxx<nx-1 && 1<=auxy<ny-1)
	  {
	    Faux[0]=qaux*(data[(auxx-1)*ny+auxy].value-data[(auxx+1)*ny+auxy].value)/(2*Delta)-gamma*aux1[0];
	    Faux[1]=qaux*(data[auxx*ny+(auxy-1)].value-data[auxx*ny+(auxy+1)].value)/(2*Delta)-gamma*aux1[1];
	    N[ii].addForce(Faux);
	  }
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
      auxx=int((nx-2)*aux[0]/l)+1;
      auxy=int((ny-2)*aux[1]/l)+1;
      if( qaux>0 && 1<=auxx<nx-1 && 1<=auxy<ny-1)
	{
	  if(data[(auxx+1)*ny+auxy].ocupation == true ||data[(auxx-1)*ny+auxy].ocupation == true || data[auxx*ny+(auxy+1)].ocupation == true ||  data[auxx*ny+(auxy-1)].ocupation == true)
	    {
	      Vaux[0]=0;
	      Vaux[1]=0;
	      N[ii].setV(Vaux);
	      N[ii].setoc(true);
	    }
	}
    }

  return;
}

void update_and_check_pos(Body * N, int nx, int ny, double l, int Nmax, data_t & data, double mu, double sigma, double dt)
{
  Vector3D Rold, Rnew, Vold, Vnew, Dr, move, Vaux;
  double qaux, move_x, move_y, b, p;
  int auxx, auxy, auxx1, auxy1;
  std::vector<int> cond (4,{0});
  CRandom rand(42);
  double ltot=l+2*DELTA;
  for(int ii=0;ii<Nmax;ii++)
    {
      if(N[ii].getoc()==false)
	{
	  Rold=N[ii].getR();
	  Vold=N[ii].getV();
	  qaux=N[ii].getQ();
	  Get_EF(N, nx, ny, l, Nmax, data, DELTA, 100);
	  N[ii].moveV(dt, 1);
	  N[ii].moveR(dt, 1);
	  move_x = rand.gauss(mu, sigma);
	  move_y = rand.gauss(mu, sigma);
	  move.load(move_x, move_y, 0);
	  N[ii].setR(N[ii].getR() + move);
	  Rnew=N[ii].getR();
	  Vnew=N[ii].getV();
	  Dr=Rnew-Rold;
	  if(Rnew[0]>=ltot) cond[0]=1;
	  else if(Rnew[1]>=ltot) cond[1]=1;
	  else if(Rnew[0]<=0) cond[2]=1;
	  else if(Rnew[1]<=0) cond[3]=1;
	  b =std::accumulate(cond.begin(), cond.end(), 0);
	  auxx=int((nx-2)*Rnew[0]/l);
          auxy=int((ny-2)*Rnew[1]/l);
	  if(auxx<0)auxx=(nx-2)-((-auxx)%(nx-2));
	  if(auxy<0)auxy=ny-((-auxy)%ny);
	  if(b>0)
	    {
	      if(b>1){
		auxx=(auxx)%(nx-2);
		auxy=(auxy)%(ny-2);
		Rnew[0]=auxx*l/(nx-2)+DELTA;
		Rnew[1]=auxy*l/(ny-2)+DELTA;
		N[ii].setV(Vnew);
		N[ii].setR(Rnew);
	      }
	      else if (qaux>0 && cond[2]==1)
		{
		  p=Dr[1]/Dr[0];
		  Vaux[0]=0;
		  Vaux[1]=Rold[1]+p*Rold[0];
		  Vaux=Vaux-Rold;
		  auxx=(auxx)%(nx-2);
		  Rnew[0]=auxx*l/(nx-2)+DELTA;
		  for (int kk=0; kk<500; kk++)
		    {
		      double daux=0.002;
		      auxx1=int(nx*(Rold[0]+kk*daux*Vaux[0])/ltot);
		      auxy1=int(ny*(Rold[1]+kk*daux*Vaux[1])/ltot);
		      if(auxx1 <nx && auxy1 < ny && data[(auxx1)*ny+auxy1].ocupation == true)
			{
			  Rold[0]=Rold[0]+(kk-1)*daux*Vaux[0];
			  Rold[1]=Rold[1]+(kk-1)*daux*Vaux[1];
			  Vold[0]=0;
			  Vold[1]=0;
			  N[ii].setR(Rold);
			  N[ii].setV(Vold);
			  N[ii].setoc(true);
			  break;
			}
		      else{
			N[ii].setV(Vnew);
			N[ii].setR(Rnew);
		      }
		    }
		}
	      
	      else if(cond[1]==1  || cond[3]==1)
		{
		  auxy=(auxy)%(ny-2);
		  Rnew[1]=auxy*l/(ny-2)+DELTA;
		  N[ii].setV(Vnew);
		  N[ii].setR(Rnew);
		}
	      else
		{
		  auxx=(auxx)%(nx-2);
		  Rnew[0]=auxx*l/(nx-2)+DELTA;
		  N[ii].setV(Vnew);
		  N[ii].setR(Rnew);
		}  
	    }
	}
    }
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
