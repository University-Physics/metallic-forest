#ifndef Vector3D
#include "encabezado.h"
#endif

/*
This function creates the initial conditions of the electric potential, it evaluates to zero the magnitude and false the variables ocupation and electrode in all the boxes of the grid.
 */

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
/*
As its name says , the boundary_conditions function sets the boundary in the potential creating the electrodes, for this purpose ones boxes in the grid are set to a potential depending if  they are part of the anode or cathode (-V_diff/2 and V_diff/2 respectively) and true for the bools electrode and ocupation
 */


void boundary_conditions(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
{
  //We have to note that three sides are part of the anode and only one of the cathode

    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {       //     #----
      data[ix*ny + iy].value = -V_diff/2;  //     #----
      data[ix*ny + iy].electrode = true;   //     #----
      data[ix*ny + iy].ocupation = true;   //     #----
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      data[ix*ny + iy].value = V_diff/2;   //     ----#
      data[ix*ny + iy].electrode = true;   //     ----#
      data[ix*ny + iy].ocupation = false;  //     ----#

    }
    // first row

    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      data[ix*ny + iy].value = 0;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    -----
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {       //    -----
      data[ix*ny + iy].value =0;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    #####
    }
}

void boundary_conditions_wn(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff, double amp, int seed)
{
  //We have to note that three sides are part of the anode and only one of the cathode
    // Set random number generator for white noise
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {       //     #----
      data[ix*ny + iy].value =-V_diff/2+amp*dis(gen);  //     #----
      data[ix*ny + iy].electrode = true;   //     #----
      data[ix*ny + iy].ocupation = true;   //     #----
	}
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      data[ix*ny + iy].value =V_diff/2+amp*dis(gen);   //     ----#
      data[ix*ny + iy].electrode = true;   //     ----#
      data[ix*ny + iy].ocupation = false;  //     ----#

    }
    // first row

    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      data[ix*ny + iy].value =V_diff/2+amp*dis(gen);   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    -----
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {       //    -----
      data[ix*ny + iy].value =V_diff/2+amp*dis(gen);   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    #####
    }
}

void boundary_conditions_wnR(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff, double amp, int seed)
{
  //We have to note that three sides are part of the anode and only one of the cathode
    // Set random number generator for white noise
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(0, 1.0);
    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy)
      {//     #----
	for(int xx=0;xx<=int(amp*dis(gen));xx++)
	 {
      data[xx*ny + iy].value =-V_diff/2;  //     #----
      data[xx*ny + iy].electrode = true;   //     #----
      data[xx*ny + iy].ocupation = true;   //     #----
         }
      }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      for(int xx=0;xx<=int(amp*dis(gen));xx++)
	 {
      data[xx*ny + iy].value =V_diff/2;  //     #----
      data[xx*ny + iy].electrode = true;   //     #----
      data[xx*ny + iy].ocupation = false;   //     #----
         }

    }
    // first row
    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      for(int xx=0;xx<=int(amp*dis(gen));xx++)
	 {
      data[xx*ny + iy].value =V_diff/2;  //     #----
      data[xx*ny + iy].electrode = true;   //     #----
      data[xx*ny + iy].ocupation = false;   //     #----
         }
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {       //    -----
      for(int xx=0;xx<=int(amp*dis(gen));xx++)
	 {
      data[xx*ny + iy].value =V_diff/2;  //     #----
      data[xx*ny + iy].electrode = true;   //     #----
      data[xx*ny + iy].ocupation = false;   //     #----
         }
    }
}



void boundary_conditions_pn(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff, double amp, int seed)
{
  //We have to note that three sides are part of the anode and only one of the cathode
    // Set the pink noise text files cause I dont know how to do it in c :(
    std::string filename="pn"+std::to_string(seed)+".txt";
    std::ifstream myfile;
    myfile.open(filename);
    double noise=0;
    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {       //     #----
      myfile >> noise;
      data[ix*ny + iy].value = -V_diff/2+amp*noise;  //     #----
      data[ix*ny + iy].electrode = true;   //     #----
      data[ix*ny + iy].ocupation = true;   //     #----
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      myfile >> noise;
      data[ix*ny + iy].value = V_diff/2+amp*noise;   //     ----#
      data[ix*ny + iy].electrode = true;   //     ----#
      data[ix*ny + iy].ocupation = false;  //     ----#

    }
    // first row

    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      myfile >> noise;
      data[ix*ny + iy].value = 0+amp*noise;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    -----
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {       //    -----
      data[ix*ny + iy].value =0+amp*noise;   //    -----
      myfile >> noise;
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    #####
    }
}

void deposited_particles(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
{
  Vector3D aux;
   for(int ii=0;ii<Nmax;ii++) // This for is concerning for all the particles, it evaluates if a particle have collided with the anode and now it's part of it.
    {
      aux=N[ii].getR();
      int auxx=int(nx*aux[0]/l);
      int auxy=int(ny*aux[1]/l);
      if(auxx>=1 && auxx<nx-1 && auxy>=1 && auxy<ny-1)
	{
      if(N[ii].getoc()==true && data[auxx*ny + auxy].electrode == false) // Collision condition
	{
	  // Set the electrode conditions in the box where the particle remains
	  double prom=0;
	  double count=0;
	  /*
	  for(int xx=0; xx<3;xx++)
	    {
	      for(int yy=0;yy<3;yy++)
		{
		  if(data[(auxx-1+xx)*ny+auxy-1+yy].electrode==true)
		    count+=1;
		  prom+=data[(auxx-1+xx)*ny+auxy-1+yy].value;
		}
	    }
	  */
	  //data[auxx*ny+auxy].value=prom/count;
	  data[auxx*ny+auxy].value=-V_diff/2;
	  data[auxx*ny+auxy].ocupation=true;
	  data[auxx*ny+auxy].electrode=true;
	}
	}
    }
}
void boundary_conditions1(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
{
  //We have to note that three sides are part of the anode and only one of the cathode

    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {       //     #----
      data[ix*ny + iy].value = -V_diff/2;  //     #----
      data[ix*ny + iy].electrode = true;   //     #----
      data[ix*ny + iy].ocupation = true;   //     #----
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      data[ix*ny + iy].value = V_diff/2;   //     ----#
      data[ix*ny + iy].electrode = true;   //     ----#
      data[ix*ny + iy].ocupation = false;  //     ----#

    }
    // first row

    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      data[ix*ny + iy].value = V_diff/2;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    -----
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {       //    -----
      data[ix*ny + iy].value =V_diff/2;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    #####
    }
}

void boundary_conditions2(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
{
  //We have to note that three sides are part of the anode and only one of the cathode

    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {       //     #----
      data[ix*ny + iy].value = -V_diff/2;  //     #----
      data[ix*ny + iy].electrode = true;   //     #----
      data[ix*ny + iy].ocupation = true;   //     #----
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      data[ix*ny + iy].value = -V_diff/2;   //     ----#
      data[ix*ny + iy].electrode = true;   //     ----#
      data[ix*ny + iy].ocupation = true;  //     ----#

    }
    // first row

    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      data[ix*ny + iy].value = -V_diff/2;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = true;  //    -----
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {       //    -----
      data[ix*ny + iy].value =-V_diff/2;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = true;  //    #####
    }
    data[(nx-1)/2*ny+ny/2].value=V_diff/2;
    data[(nx-1)/2*ny+ny/2].electrode=true;
    data[(nx-1)/2*ny+ny/2].ocupation=false;
}



void boundary_conditions4(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
{
  //We have to note that three sides are part of the anode and only one of the cathode

    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {       //     #----
      data[ix*ny + iy].value = -V_diff/2;  //     #----
      data[ix*ny + iy].electrode = true;   //     #----
      data[ix*ny + iy].ocupation = true;   //     #----
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      data[ix*ny + iy].value = V_diff/2;   //     ----#
      data[ix*ny + iy].electrode = true;   //     ----#
      data[ix*ny + iy].ocupation = false;  //     ----#

    }
    // first row

    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      data[ix*ny + iy].value = -V_diff/2*(1-(2.00*ix)/((nx-1)));   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    -----
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {                      //    -----
      data[ix*ny + iy].value =-V_diff/2*(1-(2.00*ix)/((nx-1)));//    -----
      data[ix*ny + iy].electrode = true;                  //    -----
      data[ix*ny + iy].ocupation = false;                 //    #####
    }
}

void boundary_conditions3(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
{
  //We have to note that three sides are part of the anode and only one of the cathode

    int ix, iy;
    Vector3D aux, aux1;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {       //     #----
      data[ix*ny + iy].value = V_diff/2;  //     #----
      data[ix*ny + iy].electrode = true;   //     #----
      data[ix*ny + iy].ocupation = false;   //     #----
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {       //     ----#
      data[ix*ny + iy].value = V_diff/2;   //     ----#
      data[ix*ny + iy].electrode = true;   //     ----#
      data[ix*ny + iy].ocupation = false;  //     ----#

    }
    // first row

    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {       //    #####
      data[ix*ny + iy].value = V_diff/2;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    -----
    }
    // last row
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {       //    -----
      data[ix*ny + iy].value =V_diff/2;   //    -----
      data[ix*ny + iy].electrode = true;   //    -----
      data[ix*ny + iy].ocupation = false;  //    #####
    }
    data[(nx-1)/2*ny+ny/2].value=-V_diff/2;
    data[(nx-1)/2*ny+ny/2].electrode=true;
    data[(nx-1)/2*ny+ny/2].ocupation=true;
}
double Probability_distribution(Body * molecule, int N)
{
  double count=0;
  for(int i=0;i<N ;i++)
    {
      if(molecule[i].getoc()==true)
	{
	  count+=1;
	}
    }
  return count;
}

bool evolve(data_t & data, int nx, int ny, int nsteps, int ns_est)
{
  //start_gnuplot();
  bool seguir=true;
    for(int istep = 0; istep < nsteps; ++istep) {
      //Perform usual relaxation on the big lattice
     seguir=relaxation_step(data, nx, ny);
     if(seguir==false)
       {
	 return false;
       }

      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
    for(int istep = 0; istep < nsteps; ++istep) {
      //does relaxation on the lattice with less sites
      relaxation_step_pivot(data, nx, ny);
      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
    for(int istep = 0; istep < 3*ns_est; ++istep){
      //Does relaxation over the big lattice considering the points in the
      //smaller lattice as fixed, thus the value of nearst neighbors gets
      //corrected by the last relaxation rutine
      stabilization_step(data, nx, ny);
      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
    for(int istep = 0; istep < 3*nsteps; ++istep) {
      //Perform relaxation over the big lattice to ensure convergence
      seguir=relaxation_step(data, nx, ny);
      if(seguir==false)
	{
	  return false;
	}
      //print_screen(data, nx, ny);
      //print_gnuplot(data, nx, ny);
    }
    return seguir;
}

void stabilization_step(data_t & data, int nx, int ny)
{
  for(int ix = 1; ix < nx-1; ++ix) {
        for(int iy = 1; iy < ny-1; ++iy) {
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
      if(data[(ix*10)*ny + (iy*10)].ocupation==false)
	{
	  data[(ix*10)*ny + (iy*10)].value = (data[((ix+1)*10)*ny + (iy*10)].value + data[(ix-1)*10*ny + (iy*10)].value + data[(ix*10)*ny + (iy+1)*10].value + data[(ix*10)*ny + 10*(iy-1)].value)/4.0;
	  data[(ix*10)*ny+(iy*10)].pivot=true;
	}

    }
  }

}


bool relaxation_step(data_t & data, int nx, int ny)
{
    // recorrer toda la matriz y aplicar el algoritmo,
    // teniendo cuidado con no modificar las condiciones de
    // frontera
  bool seguir=false;
  data_t aux=data;
    for(int ix = 1; ix < nx-1; ++ix)
      {
        for(int iy = 1; iy < ny-1; ++iy)
	  {
	  if(data[ix*ny + iy].ocupation==false)
	    {
	      double next_value=(aux[(ix+1)*ny + iy].value + aux[(ix-1)*ny + iy].value + aux[ix*ny + iy+1].value + aux[ix*ny + iy-1].value)/4.0;
	      if(std::fabs(data[ix*ny + iy].value-next_value)>std::fabs(next_value)*0.001)//(v(x,y,t)-v(x,y,t+dt))<0.01*v(x,y,t)
		{
		  seguir=true;
		  data[ix*ny + iy].value = next_value;
		}
	    }
         }
    }
    return seguir;
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
  std::cout<<"set hidden3d\n"<<std::endl;
  std::cout<<"set isosamples 50"<<std::endl;
  std::cout<<"set pm3d at s\n"<<std::endl;
  std::cout<<"set contour\n"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set cntrparam levels 10\n"<<std::endl;
  std::cout << "set term postscript eps enhanced color\n";
  std::cout << "set output 'color.eps'\n";
    std::cout << "splot '-' w l\n";
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
Vector3D E_field(int nx,int ny, double x, double y,double DELTA, data_t & data)
{
      int auxx=int(x/DELTA);
      int auxy=int(y/DELTA);
      double Electric_Field;
      Vector3D E_00;
      Vector3D E_01;
      Vector3D E_10;
      Vector3D E_11;
      Vector3D F_xy1;
      Vector3D F_xy2;
      E_00.load((data[(auxx-1)*ny+auxy].value-data[(auxx+1)*ny+auxy].value)/(2*DELTA),(data[auxx*ny+(auxy-1)].value-data[auxx*ny+(auxy+1)].value)/(2*DELTA),0);
      E_01.load((data[(auxx-1)*ny+auxy+1].value-data[(auxx+1)*ny+auxy+1].value)/(2*DELTA),(data[auxx*ny+(auxy)].value-data[auxx*ny+(auxy+2)].value)/(2*DELTA),0);
      E_10.load((data[(auxx)*ny+auxy].value-data[(auxx+2)*ny+auxy].value)/(2*DELTA),(data[(auxx+1)*ny+(auxy-1)].value-data[(auxx+1)*ny+(auxy+1)].value)/(2*DELTA),0);
      E_11.load((data[(auxx)*ny+auxy+1].value-data[(auxx+2)*ny+auxy+1].value)/(2*DELTA),(data[(auxx+1)*ny+(auxy)].value-data[(auxx+1)*ny+(auxy+2)].value)/(2*DELTA),0);
      F_xy1=((auxx+1)*DELTA-x)/(DELTA)*E_00+(x-auxx*DELTA)/(DELTA)*E_10;
      F_xy2=((auxx+1)*DELTA-x)/(DELTA)*E_01+(x-auxx*DELTA)/(DELTA)*E_11;
return ((auxy+1)*DELTA-y)/(DELTA)*F_xy1+(y-auxy*DELTA)/(DELTA)*F_xy2;
}

void Get_EF(Body * N, int nx, int ny, int Nmax, data_t & data, double Delta, double gamma, bool interacciones)
{
  Vector3D aux, aux1, aux2, Faux, dr;
  double q1, q2, d_aux, r_aux;
  int auxx, auxy;
  if(interacciones==true)
    {
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
	if(1<=auxx && auxx<nx-1 && 1<=auxy && auxy<ny-1)
	  {

	    for(int jj = 0; jj<ii; jj++)
	      {

		if(N[jj].getoc()==false)
		  {
		  aux2=N[jj].getR();
		  dr=aux-aux2;
		  q2=N[jj].getQ();
		  d_aux=norm(dr);
		  if(d_aux<2*N[ii].getrad())d_aux=2*N[ii].getrad();
		  Faux=q1*q2*dr/(d_aux*d_aux*d_aux);
		  N[ii].addForce(Faux);
		  N[jj].addForce((-1)*Faux);
		  }
	      }
	    //Faux[0]=q1*(data[(auxx-1)*ny+auxy].value-data[(auxx+1)*ny+auxy].value)/(2*Delta)-gamma*aux1[0];
	    //Faux[1]=q1*(data[auxx*ny+(auxy-1)].value-data[auxx*ny+(auxy+1)].value)/(2*Delta)-gamma*aux1[1];
	     Faux=q1*E_field(nx,ny,aux[0],aux[1],Delta,data)-gamma*aux1;
	    N[ii].addForce(Faux);
	  }
      }
    }
    }
  else
    {

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
	if(1<=auxx && auxx<nx-1 && 1<=auxy && auxy<ny-1)
	  {

	    Faux[0]=q1*(data[(auxx-1)*ny+auxy].value-data[(auxx+1)*ny+auxy].value)/(2*Delta)-gamma*aux1[0];
	    Faux[1]=q1*(data[auxx*ny+(auxy-1)].value-data[auxx*ny+(auxy+1)].value)/(2*Delta)-gamma*aux1[1];
	    // E_field(int nx,int ny, double x, double y,double DELTA, data_t & data)
	    // Faux=q1*E_field(nx,ny,aux[0],aux[1],Delta,data)-gamma*aux1;
	    N[ii].addForce(Faux);
	  }
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
      if(N[ii].getoc()==false)
       {
      if( qaux>0 && 1<=auxx && auxx<nx-1 && 1<=auxy && auxy<ny-1)
	{
	  if(data[(auxx+1)*ny+auxy].ocupation == true ||data[(auxx-1)*ny+auxy].ocupation == true || data[auxx*ny+(auxy+1)].ocupation == true ||  data[auxx*ny+(auxy-1)].ocupation == true || data[(auxx+1)*ny+auxy-1].ocupation == true || data[(auxx+1)*ny+auxy+1].ocupation == true ||data[(auxx-1)*ny+auxy+1].ocupation == true || data[(auxx-1)*ny+auxy-1].ocupation == true)
	    {
	      Vector3D aux2;
	      double radio=N[ii].getrad();
	      for(int jj=0;jj<Nmax;jj++)
	      {
		if(N[jj].getoc()==true)
		  {
		  aux2=N[jj].getR();
		  double norma=norm(aux2-aux);
		  if(norma<2*radio && not std::isnan(radio/norma))
		    {
		      aux+=aux*(radio/norma)-aux2*(radio/norma);
		      }
		  }
	      }
	      cond=true;
	      Vaux[0]=0;
	      Vaux[1]=0;
	      N[ii].setV(Vaux);
	      N[ii].setR(aux);
	      N[ii].setoc(true);
	    }
	}
       }
       }

  return cond;
}
void update_and_check_pos2(Body * N, int nx, int ny, int Nmax, data_t & data, double mu, double sigma, double dt, double coefx, double coefv, int seed, bool interacciones)
{
  Vector3D Rnew, Vnew, move, DR, Rold, Vold, AUXX;
  double move_x, move_y;
  CRandom rand(seed);
  Get_EF(N, nx, ny, Nmax, data, DELTA, 1,interacciones);
  for(int ii=0;ii<Nmax;ii++)
    {
      if(N[ii].getoc()==false)
	{
	  Rold=N[ii].getR();
	  Vold=N[ii].getV();
	  N[ii].moveV(dt, coefx);
	  N[ii].moveR(dt, coefv);
	  move_x = rand.gauss(mu, sigma);
	  move_y = rand.gauss(mu, sigma);
	  //std::cout<<"MOVES"<<" x"<<move_x<<" "<<move_y<<std::endl;
	  //std::cout<<"FORCE"<<N[ii].getF()[0]<<" "<<N[ii].getF()[1]<<std::endl;
	  move.load(move_x, move_y, 0);
	  AUXX=N[ii].getR();
	  AUXX+=coefx*move;
	  //std::cout<<AUXX[0]<<" "<<AUXX[1]<<std::endl;
	  N[ii].setR(N[ii].getR() + (coefx)* move);
	  Rnew=N[ii].getR();
	  Vnew=N[ii].getV();
	  DR=Rnew-Rold;
	  //std::cout<<Rold[0]<<" "<<Rold[1]<<" "<<Rnew[0]<<" "<<Rnew[1]<<" "<<coefx<<" "<<coefv<<std::endl;
	  if(Rnew[0]>1 || Rnew[0]<0)
	    {


	      if(Rnew[0]<0 && N[ii].getQ()>0)
		{
		  Rnew[0]=0;
		  Rnew[0]=Rold[0];
		  Rnew[1]=Rold[1];
		  int auxx=1;
		  int auxy=1;
		  double A=Rnew[0];
		  double B=Rnew[1];

		  for(int t=0;t<10;t++)
		    {
		      auxx=int(A/DELTA);
		      auxy=int(B/DELTA);
		      if(auxx==0 && auxy>=1 && auxy<ny-1 && data[0].ocupation)
			{

			  data[auxx*ny+auxy].ocupation=true;
			  Vnew[0]=0;
			  Vnew[1]=0;
			  N[ii].setoc(true);
			  break;
			}
		      if(auxx>=1 && auxy>=1 && auxy<ny-1 && auxx<nx-1)
			{
			  Rnew[0]=A;
			  Rnew[1]=B;
			  if(data[(auxx+1)*ny+auxy].ocupation == true || data[(auxx-1)*ny+auxy].ocupation == true || data[auxx*ny+(auxy+1)].ocupation == true ||  data[auxx*ny+(auxy-1)].ocupation == true || data[(auxx+1)*ny+auxy-1].ocupation == true || data[(auxx+1)*ny+auxy+1].ocupation == true ||data[(auxx-1)*ny+auxy+1].ocupation == true || data[(auxx-1)*ny+auxy-1].ocupation == true || data[(auxx)*ny+auxy].ocupation == true)
			{

			  Vnew[0]=0;
			  Vnew[1]=0;
			  break;

			}
			}
		      else
			{
			  Rnew[0]=Rold[0];
			  Vnew[0]*=-1;
			  break;
			}
		      A+=DR[0]/10;
		      B+=DR[1]/10;
		    }
		}
		  else
		    {
		      Rnew[0]=Rold[0];
		      Vnew[0]*=-1;
		    }
	    }

	  if(Rnew[1]>=1 || Rnew[1]<=0)
	    {

	      Rnew[1]=Rold[1];
	      Vnew[1]*=-1;
	    }
	  N[ii].setV(Vnew);
	  N[ii].setR(Rnew);
	}
    }
}

void evolve_system(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, double coefx, double coefv, int seed, double V_diff,bool interacciones)
{
  bool cond;
  //update positions and boundaries
  update_and_check_pos2(N, nx, ny, Nmax, data, mu, sigma, dt, coefx, coefv, seed,interacciones);
  cond=Update_boundary(N, nx, ny, Nmax, data);
  //Obtain potential using fast algorithm
  if(cond==true)
    {
      deposited_particles(data, nx, ny, N ,l, Nmax, V_diff);
      while(relaxation_step(data,nx,ny)==true)
	{
	}
    }
}

void PEFRL(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, int seed, double V_diff, bool interacciones)
{
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Zi, 0.0, seed, V_diff,interacciones);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Xi, coef1, seed, V_diff,interacciones);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, coef2, Lambda, seed, V_diff,interacciones);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Xi, Lambda, seed, V_diff,interacciones);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Zi, coef1, seed, V_diff,interacciones);
}

void print_fractal (int nx, int ny, data_t & data, std::string filename)
{
  std::ofstream myfile(filename);
  for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
	  myfile<< ix << "  "<<iy<<"   "<< data[ix*ny + iy].ocupation << "\n";
        }
        myfile << "\n";
    }
    myfile << "\n";
    myfile.close();
}
void print(data_q probability,std::string filename)
{
  int a=probability.size();
  std::ofstream myfile(filename);
  for(int ii=0;ii<a;ii++)
    {
      myfile<<ii<<"\t"<<probability[ii]<<"\n";
    }
  myfile<<"\n";
  myfile.close();

}
void print_potential_size (int nx, int ny, data_t & data, std::string filename, int t)
{
  std::ofstream myfile;
  int N=0;
  myfile.open (filename, std::fstream::app);
  for(int ix = 1; ix < nx-1; ++ix) {
        for(int iy = 1; iy < ny-1; ++iy) {
	  if(data[ix*ny+iy].ocupation==true) N+=1;
        }
    }
    myfile<< t << "\t"<<N<<"\n";
    myfile.close();
}

void print_potential (int nx, int ny, data_t & data, std::string filename)
{
  std::ofstream myfile;
   myfile.open (filename);
  for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
	  myfile<< ix << "\t"<<iy<<"\t"<< data[ix*ny + iy].value << "\n";
        }
    }
    myfile << "\n";
    myfile.close();
}

void print_potential_gnuplot (int nx, int ny,data_t & data)
{
  
  for(int ii=0;ii<nx;ii++)
    {
       for(int jj=0;jj<ny;jj++)
	 {
	   if(data[ii*nx+jj].electrode==true)
	     {
	       //std::cout<<", "<<ii/(201.0)<<"+"<<(ii+1)/(201.0)<<"*(t-1)/6 ,"<<(jj+1)/(201.0)<<"+t*0.0 with filledcurves above y="<<jj/(201.0)<<" lt \"gray \" ";
	       std::cout<<", "<<(ii)*(1.0)*DELTA+0.5*DELTA<<"+"<<DELTA/2<<"*cos(t) ,"<<(jj)*(1.0)*DELTA+0.5*DELTA<<"+"<<DELTA/2<<"*sin(t) with filledcurves closed  lt \"gray \" ";
	       
	     }
	 }
    }
 
   
}



void evolve_opt(data_t & data, int nx, int ny, double l, int Nmax, Body * N, double V_diff)
{
  //doesn't reset values in order to obtain faster convergence and increase performance
  boundary_conditions(data, nx, ny, N ,l, Nmax, V_diff);
  // evolve
  evolve(data, nx, ny, NSTEPS/10, NSTEPS/10);
}
bool check_fractal(Body * molecule,int nx, int Nmax, int I)
{
  bool prueba=false;
    for(int i=0;i<Nmax;i++)
    {
      double deposit=Probability_distribution(molecule,N);
      double percent=Nmax/I;
      if(deposit > percent*0.8)
	{
	  std::cout<<deposit<<" "<<percent<<std::endl;
	  prueba=true;
	}
    }
  return prueba;
}
