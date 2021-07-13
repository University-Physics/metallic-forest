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

void deposited_particles(data_t & data, int nx, int ny, Body * N, double l, int Nmax, double V_diff)
{
  Vector3D aux;
   for(int ii=0;ii<Nmax;ii++) // This for is concerning for all the particles, it evaluates if a particle have collided with the anode and now it's part of it.
    {
      aux=N[ii].getR();
      if(N[ii].getoc()==true && data[int(nx*aux[0]/l)*ny + int(ny*aux[1]/l)].electrode == false) // Collision condition
	{
	  // Set the electrode conditions in the box where the particle remains
	  
	  data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].value=-V_diff/2; 
	  data[int(nx*aux[0]/l)*ny+int(ny*aux[1]/l)].ocupation=true;  
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

void Get_EF(Body * N, int nx, int ny, int Nmax, data_t & data, double Delta, double gamma)
{
  Vector3D aux, aux1, aux2, Faux, dr;
  double q1, q2, d_aux, r_aux;
  int auxx, auxy;
  
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
		
		if(N[jj].getoc()==false){
		  aux2=N[jj].getR();
		  dr=aux2-aux1;
		  q2=N[jj].getQ();
		  d_aux=norm(dr);
		  if(d_aux<2*N[ii].getrad())d_aux=2*N[ii].getrad();
		  Faux=q1*q2*dr/(d_aux*d_aux*d_aux);
		  
		  if(norm(dr)<N[ii].getrad())
		    {
		      r_aux=std::sqrt(1/(1/N[ii].getrad()+1/N[jj].getrad()));
		      Faux-=4/3*r_aux*std::sqrt(norm(dr))*dr;
		    }
		  
		  N[ii].addForce(Faux);
		  N[jj].addForce((-1)*Faux);
		}
	      }
	    
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
      if(N[ii].getoc()==false)
       {
      if( qaux>0 && 1<=auxx && auxx<nx-1 && 1<=auxy && auxy<ny-1)
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
       }

  return cond;
}
void update_and_check_pos2(Body * N, int nx, int ny, int Nmax, data_t & data, double mu, double sigma, double dt, double coefx, double coefv, int seed)
{
  Vector3D Rnew, Vnew, move, DR, Rold, Vold;
  double move_x, move_y;
  CRandom rand(seed);
  Get_EF(N, nx, ny, Nmax, data, DELTA, 1);
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
	  move.load(move_x, move_y, 0);
	  N[ii].setR(N[ii].getR() + (coefx)* move);
	  Rnew=N[ii].getR();
	  Vnew=N[ii].getV();
	  DR=Rnew-Rold;
	  if(Rnew[0]>1 || Rnew[0]<0) 
	    {//Rnew=Rnew-DR
	      Rnew[0]=Rnew[0]-DR[0];
	      Vnew[0]*=-1;
	    }
	  if(Rnew[1]>1 || Rnew[1]<0)
	    {
	      
	      Rnew[1]=Rnew[1]-DR[1];
	      Vnew[1]*=-1;
	    }
	  N[ii].setV(Vnew);
	  N[ii].setR(Rnew);
	}
    }
}

void evolve_system(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, double coefx, double coefv, int seed, double V_diff)
{
  bool cond;
  //update positions and boundaries
  update_and_check_pos2(N, nx, ny, Nmax, data, mu, sigma, dt, coefx, coefv, seed);
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

void PEFRL(Body * N, data_t & data, int nx, int ny, double l, int Nmax, double mu, double sigma, double dt, int seed, double V_diff)
{
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Zi, 0.0, seed, V_diff);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Xi, coef1, seed, V_diff);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, coef2, Lambda, seed, V_diff);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Xi, Lambda, seed, V_diff);
  evolve_system(N, data, nx, ny, l, Nmax, mu, sigma, dt, Zi, coef1, seed, V_diff);
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

void evolve_opt(data_t & data, int nx, int ny, double l, int Nmax, Body * N, double V_diff)
{
  //doesn't reset values in order to obtain faster convergence and increase performance
  boundary_conditions(data, nx, ny, N ,l, Nmax, V_diff);
  // evolve
  evolve(data, nx, ny, NSTEPS/10, NSTEPS/10);
}
bool check_fractal(Body * molecule,int nx, int Nmax)
{
  bool prueba=false;
    for(int i=0;i<Nmax;i++)
    {
      int percent=int((nx-1)*(0.9));
      double x=molecule[i].getR()[0]/DELTA;
      bool verdad= molecule[i].getoc();
      if(int(x)>=percent && verdad==true)
	{
	  prueba=true;
	}
    }
  return prueba;
}
