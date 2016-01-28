#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* f1, cmplx* const f0, const int Nx, 
	  const double k, const double dx, const double xmin, const double dt);

//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40.0;
	const double xmax = 40.0;
	const double Tend = 10.0*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt = dx/100.0;
	double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
	const double omega = 0.2;
	const double k = omega*omega;
	const double alpha = pow(k, 1./4.);

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);


	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		  
	step(psi1, psi0, Nx, k, dx, xmin, dt);	  
		 
	h = psi0;
	
	psi0 = psi1;
	
	psi1 = h;
		  

         t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
	cout << "t = " << t << endl;
	
	delete[] psi0;
	delete[] psi1;
	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
	
    out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx;i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}

void step(cmplx* f1, cmplx* f0, const int Nx, const double k, const double dx, const double xmin, const double dt)
{	
	double x;
	cmplx* d = new cmplx[Nx]; 
	cmplx* u = new cmplx[Nx]; 
	cmplx* l = new cmplx[Nx]; 
	cmplx* r = new cmplx[Nx];
	double* V = new double[Nx];
	
	for(int i=0;i<Nx;i++){
	  x = xmin + i * dx;
	  V[i] = 0.5*k*x*x;
	  d[i] = cmplx(1.0, dt/(2*dx*dx)+dt*V[i]/2.);
	  u[i] = cmplx(0.0, -dt/(4.*dx*dx));
	  l[i] = cmplx(0.0, -dt/(4.*dx*dx));
	}
	
	//right side of equation before modified by Forward substitution
	//for the following caption it does not matter whether l or u is used
	r[0] = f0[0] - l[0]*(f0[1]-cmplx(2.,0.0)*f0[0]) - cmplx(0.0, dt/2.0)*V[0]*f0[0];
	
	for(int i=1;i<Nx-1;i++){
	  x = xmin + i * dx;
	  r[i] = f0[i] - l[i]*(f0[i+1]-cmplx(2.,0.0)*f0[i]+f0[i-1]) - cmplx(0.0, dt/2.0)*V[i]*f0[i];
	}
	
	r[Nx-1] = f0[Nx-1] - l[Nx-1]*(-cmplx(2.,0.0)*f0[Nx-1]+f0[Nx-2]) - cmplx(0.0, dt/2.0)*V[Nx-1]*f0[Nx-1];
	
	//Forward substitution
	
	for(int i=1;i<Nx;i++){
	  d[i] = d[i] - l[i]*u[i-1]/d[i-1];
	  r[i] = r[i] - l[i]*r[i-1]/d[i-1];
	}
	
	//Backward substitution
	f1[Nx-1] = r[Nx-1]/d[Nx-1];
	
	for(int i=Nx-2;i>=0;i--){
	  f1[i] = (r[i] - u[i]*f1[i+1])/d[i];
	}
  
  delete[] V;
  delete[] d;
  delete[] u;
  delete[] l;
  delete[] r;
}
