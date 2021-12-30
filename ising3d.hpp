#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#define Nbrs 6

class ising3d
{
	int N,ns,*s,**nn,teq,mts;
	double J,kB,mo,T,B,rN;
	double mM,mU,mC,mU2,rmts,M2;
	
	public:
        int nx,ny,nz;
		ising3d();
		~ising3d();
		void init(int,int,int,double,double,double);
		void init(std::string);
		void ScanMicroStates(double,double);
		void ScanMicroStates(double,double,int);
		void SMS(double,double,bool);
		double SumNN(int);
		double U();
		double u();
		double M();
		double avgM();
		double SumSpinProducts();
		double SumSpins();
		double avgU();
		double avgu();
		double avgC();
		double avgc();
		double avgX();
		void rescale();
		void initM(double&);
		void update1M(double&,int);
		void update2M(double);
		void initU(double&);
		void update1U(double&,double);
		void update2U(double);
        int getN();
		
};
#include "ising3d.cpp"
#undef Nbrs

