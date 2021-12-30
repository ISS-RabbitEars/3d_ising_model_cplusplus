#include "ran3.h"

ising3d::ising3d()
{
	nx=ny=nz=N=ns=teq=mts=0;
	J=kB=mo=T=B=rN=rmts=mM=mU=mC=mU2=M2=0;	
}

ising3d::~ising3d()
{
	if(N!=0)
	{
		delete[] s;
		s=NULL;
		for(int i=0;i<N;i++){delete[] nn[i];}delete[] nn; nn=NULL;
	}
}

void ising3d::init(std::string fn)
{
	int np=6;
	double p[np];
	std::ifstream fin;
	fin.open(fn.c_str());
	for(int i=0;i<np;i++)fin>>p[i];
	fin.close();
	init(int(p[0]),int(p[1]),int(p[2]),p[3],p[4],p[5]);
}
	

void ising3d::init(int x,int y,int z,double j,double kb,double dm)
{
    nx=x;ny=y;nz=z;J=j;kB=kb;mo=dm;
	N=nx*ny*nz;
	rN=1./double(N);
    ns=100*N;
	//ns=int(pow(double(N),3));
    teq=int(0.666*double(ns)); mts=ns-teq; rmts=1./double(mts);//pow(2.,n);
	s=new int[N];
	nn=new int*[N];
	int seed=-842135656;
	for(int i=0;i<N;i++){(ran3(&seed)>=0.5)?s[i]=-1:s[i]=1; nn[i]=new int[Nbrs];}
    
    //corners
    int nxy=nx*ny, nyz=ny*nz;
    int a=0, b=nx-1, c=nx*(ny-1), d=nxy-1, e=nxy*(nz-1), f=nx*(ny*(nz-1)+1)-1, g=nx*(nyz-1), h=N-1;
    // corner a
    nn[a][0]=c;
    nn[a][1]=b;
    nn[a][2]=a+1;
    nn[a][3]=b+1;
    nn[a][4]=e;
    nn[a][5]=nxy;
    // corner b
    nn[b][0]=d;
    nn[b][1]=b-1;
    nn[b][2]=a;
    nn[b][3]=b+nx;
    nn[b][4]=f;
    nn[b][5]=d+nx;
    // corner c
    nn[c][0]=c-nx;
    nn[c][1]=d;
    nn[c][2]=c+1;
    nn[c][3]=a;
    nn[c][4]=g;
    nn[c][5]=c+nxy;
    // corner d
    nn[d][0]=d-nx;
    nn[d][1]=d-1;
    nn[d][2]=c;
    nn[d][3]=b;
    nn[d][4]=h;
    nn[d][5]=d+nxy;
    // corner e
    nn[e][0]=g;
    nn[e][1]=f;
    nn[e][2]=e+1;
    nn[e][3]=f+1;
    nn[e][4]=e-nxy;
    nn[e][5]=a;
    // corner f
    nn[f][0]=h;
    nn[f][1]=f-1;
    nn[f][2]=e;
    nn[f][3]=f+nx;
    nn[f][4]=f-nxy;
    nn[f][5]=b;
    // corner g
    nn[g][0]=g-nx;
    nn[g][1]=h;
    nn[g][2]=g+1;
    nn[g][3]=e;
    nn[g][4]=g-nxy;
    nn[g][5]=c;
    // corner h
    nn[h][0]=h-nx;
    nn[h][1]=h-1;
    nn[h][2]=g;
    nn[h][3]=f;
    nn[h][4]=h-nxy;
    nn[h][5]=d;
    
    //edges
    for(int i=1;i<b;i++)
    {
        //edge ab
        nn[i][0]=c+i;
        nn[i][1]=i-1;
        nn[i][2]=i+1;
        nn[i][3]=i+nx;
        nn[i][4]=e+i;
        nn[i][5]=nxy+i;
        //edge cd
        int ii=c+i;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii-1;
        nn[ii][2]=ii+1;
        nn[ii][3]=i;
        nn[ii][4]=g+i;
        nn[ii][5]=ii+nxy;
        //edge ef
        ii=e+i;
        nn[ii][0]=g+i;
        nn[ii][1]=ii-1;
        nn[ii][2]=ii+1;
        nn[ii][3]=ii+nx;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=i;
        //edge gh
        ii=g+i;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii-1;
        nn[ii][2]=ii+1;
        nn[ii][3]=e+i;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=c+i;
    }
    for(int i=1;i<ny-1;i++)
    {
        //edge ac
        int ii=i*nx;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii+b;
        nn[ii][2]=ii+1;
        nn[ii][3]=ii+nx;
        nn[ii][4]=ii+e;
        nn[ii][5]=ii+nxy;
        //edge bd
        ii+=b;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii-1;
        nn[ii][2]=nn[ii][0]+1;
        nn[ii][3]=ii+nx;
        nn[ii][4]=f+i*nx;
        nn[ii][5]=ii+nxy;
        //edge eg
        ii=e+i*nx;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii+b;
        nn[ii][2]=ii+1;
        nn[ii][3]=ii+nx;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=a+i*nx;
        //edge fh
        ii+=b;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii-1;
        nn[ii][2]=ii-b;
        nn[ii][3]=ii+nx;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=b+i*nx;
    }
    for(int i=1;i<nz-1;i++)
    {
        int in=i*nxy;
        //edge ae
        int ii=a+in;
        nn[ii][0]=ii+c;
        nn[ii][1]=ii+b;
        nn[ii][2]=ii+1;
        nn[ii][3]=ii+nx;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=ii+nxy;
        //edge bf
        ii=b+in;
        nn[ii][0]=d+in;
        nn[ii][1]=ii-1;
        nn[ii][2]=ii-b;
        nn[ii][3]=ii+nx;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=ii+nxy;
        //edge cg
        ii=c+in;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii+b;
        nn[ii][2]=ii+1;
        nn[ii][3]=a+in;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=ii+nxy;
        //edge dh
        ii=d+in;
        nn[ii][0]=ii-nx;
        nn[ii][1]=ii-1;
        nn[ii][2]=ii-b;
        nn[ii][3]=b+in;
        nn[ii][4]=ii-nxy;
        nn[ii][5]=ii+nxy;
    }
    
    //faces
    for(int i=1;i<nz-1;i++)
    {
        for(int j=1;j<nx-1;j++)
        {
            //face abe
            int ii=i*nxy+j;
            nn[ii][0]=ii+c;
            nn[ii][1]=ii-1;
            nn[ii][2]=ii+1;
            nn[ii][3]=ii+nx;
            nn[ii][4]=ii-nxy;
            nn[ii][5]=ii+nxy;
            //face cdg
            ii+=nxy-nx;
            nn[ii][0]=ii-nx;
            nn[ii][1]=ii-1;
            nn[ii][2]=ii+1;
            nn[ii][3]=ii-nxy+nx;
            nn[ii][4]=ii-nxy;
            nn[ii][5]=ii+nxy;
        }
        for(int j=1;j<ny-1;j++)
        {
            //face ace
            int ii=i*nxy+j*nx;
            nn[ii][0]=ii-nx;
            nn[ii][1]=ii+b;
            nn[ii][2]=ii+1;
            nn[ii][3]=ii+nx;
            nn[ii][4]=ii-nxy;
            nn[ii][5]=ii+nxy;
            //face bdf
            ii+=b;
            nn[ii][0]=ii-nx;
            nn[ii][1]=ii-1;
            nn[ii][2]=ii-b;
            nn[ii][3]=ii+nx;
            nn[ii][4]=ii-nxy;
            nn[ii][5]=ii+nxy;
        }
    }
    for(int i=1;i<ny-1;i++)
    {
        for(int j=1;j<nx-1;j++)
        {
            //face abc
            int ii=i*nx+j;
            nn[ii][0]=ii-nx;
            nn[ii][1]=ii-1;
            nn[ii][2]=ii+1;
            nn[ii][3]=ii+nx;
            nn[ii][4]=ii+e;
            nn[ii][5]=ii+nxy;
            //face efg
            ii+=e;
            nn[ii][0]=ii-nx;
            nn[ii][1]=ii-1;
            nn[ii][2]=ii+1;
            nn[ii][3]=ii+nx;
            nn[ii][4]=ii-nxy;
            nn[ii][5]=i*nx+j;
        }
    }
    
    //bulk
    for(int k=1;k<nz-1;k++)
    {
        for(int j=1;j<ny-1;j++)
        {
            for(int i=1;i<nx-1;i++)
            {
                int ii=k*nxy+j*nx+i;
                nn[ii][0]=ii-nx;
                nn[ii][1]=ii-1;
                nn[ii][2]=ii+1;
                nn[ii][3]=ii+nx;
                nn[ii][4]=ii-nxy;
                nn[ii][5]=ii+nxy;
            }
        }
    }
}

void ising3d::ScanMicroStates(double t,double b)
{
	T=t; B=b;
	int seed=-785234291;
	double rkbt=1./(kB*T);
	double sum=0.;
	std::ofstream fout;
	fout.open("m.dat");
	
	for(int j=0;j<ns;j++)
	{
		int i=int(double(N-1)*ran3(&seed));
		double dE=2.*double(s[i])*((J*SumNN(i))+(B*mo));
		if(dE<=0.) {s[i]=-s[i];}
		else {double p=exp(-dE*rkbt); if(p>=ran3(&seed)) {s[i]=-s[i];}}
		//sum+=M();
		fout<<M()<<std::endl;
	}
	sum/=ns;
	std::cout<<sum<<" "<<M()<<std::endl;
	fout.close();
}

void ising3d::ScanMicroStates(double t,double b,int mod)
{
	T=t; B=b; mM=0.;
	int seed=-785234291;
	double rkbt=1./(kB*T);
	
	for(int j=0;j<ns;j++)
	{
		int i=int(double(N-1)*ran3(&seed));
		double dE=2.*double(s[i])*((J*SumNN(i))+(B*mo));
		if(dE<=0.) {s[i]=-s[i];}
		else {double p=exp(-dE*rkbt); if(p>=ran3(&seed)) {s[i]=-s[i];}}
		if(j>=teq){int jj=(j-teq)+1; if((jj%mod)==0){mM+=fabs(M());}}
	}
	mM*=(rmts/double(mod));
}

void ising3d::SMS(double t,double b,bool ppos)
{
	T=t; B=b; mM=0.;
	int seed=-785234291;
	double rkbt=1./(kB*T), mss=0., uss=0.;
	bool flip=false;
    int fps=30;
    double nsps=0.45;
    int fc=(ns/int(double(fps)*nsps))-1;
    //std::cout<<fc<<std::endl;
	
    std::ofstream fout;
    fout.open("./images/spin_data.dat", std::ios::app);
	for(int j=0;j<ns;j++)
	{
		int i=int(double(N-1)*ran3(&seed));
		double dE=2.*double(s[i])*((J*SumNN(i))+(B*mo));
		if(dE<=0.) {s[i]=-s[i]; flip=true;}
		else {double p=exp(-dE*rkbt); if(p>=ran3(&seed)) {s[i]=-s[i]; flip=true;}}
		if(j==teq){initM(mss); initU(uss);}
		if(j>=teq){if(flip) {update1M(mss,i); update1U(uss,dE); flip=false;}update2M(mss);update2U(uss);}
        if(ppos&&((j%fc==0)&&(j>0))) {for(int k=0;k<N;k++) fout<<s[k]<<std::endl;}
	}
	rescale();
    fout.close();
}

double ising3d::SumNN(int i)
{
	double sum=0.;
	for(int j=0; j<Nbrs;j++)sum+=double(s[nn[i][j]]);
	return sum;
}

double ising3d::SumSpinProducts()
{
	double sum=0.;
	for(int i=0;i<N;i++){sum+=(double(s[i])*SumNN(i));}
	return sum;
}

double ising3d::SumSpins()
{
	double sum=0.;
	for(int i=0;i<N;i++){sum+=double(s[i]);}
	return sum;
}

double ising3d::U()
{
	return((-0.5*J*SumSpinProducts())-(B*mo*SumSpins()));
}

double ising3d::u()
{
	return(U()*rN);
}

double ising3d::M()
{
	return(SumSpins()*rN);
}

double ising3d::avgM(){ return mM;}

double ising3d::avgU(){ return mU;}

double ising3d::avgu(){return (avgU()*rN);}

double ising3d::avgC(){return (((mU2*rmts)-(pow(avgU(),2)))/(kB*pow(T,2)));}

double ising3d::avgc(){return (avgC()*rN);}

void ising3d::rescale(){mM*=rmts;mU*=rmts;}

void ising3d::initM(double &a){mM=a=(M());/* mM=fabs(a);*/ M2=pow(a,2);}

void ising3d::update1M(double &a,int i){a+=(2*double(s[i]))*rN;}

void ising3d::update2M(double a){mM+=a;/*mM+=fabs(a);*/ M2+=pow(a,2);}

void ising3d::initU(double &a){mU=a=U(); mU2=pow(mU,2);}

void ising3d::update1U(double &a,double e){a+=e;}

void ising3d::update2U(double a){mU+=a; mU2+=pow(a,2);}

double ising3d::avgX(){return (((M2*rmts)-(pow(avgM(),2)))/(kB*T));}

int ising3d::getN() {return N;}

