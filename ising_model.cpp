#include "ising3d.hpp"
#include <sstream>



	struct itr
	{
		double o,f,i;
		itr operator=(double a){o=f=i=a;return *this;}
		friend std::istream& operator>>(std::istream &in,itr &a)
			{in>>a.o;in>>a.f;in>>a.i;return in;}
	};

	struct array2d
	{
		int r,c;
		double **a;
		array2d(){r=c=0;}
		~array2d(){for(int i=0;i<r;i++){delete[] a[i];}delete[] a;a=NULL;}
		void init(int x,int y){r=x;c=y;a=new double*[x];for(int i=0;i<r;i++){a[i]=new double[c];}for(int i=0;i<r;i++)for(int j=0;j<c;j++)a[i][j]=0.;}
	};

	void ReadInItr(itr&,std::string&);

	int main(int argc, char* argv[])
	{
		ising3d sim;
		itr T,B;T=0.;B=0.;
		std::string pfn=argv[1];
		std::string tfn=argv[2];
		std::string bfn=argv[3];
		std::string fu="u(",fM="M(",fC="C(",fX="X(";
		std::string comma =",";
		std::string delim1= ")~[(";
		std::string delim2= ");(";
		std::string end= ")].dat";
		std::string temp= "T,";
		std::string mag= ",B";
		ReadInItr(T,tfn);
		ReadInItr(B,bfn); 
		std::stringstream  to,tf,bo,bf;
		to<<T.o;tf<<T.f;bo<<B.o;bf<<B.f; 
		std::string tr= to.str()+comma+tf.str();
		std::string bru= bo.str()+comma+bf.str();
		std::string brd= bf.str()+comma+bo.str();
		array2d uU,uM,uC,uX,dU,dM,dC,dX;
		std::string up= delim1 + tr + delim2 + bru + end, down= delim1 + tr + delim2 + brd + end;
		int x,y;
		x=int ((T.f-T.o)/T.i)+1;
		y=int ((B.f-B.o)/B.i)+1;
		uU.init(x,y);uM.init(x,y);uC.init(x,y);uX.init(x,y);dU.init(x,y);dM.init(x,y);dC.init(x,y);dX.init(x,y);
		sim.init(pfn.c_str());
		int i=0;
		for(double t=T.o;t<=T.f;t+=T.i)
		{
			int j=0;
			for(double b=B.o; b<=B.f; b+=B.i)
			{
				sim.SMS(t,b,true);
				uU.a[i][j]=sim.avgu();uM.a[i][j]=sim.avgM();uC.a[i][j]=sim.avgc();uX.a[i][j]=sim.avgX();
				j++;
			}
            if(y!=1)
            {
                j=0;
                for(double b=B.f; b>=B.o; b-=B.i)
                {
                    sim.SMS(t,b,true);
                    dU.a[i][j]=sim.avgu();dM.a[i][j]=sim.avgM();dC.a[i][j]=sim.avgc();dX.a[i][j]=sim.avgX();
                    j++;
                }
            }
			i++;
		}
		
        if(x!=1)
        {
            std::ofstream foutuU,foutuM,foutuC,foutuX,foutdU,foutdM,foutdC,foutdX;
            for(i=0;i<y;i++)
            {
                std::string snum,fnuU,fnuM,fnuC,fnuX,fndU,fndM,fndC,fndX;
                std::stringstream ssnum;
                ssnum<<B.o+(i*B.i);
                snum=ssnum.str();
                ssnum.str("");
                fnuU=fu+temp+snum+up;
                fnuM=fM+temp+snum+up;
                fnuC=fC+temp+snum+up;
                fnuX=fX+temp+snum+up;
                foutuU.open(fnuU.c_str());
                foutuM.open(fnuM.c_str());
                foutuC.open(fnuC.c_str());
                foutuX.open(fnuX.c_str());
                
                if(y!=1)
                {
                    fndU=fu+temp+snum+down;
                    fndM=fM+temp+snum+down;
                    fndC=fC+temp+snum+down;
                    fndX=fX+temp+snum+down;
                    foutdU.open(fndU.c_str());
                    foutdM.open(fndM.c_str());
                    foutdC.open(fndC.c_str());
                    foutdX.open(fndX.c_str());
                }
                
                for(int j=0;j<x;j++)
                {
                    double Ti=T.o+(j*T.i);
                    foutuU<<Ti<<" "<<uU.a[j][i]<<std::endl;
                    foutuM<<Ti<<" "<<uM.a[j][i]<<std::endl;
                    foutuC<<Ti<<" "<<uC.a[j][i]<<std::endl;
                    foutuX<<Ti<<" "<<uX.a[j][i]<<std::endl;
                    if(y!=1)
                    {
                        foutdU<<Ti<<" "<<dU.a[j][y-1-i]<<std::endl;
                        foutdM<<Ti<<" "<<dM.a[j][y-1-i]<<std::endl;
                        foutdC<<Ti<<" "<<dC.a[j][y-1-i]<<std::endl;
                        foutdX<<Ti<<" "<<dX.a[j][y-1-i]<<std::endl;
                    }
                }

                foutuU.close();
                foutuM.close();
                foutuC.close();
                foutuX.close();
                if(y!=1)
                {
                    foutdU.close();
                    foutdM.close();
                    foutdC.close();
                    foutdX.close();		
                }
            }
		}
            
		if(y!=1)
		{
            std::ofstream foutuU,foutuM,foutuC,foutuX,foutdU,foutdM,foutdC,foutdX;
			for(i=0;i<x;i++)
			{
				std::string snum,fnuU,fnuM,fnuC,fnuX,fndU,fndM,fndC,fndX;
				std::stringstream ssnum;
				ssnum<<T.o+(i*T.i);
				snum=ssnum.str();
				ssnum.str("");
				fnuU=fu+snum+mag+up;
				fndU=fu+snum+mag+down;
				fnuM=fM+snum+mag+up;
				fndM=fM+snum+mag+down;
				fnuC=fC+snum+mag+up;
				fndC=fC+snum+mag+down;
				fnuX=fX+snum+mag+up;
				fndX=fX+snum+mag+down;
				
				foutuU.open(fnuU.c_str());
				foutdU.open(fndU.c_str());
				foutuM.open(fnuM.c_str());
				foutdM.open(fndM.c_str());
				foutuC.open(fnuC.c_str());
				foutdC.open(fndC.c_str());
				foutuX.open(fnuX.c_str());
				foutdX.open(fndX.c_str());
				
				for(int j=0;j<y;j++)
				{
					double Bi=B.o+(j*B.i);
					foutuU<<Bi<<" "<<uU.a[i][j]<<std::endl;
					foutdU<<Bi<<" "<<dU.a[i][y-1-j]<<std::endl;
					foutuM<<Bi<<" "<<uM.a[i][j]<<std::endl;
					foutdM<<Bi<<" "<<dM.a[i][y-1-j]<<std::endl;
					foutuC<<Bi<<" "<<uC.a[i][j]<<std::endl;
					foutdC<<Bi<<" "<<dC.a[i][y-1-j]<<std::endl;
					foutuX<<Bi<<" "<<uX.a[i][j]<<std::endl;
					foutdX<<Bi<<" "<<dX.a[i][y-1-j]<<std::endl;
				}
				foutuU.close();
				foutdU.close();
				foutuM.close();
				foutdM.close();
				foutuC.close();
				foutdC.close();
				foutuX.close();
				foutdX.close();
			}
		}
        
        int Nspf,Nf,Nspins=0,fps=30;
        double nsps=0.45;
        Nspf=sim.getN();
        Nf=int(double(fps)*nsps);
        if(x!=1) {Nf*=x;}
        else if(y!=1) {Nf*=2*y;}
        Nspins=Nf*Nspf;
        //cout<<Nspf<<std::endl;
        int *spins;
        spins=new int[Nspins];
        for(int i=0;i<Nspins;i++) spins[i]=0;
        std::ifstream fin;
        fin.open("./images/spin_data.dat");
        i=0;
        while(!fin.eof())
        {
            fin>>spins[i];
            //cout<<Nspins<<" "<<i<<std::endl;
            //if(spins[i]!=1&&spins[i]!=-1) {cout<<i<<std::endl;}
            i++;
        }
        fin.close();
        
        std::stringstream ssnum;
        std::string snum,tag,filename,dir;
        std::ofstream outfile;
        outfile.open("./images/nframes.dat");
        outfile<<Nf<<std::endl;
        outfile.close();
        outfile.open("./images/dim.dat");
        outfile<<sim.nx<<","<<std::endl;
        outfile<<sim.ny<<","<<std::endl;
        outfile<<sim.nz<<std::endl;
        outfile.close();
        dir="./images/";
        tag=".dat";
        int k=0;
        for(i=0;i<Nf;i++)
        {
            ssnum<<i;
            snum=ssnum.str();
            ssnum.str("");
            filename=dir+snum+tag;
            outfile.open(filename.c_str());
            for(int j=0;j<Nspf-1;j++) {outfile<<spins[k]<<","<<std::endl;k++;}
            outfile<<spins[k]<<std::endl;
            outfile.close();
        }
        delete[] spins;
        
		return 0;
	}

	void ReadInItr(itr &a,std::string &fn)
	{
		std::ifstream fin;
		fin.open(fn.c_str());
		fin>>a;
		fin.close();
	}

    

