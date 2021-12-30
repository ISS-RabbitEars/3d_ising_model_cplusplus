#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#define PI 3.141592653589793238462643383

int main()
{
	int N;
	double theta,delta;
    std::ifstream in;
    std::ofstream file;
	
    in.open("nframes.dat");
    in>>N;
    in.close();
    
    delta=2*PI/(N-1);
	theta=0;
    
    std::string pts="cam_";
	for(int i=0;i<N;i++)
	{
		theta=i*delta;
        std::string pu,snum;
        std::stringstream ssnum;
        ssnum<<i;
        snum=ssnum.str();
        ssnum.str("");
        pu=pts+snum;
        file.open(pu.c_str());
        file<<cos(theta)<<","<<std::endl;
        file<<sin(theta)<<std::endl;
        file.close();
	}
    
	return 0;
}

#undef PI

