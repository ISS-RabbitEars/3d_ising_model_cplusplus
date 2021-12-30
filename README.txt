3D Ising Model

Requirements:
1. g++
2. povray
3. ffmpeg


there are 3 parameter files that may be modified:
1. parameters.dat
2. T_itr.dat
3. B_itr.dat


parameters.dat
1. the first 3 lines of parameters.dat are the dimensions for the ising model:
nx
ny
nz
2. the 4th line is the value for J, the interaction energy
3. the 5th line is the scaled value for kB (Boltzmann Constant)
4. the 6th line is how strongly the system couples to the external magnetic field 
(2,3,and 4 have all been set to 1.)


T_itr.dat
1. initial temperature
2. final temperature
3. temperature increment

B_itr.dat
1. initial and final value for the external magnetic field
2. midpoint value for the external magnetic field
3. external magentic field increment

*(the magnetic field begins at some initial value, increments to the midpoint value and then returns to the initial value, unless both initial and midpoint are set to zero which effectively turns it "off")
*(one can range the Temperature or the External Magnetic Field while leaving the other constant but not both simultaneously)

 
