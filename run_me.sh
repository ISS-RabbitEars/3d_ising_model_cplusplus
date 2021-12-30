mkdir images
mkdir data
g++ ising_model.cpp -o ising3dsim
./ising3dsim parameters.dat T_itr.dat B_itr.dat
rm ising3dsim
mv *].dat data/
cp create_movie.sh images/
cp magnet.pov images/
cp parameters.dat images/
cp campath.cpp images/
cd images
./create_movie.sh

