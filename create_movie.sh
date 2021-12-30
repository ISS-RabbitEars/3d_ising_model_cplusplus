#!/bin/sh

Nframes=`cat nframes.dat`
fname="frame_data.dat"
cname="cam_temp.dat"
povfile="magnet"

g++ campath.cpp -o campath
./campath
rm campath

j=0
while [ "${j}" -lt "${Nframes}" ]
do
mv ${j}.dat ${fname}
mv cam_${j} ${cname}
povray -D +A +H1600 -J +Q9 +R5 -V +W1600 ${povfile}.pov
mv ${povfile}.png ${j}.png
rm ${fname}
rm ${cname}
#mv ${fname} ${j}.dat
#convert -resize 800x800 ${j}.png small_${j}.png
j=$(($j+1))
done

ffmpeg -framerate 30 -i %d.png ising3danim.mp4
#ffmpeg -i '%d.png' -s 800x800 -r 30 -vcodec qtrle movie.mov
#ffmpeg -f image2 -i small_%d.png movie.mpg
