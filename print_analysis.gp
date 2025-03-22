unset key; 
unset colorbox;
#set yrange[0:200]
#set xrange[0:256]

#set xtics 0,256,562 font "Arial, 13"
#set ytics 0,100,200 font "Arial, 13"

set terminal qt 1
plot './reduced/reduced56.dat' matrix with image
set terminal qt 2
plot './domains/domains56.dat' matrix with image
set terminal qt 3
plot 'centroids/centroids56.dat' matrix with image
set terminal qt 4
plot './raw/ising_dataR56T0.0.dat'  matrix with image
set yrange [0:20000]
set ytics 0, 1000, 10000
plot './changes/changes56.dat'