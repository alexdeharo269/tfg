unset key; 
unset colorbox;
set yrange[0:200]
set xrange[0:512]

set xtics 0,256,512 font "Arial, 13"
set ytics 0,100,200 font "Arial, 13"

set terminal qt 1
plot 'reduced_matrix.dat' matrix with image
set terminal qt 2
plot 'domain_labels_dbscan.dat' matrix with image
set terminal qt 3
plot 'domain_centroids_dbscan.dat' matrix with image
set terminal qt 4
set yrange [0:10000]
set ytics 0, 500, 1000
plot 'numberchanges40.txt' 