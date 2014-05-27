#!/bin/sh

cd /Users/arianasp/Desktop/Baboons/google_earth/slinky_interactions/noise_thresh_10/kmls

for f in `ls *.kml`;
do(
zip "${f%.kml}.kmz" baboon_black.png baboon_blue.png baboon_red.png "$f.kml"


)
done
