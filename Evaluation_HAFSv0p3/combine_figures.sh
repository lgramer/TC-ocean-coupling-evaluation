#!/bin/sh

##################### User input ########################
prefix1=05.2022090306.hycom.temp_trans_along_track.f
prefix2=05.2022090306.hycom.temp_diff_trans_along_track.f
prefix_comb=05.2022090306.hycom.temp_trans_along_track_combine.f

fhours_to_combine="000 006 012 018 024 030 036 042 048 054 060 066 072 078 084 090 096 102 108 114 120"

# direction=- is for figures to be combined vertically 
# direction=+ is for figures to be combined horizontally 
direction=-

##################### End user input ########################

# Combine figures
for fhour in $fhours_to_combine;
do 
    echo "$fhour"
    convert ${prefix1}${fhour}.png ${prefix2}${fhour}.png ${direction}append ${prefix_comb}${fhour}.png
done 

# Make gif video
echo "Combining png figures into gif video"
convert -delay 60 -loop 0 ${prefix_comb}*.png ${prefix_comb}.gif

