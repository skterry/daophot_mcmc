#!/bin/bash
for i in {1..2} # STEP 1
do
head -16 input_daophot.txt >> input_daophot.txt
sed -i -e "$((16*i)),$((16*(i+1)))s/image_1/image_$((i+1))/g" input_daophot.txt
done

#daophot-mcmc < input_daophot.txt  # STEP 2