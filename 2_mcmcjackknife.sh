#!/bin/bash
for i in {1..2} # STEP 1
do
head -16 input_2_mcmcdaophot.txt >> input_2_mcmcdaophot.txt
sed -i -e "$((16*i)),$((16*(i+1)))s/image_1/image_$((i+1))/g" input_2_mcmcdaophot.txt
done

#daophot-mcmc < input_2_mcmcdaophot.txt  # STEP 2