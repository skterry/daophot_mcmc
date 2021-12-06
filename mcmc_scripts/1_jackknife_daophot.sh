#for i in {1..4}; do #STEP 1
#head -16 input_1_daophot.txt >> input_1_daophot.txt
#sed -i -e "$((16*i)),$((16*(i+1)))s/mag_OB110950_1/mag_OB110950_$((i+1))/g" input_daophot.txt
#done


#daophot-mcmc < input_1_daophot.txt #STEP 2