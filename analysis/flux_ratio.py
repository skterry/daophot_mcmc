#============================================================
# Calculate flux ratio between 2 or 3 stars.
# Run script from the terminal
#============================================================
import sys
import numpy as np

starNum = float(input("How many stars ('2' or '3')? "))

if starNum == 3:
	s1 = float(input("Star 1 Mag : "))
	s2 = float(input("Star 2 Mag : "))
	s3 = float(input("Star 3 Mag : "))

	s1_s2_fratio = (10**(-0.4*s1))/(((10**(-0.4*s2))))
	s1_s3_fratio = (10**(-0.4*s1))/(((10**(-0.4*s3))))

	s2_s1_fratio = (10**(-0.4*s2))/(((10**(-0.4*s1))))
	s2_s3_fratio = (10**(-0.4*s2))/(((10**(-0.4*s3))))

	s3_s1_fratio = (10**(-0.4*s3))/(((10**(-0.4*s1))))
	s3_s2_fratio = (10**(-0.4*s3))/(((10**(-0.4*s2))))
	print("-----------------------------------")
	print("Star 1/Star2 ratio:", s1_s2_fratio)
	print("Star 1/Star3 ratio:", s1_s3_fratio)
	print("Star 2/Star1 ratio:", s2_s1_fratio)
	print("Star 2/Star3 ratio:", s2_s3_fratio)
	print("Star 3/Star1 ratio:", s3_s1_fratio)
	print("Star 3/Star2 ratio:", s3_s2_fratio)

	s1_total_fratio = (10**(-0.4*s1))/((10**(-0.4*s1))+(10**(-0.4*s2))+(10**(-0.4*s3)))
	s2_total_fratio = (10**(-0.4*s2))/((10**(-0.4*s1))+(10**(-0.4*s2))+(10**(-0.4*s3)))
	s3_total_fratio = (10**(-0.4*s3))/((10**(-0.4*s1))+(10**(-0.4*s2))+(10**(-0.4*s3)))
	print("-----------------------------------")
	print("Star 1/(Star 1+2+3) ratio:", s1_total_fratio)
	print("Star 2/(Star 1+2+3) ratio:", s2_total_fratio)
	print("Star 3/(Star 1+2+3) ratio:", s3_total_fratio)

elif starNum == 2:
	s1 = float(input("Star 1 Mag : "))
	s2 = float(input("Star 2 Mag : "))

	s1_s2_fratio = (10**(-0.4*s1))/(((10**(-0.4*s2))))
	s2_s1_fratio = (10**(-0.4*s2))/(((10**(-0.4*s1))))
	print("-----------------------------------")
	print("Star 1/Star 2 ratio:", s1_s2_fratio)
	print("Star 2/Star 1 ratio:", s2_s1_fratio)

	s1_total_fratio = (10**(-0.4*s1))/((10**(-0.4*s1))+(10**(-0.4*s2)))
	s2_total_fratio = (10**(-0.4*s2))/((10**(-0.4*s1))+(10**(-0.4*s2)))
	print("-----------------------------------")
	print("Star 1/(Star 1+2) ratio:", s1_total_fratio)
	print("Star 2/(Star 1+2) ratio:", s2_total_fratio)

else:
	print("Please type 2 or 3")