#============================================================
#Script to calculate flux ratio between 2 or 3 stars.
#============================================================
import sys
import numpy as np

starNum = float(input("How many stars ('2' or '3')? "))

if starNum == 3:
	s1 = float(input("Star 1 Mag : "))
	s2 = float(input("Star 2 Mag : "))
	s3 = float(input("Star 3 Mag : "))

	s1_fratio = (10**(-0.4*s1))/((10**(-0.4*s1))+(10**(-0.4*s2))+(10**(-0.4*s3)))
	s2_fratio = (10**(-0.4*s2))/((10**(-0.4*s1))+(10**(-0.4*s2))+(10**(-0.4*s3)))
	s3_fratio = (10**(-0.4*s3))/((10**(-0.4*s1))+(10**(-0.4*s2))+(10**(-0.4*s3)))

	print("Star 1 ratio:", s1_fratio)
	print("Star 2 ratio:", s2_fratio)
	print("Star 3 ratio:", s3_fratio)

elif starNum == 2:
	s1 = float(input("Star 1 Mag : "))
	s2 = float(input("Star 2 Mag : "))

	s1_fratio = (10**(-0.4*s1))/((10**(-0.4*s1))+(10**(-0.4*s2)))
	s2_fratio = (10**(-0.4*s2))/((10**(-0.4*s1))+(10**(-0.4*s2)))

	print("Star 1 ratio:", s1_fratio)
	print("Star 2 ratio:", s2_fratio)

else:
	print("Please type 2 or 3")