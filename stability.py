# This program tests stability of computing power of golden mean
# phi**n, phi =(sqrt(5.)-1.0)/2.0
# See Press et al "Numerical Recipes in C", p.30
# Since default setting for float in Python is double precision,
# we must work harder to show this instability (i.e. run upto i = 59).
# So, we use numpy half precision float16 instead.
#
import math
import numpy as np

p0 = 0.5*(math.sqrt(5.0)-1.0)

p = np.float16(1.0)
f0 = np.float16(1.0)
f1 = np.float16(p0)

for i in range(0,60):
    print("i=",i, "phi**i by mul=", p, "by **=", p0**i, "by sub=", f0)
    p = p*np.float16(p0)
    f2 = f0-f1
    f0 = f1
    f1 = f2

print("the end")
