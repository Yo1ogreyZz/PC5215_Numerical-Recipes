# polynomial interpolation according to Neville's algorithm,
# i.e.  "Numerical Recipes" polint() on page 109 translated into Python.
# polint() takes two lists, xa, ya, size n, and input x, and pass out
# two lists of size 1 to emulate C for y and dy.  y is interpolated value
# and dy is its error. Unlike the C code, here list index starts from 0.
#
import math


def polint(xa, ya, n, x, y, dy):
    # initialize c and d to be equal to ya

    c = ya.copy()
    d = ya.copy()
    # then find an index ns which is closest to the value x in xa()
    ns = 0
    dif = math.fabs(x - xa[0])
    for i in range(1, n):
        dift = math.fabs(x - xa[i])
        if (dift < dif):
            ns = i
            dif = dift
    y[0] = ya[ns]
    ns -= 1

    # do double loop over column m and row i on the triangular Tableau.
    # we don't try to catch the dividing by 0 error, let Python itself do it
    for m in range(1, n):
        for i in range(0, n - m):
            ho = xa[i] - x
            hp = xa[i + m] - x
            w = c[i + 1] - d[i]
            den = ho - hp
            den = w / den
            d[i] = hp * den
            c[i] = ho * den
        # end for i
        # tricking coding here in C
        if (2 * (ns + 1) < (n - m)):
            dy[0] = c[ns + 1]
        else:
            dy[0] = d[ns]
            ns -= 1
        # end if else
        y[0] += dy[0]
    # end for m
    print(c)
    print(d)


# end def polint(...)

# a test run
xa = [0, 1, 2, 4]
ya = [1, 2, 3, 0]
x = 3.0
y = [0]
dy = [0]
n = 4
polint(xa, ya, n, x, y, dy)
print("y=", y, "dy=", dy)