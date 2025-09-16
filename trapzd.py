# "Numerical Recipes in C", 2nd ed, trapzd() on page 137
# We need to pass a function with one argument to it.
# We define a sample square function first.

import math

# The square function
def FUNC(x):
    return (x*x)

# trapzd takes a function as input, with a limit (a,b), at n-th iteration.
# C uses static variable s. But in Python, I better pass it as argument.
# the last argument s is the value calculated in previous iteration.
def trapzd(FUNC, a, b, n, s):
    if n==1:
        s = 0.5*(b-a)*(FUNC(a)+FUNC(b))
        return s
    else:
        it = 2**(n-2)
        d = (b-a)/it
        x = a + 0.5*d
        sum = 0.0
        for j in range(0,it):
            sum += FUNC(x)
            x += d
        s = 0.5*(s + (b-a)*sum/it)
        return s
    # end if
# end def trapzd

# Call trapzd n time to get to the desired accuracy
# qtrap() returns the integral of the function FUNC from a to b.
# The parameter EPS can be set to the desired relative accuracy
# JMAX is the maximum allowed number of steps, refining the interval
# to (b-a)/2**(JMAX-1).  Integration is performed by the trapzoidal rule.
def qtrap(FUNC, a, b):
    EPS = 1.0e-8
    JMAX = 20
    olds = -1.0e300
    s = 0.0
    for j in range(1, JMAX):
        s = trapzd(FUNC,a,b,j,s)
        if (math.fabs(s-olds) < EPS*math.fabs(olds)):
            return s
        if(s == 0.0 and olds == 0.0 and j > 6):
            return s
        olds = s
# end def qtrap

# testing int_0^1 x^2 dx:
res = qtrap(FUNC,0.0,1.0)
print(res)
