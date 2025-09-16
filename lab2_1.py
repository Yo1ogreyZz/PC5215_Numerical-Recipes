import numpy as np
import matplotlib.pyplot as plt
from polint import polint

xa = [-1.0, 1.0, 2.0, 4.0]
ya = [ 1.25, 2.0, 3.0, 0.0]
n = len(xa)

def y_at(x):
    y, dy = [0.0], [0.0]
    polint(xa, ya, n, x, y, dy)
    return y[0]

xs = np.linspace(-1.0, 4.0, 600)
ys = [y_at(x) for x in xs]

plt.figure()
plt.plot(xs, ys, '-', label='Neville cubic interpolation')
plt.plot(xa, ya, 'o', label='data points')
plt.xlabel('x'); plt.ylabel('y'); plt.legend(); plt.tight_layout()
plt.show()