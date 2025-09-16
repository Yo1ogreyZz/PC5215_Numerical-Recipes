import math
from polint import polint
from trapzd import trapzd

def qromb(FUNC, a, b, EPS=1e-10, JMAX=25, K=5):
    s_vals, x_vals = [], []
    s_prev = 0.0
    y_best, dy_best = None, None

    for j in range(1, JMAX + 1):
        s_prev = trapzd(FUNC, a, b, j, s_prev)
        s_vals.append(s_prev)

        h = (b - a) / (1 << (j - 1))
        x_vals.append(h * h)

        if j >= K:
            xa = x_vals[j - K: j]
            ya = s_vals[j - K: j]
            y, dy = [0.0], [0.0]
            polint(xa, ya, K, 0.0, y, dy)
            y_best, dy_best = y[0], dy[0]

            if abs(dy_best) <= EPS * abs(y_best):
                return y_best, dy_best, j

    return y_best, dy_best, JMAX

def f1(x):  # ln(1 + x + x^2)
    return math.log(1.0 + x + x*x)

def f2(x):  # [ln(1 + x + x^2)]^2
    v = math.log(1.0 + x + x*x)
    return v*v

if __name__ == "__main__":
    a, b = 0.0, 4.0

    for name, f in [("ln(1+x+x^2)", f1), ("[ln(1+x+x^2)]^2", f2)]:
        val, err, iters = qromb(f, a, b, EPS=1e-10, JMAX=25, K=5)
        print(f"{name}: Romberg ≈ {val:.16g}, est.err ≈ {err:.2e}, iterations={iters}")