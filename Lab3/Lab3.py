from math import *
from pyLabOn import Report
from time import process_time

default_eps = 1e-6


def romberg_integral(integral_range, N,  eps, f):
    a, b = integral_range
    T = []
    S = []
    C = []
    R = []
    h = b - a
    T.append(h * (f(a) + f(b)) / 2)
    ep = eps+1
    m = 1
    while(m < N):
        t = 0
        pow_2_m = 1 << m
        dh = h/pow_2_m
        for i in range(1, pow_2_m, 2):
            t += f(a+i*dh)*dh
        t += T[-1]/2
        T.append(t)
        pow_4_m = 1 << (m << 1)
        if m >= 1:
            S.append((pow_4_m*T[-1]-T[-2])/(pow_4_m-1))
        if m >= 2:
            C.append((pow_4_m*S[-1]-S[-2])/(pow_4_m-1))
        if m >= 3:
            R.append((pow_4_m*C[-1]-C[-2])/(pow_4_m-1))
        if m > 4 and abs((R[-1]-R[-2])) < eps:
            break
        m += 1

    return R, m-1


def f_test(x):
    return x


def f1(x):
    return x*x*exp(x)


def f2(x):
    return exp(x)*sin(x)


def f3(x):
    return 4/(1+x*x)


def f4(x):
    return 1/(x+1)


def lab_log(index, r_result):
    result, m = r_result
    return [str(index), str(result[-1]), str(m)]


report = Report("Lab3 Results", "Lab3")

header = ["Index", "Result", "Iteration Count"]

data = []

data.append(
    lab_log(1,
            romberg_integral(
                (0, 1), 100, default_eps, f1)
            )
)

data.append(
    lab_log(2,
            romberg_integral(
                (1,3), 100, default_eps, f2)
            )
)

data.append(
    lab_log(3,
            romberg_integral(
                (0, 1), 100, default_eps, f3)
            )
)

data.append(
    lab_log(4,
            romberg_integral(
                (0, 1), 100, default_eps, f4)
            )
)

report.add_table(header,data)

report.compile()