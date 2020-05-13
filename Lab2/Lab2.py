import numpy as np
from math import *
from fractions import Fraction
import fractions

default_eps1 = 1e-6
default_eps2 = 1e-4

lab_result = open("Lab2\\lab_result.md", "w+")

lab_result.write("# lab2 results\n")
lab_result.write("\n")

def log_format(log_value):
    if isinstance(log_value, list):
        r = ""
        for value in log_value:
            r += log_format(value)
        return r
    else:
        b, v = log_value
        if b:
            return str(v)+'\n'
        else:
            return 'NaN\n'



def lab_log(key, value):
    global lab_result

    lab_result.write("## "+key+"\n")
    lab_result.write("\n")
    lab_result.write(log_format(value)+"\n")


def newton_iteration_method(alpha, epsilon1, epsilon2, N, f, df):
    x_this = alpha
    x_next = 0.0
    for n in range(0, N):
        F = f(x_this)
        DF = df(x_this)
        if abs(F) < epsilon1:
            return True, x_this
        if abs(DF) < epsilon2:
            return False, 0.0
        x_next = x_this-F/DF
        if abs(x_next-x_this) < epsilon1:
            return True, x_next
        x_this = x_next
    return False, 0.0


def f11(x): return cos(x)-x


def df11(x): return -sin(x)-1


alpha11 = pi/4


def f12(x): return exp(-x)-sin(x)


def df12(x): return -exp(-x)-cos(x)


alpha12 = 0.6


def f21(x): return x-exp(-x)


def df21(x): return 1+exp(-x)


alpha21 = 0.5


def f22(x): return x*x-2*x*exp(-x)+exp(-2*x)


def df22(x): return 2*x-(2*x-2)*exp(-x)-2*exp(-2*x)


alpha22 = 0.5


class polynomial:
    def __init__(self):
        super().__init__()
        self.coefficients = np.asarray(Fraction(0), dtype=Fraction)

    def __init__(self, coefficients):
        super().__init__()
        self.coefficients = coefficients

    def __zero_extend(self, other):
        coe1 = self.coefficients
        coe2 = other.coefficients
        length = max(len(coe1), len(coe2))

        coe1 = np.pad(coe1, (0, length-len(coe1)),
                      'constant', constant_values=Fraction(0))
        coe2 = np.pad(coe2, (0, length-len(coe2)),
                      'constant', constant_values=Fraction(0))
        return coe1, coe2

    def __add__(self, other):
        if not isinstance(other, polynomial):
            return self
        coe1, coe2 = self.__zero_extend(other)
        return polynomial(coe1+coe2)

    def __sub__(self, other):
        if not isinstance(other, polynomial):
            return self
        coe1, coe2 = self.__zero_extend(other)
        return polynomial(coe1-coe2)

    def mulx(self):
        return polynomial(np.pad(self.coefficients, (1, 0), 'constant', constant_values=Fraction(0)))

    def __mul__(self, other):
        if(isinstance(other, polynomial)):
            r = polynomial()
            self_base = self
            for c in other.coefficients:
                r += self_base*c
                self_base = self_base.mulx
            return r
        else:
            return polynomial(self.coefficients*other)

    def __truediv__(self, other):
        return polynomial(self.coefficients/other)

    def __call__(self, x):
        r = 0.0
        x_base = 1.0
        for c in self.coefficients:
            r += c*x_base
            x_base *= x
        return r

    def __repr__(self):
        index = 0
        r = ""
        flag = False
        for c in self.coefficients:
            if c == 0:
                index += 1
                continue
            if flag:
                r += '+' if c >= 0 else '-'
            if c != Fraction(1) or index == 0:
                r += str(abs(c.numerator)) if c.denominator == 1 else'\\frac{'+str(
                    abs(c.numerator))+"}{"+str(c.denominator)+"}"
            if index == 0:
                pass
            else:
                r += 'x'
                if index != 1:
                    r += "^{"+str(index)+"}"
            flag = True
            index += 1
        if r == "":
            return "0"
        return r

    def get_derivative(self):
        length = len(self.coefficients)
        if length <= 1:
            return polynomial(np.asarray([Fraction(0)], dtype=Fraction))
        derivative_coefficients = []
        for index in range(1, length):
            derivative_coefficients.append(
                self.coefficients[index]*Fraction(index)
            )
        return polynomial(np.asarray(
            derivative_coefficients, dtype=Fraction))

    def __str__(self):
        index = 0
        r = ""
        for c in self.coefficients:
            if c == 0:
                index += 1
                continue
            if index != 0 and c >= 0:
                r += '+'
            r += str(c)
            if index == 0:
                pass
            else:
                r += 'x'
                if index != 1:
                    r += str(index)
            index += 1

        return r


def generate_polynomial(size, gen_coe0, gen_coe1, gen_coe2):
    P = [
        polynomial(np.asarray([Fraction(1)], dtype=Fraction)),
        polynomial(np.asarray([Fraction(0), Fraction(1)], dtype=Fraction))
    ]
    for i in range(0, size-1):
        P.append(P[i+1]*gen_coe0(i)+P[i+1].mulx()
                 * gen_coe1(i)-(P[i]*gen_coe2(i)))
    return P


def generate_polynomial_md(name, prefix, P):
    f = open('lab2\\'+name+'.md', "w+")
    f.write("# "+name+"\n\n")
    index = 0
    for p in P:
        f.write("$$\n"+prefix+"_{"+str(index)+"}(x)="+repr(p)+"\n$$\n\n")
        index += 1
    f.close()

# 根据给定列表生成关于原点对称的区间


def generate_luminus_element(origin_list):
    origin_list.extend(-np.asarray(origin_list))
    return origin_list


# Legendre
P = generate_polynomial(
    6,
    lambda i: Fraction(0),
    lambda i: Fraction((2*i+3), (i+2)),
    lambda i: Fraction((i+1), (i+2))
)

alpha_legendre = generate_luminus_element(
    [0.9324695142, 0.6612093865, 0.2386191861]
)

# Chebyshev

T = generate_polynomial(
    6,
    lambda i: Fraction(0),
    lambda i: Fraction(2),
    lambda i: Fraction(1)
)

alpha_chebyshev = [cos((2*j+1)/(14)*pi) for j in range(0, 7)]

# Laguerre

L = generate_polynomial(
    6,
    lambda i: Fraction(2*i+3),
    lambda i: Fraction(-1),
    lambda i: Fraction((i+1)**2)
)

alpha_laguerre = [0.2635603197, 1.4134030591,
                  3.5964257710, 7.0858100059, 12.6408008443]

# Hermite

H = generate_polynomial(
    6,
    lambda i: Fraction(0),
    lambda i: Fraction(2),
    lambda i: Fraction(2*(i+1))
)

alpha_hermite = generate_luminus_element(
    [2.3506049737, 1.3358490740, 0.4360774119]
)

# Part 1

lab_log(
    "Part 1.1",
    newton_iteration_method(alpha11, default_eps1,
                            default_eps2, 10, f11, df11)
)

lab_log(
    "Part 1.2",
    newton_iteration_method(alpha12, default_eps1,
                            default_eps2, 10, f12, df12)
)

# Part 2

lab_log(
    "Part 2.1",
    newton_iteration_method(alpha21, default_eps1,
                            default_eps2, 10, f21, df21)
)

lab_log(
    "Part 2.2",
    newton_iteration_method(alpha22, default_eps1,
                            default_eps2, 20, f22, df22)
)

# Part 3

generate_polynomial_md("Legendre", 'P', P)

lab_log(
    "Part 3.1",
    [
        newton_iteration_method(
            alpha, default_eps1, default_eps1,
            6, P[6], P[6].get_derivative()
        )
        for alpha in alpha_legendre
    ]
)

generate_polynomial_md("Chebyshev", 'T', T)

lab_log(
    "Part 3.2",
    [
        newton_iteration_method(
            alpha, default_eps1, default_eps1,
            6, T[6], T[6].get_derivative()
        )
        for alpha in alpha_chebyshev
    ]
)

generate_polynomial_md("Laguerre", 'L', L)

lab_log(
    "Part 3.3",
    [
        newton_iteration_method(
            alpha, default_eps1, default_eps1,
            6, L[6], L[6].get_derivative()
        )
        for alpha in alpha_laguerre
    ]
)

generate_polynomial_md("Hermite", 'H', H)

lab_log(
    "Part 3.4",
    [
        newton_iteration_method(
            alpha, default_eps1, default_eps1,
            6, H[6], H[6].get_derivative()
        )
        for alpha in alpha_hermite
    ]
)

lab_result.close()
