import math


def lagrange(data, n, x):

    if len(data) != n+1:
        return
    y = 0.0
    for k in range(0, n+1):
        xk, fxk = data[k]
        l = 1.0
        for j in range(0, n+1):
            if j != k:
                xj, __ = data[j]
                l = (l*(x-xj))/(xk-xj)
        y += l*fxk
    return x, y


def get_data(xk_array, f):
    ret = []
    for xk in xk_array:
        ret.append((xk, f(xk)))
    return ret


def get_xk_linear(f, n, range_length):
    ret = []
    x0 = -range_length/2
    h = range_length/n
    for k in range(0, n+1):
        ret.append(x0+k*h)
    return ret


def get_xk_cos(f, n, range_length):
    ret = []
    x0 = -range_length/2
    for k in range(0, n+1):
        ret.append(math.cos((2*k+1)*math.pi/(2*(n+1))))
    return ret


def run_lagrange(name, x_array, n_array, data_generator, f):
    txt = open("Lab1\\"+name+".txt", "w+")
    txt.write("n/x")
    for x in x_array:
        txt.write('\t'+str(x))
    txt.write('\n')

    for n in n_array:
        txt.write(str(n))
        for x in x_array:
            _, y = lagrange(get_data(data_generator(n),f), n, x)
            txt.write('\t'+str(y))
        txt.write('\n')

    txt.write('inf')
    for x in x_array:
        y = f(x)
        txt.write('\t'+str(y))
    txt.write('\n')
    print("Lagrange "+name+ " Done.")


def f1(x):
    return 1.0/(1.0+x*x)


x_range5 = [0.75, 1.75, 2.75, 3.75, 4.75]
x_range1 = [-0.95, -0.05, 0.05, 0.95]
x_range_part4 = [5, 50, 115, 185]
n_range = [5, 10, 20]

run_lagrange("part1.1", x_range5, n_range,
             lambda n: get_xk_linear(f1, n, 10.0), f1)
run_lagrange("part1.2", x_range1, n_range,
             lambda n: get_xk_linear(math.exp, n, 2.0), math.exp)
run_lagrange("part2.1", x_range1, n_range,
             lambda n: get_xk_linear(f1, n, 2.0), f1)
run_lagrange("part2.2", x_range5, n_range,
             lambda n: get_xk_linear(math.exp, n, 10.0), math.exp)
run_lagrange("part3.1", x_range1, n_range,
             lambda n: get_xk_cos(f1, n, 2.0), f1)
run_lagrange("part3.2", x_range1, n_range,
             lambda n: get_xk_cos(math.exp, n, 2.0), math.exp)
run_lagrange("part4.1", x_range_part4, [2], lambda n: [1, 4, 9], math.sqrt)
run_lagrange("part4.2", x_range_part4, [2], lambda n: [36, 49, 64], math.sqrt)
run_lagrange("part4.3", x_range_part4, [2], lambda n: [100, 121, 144], math.sqrt)
run_lagrange("part4.4", x_range_part4, [2], lambda n: [169,196, 225], math.sqrt)

print("Done.")
