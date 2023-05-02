'''Auxiliary functions for solving subtasks'''

from math import comb
from re import compile, match
from sympy import init_printing, Expr, symbols, simplify, degree, solve
from sympy import acos, sqrt, ln, sin, cos
# from sympy.plotting.plot import plot3d


init_printing()


RING = 'Ring'
INT = 'Interior'
EXT = 'Exterior'

EXPL = 'Explicit'
IMPL = 'Implicit'


def circle_input(circle: str) -> tuple:
    '''Calculation of center and radius from the string with the circle equation
    like x^2 + a*x + y^2 + b*y + c = d'''

    pattern = r'^x\^2[+-]\d+(?:\.\d+)?\*?x\+y\^2[+-]?\d' \
              + r'+(?:\.\d+)?\*?y[+-]?\d+(?:\.\d+)?=\d+(?:\.\d+)?$'
    pattern = compile(pattern)

    match_flg = match(pattern, circle)
    if not match_flg:
        raise ValueError('Invalid circle equation: ' + circle)


    x, y = symbols('x y')
    lhs, rhs = circle.split('=')
    equation = simplify((simplify(lhs) - simplify(rhs)))

    a = equation.coeff(x, 1) / 2
    b = equation.coeff(y, 1) / 2
    c = equation.coeff(x, 0).coeff(y, 0)

    x_0, y_0 = a, b
    radius = sqrt(a**2 + b**2 - c)

    return x_0, y_0, radius


def function_input(x_0: float, y_0: float, fun_entr: str, rad=None) -> Expr:
    '''Simplifying a function obtained in a string'''

    x, y, rho, phi = symbols('x, y, rho, phi')

    function = simplify(fun_entr)
    function = simplify(function.subs([(x, x_0 + rho * cos(phi)),
                                       (y, y_0 + rho * sin(phi))]))


    max_deg = max(degree(function, sin(phi)), degree(function, cos(phi)))

    for n in range(max_deg, 1, -1):
        temp_sin = temp_cos = 0
        if n % 2 == 0:
            for k in range(n//2):
                temp_sin += ((-1)**(n//2-k) * comb(n, k) * cos((n-2*k)*phi))
                temp_cos += (comb(n, k) * cos((n-2*k)*phi))

            temp_sin += (comb(n, n//2) / 2)
            temp_cos += (comb(n, n//2) / 2)
        else:
            for k in range((n+1)//2):
                temp_sin += ((-1)**((n-1)//2-k) * comb(n, k) * sin((n-2*k)*phi))
                temp_cos += (comb(n, k) * cos((n-2*k)*phi))

        temp_sin /= 2**(n-1)
        temp_cos /= 2**(n-1)

        function = function.subs(sin(phi)**n, temp_sin)
        function = function.subs(cos(phi)**n, temp_cos)


    if rad is not None:
        function = simplify(function.subs(rho, rad))

    return function


def runner(task: str, expl: str, cent, rad, circ, fun, inh) -> None:
    '''Running side functions with the task in mind'''

    cent = cent.get().replace(' ', '')
    rad = rad.get().replace(' ', '')
    circ = circ.get().replace(' ', '')
    fun = fun.get().replace(' ', '')
    inh = inh.get().replace(' ', '')

    if task == RING:
        if expl == EXPL:
            x_0, y_0 = map(float, cent.split(';'))
            rad_1, rad_2 = map(float, rad.split(';'))
        else:
            circ_1, circ_2 = circ.split(';')
            x_0, y_0, rad_1 = map(float, circle_input(circ_1))
            x_0, y_0, rad_2 = map(float, circle_input(circ_2))

        fun_1, fun_2 = fun.split(';')
        fun_1 = function_input(x_0, y_0, fun_1, rad_1)
        fun_2 = function_input(x_0, y_0, fun_2, rad_2)
        inh_ex = function_input(x_0, y_0, inh)

        u = solver(task, rad_1, fun_1, inh_ex,
                   rad_2=rad_2, fun_2=fun_2)
    else:
        if expl == EXPL:
            x_0, y_0 = map(float, cent.split(';'))
            rad_ex = float(rad.replace(' ', ''))
        else:
            x_0, y_0, rad_ex = map(float, circle_input(circ))

        fun_ex = function_input(x_0, y_0, fun, rad)
        inh_ex = function_input(x_0, y_0, inh)

        u = solver(task, rad_ex, fun_ex, inh_ex)

    x, y, rho, phi = symbols('x y rho phi')
    rho_new = sqrt((x - x_0)**2 + (y - y_0**2))
    phi_new = acos((x - x_0) / rho_new)

    u = simplify(u.subs((rho, rho_new), (phi, phi_new)))

    # plot3d(u)


def solver(task: str, rad_1: float, fun_1: Expr, inh,
           rad_2=None, fun_2=None) -> Expr:
    '''Solving the problem according to all obtained conditions'''

    rho, phi = symbols('rho phi')
    max_deg = int(max(degree(fun_1, rho), degree(fun_2, rho),
                      degree(fun_1, rho**-1), degree(fun_2, rho**-1)))

    if task == RING:
        c_00, c_01 = symbols('c_00 c_01')
        u = c_00*rho + c_01*ln(rho)
        c_arr = [c_00, c_01]

        for i in range(1, max_deg+1):
            c_0, c_1, c_2, c_3 = symbols(f'c_{i}0 c_{i}1 c_{i}2 c_{i}3')
            c_arr.extend([c_0, c_1, c_2, c_3])

            u += c_0*(rho**i) * cos(phi) + c_1*(rho**i) * sin(phi) \
                 + c_2*(rho**(-i)) * cos(phi) + c_3*(rho**(-i)) * sin(phi)

        u_1 = simplify(u.subs(rho, rad_1) - fun_1)
        u_2 = simplify(u.subs(rho, rad_2) - fun_2)

        eq_arr = [u_1.coeff(cos(phi), 0).coeff(sin(phi), 0),
                  u_2.coeff(cos(phi), 0).coeff(sin(phi), 0)]

        for i in range(1, max_deg+1):
            eq_arr.extend([u_1.coeff(cos(phi), i), u_1.coeff(sin(phi), i),
                           u_2.coeff(cos(phi), i), u_2.coeff(sin(phi), i)])

    elif task == INT:
        c_0 = symbols('c_00')
        u = c_0*rho
        c_arr = [c_0]

        for i in range(1, max_deg+1):
            c_0, c_1 = symbols(f'c_{i}0 c_{i}1')
            c_arr.extend([c_0, c_1])

            u += c_0*(rho**i) * cos(phi) + c_1*(rho**i) * sin(phi)

        u_1 = simplify(u.subs(rho, rad_1) - fun_1)

        eq_arr = [u_1.coeff(cos(phi), 0).coeff(sin(phi), 0)]

        for i in range(1, max_deg+1):
            eq_arr.extend([u_1.coeff(cos(phi), i), u_1.coeff(sin(phi), i)])

    else:
        c_0 = symbols('c_00')
        u = c_0*rho
        c_arr = [c_0]

        for i in range(1, max_deg+1):
            c_0, c_1 = symbols(f'c_{i}0 c_{i}1')
            c_arr.extend([c_0, c_1])

            u += c_0*(rho**(-i)) * cos(phi) + c_1*(rho**(-i)) * sin(phi)

        u_1 = simplify(u.subs(rho, rad_1) - fun_1)

        eq_arr = [u_1.coeff(cos(phi), 0).coeff(sin(phi), 0)]

        for i in range(1, max_deg+1):
            eq_arr.extend([u_1.coeff(cos(phi), i), u_1.coeff(sin(phi), i)])


    u = simplify(u.subs(c_arr, solve(eq_arr, c_arr)))

    return u
