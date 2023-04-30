'''Auxiliary functions for solving subtasks'''

from math import comb
from re import compile, match
from sympy import Expr, symbols, simplify, degree, solve
from sympy import ln, sqrt, sin, cos


RING = 'Ring'
INTERIOR = 'Interior'


def circle_input(circle: str) -> tuple:
    '''Calculation of center and radius from the string with the circle equation
    like x^2 + a*x + y^2 + b*y + c = 0'''

    pattern = r'^x\^2[+-]\d+(?:\.\d+)?\*?x\+y\^2[+-]?\d+(?:\.\d+)?\*?y[+-]?\d+(?:\.\d+)?=\d+(?:\.\d+)?$'
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


def solver(task: str, rad_1: float, fun_1: Expr,
           rad_2=None, fun_2=None, inh=None) -> Expr:
    '''Solving the problem according to all obtained conditions'''

    rho, phi = symbols('rho phi')
    deg = int(max(degree(fun_1, rho), degree(fun_2, rho),
                  degree(fun_1, rho**-1), degree(fun_2, rho**-1)))

    if task == RING:
        c_00, c_01 = symbols('c_00 c_01')
        u = c_00*rho + c_01*ln(rho)
        c_arr = [c_00, c_01]

        for i in range(1, deg+1):
            c_0, c_1, c_2, c_3 = symbols(f'c_{i}0 c_{i}1 c_{i}2 c_{i}3')
            c_arr.extend([c_0, c_1, c_2, c_3])

            u += c_0*(rho**i) * cos(phi) + c_1*(rho**i) * sin(phi) \
                 + c_2*(rho**(-i)) * cos(phi) + c_3*(rho**(-i)) * sin(phi)

        u_1 = simplify(u.subs(rho, rad_1) - fun_1)
        u_2 = simplify(u.subs(rho, rad_2) - fun_2)

        eq_arr = [u_1.coeff(cos(phi), 0).coeff(sin(phi), 0),
                  u_2.coeff(cos(phi), 0).coeff(sin(phi), 0)]

        for i in range(1, deg+1):
            eq_arr.extend([u_1.coeff(cos(phi), i), u_1.coeff(sin(phi), i),
                           u_2.coeff(cos(phi), i), u_2.coeff(sin(phi), i)])
    elif task == INTERIOR:
        c_0 = symbols('c_00')
        u = c_0*rho
        c_arr = [c_0]

        for i in range(1, deg+1):
            c_0, c_1 = symbols(f'c_{i}0 c_{i}1')
            c_arr.extend([c_0, c_1])

            u += c_0*(rho**i) * cos(phi) + c_1*(rho**i) * sin(phi)

        u_1 = simplify(u.subs(rho, rad_1) - fun_1)

        eq_arr = [u_1.coeff(cos(phi), 0).coeff(sin(phi), 0)]

        for i in range(1, deg+1):
            eq_arr.extend([u_1.coeff(cos(phi), i), u_1.coeff(sin(phi), i)])
    else:
        c_0 = symbols('c_00')
        u = c_0*rho
        c_arr = [c_0]

        for i in range(1, deg+1):
            c_0, c_1 = symbols(f'c_{i}0 c_{i}1')
            c_arr.extend([c_0, c_1])

            u += c_0*(rho**(-i)) * cos(phi) + c_1*(rho**(-i)) * sin(phi)

        u_1 = simplify(u.subs(rho, rad_1) - fun_1)

        eq_arr = [u_1.coeff(cos(phi), 0).coeff(sin(phi), 0)]

        for i in range(1, deg+1):
            eq_arr.extend([u_1.coeff(cos(phi), i), u_1.coeff(sin(phi), i)])


    u = simplify(u.subs(c_arr, solve(eq_arr, c_arr)))

    return u
