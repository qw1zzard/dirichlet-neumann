'''Auxiliary functions for solving subtasks'''

from math import comb
from re import compile, match
from sympy import Expr, Float, symbols, simplify, degree, solve
from sympy import acos, sqrt, ln, sin, cos, preorder_traversal


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


def function_input(x_0: float, y_0: float, fun_entr: str) -> Expr:
    '''Simplifying a function obtained in a string'''

    x, y, rho, phi = symbols('x, y, rho, phi')

    function = simplify(fun_entr)
    function = simplify(function.subs([(x, x_0 + rho * cos(phi)),
                                       (y, y_0 + rho * sin(phi))]))


    max_deg = int(max(0, degree(function, sin(phi)),
                      degree(function, cos(phi))))

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

    return function


def solver(task: str, rad_1: float, fun_1: Expr, inh: Expr,
           rad_2=0.0, fun_2=Expr()) -> Expr:
    '''Solving the problem according to all obtained conditions'''

    rho, phi, c_0, c_ln  = symbols('rho phi c_0 c_ln')
    u = c_0
    c_arr = [c_0]
    if task == RING:
        u += c_ln * ln(rho)
        c_arr.append(c_ln)

    cos_args, sin_args = set(), set()
    for expr in preorder_traversal(fun_1 + fun_2):
        if isinstance(expr, cos):
            cos_args.add(1)
            if expr.args[0].args:
                cos_args.add(expr.args[0].args[0])

        elif isinstance(expr, sin):
            sin_args.add(1)
            if expr.args[0].args:
                sin_args.add(expr.args[0].args[0])

    for i in cos_args:
        c_c0, c_c1 = symbols(f'c_c0_{i} c_c1_{i}')
        if task == RING:
            c_arr.extend([c_c0,  c_c1])
            u += c_c0*(rho**i) * cos(i*phi) \
                + c_c1*(rho**(-i)) * cos(i*phi)

        elif task == INT:
            c_arr.extend([c_c0])
            u += c_c0*(rho**i) * cos(i*phi)

        else:
            c_arr.extend([c_c1])
            u += c_c1*(rho**(-i)) * cos(i*phi)

    for i in sin_args:
        c_s0, c_s1 = symbols(f'c_s0_{i} c_s1_{i}')
        if task == RING:
            c_arr.extend([c_s0, c_s1])
            u += c_s0*(rho**i) * sin(i*phi) \
                + c_s1*(rho**(-i)) * sin(i*phi)

        elif task == INT:
            c_arr.extend([c_s0])
            u += c_s0*(rho**i) * sin(i*phi)

        else:
            c_arr.extend([c_s1])
            u += c_s1*(rho**(-i)) * sin(i*phi)

    u_1 = simplify(u.subs(rho, rad_1) - fun_1.subs(rho, rad_1))
    if task == RING:
        u_2 = simplify(u.subs(rho, rad_2) - fun_2.subs(rho, rad_2))

    eq_arr = []
    for i in cos_args:
        eq_arr.extend([u_1.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0),
                        u_1.coeff(cos(i*phi), i)])
        if task == RING:
            eq_arr.extend([u_2.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0),
                           u_2.coeff(cos(i*phi), i)])

    for i in sin_args:
        eq_arr.extend([u_1.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0),
                        u_1.coeff(sin(i*phi), i)])
        if task == RING:
            eq_arr.extend([u_2.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0),
                           u_2.coeff(sin(i*phi), i)])

    eq_arr = set(eq_arr)
    c_sol = solve(eq_arr, c_arr)
    u = simplify(u.subs(c_sol))

    return u


def runner(task: str, expl: str, cent: str, rad: str,
           circ: str, fun: str, inh: str) -> tuple:
    '''Running side functions with the task in mind'''

    if task == RING:
        if expl == EXPL:
            x_0, y_0 = map(float, cent.split(';'))
            rad_1, rad_2 = map(float, rad.split(';'))
        else:
            circ_1, circ_2 = circ.split(';')
            x_0, y_0, rad_1 = map(float, circle_input(circ_1))
            x_0, y_0, rad_2 = map(float, circle_input(circ_2))

        fun_1, fun_2 = fun.split(';')
        fun_1 = function_input(x_0, y_0, fun_1)
        fun_2 = function_input(x_0, y_0, fun_2)
        inh_x = function_input(x_0, y_0, inh)

        u_rp = solver(task, rad_1, fun_1, inh_x,
                      rad_2=rad_2, fun_2=fun_2)
    else:
        if expl == EXPL:
            x_0, y_0 = map(float, cent.split(';'))
            rad_ex = float(rad.replace(' ', ''))
        else:
            x_0, y_0, rad_ex = map(float, circle_input(circ))

        fun_x = function_input(x_0, y_0, fun)
        inh_x = function_input(x_0, y_0, inh)

        u_rp = solver(task, rad_ex, fun_x, inh_x)

    x, y, rho, phi = symbols('x y rho phi')
    rho_new = sqrt((x - x_0)**2 + (y - y_0**2))
    phi_new = acos((x - x_0) / rho_new)

    u_xy = simplify(u_rp.subs([(rho, rho_new), (phi, phi_new)]))


    for arg in preorder_traversal(u_rp):
        if isinstance(arg, Float):
            u_rp = u_rp.subs(arg, round(arg, 4))

    for arg in preorder_traversal(u_xy):
        if isinstance(arg, Float):
            u_xy = u_xy.subs(arg, round(arg, 4))

    return u_rp, u_xy
