'''Auxiliary functions for solving subtasks'''

from math import comb
from re import compile, match

from sympy import (Add, Expr, Float, Mul, Pow, acos, cos, exp, expand,
                   integrate, ln, preorder_traversal, simplify, sin, solve,
                   sqrt, symbols)
from sympy.abc import phi, rho, t, x, y

RING = 'Ring'
INT = 'Interior'
EXT = 'Exterior'

EXPL = 'Explicit'
IMPL = 'Implicit'


def circle_input(circle: str) -> tuple:
    '''Calculate a center and radius from the string with
    the circle equation like x^2 + a*x + y^2 + b*y = c'''

    pattern = r'^x\^2([+-](\d+(?:\.\d+)?\*?)?x)?\+y\^2([+-](\d+(?:\.\d+)?\*?)?y)?=[+-]?\d+(?:\.\d+)?$'
    pattern = compile(pattern)
    match_flg = match(pattern, circle)
    if not match_flg:
        raise ValueError('Invalid circle equation: ' + circle)

    lhs, rhs = circle.split('=')
    equation = simplify((simplify(lhs) - simplify(rhs)))

    a = equation.coeff(x, 1) / 2.0
    b = equation.coeff(y, 1) / 2.0
    c = - equation.coeff(x, 0).coeff(y, 0)

    x_0, y_0 = -a, -b
    radius = sqrt(c + a**2 + b**2)

    return x_0, y_0, radius


def humbler(function: str|Expr, x_0=0.0, y_0=0.0) -> Expr:
    '''Simplify a function and decrease the degree of cosines and sines'''

    if isinstance(function, str):
        function = expand(simplify(function).subs([(x, x_0 + rho * cos(phi)),
                                                   (y, y_0 + rho * sin(phi))]))

    cos_deg = {expr.args[1] for expr in preorder_traversal(function) \
                if isinstance(expr, Pow) and isinstance(expr.args[0], cos)}
    sin_deg = {expr.args[1] for expr in preorder_traversal(function) \
                if isinstance(expr, Pow) and isinstance(expr.args[0], sin)}

    for n in cos_deg:
        sum_cos = 0
        if n % 2 == 0:
            for k in range(n//2):
                sum_cos += (comb(n, k) * cos((n-2*k)*phi))

            sum_cos += (comb(n, n//2) / 2)
        else:
            for k in range((n+1)//2):
                sum_cos += (comb(n, k) * cos((n-2*k)*phi))

        sum_cos /= 2**(n-1)
        function = function.subs({cos(phi)**n: sum_cos})

    for n in sin_deg:
        sum_sin = 0
        if n % 2 == 0:
            for k in range(n//2):
                sum_sin += ((-1)**(n//2-k) * comb(n, k) * cos((n-2*k)*phi))

            sum_sin += (comb(n, n//2) / 2)
        else:
            for k in range((n+1)//2):
                sum_sin += ((-1)**((n-1)//2-k) * comb(n, k) * sin((n-2*k)*phi))

        sum_sin /= 2**(n-1)
        function = function.subs(sin(phi)**n, sum_sin)

    return expand(function).rewrite(exp).expand().rewrite(cos).expand()


def inhomogeneity(x_0: float, y_0: float, inh_s: str) -> Expr | int:
    '''Solve the inhomogeneity in the equation, give a partial solution'''

    if inh_s == '0':
        return 0

    inh = humbler(inh_s, x_0, y_0)

    if isinstance(inh, Add):
        u_arr = list(inh.args)
    else:
        u_arr = [inh]

    for i, arg in enumerate(u_arr):
        mul = arg.args
        f_r = Mul(mul[0], mul[1]) if len(mul) == 3 else mul[0]
        g_p = mul[-1]

        if g_p.args[0].args:
            k = g_p.args[0].args[0]
        else:
            k = 1

        f_r = f_r.subs(rho, exp(rho))

        u = (exp(k*t) * integrate(exp(rho * (2-k)) * f_r, (rho, 1, t))
             - exp(-k*t) * integrate(exp(rho * (2+k)) * f_r, (rho, 1, t))) \
            / (2*k)

        u = u.expand()
        u -= (exp(k*t) * u.coeff(exp(k*t)) + exp(-k*t) * u.coeff(exp(-k*t)))

        u = humbler(u.subs(t, ln(rho)))
        u_arr[i] = u * g_p

    return Add(*u_arr)


def solver(task: str, rad_1: float, fun_1: Expr,
           rad_2=0.0, fun_2=Expr()) -> Expr:
    '''Solve the problem according to all obtained conditions'''

    c_0, c_ln = symbols('c_0 c_ln')
    u = c_0
    c_arr = {c_0}
    if task == RING:
        u += c_ln * ln(rho)
        c_arr.add(c_ln)

    cos_args = {expr.args[0].args[0] if expr.args[0].args else 1 \
                for expr in (fun_1 + fun_2).atoms(cos)}
    sin_args = {expr.args[0].args[0] if expr.args[0].args else 1 \
                for expr in (fun_1 + fun_2).atoms(sin)}

    for i in cos_args:
        c_c0, c_c1 = symbols(f'c_c0_{i} c_c1_{i}')
        if task == RING:
            c_arr.update({c_c0, c_c1})
            u += c_c0*(rho**i) * cos(i*phi) \
                + c_c1*(rho**-i) * cos(i*phi)
        elif task == INT:
            c_arr.add(c_c0)
            u += c_c0*(rho**i) * cos(i*phi)
        else:
            c_arr.add(c_c1)
            u += c_c1*(rho**-i) * cos(i*phi)

    for i in sin_args:
        c_s0, c_s1 = symbols(f'c_s0_{i} c_s1_{i}')
        if task == RING:
            c_arr.update({c_s0, c_s1})
            u += c_s0*(rho**i) * sin(i*phi) \
                + c_s1*(rho**-i) * sin(i*phi)
        elif task == INT:
            c_arr.add(c_s0)
            u += c_s0*(rho**i) * sin(i*phi)
        else:
            c_arr.add(c_s1)
            u += c_s1*(rho**-i) * sin(i*phi)

    u_1 = humbler(u.subs(rho, rad_1) - fun_1.subs(rho, rad_1))
    if task == RING:
        u_2 = humbler(u.subs(rho, rad_2) - fun_2.subs(rho, rad_2))

    u_10 = u_1.coeff(cos(phi), 0).coeff(sin(phi), 0)
    if task == RING:
        u_20 = u_2.coeff(cos(phi), 0).coeff(sin(phi), 0)

    for i in cos_args | sin_args:
        u_10 = u_10.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0)
        if task == RING:
            u_20 = u_20.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0)

    eq_arr = {u_10} | {u_1.coeff(cos(i*phi)) for i in cos_args} | \
        {u_1.coeff(sin(i*phi)) for i in sin_args}
    if task == RING:
        eq_arr |= {u_20} | {u_2.coeff(cos(i*phi)) for i in cos_args} | \
            {u_2.coeff(sin(i*phi)) for i in sin_args}

    c_sol = solve(eq_arr, c_arr)
    u = simplify(u.subs(c_sol))

    return u


def runner(task: str, expl: str, cent: str, rad_s: str,
           circ: str, fun_s: str, inh_s: str) -> tuple:
    '''Run side functions with the task in mind'''

    if task == RING:
        if expl == EXPL:
            x_0, y_0 = map(float, cent.split(';'))
            rad_1, rad_2 = map(float, rad_s.split(';'))
        else:
            circ_1, circ_2 = circ.split(';')
            x_0, y_0, rad_1 = map(float, circle_input(circ_1))
            x_0, y_0, rad_2 = map(float, circle_input(circ_2))

        fun_1, fun_2 = fun_s.split(';')
        fun_1 = humbler(fun_1, x_0, y_0)
        fun_2 = humbler(fun_2, x_0, y_0)
        u_par = inhomogeneity(x_0, y_0, inh_s)

        u_rp = solver(task, rad_1, fun_1-u_par,
                      rad_2=rad_2, fun_2=fun_2-u_par) + u_par
    else:
        if expl == EXPL:
            x_0, y_0 = map(float, cent.split(';'))
            rad = float(rad_s.replace(' ', ''))
        else:
            x_0, y_0, rad = map(float, circle_input(circ))

        fun = humbler(fun_s, x_0, y_0)
        u_par = inhomogeneity(x_0, y_0, inh_s)

        u_rp = solver(task, rad, fun-u_par) + u_par

    rho_new = sqrt((x - x_0)**2 + (y - y_0**2))
    phi_new = acos((x - x_0) / rho_new)

    u_xy = simplify(u_rp.subs([(rho, rho_new), (phi, phi_new)]))

    u_rp_args = [arg for arg in preorder_traversal(u_rp) \
                 if isinstance(arg, Float)]
    u_rp = u_rp.subs([(arg, round(arg, 4)) for arg in u_rp_args])

    u_xy_args = [arg for arg in preorder_traversal(u_xy) \
                 if isinstance(arg, Float)]
    u_xy = u_xy.subs([(arg, round(arg, 4)) for arg in u_xy_args])

    u_rp, u_xy = map(humbler, [u_rp, u_xy])

    return u_rp, u_xy
