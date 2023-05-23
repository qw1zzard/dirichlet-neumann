"""Auxiliary functions to solve subtasks"""

from math import comb
from re import compile as compile_pattern
from re import match

from sympy import (Add, Expr, Float, Mul, Pow, acos, cos, exp, expand,
                   integrate, ln, simplify, sin, solve, sqrt, symbols)
from sympy.abc import phi, rho, t, x, y

RING = 'Ring'
INT = 'Interior'
EXT = 'Exterior'

EXPL = 'Explicit'
IMPL = 'Implicit'


def circle_input(circle: str) -> tuple:
    """ Calculate the center and radius from the line with the equation
    of a circle of the form x^2 + a*x + y^2 + b*y = c"""

    pattern = r'^x\^2([+-](\d+(?:\.\d+)?\*?)?x)?\+y\^2([+-](\d+(?:\.\d+)?\*?)?y)?=[+-]?\d+(?:\.\d+)?$'
    pattern = compile_pattern(pattern)
    match_flg = match(pattern, circle)
    if not match_flg:
        raise ValueError('Invalid circle equation: ' + circle)

    lhs, rhs = circle.split('=')
    equation = simplify((simplify(lhs) - simplify(rhs)))

    a_coef = equation.coeff(x, 1) / 2.0
    b_coef = equation.coeff(y, 1) / 2.0
    c_coef = - equation.coeff(x, 0).coeff(y, 0)

    x_0, y_0 = -a_coef, -b_coef
    radius = sqrt(c_coef + a_coef**2 + b_coef**2)

    return x_0, y_0, radius


def humbler(function: str | Expr, x_0=0.0, y_0=0.0) -> Expr:
    """Simplify to the power of cosines and sines"""

    if isinstance(function, str):
        function = expand(simplify(function).subs([(x, x_0 + rho * cos(phi)),
                                                   (y, y_0 + rho * sin(phi))]))

    cos_deg = {int(expr.args[1]) for expr in function.atoms(Pow)
               if isinstance(expr.args[0], cos)}
    sin_deg = {int(expr.args[1]) for expr in function.atoms(Pow)
               if isinstance(expr.args[0], sin)}

    for deg in cos_deg:
        sum_cos = 0
        if deg % 2 == 0:
            for k in range(deg//2):
                sum_cos += Mul(comb(deg, k), cos((deg-2*k)*phi))

            sum_cos += (comb(deg, deg//2) / 2)
        else:
            for k in range((deg + 1)//2):
                sum_cos += Mul(comb(deg, k), cos((deg-2*k)*phi))

        sum_cos /= 2**(deg-1)
        function = function.subs({Pow(cos(phi), deg): sum_cos})

    for deg in sin_deg:
        sum_sin = 0
        if deg % 2 == 0:
            for k in range(deg//2):
                sum_sin += ((-1)**(deg//2-k) * comb(deg, k)
                            * cos((deg-2*k)*phi))

            sum_sin += (comb(deg, deg//2) / 2)
        else:
            for k in range((deg + 1)//2):
                sum_sin += ((-1)**((deg-1)//2-k) *
                            comb(deg, k) * sin((deg-2*k)*phi))

        sum_sin /= 2**(deg-1)
        function = function.subs(Pow(sin(phi), deg), sum_sin)

    return expand(function).rewrite(exp).expand().rewrite(cos).expand()


def inhomogeneity(x_0: float, y_0: float, inh_s: str) -> Expr | int:
    """Solve the inhomogeneity in the equation, give a partial solution"""

    if inh_s == '0':
        return 0

    inh = humbler(inh_s, x_0, y_0)

    if isinstance(inh, Add):
        u_arr = list(inh.args)
    else:
        u_arr = [inh]

    for i, arg in enumerate(u_arr):
        mul = arg.args
        f_rho = Mul(mul[0], mul[1]) if len(mul) == 3 else mul[0]
        g_phi = mul[-1]

        if g_phi.args[0].args:
            k = g_phi.args[0].args[0]
        else:
            k = 1

        f_rho = f_rho.subs(rho, exp(rho))

        u_par = (
            Mul(exp(k*t),  integrate(Mul(exp(rho * (2-k)), f_rho), (rho, 1, t)))
            - Mul(exp(-k*t), integrate(Mul(exp(rho * (2+k)), f_rho), (rho, 1, t)))
        ) / Mul(2, k)

        u_par = u_par.expand()
        u_par -= (exp(k*t) * u_par.coeff(exp(k*t)) +
                  exp(-k*t) * u_par.coeff(exp(-k*t)))

        u_par = humbler(u_par.subs(t, ln(rho)))
        u_arr[i] = u_par * g_phi

    return Add(*u_arr)


def dirichlet(task: str, rad_1: float, fun_1: Expr,
              rad_2=0.0, fun_2=Expr()) -> Expr:
    """Solve the Dirichlet problem according to all obtained conditions"""

    c_0, c_ln = symbols('c_0 c_ln')
    u_gen = c_0
    c_arr = {c_0}
    if task == RING:
        u_gen += c_ln * ln(rho)
        c_arr.add(c_ln)

    cos_args = {expr.args[0].args[0] if expr.args[0].args else 1
                for expr in (fun_1 + fun_2).atoms(cos)}
    sin_args = {expr.args[0].args[0] if expr.args[0].args else 1
                for expr in (fun_1 + fun_2).atoms(sin)}

    for i in cos_args:
        c_c0, c_c1 = symbols(f'c_c0_{i} c_c1_{i}')
        if task == RING:
            c_arr.update({c_c0, c_c1})
            u_gen += (c_c0*(rho**i) * cos(i*phi) + c_c1*(rho**-i) * cos(i*phi))
        elif task == INT:
            c_arr.add(c_c0)
            u_gen += c_c0*(rho**i) * cos(i*phi)
        else:
            c_arr.add(c_c1)
            u_gen += c_c1*(rho**-i) * cos(i*phi)

    for i in sin_args:
        c_s0, c_s1 = symbols(f'c_s0_{i} c_s1_{i}')
        if task == RING:
            c_arr.update({c_s0, c_s1})
            u_gen += (c_s0*(rho**i) * sin(i*phi) + c_s1*(rho**-i) * sin(i*phi))
        elif task == INT:
            c_arr.add(c_s0)
            u_gen += c_s0*(rho**i) * sin(i*phi)
        else:
            c_arr.add(c_s1)
            u_gen += c_s1*(rho**-i) * sin(i*phi)

    u_1 = humbler(u_gen.subs(rho, rad_1) - fun_1.subs(rho, rad_1))
    u_10 = u_1.coeff(cos(phi), 0).coeff(sin(phi), 0)
    if task == RING:
        u_2 = humbler(u_gen.subs(rho, rad_2) - fun_2.subs(rho, rad_2))
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
    u_gen = simplify(u_gen.subs(c_sol))

    return u_gen


def neumann(task: str, rad_1: float, fun_1: Expr,
            rad_2=0.0, fun_2=Expr()) -> Expr | str:
    """Solve the Neumann problem according to all obtained conditions"""

    fun_1 = humbler(fun_1.subs(rho, rad_1))

    cos_1_args = {expr.args[0].args[0] if expr.args[0].args else 1
                  for expr in fun_1.atoms(cos)}
    sin_1_args = {expr.args[0].args[0] if expr.args[0].args else 1
                  for expr in fun_1.atoms(sin)}

    cos_1_coef = [fun_1.coeff(cos(arg*phi)) for arg in cos_1_args]
    sin_1_coef = [fun_1.coeff(sin(arg*phi)) for arg in sin_1_args]

    u_10 = fun_1.coeff(cos(phi), 0).coeff(sin(phi), 0)
    for i in cos_1_args | sin_1_args:
        u_10 = u_10.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0)

    if task == RING:
        fun_2 = humbler(fun_2.subs(rho, rad_2))

        cos_2_args = {expr.args[0].args[0] if expr.args[0].args else 1
                      for expr in fun_2.atoms(cos)}
        sin_2_args = {expr.args[0].args[0] if expr.args[0].args else 1
                      for expr in fun_2.atoms(sin)}

        cos_2_coef = [fun_2.coeff(cos(arg*phi)) for arg in cos_2_args]
        sin_2_coef = [fun_2.coeff(sin(arg*phi)) for arg in sin_2_args]

        u_20 = fun_2.coeff(cos(phi), 0).coeff(sin(phi), 0)
        for i in cos_2_coef + sin_2_coef:
            u_20 = u_20.coeff(cos(i*phi), 0).coeff(sin(i*phi), 0)

        if u_10 != u_20:
            return 'не выполнено условие'

        rad = rad_2 - rad_1

        u_gen = u_10 * ((rad_2 - rho) / rad) ** 2 / 2.0 \
            + u_20 * ((rho - rad_1) / rad) ** 2 / 2.0

        for coef, arg in zip(cos_1_coef, cos_1_args):
            u_gen += (((rad_2 - rho) / rad) ** (arg + 1)
                      * (rad * coef / (arg + 1)) * cos(arg * phi))

        for coef, arg in zip(cos_2_coef, cos_2_args):
            u_gen += (((rho - rad_1) / rad) ** (arg + 1)
                      * (rad * coef / (arg + 1)) * cos(arg * phi))

        for coef, arg in zip(sin_1_coef, sin_1_args):
            u_gen += (((rad_2 - rho) / rad) ** (arg + 1)
                      * (rad * coef / (arg + 1)) * sin(arg * phi))

        for coef, arg in zip(sin_2_coef, sin_2_args):
            u_gen += (((rho - rad_1) / rad) ** (arg + 1)
                      * (rad * coef / (arg + 1)) * sin(arg * phi))
    else:
        if u_10:
            return 'не выполнено условие'

        delta = 1
        if task == EXT:
            delta = -1

        u_gen = 0
        for coef, arg in zip(cos_1_coef, cos_1_args):
            u_gen += ((rho / rad_1) ** (delta * arg)
                      * (coef * rad_1 / arg) * cos(arg * phi))

        for coef, arg in zip(sin_1_coef, sin_1_args):
            u_gen += ((rho / rad_1) ** (delta * arg)
                      * (coef * rad_1 / arg) * sin(arg * phi))

        u_gen = humbler(u_gen)

    return u_gen


def runner(task: str, expl: str, cent: str, rad_s: str,
           circ: str, fun_s: str, inh_s: str) -> tuple:
    """Run side functions with the task in mind"""

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

        u_dir_rp = Add(dirichlet(task, rad_1, fun_1-u_par,
                                 rad_2=rad_2, fun_2=fun_2-u_par), u_par)
        u_neu_rp = neumann(task, rad_1, fun_1, rad_2=rad_2, fun_2=fun_2)
    else:
        if expl == EXPL:
            x_0, y_0 = map(float, cent.split(';'))
            rad = float(rad_s.replace(' ', ''))
        else:
            x_0, y_0, rad = map(float, circle_input(circ))

        fun = humbler(fun_s, x_0, y_0)
        u_par = inhomogeneity(x_0, y_0, inh_s)

        u_dir_rp = Add(dirichlet(task, rad, fun-u_par), u_par)
        u_neu_rp = neumann(task, rad, fun)

    u_dir_rp_args = [arg for arg in u_dir_rp.atoms(Float)]
    u_dir_rp = u_dir_rp.subs([(arg, round(arg, 4)) for arg in u_dir_rp_args])

    if isinstance(u_neu_rp, Expr):
        u_neu_rp_args = [arg for arg in u_neu_rp.atoms(Float)]
        u_neu_rp = u_neu_rp.subs([(arg, round(arg, 4))
                                 for arg in u_neu_rp_args])

    rho_new = sqrt((x - x_0)**2 + (y - y_0**2))
    phi_new = acos((x - x_0) / rho_new)
    u_dir_xy = u_dir_rp.subs([(rho, rho_new), (phi, phi_new)])

    u_dir_xy_args = [arg for arg in u_dir_xy.atoms(Float)]
    u_dir_xy = u_dir_xy.subs([(arg, round(arg, 4)) for arg in u_dir_xy_args])

    u_dir_rp, u_dir_xy = map(humbler, [u_dir_rp, u_dir_xy])

    return u_dir_rp, u_neu_rp
