'''Документация к модулю'''
from sympy import Expr, symbols, ln, sin, cos, solve


def solver(task: int, radius,
           function, heterogeneity: Expr) -> Expr:
    '''Решение полученной задачи'''

    u, rho, phi = symbols('u, rho, phi')
    c1, c2, c3, c4, c5, c6 = symbols('c1, c2, c3, c4, c5, c6')
    u = c1*rho + c2*ln(rho) + c3*rho*cos(phi) + c4*rho*sin(phi) \
        + c5*(1/rho)*sin(phi) + c6*(1/rho)*cos(phi)

    if task == 1:
        function = tuple(map(Expr, function))
        function_in, function_out = (function[0].subs(radius[0]),
                                     function[1].subs(radius[1]))

        u = solve([function_in, function_out],
                      c1, c2, c3, c4, c5, c6)
    else:
        function = Expr(function)
        function = function.subs(radius)
        u = solve(function,
                      [c1, c2, c3, c4, c5, c6])

    return u
