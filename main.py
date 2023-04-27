'''Документация к модулю'''

import re
from tkinter import *
from tkinter import ttk
from sympy import init_printing, Expr, symbols, solve, simplify
from sympy import ln, sqrt, sin, cos
from sympy.plotting.plot import plot3d

init_printing()


def circle_input(circle: str) -> tuple:
    '''Обработка строки с уравнением окружности'''

    circle = circle.replace(' ', '')
    a = b = c = ''

    i = 3 # x^2 пропускаем
    while circle[i] != '*':
        a += circle[i]
        i += 1
    a = float(a) / 2

    i += 6 # *x+y^2 пропускаем
    while circle[i] != '*':
        b += circle[i]
        i += 1
    b = float(b) / 2

    i += 2 # *y пропускаем
    while circle[i] != '=':
        c += circle[i]
        i += 1
    c = float(c)

    x_0, y_0 = a, b
    radius = sqrt(a**2 + b**2 - c)

    return x_0, y_0, radius


def function_input(x_0: float, y_0: float, fun_entr: str) -> Expr:
    '''Обработка введённой функции'''

    x, y, rho, phi = symbols('x, y, rho, phi')

    function = simplify(fun_entr)
    function = simplify(function.subs([(x, x_0 + rho * cos(phi)),
                                       (y, y_0 + rho * sin(phi))]))

    return function


def solver(rad_1: float, fun_1: Expr,
           rad_2=None, fun_2=None, inh=None) -> Expr:
    '''Решение полученной задачи'''

    u, rho, phi = symbols('u rho phi')
    c1, c2, c3, c4, c5, c6 = symbols('c1 c2 c3 c4 c5 c6')

    # Define the solution function
    u_expr = c1*rho + c2*ln(rho) + c3*rho*cos(phi) + c4*rho*sin(phi) \
        + c5*(1/rho)*sin(phi) + c6*(1/rho)*cos(phi)

    # Substitute the radii
    u_expr = u_expr.subs({rad_1: rho, rad_2: rho})

    # Substitute the functions on the circles
    if fun_1:
        u_expr = u_expr.subs(u.subs(rho, rad_1), fun_1)
    if fun_2:
        u_expr = u_expr.subs(u.subs(rho, rad_2), fun_2)

    # Substitute the inhomogeneity function
    if inh:
        u_expr = u_expr.subs(c1, inh)

    # Solve for the coefficients
    coeffs = solve(u_expr, [c1, c2, c3, c4, c5, c6])

    return u_expr.subs(coeffs)


def main() -> None:
    '''Build the GUI and select the problem'''

    def runner():
        # implementation of the solver function

        x_0, y_0 = map(float, cent_entr.get().split())

        if task.get() == type_ring:
            if expl.get() == expl_yes:
                rad_1, rad_2 = map(float, rad_entr.get().split())
            else:
                circ_1, circ_2 = circ_entr.get().split()
                x_0, y_0, rad_1 = map(float, circle_input(circ_1))
                x_0, y_0, rad_2 = map(float, circle_input(circ_2))

            fun_1, fun_2 = fun_entr.get().split()
            fun_1 = function_input(x_0, y_0, fun_1)
            fun_2 = function_input(x_0, y_0, fun_2)

            u = solver(rad_1, fun_1, rad_2=rad_2,
                       fun_2=fun_2, inh=inh_entr.get())
        else:
            if expl.get() == expl_yes:
                rad = float(rad_entr.get())
            else:
                circ = circ_entr.get()
                x_0, y_0, rad = map(float, circle_input(circ))

            fun = function_input(x_0, y_0, fun_entr.get())

            u = solver(rad, fun, inh=inh_entr.get())

        plot3d(u)


    root = Tk()
    root.title('Solution of the Dirichlet-Neumann problem')
    root.geometry('300x340')


    task_lbl = Label(text='Problem type:')
    task_lbl.grid(row=0, column=1)

    type_ring = 'Ring'
    type_in = 'Interior'
    type_out = 'Exterior'

    task = StringVar()

    ring_rdbtn = ttk.Radiobutton(text=type_ring, value=type_ring, variable=task)
    ring_rdbtn.grid(row=1, column=0)
    in_rdbtn = ttk.Radiobutton(text=type_in, value=type_in, variable=task)
    in_rdbtn.grid(row=1, column=1)
    out_rdbtn = ttk.Radiobutton(text=type_out, value=type_out, variable=task)
    out_rdbtn.grid(row=1, column=2)


    expl_lbl = Label(text='Circle parameters are given:')
    expl_lbl.grid(row=2, column=1)

    expl_yes = 'Explicit'
    expl_no = 'Implicit'

    expl = StringVar()

    expl_yes_rdbtn = ttk.Radiobutton(text=expl_yes, value=expl_yes, variable=expl)
    expl_yes_rdbtn.grid(row=3, column=0)
    expl_no_rdbtn = ttk.Radiobutton(text=expl_no, value=expl_no, variable=expl)
    expl_no_rdbtn.grid(row=3, column=1)


    cent_lbl = Label(text='Circle center:')
    cent_lbl.grid(row=4, column=1)

    cent_entr = ttk.Entry()
    cent_entr.grid(row=5, column=1)


    rad_lbl = Label(text='Circle radius:')
    rad_lbl.grid(row=6, column=1)

    rad_entr = ttk.Entry()
    rad_entr.grid(row=7, column=1)


    circ_lbl = Label(text='Circle equation:')
    circ_lbl.grid(row=8, column=1)

    circ_entr = ttk.Entry()
    circ_entr.grid(row=9, column=1)


    fun_lbl = Label(text='Function on the boundary:')
    fun_lbl.grid(row=10, column=1)

    fun_entr = ttk.Entry()
    fun_entr.grid(row=11, column=1)


    inh_lbl = Label(text='Non-homogeneous function:')
    inh_lbl.grid(row=12, column=1)

    inh_entr = ttk.Entry()
    inh_entr.grid(row=13, column=1)


    solve_btn = ttk.Button(text='Solve problem', command=runner)
    solve_btn.grid(row=14, column=1)

    root.mainloop()


if __name__ == '__main__':
    main()
