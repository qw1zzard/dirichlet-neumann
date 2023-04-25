''''''

from tkinter import ttk, Tk, Label, StringVar
from sympy import init_printing, Expr, symbols, solve, simplify
from sympy import ln, sqrt, sin, cos
# from sympy.plotting.plot import plot3d

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

    center = [a, b]
    radius = sqrt(a**2 + b**2 - c)

    return center, radius


def function_input(x_0, y_0, fun_entr: ttk.Entry) -> Expr:
    '''Обработка введённой функции'''

    x, y, rho, phi = symbols('x, y, rho, phi')

    function = simplify(fun_entr.get())
    function = simplify(function.subs([(x, x_0 + rho * cos(phi)),
                                       (y, y_0 + rho * sin(phi))]))

    return function


def solver(rad_1, rad_2=None,
           fun_1, fun_2=None, inh=None) -> Expr:
    '''Решение полученной задачи'''

    u, rho, phi = symbols('u, rho, phi')
    c1, c2, c3, c4, c5, c6 = symbols('c1, c2, c3, c4, c5, c6')
    u = c1*rho + c2*ln(rho) + c3*rho*cos(phi) + c4*rho*sin(phi) \
        + c5*(1/rho)*sin(phi) + c6*(1/rho)*cos(phi)

    function = Expr(function)
    function = function.subs(radius)
    u = solve(function, [c1, c2, c3, c4, c5, c6])

    return u


def runner():
    '''Запуск всех побочных функций'''



def main() -> None:
    '''Построение интерфейса и выбор задачи'''

    root = Tk()
    root.title('Решение задачи Дирихле и Неймана')
    root.geometry('320x600')


    task_lbl = Label(text='Тип задачи:')
    task_lbl.grid(row=0, column=1)

    type_ring = 'Кольцо'
    type_in = 'Внутренняя'
    type_out = 'Внешняя'

    global task
    task = StringVar()

    ring_rdbtn = ttk.Radiobutton(text=type_ring, value=type_ring, variable=task)
    ring_rdbtn.grid(row=1, column=0)
    in_rdbtn = ttk.Radiobutton(text=type_in, value=type_in, variable=task)
    in_rdbtn.grid(row=1, column=1)
    out_rdbtn = ttk.Radiobutton(text=type_out, value=type_out, variable=task)
    out_rdbtn.grid(row=1, column=2)


    expl_lbl = Label(text='Параметры окружности заданы:')
    expl_lbl.grid(row=2, column=1)

    expl_yes = 'Явно'
    expl_no = 'Неявно'

    global expl
    expl = StringVar()

    expl_yes_rdbtn = ttk.Radiobutton(text=expl_yes, value=expl_yes, variable=expl)
    expl_yes_rdbtn.grid(row=3, column=0)
    expl_no_rdbtn = ttk.Radiobutton(text=expl_no, value=expl_no, variable=expl)
    expl_no_rdbtn.grid(row=3, column=1)



    # if expl == expl_yes:
    cent_lbl = Label(text='Центр окружностей:')
    cent_lbl.grid(row=4, column=1)

    cent_entr = ttk.Entry()
    cent_entr.grid(row=5, column=1)

    rad_lbl = Label(text='Радиусы окружностей:')
    rad_lbl.grid(row=6, column=1)

    rad_1_entr = ttk.Entry()
    rad_1_entr.grid(row=7, column=1)
    rad_2_entr = ttk.Entry()
    rad_2_entr.grid(row=8, column=1)
    # else:
    circ_lbl = Label(text='Уравнения окружностей:')
    circ_lbl.grid(row=9, column=1)

    circ_1_entr = ttk.Entry()
    circ_1_entr.grid(row=10, column=1)
    circ_2_entr = ttk.Entry()
    circ_2_entr.grid(row=11, column=1)


    fun_lbl = Label(text='Функции на границе:')
    fun_lbl.grid(row=12, column=1)

    fun_1_entr = ttk.Entry()
    fun_1_entr.grid(row=13, column=1)
    fun_2_entr = ttk.Entry()
    fun_2_entr.grid(row=14, column=1)


    inh_lbl = Label(text='Функция неоднородности:')
    inh_lbl.grid(row=15, column=1)

    inh_entr = ttk.Entry()
    inh_entr.grid(row=16, column=1)


    solve_btn = ttk.Button(text='Решить задачу', command=runner)
    solve_btn.grid(row=17, column=1)



    # fun_1 = function_input(fun_1_entr, x_0, y_0)
    # fun_2 = function_input(fun_2_entr, x_0, y_0)
    # inh = function_input(inh_entr, x_0, y_0)


    # u = solver(task, *condition)

    # plot3d(u)


    root.mainloop()


if __name__ == '__main__':
    main()
