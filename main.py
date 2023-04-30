'''There's a lot of spaghetti code. Sorry.'''

from tkinter import ttk, Tk, Label, StringVar
from sympy import init_printing, simplify, symbols, acos, sqrt
from sympy.plotting.plot import plot3d

from stuff import circle_input, function_input, solver


init_printing()


RING = 'Ring'
INTERIOR = 'Interior'
EXTERIOR = 'Exterior'

EXPLICIT = 'Explicit'
IMPLICIT = 'Implicit'


def main() -> None:
    '''Build the GUI and select the problem'''

    def runner():
        '''Running side functions with the task in mind'''

        if task.get() == RING:
            if expl.get() == EXPLICIT:
                cent = cent_entr.get().replace(' ', '')
                x_0, y_0 = map(float, cent.split(';'))
                rad = rad_entr.get().replace(' ', '')
                rad_1, rad_2 = map(float, rad.split(';'))
            else:
                circ = circ_entr.get().replace(' ', '')
                circ_1, circ_2 = circ.split(';')
                x_0, y_0, rad_1 = map(float, circle_input(circ_1))
                x_0, y_0, rad_2 = map(float, circle_input(circ_2))

            # inh = inh_entr.get().replace(' ', '')
            # inh = function_input(x_0, y_0, inh)

            fun = fun_entr.get().replace(' ', '')
            fun_1, fun_2 = fun.split(';')
            fun_1 = function_input(x_0, y_0, fun_1, rad_1)
            fun_2 = function_input(x_0, y_0, fun_2, rad_2)

            u = solver(task.get(), rad_1, fun_1, rad_2=rad_2,
                       fun_2=fun_2) # , inh=inh)
        else:
            if expl.get() == EXPLICIT:
                cent = cent_entr.get().replace(' ', '')
                x_0, y_0 = map(float, cent.split(';'))
                rad = float(rad_entr.get().replace(' ', ''))
            else:
                circ = circ_entr.get().replace(' ', '')
                x_0, y_0, rad = map(float, circle_input(circ))

            # inh = inh_entr.get().replace(' ', '')
            # inh = function_input(x_0, y_0, inh)

            fun = function_input(x_0, y_0, fun_entr.get().replace(' ', ''), rad)

            u = solver(task.get(), rad, fun) # , inh=inh)
        
        x, y, rho, phi = symbols('x y rho phi')
        rho_new = sqrt((x - x_0)**2 + (y - y_0**2))
        phi_new = acos((x - x_0) / rho_new)
        u = simplify(u.subs((rho, rho_new),
                            (phi, phi_new)))

        plot3d(u)


    root = Tk()
    root.title('Solution of the Dirichlet-Neumann problem')
    root.geometry('300x340')


    task_lbl = Label(text='Problem type:')
    task_lbl.grid(row=0, column=1)

    task = StringVar(value=RING)

    ring_rdbtn = ttk.Radiobutton(text=RING, value=RING, variable=task)
    ring_rdbtn.grid(row=1, column=0)
    in_rdbtn = ttk.Radiobutton(text=INTERIOR, value=INTERIOR, variable=task)
    in_rdbtn.grid(row=1, column=1)
    out_rdbtn = ttk.Radiobutton(text=EXTERIOR, value=EXTERIOR, variable=task)
    out_rdbtn.grid(row=1, column=2)


    expl_lbl = Label(text='Circle parameters are given:')
    expl_lbl.grid(row=2, column=1)

    expl = StringVar(value=EXPLICIT)

    expl_yes_rdbtn = ttk.Radiobutton(text=EXPLICIT, value=EXPLICIT, variable=expl)
    expl_yes_rdbtn.grid(row=3, column=0)
    expl_no_rdbtn = ttk.Radiobutton(text=IMPLICIT, value=IMPLICIT, variable=expl)
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
