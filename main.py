'''There's a lot of spaghetti code. Sorry.'''

from functools import partial
from tkinter import Tk, Label, StringVar
from tkinter.ttk import Radiobutton, Entry, Button
from sympy.plotting.plot import plot3d
from sympy import init_printing

from stuff import runner


init_printing()


RING = 'Ring'
INT = 'Interior'
EXT = 'Exterior'

EXPL = 'Explicit'
IMPL = 'Implicit'

RELY = 1/19.0


def informer(info: Label, task: StringVar, expl: StringVar, cent: Entry,
             rad: Entry, circ: Entry, fun: Entry, inh: Entry) -> None:
    '''Collect the information, transfer it to the solver and build a graph'''

    task_s = task.get()
    expl_s = expl.get()
    cent_s = cent.get().replace(' ', '')
    rad_s = rad.get().replace(' ', '')
    circ_s = circ.get().replace(' ', '')
    fun_s = fun.get().replace(' ', '')
    inh_s = inh.get().replace(' ', '')

    u_rp, u_xy = runner(task_s, expl_s, cent_s, rad_s, circ_s, fun_s, inh_s)

    u_rp_s = 'u(rho, phi) = ' + str(u_rp)
    Label(text=u_rp_s).place(relx=.5, rely=17*RELY, anchor='center')

    u_xy_s = 'u(x, y) = ' + str(u_xy)
    Label(text=u_xy_s).place(relx=.5, rely=18*RELY, anchor='center')

    plot3d(u_xy)


def main() -> None:
    '''Build the GUI and select the problem'''

    root = Tk()
    root.title('Solution of the Dirichlet-Neumann problem')
    root.geometry('600x400')


    info = Label(text='Use  ;  to separate input.')
    info.place(relx=.5, rely=RELY, anchor='center')


    Label(text='Problem type:').place(relx=.5, rely=2*RELY, anchor='center')
    task = StringVar(value=RING)

    Radiobutton(text=RING, value=RING, variable=task) \
        .place(relx=.25, rely=3*RELY, anchor='center')
    Radiobutton(text=INT, value=INT, variable=task) \
        .place(relx=.5, rely=3*RELY, anchor='center')
    Radiobutton(text=EXT, value=EXT, variable=task) \
        .place(relx=.75, rely=3*RELY, anchor='center')


    Label(text='Circle parameters are given:') \
        .place(relx=.5, rely=4*RELY, anchor='center')
    expl = StringVar(value=EXPL)

    Radiobutton(text=EXPL, value=EXPL, variable=expl) \
        .place(relx=.33, rely=5*RELY, anchor='center')
    Radiobutton(text=IMPL, value=IMPL, variable=expl) \
        .place(relx=.66, rely=5*RELY, anchor='center')


    Label(text='Circle center:').place(relx=.5, rely=6*RELY, anchor='center')
    cent = Entry()
    cent.place(relx=.5, rely=7*RELY, anchor='center')
    cent.insert(0, '0 ; 0')


    Label(text='Radius:').place(relx=.5, rely=8*RELY, anchor='center')
    rad = Entry()
    rad.place(relx=.5, rely=9*RELY, anchor='center')
    rad.insert(0, '1 ; 2')


    Label(text='Circle equation:') \
        .place(relx=.5, rely=10*RELY, anchor='center')
    circ = Entry()
    circ.place(relx=.5, rely=11*RELY, anchor='center', relwidth=.4)
    circ.insert(0, 'x^2+a*x+y^2+b*y=c ; x^2+a*x+y^2+b*y=c')


    Label(text='Function on boundary:') \
        .place(relx=.5, rely=12*RELY, anchor='center')
    fun = Entry()
    fun.place(relx=.5, rely=13*RELY, anchor='center', relwidth=.4)
    fun.insert(0, 'x^2 ; 1/32*x^2+1/8')


    Label(text='Non-homogeneous:') \
        .place(relx=.5, rely=14*RELY, anchor='center')
    inh = Entry()
    inh.place(relx=.5, rely=15*RELY, anchor='center', relwidth=.4)
    inh.insert(0, '0')


    solve_btn = Button(text='Solve problem',
                       command=partial(informer, info, task, expl, cent,
                                       rad, circ, fun, inh))
    solve_btn.place(relx=.5, rely=16*RELY, anchor='center')


    root.mainloop()


if __name__ == '__main__':
    main()
