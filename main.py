'''There's a lot of spaghetti code. Sorry.'''

from functools import partial
from tkinter import Tk, Label, StringVar
from tkinter.ttk import Radiobutton, Entry, Button
from sympy.plotting.plot import plot3d

from stuff import runner


RING = 'Ring'
INT = 'Interior'
EXT = 'Exterior'

EXPL = 'Explicit'
IMPL = 'Implicit'


def informer(task: StringVar, expl: StringVar, cent: Entry, rad: Entry,
             circ: Entry, fun: Entry, inh: Entry) -> None:
    '''Collect the information, transfer it to the solver and build a graph'''

    task_s = task.get()
    expl_s = expl.get()
    cent_s = cent.get().replace(' ', '')
    rad_s = rad.get().replace(' ', '')
    circ_s = circ.get().replace(' ', '')
    fun_s = fun.get().replace(' ', '')
    inh_s = inh.get().replace(' ', '')

    u_rp, u_xy = runner(task_s, expl_s, cent_s, rad_s, circ_s, fun_s, inh_s)

    Label(text=str(u_rp)).grid(row=16, column=1)
    Label(text=str(u_xy)).grid(row=17, column=1)

    plot3d(u_xy)


def main() -> None:
    '''Build the GUI and select the problem'''

    root = Tk()
    root.title('Solution of the Dirichlet-Neumann problem')
    root.geometry('280x350')


    Label(text='Problem type:').grid(row=1, column=1)
    task = StringVar(value=RING)

    Radiobutton(text=RING, value=RING, variable=task).grid(row=2, column=0)
    Radiobutton(text=INT, value=INT, variable=task).grid(row=2, column=1)
    Radiobutton(text=EXT, value=EXT, variable=task).grid(row=2, column=2)


    Label(text='Circle parameters are given:').grid(row=3, column=1)
    expl = StringVar(value=EXPL)

    Radiobutton(text=EXPL, value=EXPL, variable=expl).grid(row=4, column=0)
    Radiobutton(text=IMPL, value=IMPL, variable=expl).grid(row=4, column=1)


    Label(text='Circle center:').grid(row=5, column=1)
    cent = Entry()
    cent.grid(row=6, column=1)
    cent.insert(0, '0 ; 0')


    Label(text='Radius:').grid(row=7, column=1)
    rad = Entry()
    rad.grid(row=8, column=1)
    rad.insert(0, '1 ; 2')


    Label(text='Circle equation:').grid(row=9, column=1)
    circ = Entry()
    circ.grid(row=10, column=1)
    circ.insert(0, 'x^2+a*x+y^2+b*y=c ; x^2+a*x+y^2+b*y=c')


    Label(text='Function on boundary:').grid(row=11, column=1)
    fun = Entry()
    fun.grid(row=12, column=1)
    fun.insert(0, 'x^2 ; 1/32*x^2+1/8')


    Label(text='Non-homogeneous:').grid(row=13, column=1)
    inh = Entry()
    inh.grid(row=14, column=1)
    inh.insert(0, '0')


    solve_btn = Button(text='Solve problem',
                       command=partial(informer, task, expl, cent,
                                       rad, circ, fun, inh))
    solve_btn.grid(row=15, column=1)


    root.mainloop()


if __name__ == '__main__':
    main()
