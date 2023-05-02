'''There's a lot of spaghetti code. Sorry.'''

from functools import partial
from tkinter import Tk, Label, StringVar
from tkinter.ttk import Radiobutton, Entry, Button

from stuff import runner


RING = 'Ring'
INT = 'Interior'
EXT = 'Exterior'

EXPL = 'Explicit'
IMPL = 'Implicit'


def main() -> None:
    '''Build the GUI and select the problem'''

    root = Tk()
    root.title('Solution of the Dirichlet-Neumann problem')
    root.geometry('300x350')


    task_lbl = Label(text='Problem type:')
    task_lbl.grid(row=0, column=1)

    task = StringVar(value=RING)

    ring_rdbtn = Radiobutton(text=RING, value=RING, variable=task)
    ring_rdbtn.grid(row=1, column=0)
    int_rdbtn = Radiobutton(text=INT, value=INT, variable=task)
    int_rdbtn.grid(row=1, column=1)
    ext_rdbtn = Radiobutton(text=EXT, value=EXT, variable=task)
    ext_rdbtn.grid(row=1, column=2)


    expl_lbl = Label(text='Circle parameters are given:')
    expl_lbl.grid(row=2, column=1)

    expl = StringVar(value=EXPL)

    expl_yes_rdbtn = Radiobutton(text=EXPL, value=EXPL, variable=expl)
    expl_yes_rdbtn.grid(row=3, column=0)
    expl_no_rdbtn = Radiobutton(text=IMPL, value=IMPL, variable=expl)
    expl_no_rdbtn.grid(row=3, column=1)


    cent_lbl = Label(text='Circle center:')
    cent_lbl.grid(row=4, column=1)

    cent = Entry()
    cent.grid(row=5, column=1)


    rad_lbl = Label(text='Radius:')
    rad_lbl.grid(row=6, column=1)

    rad = Entry()
    rad.grid(row=7, column=1)


    circ_lbl = Label(text='Circle equation:')
    circ_lbl.grid(row=8, column=1)

    circ = Entry()
    circ.grid(row=9, column=1)


    fun_lbl = Label(text='Function on boundary:')
    fun_lbl.grid(row=10, column=1)

    fun = Entry()
    fun.grid(row=11, column=1)


    inh_lbl = Label(text='Non-homogeneous:')
    inh_lbl.grid(row=12, column=1)

    inh = Entry()
    inh.grid(row=13, column=1)


    solve_btn = Button(text='Solve problem',
                           command=partial(runner, task.get(), expl.get(),
                                           cent, rad, circ, fun, inh))
    solve_btn.grid(row=14, column=1)


    root.mainloop()


if __name__ == '__main__':
    main()
