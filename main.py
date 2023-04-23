'''Документация к модулю'''
# from tkinter import ttk, Tk, Label
from sympy import init_printing
from sympy.plotting.plot import plot3d

from handler import handler
from solver import solver

init_printing()


def main() -> None:
    '''Основная функция программы, выбор задачи'''

    task = None
    while task is None:
        print('Выберите тип задачи:',
              '1. Кольцо',
              '2. Круг, внутренняя задача',
              '3. Круг, внешняя задача', sep='\n')
        task = int(input())
        if not 0 < task < 4:
            task = None

    condition = handler(task)
    u = solver(task, *condition)


    # root = Tk()
    # root.title('Решение задачи Дирихле и Неймана')
    # root.geometry('600x400')

    # label = Label(text='Выберите тип задачи:')
    # label.pack()


    # btn = ttk.Button(text='Click')
    # btn.pack()



    plot3d(u)


    # root.mainloop()


if __name__ == '__main__':
    main()
