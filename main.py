import sympy as sym
from sympy.plotting.plot import plot3d

from handler import handler
from solver import solver


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
    u_function = solver(task, *condition)

    plot3d(u_function)


if __name__ == '__main__':
    main()
