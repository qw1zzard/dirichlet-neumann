import sympy as sym

from stuff import circle_input, function_input


def handler(task) -> tuple:
    '''Получение условий задачи от пользователя'''

    sym.init_printing()

    if task == 1:
        print('Введите информацию о внутренней и внешней окружностях последовательно:')
        circle = (circle_input(), circle_input())

        print('И функциях на них соответственно:')
        function = (function_input(), function_input())
    else:
        print('Введите информацию об окружности:')
        circle = circle_input()

        print('И функции на ней:')
        function = function_input()

    print('Если есть неоднородность, напишите её, иначе 0:')
    heterogeneity = function_input()

    return circle, function, heterogeneity
