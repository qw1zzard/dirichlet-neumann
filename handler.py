'''Документация к модулю'''
from sympy import symbols, sin, cos, simplify

from stuff import circle_input, function_input


def handler(task) -> tuple:
    '''Получение условий задачи от пользователя'''

    x, y, rho, phi = symbols('x, y, rho, phi')

    if task == 1:
        print('Введите информацию о внутренней и внешней ' +
              'окружностях последовательно:')
        circle = (circle_input(), circle_input())

        (x_0, y_0), r_in = circle[0]
        (x_0, y_0), r_out = circle[1]

        print('И функциях на них соответственно:')
        function_in, function_out = (function_input(), function_input())
        function_in = simplify(function_in.subs([(x, x_0 + rho * cos(phi)),
                                                 (y, y_0 + rho * sin(phi))]))
        function_out = simplify(function_out.subs([(x, x_0 + rho * cos(phi)),
                                                   (y, y_0 + rho * sin(phi))]))

        radius = (r_in, r_out)
        function = (function_in, function_out)
    else:
        print('Введите информацию об окружности:')
        (x_0, y_0), radius = circle_input()

        print('И функции на ней:')
        function = function_input()
        function = simplify(function.subs([(x, x_0 + rho * cos(phi)),
                                           (y, y_0 + rho * sin(phi))]))

    print('Если есть функция неоднородности, напишите её, иначе Enter:')
    heterogeneity = function_input() or 0
    heterogeneity = simplify(heterogeneity.subs([(x, x_0 + rho * cos(phi)),
                                                 (y, y_0 + rho * sin(phi))]))

    return radius, function, heterogeneity
