import sympy as sym


def circle_input() -> tuple:
    '''Обработка ввода про окружность'''

    explicitly = None
    while explicitly is None:
        print('Параметры окружности заданы:',
              '1. Явно',
              '2. Неявно', sep='\n')
        explicitly = int(input())
        if not 0 < explicitly < 3:
            explicitly = None

    if explicitly == 1:
        print('Введите центр окружности через запятую или пробел:')
        center = tuple(map(int, (input().split(','))))

        print('Введите радиус окружности:')
        radius = int(input())
    else:
        print('Введите уравнение окружности в виде',
              'x^2 + ax + y^2 + by + c = 0:')

        string = input().replace(' ', '')
        a = b = c = ''

        i = 3 # x^2 пропускаем
        while string[i] != 'x':
            print('a', i, string[i])
            a += string[i]
            i += 1
        a = float(a) / 2

        i += 5 # x+y^2 пропускаем
        while string[i] != 'y':
            b += string[i]
            i += 1
        b = float(b) / 2

        i += 1 # y пропускаем
        while string[i] != '=':
            c += string[i]
            i += 1
        c = float(c)

        center = [a, b]
        radius = sym.sqrt(a**2 + b**2 - c)

    return center, radius

def function_input() -> None:
    '''Обработка ввода про функцию'''

    return None
