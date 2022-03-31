import time
import math


def main(number):
    choice = "0"
    algorithms = [pollard_rho, elliptic_curve, pollard_p1, rational_sieve, dixon]
    algorithm_names = ["Pollard's Rho", "Lenstra Elliptic Curve", "Pollard's p-1", "Rational Sieve", "Dixon's Method"]
    menu = "\nChoose the algorithm:\n\n"
    for i in range(len(algorithm_names)):
        menu += str(i+1) + ". " + algorithm_names[i] + "\n"

    menu += str(len(algorithm_names)+1) + ". Exit\n\nChoice: "

    while choice != str(len(algorithms)+1):
        choice = input(menu)

        if choice != str(len(algorithms)+1):
            print("\nStarting Factorisation...")
            time_start = time.time()
            f = algorithms[int(choice)-1](number)
            time_end = time.time()
            t = time_end - time_start
            print("\n********************************")
            print("Number Factorised")
            print("Time taken: " + str(t) + " seconds")
            print("Factors: " + str(f) + ", " + str(number // f))
            print("********************************")
            time.sleep(3)


def pollard_rho(n):
    x = 2
    y = 2
    d = 1

    while d == 1:
        x = (x**2 + 1) % n
        y = (((y**2 + 1) % n)**2 + 1) % n
        d = math.gcd(abs(x - y), n)

    return d


def elliptic_curve(n):
    limit = 6
    x = 0
    y = 1
    h = 0
    found = False
    output = []
    a = 1

    while not found:
        h += 1
        if h > limit:
            h = 1
            a += 1
            x = 0
            y = a

        output = multiply(n, h, [x, y], a)

        if len(output) != 2:
            found = True
        else:
            x = output[0]
            y = output[1]

    print(h)
    print(a)
    factor = math.gcd(output[0], n)
    return factor


def double(n, point, a):
    x = point[0]
    y = point[1]

    numerator = (3 * (x**2) + a) % n
    denominator = modular_inverse((2 * y), n)

    if denominator != 0:
        slope = (numerator * denominator) % n
        x2 = (slope**2 - (2 * x)) % n
        y2 = (slope * (x - x2) - y) % n

        return [x2, y2]

    else:
        return [2*y % n, True, 1]


def add(n, point1, point2):
    x1 = point1[0]
    y1 = point1[1]
    x2 = point2[0]
    y2 = point2[1]

    denominator = ((y1 - y2) * modular_inverse(x1 - x2, n)) % n

    if denominator != 0:
        numerator = (y1 - y2)
        slope = numerator * denominator
        x3 = (slope**2 - x1 - x2) % n
        y3 = (slope * (x1 - x3) - y1) % n

        return [x3, y3]

    else:
        return [(x1 - x2) % n, True, 1]


def egcd(a, b):
    if a == 0:
        return b, 0, 1
    else:
        g, y, x = egcd(b % a, a)
        return g, x - (b // a) * y, y


def modular_inverse(a, m):
    a = a % m
    g, x, y = egcd(a, m)
    if g != 1:
        return 0
    else:
        return x % m


def multiply(n, e, point, a):
    flipbin = format(e, 'b')[::-1]
    points = []

    for i in range(len(flipbin)):
        point_temp = point
        value = int(flipbin[i])
        if value == 1:
            for j in range(i):
                point_temp = double(n, point_temp, a)

                if len(point_temp) != 2:
                    return point_temp

            points.append(point_temp)

    new_point = points[0]

    for k in range(len(points) - 1):
        new_point = add(n, new_point, points[k + 1])

        if len(new_point) != 2:
            return new_point

    return new_point


def pollard_p1(n):
    if n % 2 == 0:
        return 2

    a = 2
    k = 1
    limit = 15

    while True:
        k += 1
        if k > limit:
            k = 1
            a += 1
            if n % a == 0:
                gcd = a
                break

        power = math.factorial(k)
        gcd = math.gcd((a**power - 1) % n, n)
        if gcd != 1 and gcd != n:
            break

    return gcd


def rational_sieve(n):
    if n < 3:
        print('Number must be greater than 2')
        return

    smoothness_bound = math.floor(
        math.e ** (1 / 2 * math.sqrt(math.log(n, math.e) * math.log(math.log(n, math.e), math.e))))

    if smoothness_bound > primes[-1]:
        choice = input('WARNING: Insufficient primes supplied for large number.\nContinue? (y/n)')
        if choice == 'n':
            return

    factor_base = primes[:list(filter(lambda k: k > smoothness_bound, primes))[0]-1]

    count = 0
    z = 1
    smooth_number_pool = []
    factorisation_pool = []

    while count <= smoothness_bound+5 and z < n:
        if smooth(z, factor_base):
            if smooth(z+n, factor_base):
                smooth_number_pool.append(z)
                count += 1

        z += 1

    for num in smooth_number_pool:
        factorisation_pool.append([factorise(num, factor_base), factorise(num+n, factor_base)])

    print(factorisation_pool)

    return


def smooth(n, factor_base):
    k = math.prod(factor_base)
    g = 0

    while g != 1:
        g = math.gcd(n, k)
        n = n // g
        if n == 1:
            return True

    return False


def factorise(n, factor_base):
    factorisation = []
    
    for factor in factor_base:
        count = 0
        t = n
        while t % factor == 0:
            count += 1
            t = t // factor

        factorisation.append(count % 2)

    return factorisation


def dixon(n):
    if n < 3:
        print('Number must be greater than 2')
        return

    smoothness_bound = math.floor(
        math.e ** (1 / 2 * math.sqrt(math.log(n, math.e) * math.log(math.log(n, math.e), math.e))))

    if smoothness_bound > primes[-1]:
        choice = input('WARNING: Insufficient primes supplied for large number.\nContinue? (y/n)')
        if choice == 'n':
            return

    factor_base = primes[:list(filter(lambda k: k > smoothness_bound, primes))[0]-1]

    for f in factor_base:
        if n % f == 0:
            return f

    count = 0
    z = math.ceil(math.sqrt(n))
    factorisation_pool = []
    zero_vector = [0 for _ in range(len(factor_base))]

    while count <= smoothness_bound+1 and z < n:
        if smooth(z**2 % n, factor_base):
            r = factorise(z**2 % n, factor_base)
            if r == zero_vector:
                x = int(math.sqrt(z**2 % n))
                if (z - x) % n != 0:
                    a = math.gcd(z - x, n)
                    b = math.gcd(z + x, n)
                    if a == 1:
                        if b != 1:
                            return b
                    else:
                        return a
            else:
                count += 1
                factorisation_pool.append([r, z])

        z += 1

    matrix = factorisation_pool
    multipliers = []

    while matrix != []:
        next_row = get_top(matrix)
        matrix, used = eliminate_row(matrix, next_row)
        s = len(matrix)
        k = 0
        if used:
            multipliers.append(next_row[1])
        while k < s:
            vector = matrix[k]
            if vector[0] == zero_vector:
                try:
                    product_y = math.sqrt(vector[1]**2 % n)
                    for multiplier in multipliers:
                        product_y = product_y * math.sqrt(multiplier**2 % n)

                    x = int(math.prod(multipliers))
                    y = int(product_y)
                    a = math.gcd(x+y, n)
                    b = math.gcd(x-y, n)
                    if a == 1 or a == n:
                        del matrix[k]
                        s -= 1
                    else:
                        return a
                except:
                    pass
            else:
                k += 1
    print('Failed')
    return 1


def add_rows(a, b):
    new_row = [[], a[1]]
    for i in range(len(a[0])):
        new_row[0].append((a[0][i] + b[0][i]) % 2)

    return new_row


def eliminate_row(matrix, row):
    flag = False
    matrix.remove(row)
    l = next((i for i, x in enumerate(row[0]) if x), None)
    for j in range(len(matrix)):
        if matrix[j][0][l]:
            matrix[j] = add_rows(matrix[j], row)
            flag = True

    return matrix, flag


def get_top(matrix):
    for i in range(len(matrix[0][0])):
        for j in range(len(matrix)):
            if matrix[j][0][i]:
                return matrix[j]


if __name__ == '__main__':
    number1 = input("Enter the number: ")
    #number1 = str(1073676287*16769023)

    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
              107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
              227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
              349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
              467, 479, 487, 491, 499, 503, 509, 521, 523, 541]

    while (not number1.isdigit()) or int(number1) < 2:
        print("Number must be an integer greater than 1")
        number1 = input("Enter the number: ")

    number1 = int(number1)

    main(number1)
