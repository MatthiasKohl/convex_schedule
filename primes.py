import numpy as np
import timeit

def prime_factor_table(n):
    if n < 1: return np.zeros(1)
    table = np.zeros(n+1) # index table by the number itself
    table[1] = 1
    for i in range(2, n+1):
        # trial division to fill the table
        factor = 2
        while factor <= i:
            if (i % factor == 0):
                table[i] = factor
                break
            factor += 1

    return table

def prime_factors(table, n):
    def factors(table, n, l):
        if (n <= 1): return l
        factor = table[n]
        l.append(factor)
        return factors(table, n // factor, l)
    return factors(table, n, [])


if __name__ == '__main__':
    import timeit
    for n in map(lambda x: x*1000, range(1, 11)):
        print(str(n) + ': ' + str(timeit.timeit('prime_factor_table(n)',
                                            setup='from __main__ import prime_factor_table\n' +
                                            'n=' + str(n),
                                            number=10)))

        print(str(n) + ': ' + str(timeit.timeit('prime_factors(t, n)',
                                            setup='import random\n' +
                                            'from __main__ import prime_factor_table, prime_factors\n' +
                                            't=prime_factor_table(' + str(n) + ')\n' +
                                            'n=random.randint(2, ' + str(n) +')',
                                            number=1000)))
