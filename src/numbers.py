"""
Some functions for numbers theory.
"""


def gcd(n, m):
    """
    Greatest common divisor.

    Parameters
    ----------
    n : int
        First number.
    m : int
        Second number.

    Returns
    -------
    int
        Greatest common divisor.
    """

    assert (n > 0) and (m > 0)

    if m > n:
        return gcd(m, n)
    elif n % m == 0:
        return m
    else:
        return gcd(m, n % m)


def factorization(n):
    """
    Get factorization of a number.

    Parameters
    ----------
    n : int
        Number.

    Returns
    -------
    [int]
        List of factors.
    """

    assert n > 0

    factors = []
    cn = n
    d = 2

    while True:
        if (cn == d) or (d * d > cn):
            factors.append(cn)
            break
        if cn % d == 0:
            factors.append(d)
            cn //= d
        else:
            d += 1

    return factors


if __name__ == '__main__':

    # gcd
    assert gcd(10, 1) == 1
    assert gcd(4, 6) == 2
    assert gcd(5, 25) == 5

    # factorization
    assert factorization(1) == [1]
    assert factorization(37) == [37]
    assert factorization(10) == [2, 5]
    assert factorization(27) == [3, 3, 3]
    assert factorization(42) == [2, 3, 7]
