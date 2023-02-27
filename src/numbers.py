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

    if m > n:
        return gcd(m, n)
    elif n % m == 0:
        return m
    else:
        return gcd(m, n % m)


if __name__ == '__main__':

    # gcd
    assert gcd(10, 1) == 1
    assert gcd(4, 6) == 2
    assert gcd(5, 25) == 5
