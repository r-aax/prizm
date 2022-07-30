def det(a1, b1, a2, b2):
    return a1 * b2 - a2 * b1
def solve(a1, b1, c1, a2, b2, c2):
    d = det(a1, b1, a2, b2)
    return [det(c1, b1, c2, b2) / d, det(a1, c1, a2, c2) / d]

if __name__ == '__main__':
    S = [12., 12., 12., 12., 12., 12.]
    m = [7., 5., 2., 0., 0., 0.]
    ii = range(1, len(m) + 1)
    print([i for i in ii])
    # R, A
    a1 = -len(m)
    b1 = sum([1. / i for i in ii])
    c1 = -sum([((S[i - 1] - m[i - 1]) / S[i - 1]) for i in ii])
    a2 = -sum([(1. / i) for i in ii])
    b2 = sum([1. / (i * i) for i in ii])
    c2 = -sum([((S[i - 1] - m[i - 1]) / (S[i - 1] * i)) for i in ii])
    [P, A] = solve(a1, b1, c1, a2, b2, c2)
    print(P - (A / len(m)))
