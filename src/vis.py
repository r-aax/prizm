# -*- coding: utf-8 -*-
"""
Visualization functions.

Created on Tue Mar 26 15:34:08 2019

@author: Rybakov
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi

#---------------------------------------------------------------------------------------------------
# Graphics construction.
#---------------------------------------------------------------------------------------------------

def simple_graphic(xs, ys, title):
    """
    Draw simple graphic.

    Arguments:
        xs -- array of X coordinates,
        ys -- array of Y coordinates,
        title -- title.
    """

    plt.figure(num=1, figsize=(10, 6))
    plt.title(title, size=14)
    plt.xlabel('t (градусы Цельсия)', size=14)
    plt.ylabel('L_ev (МДж/кг)', size=14)

    # Draw.
    plt.plot(xs, [2.5] * len(xs), color='b', linestyle='-', label='постоянное значение')
    f = spi.interp1d(xs, ys, kind='nearest')
    plt.plot(range(0, 374), f(range(0, 374)), color='y', linestyle='-', label='кусочно-постоянная интерполяция')
    plt.plot(xs, ys, color='r', linestyle='-', marker='o', label='кусочно-линейная интерполяция')
    tck = spi.splrep(xs, ys, s=0)
    r = range(0, 375)
    nys = spi.splev(r, tck, der=0)
    plt.plot(r, nys, color='g', linestyle='-', label='использование гладкой интерполяции')

    # Legend.
    plt.legend(loc='lower left')

    # Save figure.
    plt.savefig('simple_graphic.png', format='png')

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    xs = [-50.0, 0.0, 10.0, 20.0, 30.0, 50.0, 70.0, 90.0, 100.0, 120.0, 150.0, 180.0, 200.0, 220.0, 250.0, 300.0, 350.0, 370.0, 374.0, 374.15, 500.0]
    ys = [2.5, 2.5, 2.47, 2.45, 2.4, 2.38, 2.32, 2.28, 2.26, 2.2, 2.11, 2.01, 1.94, 1.86, 1.7, 1.4, 0.89, 0.44, 0.11, 0.0, 0.0]
    simple_graphic(xs, ys, "Удельная теплота парообразования")

#---------------------------------------------------------------------------------------------------
