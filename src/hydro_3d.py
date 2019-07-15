# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:26:02 2019

@author: Rybakov
"""

import draw
import aggdraw
import mth
import math
import vis

#---------------------------------------------------------------------------------------------------
# Constants.
#---------------------------------------------------------------------------------------------------

Gamma = 1.4

#---------------------------------------------------------------------------------------------------
# Utilitiees.
#---------------------------------------------------------------------------------------------------

def array_3d(x, y, z):
    arr = [None] * x
    for i in range(x):
        arr[i] = [None] * y
        for j in range(y):
            arr[i][j] = [None] * z
    return arr

#---------------------------------------------------------------------------------------------------

def convert_data_d_to_data_u(cell):
    cell.U = cell.D.CreateDataU()

#---------------------------------------------------------------------------------------------------

def convert_data_u_to_data_d(cell):
    cell.D = cell.U.CreateDataD()

#---------------------------------------------------------------------------------------------------

def l1(u, a):
    return u - a

#---------------------------------------------------------------------------------------------------
def l2(u):
    return u

#---------------------------------------------------------------------------------------------------

def l5(u, a):
    return u + a

#---------------------------------------------------------------------------------------------------

def lp(l):
    return 0.5 * (l + abs(l))

#---------------------------------------------------------------------------------------------------

def lm(l):
    return 0.5 * (l - abs(l))

#---------------------------------------------------------------------------------------------------
# Class DataD.
#---------------------------------------------------------------------------------------------------

class DataD:

#---------------------------------------------------------------------------------------------------

    def __init__(self,
                 r = 0.0,
                 u = 0.0, v = 0.0, w = 0.0,
                 p = 0.0):
        self.r = r
        self.u = u
        self.v = v
        self.w = w
        self.p = p

#---------------------------------------------------------------------------------------------------

    def __add__(self, d):
        return DataD(self.r + d.r, self.u + d.u, self.v + d.v, self.w + d.w, self.p + d.p)

#---------------------------------------------------------------------------------------------------

    def __mul__(self, k):
        return DataD(self.r * k, self.u * k, self.v * k, self.w * k, self.p * k)

#---------------------------------------------------------------------------------------------------

    def __rmul__(self, k):
        return self * k

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        return '%f / %f / %f / %f / %f' % (self.r, self.u, self.v, self.w, self.p)

#---------------------------------------------------------------------------------------------------

    def V2(self):
        return self.u * self.u + self.v * self.v + self.w * self.w

#---------------------------------------------------------------------------------------------------

    def e(self):
        return self.p / ((Gamma - 1.0) * self.r)

#---------------------------------------------------------------------------------------------------

    def E(self):
        return self.r * (0.5 * self.V2() + self.e())

#---------------------------------------------------------------------------------------------------

    def H(self):
        return (self.E() + self.p) / self.r

#---------------------------------------------------------------------------------------------------

    def a(self):
        return math.sqrt((Gamma - 1.0) * (self.H() - 0.5 * self.V2()))

#---------------------------------------------------------------------------------------------------

    def CreateDataU(self):
        return DataU(self.r,
                     self.r * self.u,
                     self.r * self.v,
                     self.r * self.w,
                     self.E())

#---------------------------------------------------------------------------------------------------

    def CreateFlowF(self):
        return DataU(self.r * self.u,
                     self.r * self.u * self.u + self.p,
                     self.r * self.u * self.v,
                     self.r * self.u * self.w,
                     self.u * (self.E() + self.p))

#---------------------------------------------------------------------------------------------------

    def CreateFlowG(self):
        return DataU(self.r * self.v,
                     self.r * self.u * self.v,
                     self.r * self.v * self.v + self.p,
                     self.r * self.v * self.w,
                     self.v * (self.E() + self.p))

#---------------------------------------------------------------------------------------------------

    def CreateFlowH(self):
        return DataU(self.r * self.w,
                     self.r * self.u * self.w,
                     self.r * self.v * self.w,
                     self.r * self.w * self.w + self.p,
                     self.w * (self.E() + self.p))

#---------------------------------------------------------------------------------------------------

    def CreateFlowF_StegerWarming(self, l1, l2, l5):
        a = self.a()
        g = Gamma
        g1 = g - 1.0
        f = DataU(l1 + 2.0 * g1 * l2 + l5,
                  (self.u - a) * l1 + 2.0 * g1 * self.u * l2 + (self.u + a) * l5,
                  self.v * l1 + 2.0 * g1 * self.v * l2 + self.v * l5,
                  self.w * l1 + 2.0 * g1 * self.w * l2 + self.w * l5,
                  (self.H() - self.u * a) * l1 + g1 * self.V2() * l2 + (self.H() + self.u * a) * l5)
        return (self.r / (2.0 * g)) * f

#---------------------------------------------------------------------------------------------------

    def CreateFlowG_StegerWarming(self, l1, l2, l5):
        a = self.a()
        g = Gamma
        g1 = g - 1.0
        f = DataU(l1 + 2.0 * g1 * l2 + l5,
                  self.u * l1 + 2.0 * g1 * self.u * l2 + self.v * l5,
                  (self.v - a) * l1 + 2.0 * g1 * self.v * l2 + (self.v + a) * l5,
                  self.w * l1 + 2.0 * g1 * self.w * l2 + self.w * l5,
                  (self.H() - self.v * a) * l1 + g1 * self.V2() * l2 + (self.H() + self.v * a) * l5)
        return (self.r / (2.0 * g)) * f

#---------------------------------------------------------------------------------------------------

    def CreateFlowH_StegerWarming(self, l1, l2, l5):
        a = self.a()
        g = Gamma
        g1 = g - 1.0
        f = DataU(l1 + 2.0 * g1 * l2 + l5,
                  self.u * l1 + 2.0 * g1 * self.u * l2 + self.u * l5,
                  self.v * l1 + 2.0 * g1 * self.v * l2 + self.v * l5,
                  (self.w - a) * l1 + 2.0 * g1 * self.w * l2 + (self.w + a) * l5,
                  (self.H() - self.w * a) * l1 + g1 * self.V2() * l2 + (self.H() + self.w * a) * l5)
        return (self.r / (2.0 * g)) * f

#---------------------------------------------------------------------------------------------------
# Class DataU.
#---------------------------------------------------------------------------------------------------

class DataU:

#---------------------------------------------------------------------------------------------------

    def __init__(self,
                 r = 0.0,
                 ru = 0.0, rv = 0.0, rw = 0.0,
                 E = 0.0):
        self.r = r
        self.ru = ru
        self.rv = rv
        self.rw = rw
        self.E = E

#---------------------------------------------------------------------------------------------------

    def __add__(self, u):
        return DataU(self.r + u.r, self.ru + u.ru, self.rv + u.rv, self.rw + u.rw, self.E + u.E)

#---------------------------------------------------------------------------------------------------

    def __sub__(self, u):
        return DataU(self.r - u.r, self.ru - u.ru, self.rv - u.rv, self.rw - u.rw, self.E - u.E)

#---------------------------------------------------------------------------------------------------

    def __mul__(self, k):
        return DataU(self.r * k, self.ru * k, self.rv * k, self.rw * k, self.E * k)

#---------------------------------------------------------------------------------------------------

    def __rmul__(self, k):
        return self * k

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        return '%f / %f / %f / %f / %f' % (self.r, self.ru, self.rv, self.rw, self.E)

#---------------------------------------------------------------------------------------------------

    def CreateDataD(self):
        data_d = DataD(self.r,
                       self.ru / self.r,
                       self.rv / self.r,
                       self.rw / self.r,
                       0.0)
        data_d.p = (self.E / self.r - 0.5 * data_d.V2()) * ((Gamma - 1.0) * self.r)
        return data_d

#---------------------------------------------------------------------------------------------------
# Class Cell.
#---------------------------------------------------------------------------------------------------

class Cell:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        self.D = DataD()
        self.U = DataU()

#---------------------------------------------------------------------------------------------------
# Class Face.
#---------------------------------------------------------------------------------------------------

class Face:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        F = None
        G = None
        H = None

#---------------------------------------------------------------------------------------------------
# Class Grid.
#---------------------------------------------------------------------------------------------------

class Grid:

#---------------------------------------------------------------------------------------------------

    def __init__(self,
                 size_x, size_y, size_z,
                 cells_x, cells_y, cells_z):

        # Sizes.
        self.SizeX, self.SizeY, self.SizeZ = size_x, size_y, size_z
        self.CellsX, self.CellsY, self.CellsZ = cells_x, cells_y, cells_z
        self.dx, self.dy, self.dz = size_x / cells_x, size_y / cells_y, size_z / cells_z

        # Cells.
        self.Cells = array_3d(cells_x, cells_y, cells_z)
        for i in range(cells_x):
            for j in range(cells_y):
                for k in range(cells_z):
                    self.Cells[i][j][k] = Cell()

        # Faces.
        self.FacesX = array_3d(cells_x + 1, cells_y, cells_z)
        for i in range(cells_x + 1):
            for j in range(cells_y):
                for k in range(cells_z):
                    self.FacesX[i][j][k] = Face()
        self.FacesY = array_3d(cells_x, cells_y + 1, cells_z)
        for i in range(cells_x):
            for j in range(cells_y + 1):
                for k in range(cells_z):
                    self.FacesY[i][j][k] = Face()
        self.FacesZ = array_3d(cells_x, cells_y, cells_z + 1)
        for i in range(cells_x):
            for j in range(cells_y):
                for k in range(cells_z + 1):
                    self.FacesZ[i][j][k] = Face()

#---------------------------------------------------------------------------------------------------

    def CellsCount(self):
        return self.CellsX * self.CellsY * self.CellsZ

#---------------------------------------------------------------------------------------------------

    def FacesCount(self):
        cx, cy, cz = self.CellsX, self.CellsY, self.CellsZ
        return (cx + 1) * cy * cz + cx * (cy + 1) * cz + cx * cy * (cz + 1)

#---------------------------------------------------------------------------------------------------

    def Apply(self, fun):
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                for k in range(self.CellsZ):
                    fun(self.Cells[i][j][k])

#---------------------------------------------------------------------------------------------------

    def Init(self, case):
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                for k in range(self.CellsZ):
                    cell = self.Cells[i][j][k]
                    x = (i + 0.5) * self.dx
                    y = (j + 0.5) * self.dy

                    if case == 1:
                        # case 1
                        if x < 0.5:
                            cell.D = DataD(10.0, 0.0, 0.0, 0.0, 10.0)
                        else:
                            cell.D = DataD(1.0, 0.0, 0.0, 0.0, 1.0)
                    elif case == 2:
                        # case 1
                        if y < 0.5:
                            cell.D = DataD(10.0, 0.0, 0.0, 0.0, 10.0)
                        else:
                            cell.D = DataD(1.0, 0.0, 0.0, 0.0, 1.0)
                    elif case == 3:
                        # case 3
                        if x < 0.5:
                            cell.D = DataD(1.0, 0.75, 0.0, 0.0, 1.0)
                        else:
                            cell.D = DataD(0.125, 0.0, 0.0, 0.0, 0.1)
                    elif case == 4:
                        # case 4
                        if x * x + y * y < 0.64:
                            cell.D = DataD(10.0, 0.0, 0.0, 0.0, 10.0)
                        else:
                            cell.D = DataD(1.0, 0.0, 0.0, 0.0, 1.0)
                    else:
                        raise Exception('unknown case number')

#---------------------------------------------------------------------------------------------------

    def DtoU(self):
        self.Apply(convert_data_d_to_data_u)

#---------------------------------------------------------------------------------------------------

    def UtoD(self):
        self.Apply(convert_data_u_to_data_d)

#---------------------------------------------------------------------------------------------------

    def CalcFlowsTrivial(self):

        cs = self.Cells

        # X
        for i in range(self.CellsX + 1):
            for j in range(self.CellsY):
                for k in range(self.CellsZ):
                    if i == 0:
                        # L bc
                        d = cs[i][j][k].D
                    elif i == self.CellsX:
                        # R bc
                        d = cs[i - 1][j][k].D
                    else:
                        # LR
                        d = 0.5 * (cs[i - 1][j][k].D + cs[i][j][k].D)
                    self.FacesX[i][j][k].F = d.CreateFlowF()

        # Y
        for i in range(self.CellsX):
            for j in range(self.CellsY + 1):
                for k in range(self.CellsZ):
                    if j == 0:
                        # D bc
                        d = cs[i][j][k].D
                    elif j == self.CellsY:
                        # U bc
                        d = cs[i][j - 1][k].D
                    else:
                        # DU
                        d = 0.5 * (cs[i][j - 1][k].D + cs[i][j][k].D)
                    self.FacesY[i][j][k].G = d.CreateFlowG()

        # Z
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                for k in range(self.CellsZ + 1):
                    if k == 0:
                        # B bc
                        d = cs[i][j][k].D
                    elif k == self.CellsZ:
                        # F bc
                        d = cs[i][j][k - 1].D
                    else:
                        # BF
                        d = 0.5 * (cs[i][j][k - 1].D + cs[i][j][k].D)
                    self.FacesZ[i][j][k].H = d.CreateFlowH()

#---------------------------------------------------------------------------------------------------

    def CalcFlowsStegerWarming(self):

        cs = self.Cells

        # X
        for i in range(self.CellsX + 1):
            for j in range(self.CellsY):
                for k in range(self.CellsZ):
                    if i == 0:
                        # L bc
                        d = cs[i][j][k].D
                        self.FacesX[i][j][k].F = d.CreateFlowF()
                    elif i == self.CellsX:
                        # R bc
                        d = cs[i - 1][j][k].D
                        self.FacesX[i][j][k].F = d.CreateFlowF()
                    else:
                        # LR
                        dp = cs[i - 1][j][k].D
                        dm = cs[i][j][k].D
                        ap = dp.a()
                        am = dm.a()
                        lp1, lp2, lp5 = lp(l1(dp.u, ap)), lp(l2(dp.u)), lp(l5(dp.u, ap))
                        lm1, lm2, lm5 = lm(l1(dm.u, am)), lm(l2(dm.u)), lm(l5(dm.u, am))
                        fp = dp.CreateFlowF_StegerWarming(lp1, lp2, lp5)
                        fm = dm.CreateFlowF_StegerWarming(lm1, lm2, lm5)
                        self.FacesX[i][j][k].F = fp + fm

        # Y
        for i in range(self.CellsX):
            for j in range(self.CellsY + 1):
                for k in range(self.CellsZ):
                    if j == 0:
                        # D bc
                        d = cs[i][j][k].D
                        self.FacesY[i][j][k].G = d.CreateFlowG()
                    elif j == self.CellsY:
                        # U bc
                        d = cs[i][j - 1][k].D
                        self.FacesY[i][j][k].G = d.CreateFlowG()
                    else:
                        # DU
                        dp = cs[i][j - 1][k].D
                        dm = cs[i][j][k].D
                        ap = dp.a()
                        am = dm.a()
                        lp1, lp2, lp5 = lp(l1(dp.u, ap)), lp(l2(dp.u)), lp(l5(dp.u, ap))
                        lm1, lm2, lm5 = lm(l1(dm.u, am)), lm(l2(dm.u)), lm(l5(dm.u, am))
                        fp = dp.CreateFlowG_StegerWarming(lp1, lp2, lp5)
                        fm = dm.CreateFlowG_StegerWarming(lm1, lm2, lm5)
                        self.FacesY[i][j][k].G = fp + fm

        # Z
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                for k in range(self.CellsZ + 1):
                    if k == 0:
                        # B bc
                        d = cs[i][j][k].D
                        self.FacesZ[i][j][k].H = d.CreateFlowH()
                    elif k == self.CellsZ:
                        # F bc
                        d = cs[i][j][k - 1].D
                        self.FacesZ[i][j][k].H = d.CreateFlowH()
                    else:
                        # BF
                        dp = cs[i][j][k - 1].D
                        dm = cs[i][j][k].D
                        ap = dp.a()
                        am = dm.a()
                        lp1, lp2, lp5 = lp(l1(dp.u, ap)), lp(l2(dp.u)), lp(l5(dp.u, ap))
                        lm1, lm2, lm5 = lm(l1(dm.u, am)), lm(l2(dm.u)), lm(l5(dm.u, am))
                        fp = dp.CreateFlowH_StegerWarming(lp1, lp2, lp5)
                        fm = dm.CreateFlowH_StegerWarming(lm1, lm2, lm5)
                        self.FacesZ[i][j][k].H = fp + fm

#---------------------------------------------------------------------------------------------------

    def Renew(self, dt):
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                for k in range(self.CellsZ):
                    c = self.Cells[i][j][k]
                    u = c.U
                    nu = u \
                         - (dt / self.dx) * (self.FacesX[i + 1][j][k].F - self.FacesX[i][j][k].F) \
                         - (dt / self.dy) * (self.FacesY[i][j + 1][k].G - self.FacesY[i][j][k].G) \
                         - (dt / self.dz) * (self.FacesZ[i][j][k + 1].H - self.FacesZ[i][j][k].H)
                    c.U = nu

#---------------------------------------------------------------------------------------------------

    def Step(self, dt):
        self.DtoU()
        self.CalcFlowsStegerWarming()
        self.Renew(dt)
        self.UtoD()

#---------------------------------------------------------------------------------------------------

    def Steps(self, n, dt):
        for _ in range(n):
            self.Step(dt)

#---------------------------------------------------------------------------------------------------


    def Print(self):
        print('Grid : %d cells, %d faces' % (self.CellsCount(), self.FacesCount()))

#---------------------------------------------------------------------------------------------------

    def FunInterval(self, fun):
        inter = (0.0, -1.0)
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                for k in range(self.CellsZ):
                    v = fun(self.Cells[i][j][k])
                    (mn, mx) = inter
                    if mx < mn:
                        inter = (v, v)
                    else:
                        if v < mn:
                            inter = (v, mx)
                        elif v > mx:
                            inter = (mn, v)
                        else:
                            pass
        return inter

#---------------------------------------------------------------------------------------------------

    def XValues(self, fun):
        return [fun(self.Cells[i][0][0]) for i in range(self.CellsX)]

#---------------------------------------------------------------------------------------------------

    def YValues(self, fun):
        return [fun(self.Cells[0][j][0]) for j in range(self.CellsY)]

#---------------------------------------------------------------------------------------------------

    def Draw(self, fun):
        d = draw.Drawer(draw_area = (0.0, 0.0, self.SizeX, self.SizeY),
                        pic_size = (500, 500))
        (mn, mx) = self.FunInterval(fun)
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                cell = self.Cells[i][j][0]
                factor = 255 - int((fun(cell) - mn) / (mx - mn) * 255)
                color = (factor, factor, factor)
                d.Rect((i * self.dx, j * self.dy), ((i + 1) * self.dx, (j + 1) * self.dy),
                       aggdraw.Pen(color, 1.0), aggdraw.Brush(color))
            d.Rect((0.0, 0.0), (self.SizeX, self.SizeY), aggdraw.Pen('blue', 1.0))
        d.FSS()

#---------------------------------------------------------------------------------------------------
# Main.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    print('HYDRO_3D')
    case = 1
    if case == 1:
        g = Grid(1.0, 1.0, 1.0, 100, 1, 1)
    elif case == 2:
        g = Grid(1.0, 1.0, 1.0, 1, 100, 1)
    elif case == 3:
        g = Grid(1.0, 1.0, 1.0, 100, 1, 1)
    elif case == 4:
        g = Grid(1.0, 1.0, 1.0, 100, 100, 1)
    else:
        raise Exception('unknown case number')
    g.Init(case)

    pics, n, dt = 10, 10, 0.001
    fun = lambda cell: cell.D.p

    for _ in range(pics):
        g.Steps(n, dt)
        #g.Draw(fun)
        vis.simple_graphic_ys(g.YValues(fun))

#---------------------------------------------------------------------------------------------------