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
import time
import geom

#---------------------------------------------------------------------------------------------------
# Constants.
#---------------------------------------------------------------------------------------------------

Gamma = 1.4
Sph = geom.Sphere(geom.Vector(0.6, 0.5, 0.0), 0.3)

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
                 r = 1.0,
                 u = 0.0, v = 0.0, w = 0.0,
                 p = 1.0):

        # Basic values.
        self.r = r
        self.u = u
        self.v = v
        self.w = w
        self.p = p

        self.CalculateDerivedVariables()

#---------------------------------------------------------------------------------------------------

    def CalculateDerivedVariables(self):
        self.V2 = self.u * self.u + self.v * self.v + self.w * self.w
        self.e = self.p / ((Gamma - 1.0) * self.r)
        self.E = self.r * (0.5 * self.V2 + self.e)
        self.H = (self.E + self.p) / self.r
        self.a = math.sqrt((Gamma - 1.0) * (self.H - 0.5 * self.V2))

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

    def CreateDataU(self):
        return DataU(self.r,
                     self.r * self.u,
                     self.r * self.v,
                     self.r * self.w,
                     self.E)

#---------------------------------------------------------------------------------------------------

    def CreateFlowF(self):
        return DataU(self.r * self.u,
                     self.r * self.u * self.u + self.p,
                     self.r * self.u * self.v,
                     self.r * self.u * self.w,
                     self.u * (self.E + self.p))

#---------------------------------------------------------------------------------------------------

    def CreateFlowG(self):
        return DataU(self.r * self.v,
                     self.r * self.u * self.v,
                     self.r * self.v * self.v + self.p,
                     self.r * self.v * self.w,
                     self.v * (self.E + self.p))

#---------------------------------------------------------------------------------------------------

    def CreateFlowH(self):
        return DataU(self.r * self.w,
                     self.r * self.u * self.w,
                     self.r * self.v * self.w,
                     self.r * self.w * self.w + self.p,
                     self.w * (self.E + self.p))

#---------------------------------------------------------------------------------------------------

    def CreateFlowF_Zero(self):
        return DataU(0.0, self.p, 0.0, 0.0, 0.0)

#---------------------------------------------------------------------------------------------------

    def CreateFlowG_Zero(self):
        return DataU(0.0, 0.0, self.p, 0.0, 0.0)

#---------------------------------------------------------------------------------------------------

    def CreateFlowH_Zero(self):
        return DataU(0.0, 0.0, 0.0, self.p, 0.0)

#---------------------------------------------------------------------------------------------------

    def CreateFlowF_StegerWarming(self, l1, l2, l5):
        a = self.a
        g = Gamma
        g1 = g - 1.0
        f = DataU(l1 + 2.0 * g1 * l2 + l5,
                  (self.u - a) * l1 + 2.0 * g1 * self.u * l2 + (self.u + a) * l5,
                  self.v * l1 + 2.0 * g1 * self.v * l2 + self.v * l5,
                  self.w * l1 + 2.0 * g1 * self.w * l2 + self.w * l5,
                  (self.H - self.u * a) * l1 + g1 * self.V2 * l2 + (self.H + self.u * a) * l5)
        return (self.r / (2.0 * g)) * f

#---------------------------------------------------------------------------------------------------

    def CreateFlowG_StegerWarming(self, l1, l2, l5):
        a = self.a
        g = Gamma
        g1 = g - 1.0
        f = DataU(l1 + 2.0 * g1 * l2 + l5,
                  self.u * l1 + 2.0 * g1 * self.u * l2 + self.u * l5,
                  (self.v - a) * l1 + 2.0 * g1 * self.v * l2 + (self.v + a) * l5,
                  self.w * l1 + 2.0 * g1 * self.w * l2 + self.w * l5,
                  (self.H - self.v * a) * l1 + g1 * self.V2 * l2 + (self.H + self.v * a) * l5)
        return (self.r / (2.0 * g)) * f

#---------------------------------------------------------------------------------------------------

    def CreateFlowH_StegerWarming(self, l1, l2, l5):
        a = self.a
        g = Gamma
        g1 = g - 1.0
        f = DataU(l1 + 2.0 * g1 * l2 + l5,
                  self.u * l1 + 2.0 * g1 * self.u * l2 + self.u * l5,
                  self.v * l1 + 2.0 * g1 * self.v * l2 + self.v * l5,
                  (self.w - a) * l1 + 2.0 * g1 * self.w * l2 + (self.w + a) * l5,
                  (self.H - self.w * a) * l1 + g1 * self.V2 * l2 + (self.H + self.w * a) * l5)
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
        u, v, w, = self.ru / self.r, self.rv / self.r, self.rw / self.r
        v2 = u * u + v * v + w * w
        p = (self.E / self.r - 0.5 * v2) * ((Gamma - 1.0) * self.r)
        return  DataD(self.r, u, v, w, p)

#---------------------------------------------------------------------------------------------------
# Class Cell.
#---------------------------------------------------------------------------------------------------

BorderL = 0
BorderR = 1
BorderD = 2
BorderU = 3
BorderB = 4
BorderF = 5
BorderNone = 0
BorderSoft = 1
BorderHard = 2
TypeCommon = 0
TypeBorder = 1
TypeGhost = 2
TypePhantom = 3
TypeInner = 4

class Cell:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        self.D = DataD()
        self.U = DataU()
        self.Type = TypeCommon
        self.Borders = [BorderNone] * 6

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

        # Borders.
        for i in range(cells_x):
            for j in range(cells_y):
                self.Cells[i][j][0].Borders[BorderB] = BorderSoft
                self.Cells[i][j][cells_z - 1].Borders[BorderF] = BorderSoft
        for i in range(cells_x):
            for k in range(cells_z):
                self.Cells[i][0][k].Borders[BorderD] = BorderSoft
                self.Cells[i][cells_y - 1][k].Borders[BorderU] = BorderSoft
        for j in range(cells_y):
            for k in range(cells_z):
                self.Cells[0][j][k].Borders[BorderL] = BorderSoft
                self.Cells[cells_x - 1][j][k].Borders[BorderR] = BorderSoft

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

    def DtoU(self):
        self.Apply(convert_data_d_to_data_u)

#---------------------------------------------------------------------------------------------------

    def UtoD(self):
        self.Apply(convert_data_u_to_data_d)

#---------------------------------------------------------------------------------------------------

    def CalcFlowF_StegerWarming(self, face, cp, cm):
        dp = cp.D
        dm = cm.D
        up, um = dp.u, dm.u
        ap, am = dp.a, dm.a
        lp1, lp2, lp5 = lp(up - ap), lp(up), lp(up + ap)
        lm1, lm2, lm5 = lm(um - am), lm(um), lm(um + am)
        fp = dp.CreateFlowF_StegerWarming(lp1, lp2, lp5)
        fm = dm.CreateFlowF_StegerWarming(lm1, lm2, lm5)
        face.F = fp + fm

#---------------------------------------------------------------------------------------------------

    def CalcFlowG_StegerWarming(self, face, cp, cm):
        dp = cp.D
        dm = cm.D
        vp, vm = dp.v, dm.v
        ap, am = dp.a, dm.a
        lp1, lp2, lp5 = lp(vp - ap), lp(vp), lp(vp + ap)
        lm1, lm2, lm5 = lm(vm - am), lm(vm), lm(vm + am)
        fp = dp.CreateFlowG_StegerWarming(lp1, lp2, lp5)
        fm = dm.CreateFlowG_StegerWarming(lm1, lm2, lm5)
        face.G = fp + fm

#---------------------------------------------------------------------------------------------------

    def CalcFlowH_StegerWarming(self, face, cp, cm):
        dp = cp.D
        dm = cm.D
        wp, wm = dp.w, dm.w
        ap, am = dp.a, dm.a
        lp1, lp2, lp5 = lp(wp - ap), lp(wp), lp(wp + ap)
        lm1, lm2, lm5 = lm(wm - am), lm(wm), lm(wm + am)
        fp = dp.CreateFlowH_StegerWarming(lp1, lp2, lp5)
        fm = dm.CreateFlowH_StegerWarming(lm1, lm2, lm5)
        face.H = fp + fm

#---------------------------------------------------------------------------------------------------

    def CalcFlowsStegerWarming(self):

        cs = self.Cells

        # Cells.
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                for k in range(self.CellsZ):

                    c = cs[i][j][k]
                    bs = c.Borders

                    # L
                    if bs[BorderL] == BorderSoft:
                        self.FacesX[i][j][k].F = c.D.CreateFlowF()
                    elif bs[BorderL] == BorderHard:
                        self.FacesX[i][j][k].F = c.D.CreateFlowF_Zero()
                    else:
                        pass

                    # R
                    if bs[BorderR] == BorderSoft:
                        self.FacesX[i + 1][j][k].F = c.D.CreateFlowF()
                    elif bs[BorderR] == BorderHard:
                        self.FacesX[i + 1][j][k].F = c.D.CreateFlowF_Zero()
                    else:
                        self.CalcFlowF_StegerWarming(self.FacesX[i + 1][j][k], c, cs[i + 1][j][k])

                    # D
                    if bs[BorderD] == BorderSoft:
                        self.FacesY[i][j][k].G = c.D.CreateFlowG()
                    elif bs[BorderD] == BorderHard:
                        self.FacesY[i][j][k].G = c.D.CreateFlowG_Zero()
                    else:
                        pass

                    # U
                    if bs[BorderU] == BorderSoft:
                        self.FacesY[i][j + 1][k].G = c.D.CreateFlowG()
                    elif bs[BorderU] == BorderHard:
                        self.FacesY[i][j + 1][k].G = c.D.CreateFlowG_Zero()
                    else:
                        self.CalcFlowG_StegerWarming(self.FacesY[i][j + 1][k], c, cs[i][j + 1][k])

                    # B
                    if bs[BorderB] == BorderSoft:
                        self.FacesZ[i][j][k].H = c.D.CreateFlowH()
                    elif bs[BorderB] == BorderHard:
                        self.FacesZ[i][j][k].H = c.D.CreateFlowH_Zero()
                    else:
                        pass

                    # F
                    if bs[BorderF] == BorderSoft:
                        self.FacesZ[i][j][k + 1].H = c.D.CreateFlowH()
                    elif bs[BorderF] == BorderHard:
                        self.FacesZ[i][j][k + 1].H = c.D.CreateFlowH_Zero()
                    else:
                        self.CalcFlowH_StegerWarming(self.FacesZ[i][j][k + 1], c, cs[i][j][k + 1])

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

    def ZValues(self, fun):
        return [fun(self.Cells[0][0][k]) for k in range(self.CellsZ)]

#---------------------------------------------------------------------------------------------------

    def Draw(self, fun):
        d = draw.Drawer(draw_area = (0.0, 0.0, self.SizeX, self.SizeY),
                        pic_size = (800, 800))
        (mn, mx) = self.FunInterval(fun)
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                cell = self.Cells[i][j][0]
                factor = 255 - int((fun(cell) - mn) / (mx - mn) * 255)
                color = (factor, factor, factor)
                d.Rect((i * self.dx, j * self.dy), ((i + 1) * self.dx, (j + 1) * self.dy),
                       aggdraw.Pen(color, 1.0), aggdraw.Brush(color))
            d.Rect((0.0, 0.0), (self.SizeX, self.SizeY), aggdraw.Pen('blue', 1.0))

        # Draw ellipse and cells.
        d.Ellipse(Sph.LoPoint().Tuple(2), Sph.HiPoint().Tuple(2), pen = aggdraw.Pen('black', 2.0))
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                cell = self.Cells[i][j][0]
                x1, x2 = i * g.dx, (i + 1) * g.dx
                y1, y2 = j * g.dy, (j + 1) * g.dy
                if cell.Type == TypeInner:
                    d.Rect((x1, y1), (x2, y2), pen = aggdraw.Pen('orange', 2.0))
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                cell = self.Cells[i][j][0]
                x1, x2 = i * g.dx, (i + 1) * g.dx
                y1, y2 = j * g.dy, (j + 1) * g.dy
                if cell.Type == TypeGhost:
                    d.Rect((x1, y1), (x2, y2), pen = aggdraw.Pen('red', 2.0))
                    d.Point((x1 + 0.5 * self.dx, y1 + 0.5 * self.dy), 1.0,
                            pen = aggdraw.Pen('red', 2.0))
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                cell = self.Cells[i][j][0]
                x1, x2 = i * g.dx, (i + 1) * g.dx
                y1, y2 = j * g.dy, (j + 1) * g.dy
                if cell.Type == TypePhantom:
                    d.Rect((x1, y1), (x2, y2), pen = aggdraw.Pen('blue', 2.0))
                    d.Point((x1 + 0.5 * self.dx, y1 + 0.5 * self.dy), 1.0,
                            pen = aggdraw.Pen('blue', 2.0))
        for i in range(self.CellsX):
            for j in range(self.CellsY):
                cell = self.Cells[i][j][0]
                x1, x2 = i * g.dx, (i + 1) * g.dx
                y1, y2 = j * g.dy, (j + 1) * g.dy
                if cell.Type == TypeBorder:
                    d.Rect((x1, y1), (x2, y2), pen = aggdraw.Pen('green', 2.0))
                    d.Point((x1 + 0.5 * self.dx, y1 + 0.5 * self.dy), 1.0,
                            pen = aggdraw.Pen('green', 2.0))

        d.FSS()

#---------------------------------------------------------------------------------------------------
# Main.
#---------------------------------------------------------------------------------------------------

Case_1D_X = 1
Case_1D_Y = 2
Case_1D_Z = 3
Case_2D_XY = 4
Case_1D_SodMod = 5

def create_and_init_grid(case):

    if case == Case_1D_X:
        g = Grid(1.0, 1.0, 1.0, 100, 1, 1)
    elif case == Case_1D_Y:
        g = Grid(1.0, 1.0, 1.0, 1, 100, 1)
    elif case == Case_1D_Z:
        g = Grid(1.0, 1.0, 1.0, 1, 1, 100)
    elif case == Case_2D_XY:
        g = Grid(1.0, 1.0, 1.0, 40, 40, 1)
    else:
        raise Exception('unknown case number')

    for i in range(g.CellsX):
        for j in range(g.CellsY):
            for k in range(g.CellsZ):
                c = g.Cells[i][j][k]
                x, y, z = (i + 0.5) * g.dx, (j + 0.5) * g.dy, (k + 0.5) * g.dz

                if case == Case_1D_X:
                    if x < 0.5:
                        c.D = DataD(10.0, 0.0, 0.0, 0.0, 10.0)
                    else:
                        c.D = DataD(1.0, 0.0, 0.0, 0.0, 1.0)
                    if i == g.CellsX - 1:
                        c.Borders[BorderR] = BorderHard
                elif case == Case_1D_Y:
                    if y < 0.5:
                        c.D = DataD(10.0, 0.0, 0.0, 0.0, 10.0)
                    else:
                        c.D = DataD(1.0, 0.0, 0.0, 0.0, 1.0)
                elif case == Case_1D_Z:
                    if z < 0.5:
                        c.D = DataD(10.0, 0.0, 0.0, 0.0, 10.0)
                    else:
                        c.D = DataD(1.0, 0.0, 0.0, 0.0, 1.0)
                elif case == Case_2D_XY:
                    if x < 0.1:
                        c.D = DataD(10.0, 0.0, 0.0, 0.0, 10.0)
                    else:
                        c.D = DataD(1.0, 0.0, 0.0, 0.0, 1.0)
                else:
                    raise Exception('unknown case number')

    if case == Case_2D_XY:
        # Ellipse equation.
        # (x - 0.6)^2 + (y - 0.5)^2 - 0.3^2 = 0
        elfun = lambda x, y: (x - Sph.C.X) ** 2 + (y - Sph.C.Y) ** 2 - Sph.R ** 2
        for i in range(g.CellsX):
            for j in range(g.CellsY):
                cell = g.Cells[i][j][0]
                x1, x2 = i * g.dx, (i + 1) * g.dx
                y1, y2 = j * g.dy, (j + 1) * g.dy
                f11, f12, f21, f22 = elfun(x1, y1), elfun(x1, y2), elfun(x2, y1), elfun(x2, y2)
                if abs(f11) < 0.001:
                    f11 = 0.0
                if abs(f12) < 0.001:
                    f12 = 0.0
                if abs(f21) < 0.001:
                    f21 = 0.0
                if abs(f22) < 0.001:
                    f22 = 0.0
                s11, s12, s21, s22 = mth.sign(f11), mth.sign(f12), mth.sign(f21), mth.sign(f22)
                ss = int(s11 + s12 + s21 + s22)
                if ss == -4:
                    cell.Type = TypeInner
                elif ss == 4:
                    cell.Type = TypeCommon
                else:
                    # Ghost or border.
                    xc, yc = (i + 0.5) * g.dx, (j + 0.5) * g.dy
                    if elfun(xc, yc) < 0.0:
                        cell.Type = TypeGhost
                    else:
                        cell.Type = TypeBorder
        for i in range(g.CellsX):
            for j in range(g.CellsY):
                for k in range(g.CellsZ):
                    cell = g.Cells[i][j][k]
                    if cell.Type == TypeBorder:
                        if g.Cells[i - 1][j][k].Type == TypeInner:
                            g.Cells[i - 1][j][k].Type = TypePhantom
                        if g.Cells[i + 1][j][k].Type == TypeInner:
                            g.Cells[i + 1][j][k].Type = TypePhantom
                        if g.Cells[i][j - 1][k].Type == TypeInner:
                            g.Cells[i][j - 1][k].Type = TypePhantom
                        if g.Cells[i][j + 1][k].Type == TypeInner:
                            g.Cells[i][j + 1][k].Type = TypePhantom
                        #if g.Cells[i][j][k - 1].Type == TypeInner:
                        #    g.Cells[i][j][k - 1].Type = TypePhantom
                        #if g.Cells[i][j][k + 1].Type == TypeInner:
                        #    g.Cells[i][j][k + 1].Type = TypePhantom

    return g

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    print('HYDRO_3D')
    g = create_and_init_grid(case = Case_2D_XY)

    pics, n, dt = 1, 10, 0.001
    fun = lambda cell: cell.D.p

    ts = time.time()
    for _ in range(pics):
        g.Steps(n, dt)
        g.Draw(fun)
        #vis.simple_graphic_ys(g.ZValues(fun))
    tf = time.time()
    print('total time : %f' % (tf - ts))

#---------------------------------------------------------------------------------------------------
