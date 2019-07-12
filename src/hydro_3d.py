# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:26:02 2019

@author: Rybakov
"""

import draw
import aggdraw

#---------------------------------------------------------------------------------------------------
# Constants.
#---------------------------------------------------------------------------------------------------

# Gamma for GAS condition.
Gamma = 1.4

# Count of rectangle faces.
RectangleFacesCount = 6

# Directions.
DirLR = 0
DirDU = 1
DirBF = 2

# Border conditons types.
BorderConditionHard = 0
BorderConditionFree = 1

#---------------------------------------------------------------------------------------------------
# Class DataD.
#---------------------------------------------------------------------------------------------------

class D:

#---------------------------------------------------------------------------------------------------

    def __init__(self, r = 0.0, u = 0.0, v = 0.0, w = 0.0, p = 0.0):
        """
        Constructor.

        Arguments:
            r -- density,
            u -- velosity X,
            v -- velosity Y,
            w -- velosity Z,
            p -- pressure.
        """

        self.r = r
        self.u = u
        self.v = v
        self.w = w
        self.p = p

#---------------------------------------------------------------------------------------------------

    def e(self):
        """
        Inner energy.

        Result:
            Inner energy.
        """

        return self.p / ((Gamma - 1.0) * self.r)

#---------------------------------------------------------------------------------------------------

    def V2(self):
        """
        Mod2 of velosity.

        Result:
            Mod2 of velosity.
        """

        return self.u * self.u + self.v * self.v + self.w * self.w

#---------------------------------------------------------------------------------------------------

    def FromU(self, U):
        """
        Set from U data.

        Arguments:
            U -- U data.
        """

        self.r = U.r
        self.u = U.ru / U.r
        self.v = U.rv / U.r
        self.w = U.rw / U.r
        self.p = (U.E / U.r - 0.5 * self.V2()) * ((Gamma - 1.0) * U.r)

#---------------------------------------------------------------------------------------------------
# Class U.
#---------------------------------------------------------------------------------------------------

class U:

#---------------------------------------------------------------------------------------------------

    def __init__(self, r = 0.0, ru = 0.0, rv = 0.0, rw = 0.0, E = 0.0):
        """
        Constructor.

        Arguments:
            r -- density,
            ru -- density * velosity X,
            rv -- density * velosity Y,
            rw -- density * velosity Z,
            E -- full energy.
        """

        self.r = r
        self.ru = ru
        self.rv = rv
        self.rw = rw
        self.E = E

#---------------------------------------------------------------------------------------------------

    def FromD(self, D):
        """
        Set from D.

        Arguments:
            D -- D data.
        """

        self.r = D.r
        self.ru = D.r * D.u
        self.rv = D.r * D.v
        self.rw = D.r * D.w
        self.E = D.r * (0.5 * D.V2() + D.e())

#---------------------------------------------------------------------------------------------------
# Cell class.
#---------------------------------------------------------------------------------------------------

class Cell:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor:
        """

        # Physical data.
        self.D = D()
        self.U = U()

        self.Faces = []

#---------------------------------------------------------------------------------------------------

    def SetSizes(self, dx, dy, dz):
        """
        Set cell sizes.

        Arguments:
            dx -- size in X direction,
            dy -- size in Y direction,
            dz -- size in Z direction.
        """

        self.dx = dx
        self.dy = dy
        self.dz = dz

#---------------------------------------------------------------------------------------------------

    def SetCoordinates(self, left, down, back):
        """
        Set coordinates.

        Arguments:
            left -- left X,
            down -- down Y,
            back -- back Z.
        """

        self.Left = left
        self.Right = self.Left + self.dx
        self.Down = down
        self.Up = self.Down + self.dy
        self.Back = back
        self.Front = self.Back + self.dz

#---------------------------------------------------------------------------------------------------

    def Center(self):
        """
        Get center.

        Result:
            Center point.
        """

        return (0.5 * (self.Left + self.Right),
                0.5 * (self.Down + self.Up),
                0.5 * (self.Back + self.Front))

#---------------------------------------------------------------------------------------------------
# Face class.
#---------------------------------------------------------------------------------------------------

class Face:

#---------------------------------------------------------------------------------------------------

    def __init__(self, d, bc_type = BorderConditionFree):
        """
        Constructor.

        Arguments:
            d -- direction.
        """

        self.Dir = d
        self.BorderConditionType = bc_type

        # Flow.
        self.FGH = U()

        self.Cells = []

#---------------------------------------------------------------------------------------------------
# Grid class.
#---------------------------------------------------------------------------------------------------

class Grid:

#---------------------------------------------------------------------------------------------------

    def __init__(self,
                 size_x, size_y, size_z):
        """
        Constructor:

        Arguments:
            size_x -- size in X direction,
            size_y -- size in Y direction,
            size_z -- size in Z direction.
        """

        self.SizeX = size_x
        self.SizeY = size_y
        self.SizeZ = size_z

        self.Cells = []
        self.Faces = []

#---------------------------------------------------------------------------------------------------

    def Link(self, d,
             cell1, cell1_face_index,
             cell2, cell2_face_index):
        """
        Link two cells with face.

        Arguments:
            d -- direction,
            cell1 -- first cell,
            cell1_face_index -- index of face in cell1 faces list,
            cell2 -- second face,
            cell2_face_index -- index of face in cell2 faces list.
        """

        face = Face(d)
        face.Cells = [cell1, cell2]
        self.Faces.append(face)
        if cell1_face_index >= 0:
            cell1.Faces[cell1_face_index] = face
        if cell2_face_index >= 0:
            cell2.Faces[cell2_face_index] = face

#---------------------------------------------------------------------------------------------------

    def LinkLR(self, cell_l, cell_r):
        """
        Link L-R cells.

        Arguments:
            cell_l -- left cell,
            cell_r -- right cell.
        """

        self.Link(DirLR, cell_l, 1, cell_r, 0)

#---------------------------------------------------------------------------------------------------

    def LinkDU(self, cell_d, cell_u):
        """
        Link D-U cells.

        Arguments:
            cell_d -- down cell,
            cell_u -- upper cell.
        """

        self.Link(DirDU, cell_d, 3, cell_u, 2)

#---------------------------------------------------------------------------------------------------

    def LinkBF(self, cell_b, cell_f):
        """
        Link B-F cells.

        Arguments:
            cell_b -- back cell,
            cell_f -- front cell.
        """

        self.Link(DirBF, cell_b, 5, cell_f, 4)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionL(self, cell):
        """
        Set left boundary condition.

        Arguments:
            cell -- cell.
        """

        self.Link(DirLR, None, -1, cell, 0)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionR(self, cell):
        """
        Set right boundary condition.

        Arguments:
            cell -- cell.
        """

        self.Link(DirLR, cell, 1, None, -1)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionD(self, cell):
        """
        Set down boundary condition.

        Arguments:
            cell -- cell.
        """

        self.Link(DirDU, None, -1, cell, 2)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionU(self, cell):
        """
        Set up boundary condition.

        Arguments:
            cell -- cell.
        """

        self.Link(DirDU, cell, 3, None, -1)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionB(self, cell):
        """
        Set back boundary condition.

        Arguments:
            cell -- cell.
        """

        self.Link(DirBF, None, -1, cell, 4)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionF(self, cell):
        """
        Set front boundary condition.

        Arguments:
            cell -- cell.
        """

        self.Link(DirBF, cell, 5, None, -1)

#---------------------------------------------------------------------------------------------------

    def CreateUniformGrid(self,
                          cells_count_x, cells_count_y, cells_count_z):
        """
        Create uniform grid (cells with same sizes).

        Arguments:
            cells_count_x -- count of cells in X direction,
            cells_count_y -- count of cells in Y direction,
            cells_count_Z -- count of cells in Z direction.
        """

        # Create cells array.
        cells_count = cells_count_x * cells_count_y * cells_count_z
        self.Cells = [None] * cells_count

        # Cell size.
        dx = self.SizeX / cells_count_x
        dy = self.SizeY / cells_count_y
        dz = self.SizeZ / cells_count_z

        # Create cells.
        cur = 0
        for k in range(cells_count_z):
            for j in range(cells_count_y):
                for i in range(cells_count_x):
                    cell = Cell()
                    cell.SetSizes(dx, dy, dz)
                    cell.SetCoordinates(i * dx, j * dy, k * dz)
                    cell.Faces = [None] * RectangleFacesCount
                    self.Cells[cur] = cell
                    cur += 1

        # Link cells with faces.
        cur = 0
        for k in range(cells_count_z):
            for j in range(cells_count_y):
                for i in range(cells_count_x):

                    if k == 0:
                        self.SetBoundaryConditionB(self.Cells[cur])
                    if k != cells_count_z - 1:
                        self.LinkBF(self.Cells[cur],
                                    self.Cells[cur + cells_count_x * cells_count_y])
                    else:
                        self.SetBoundaryConditionF(self.Cells[cur])

                    if j == 0:
                        self.SetBoundaryConditionD(self.Cells[cur])
                    if j != cells_count_y - 1:
                        self.LinkLR(self.Cells[cur], self.Cells[cur + cells_count_x])
                    else:
                        self.SetBoundaryConditionU(self.Cells[cur])

                    if i == 0:
                        self.SetBoundaryConditionL(self.Cells[cur])
                    if i != cells_count_x - 1:
                        self.LinkLR(self.Cells[cur], self.Cells[cur + 1])
                    else:
                        self.SetBoundaryConditionR(self.Cells[cur])

                    cur += 1

#---------------------------------------------------------------------------------------------------

    def InitD(self):
        """
        Init D data.
        """

        for cell in self.Cells:
            (x, y, z) = cell.Center()
            D = cell.D
            D.r = 200.0 / (x + y)
            D.u = 0.0
            D.v = 0.0
            D.w = 0.0
            D.p = 1.0

#---------------------------------------------------------------------------------------------------

    def DtoU(self):
        """
        Convert D to U.
        """

        for cell in self.Cells:
            cell.U.FromD(cell.D)

#---------------------------------------------------------------------------------------------------

    def UtoD(self):
        """
        Convert U to D.
        """

        for cell in self.Cells:
            cell.D.FromU(cell.U)

#---------------------------------------------------------------------------------------------------

    def Step(self, dt):
        """
        Step.

        Arguments:
            dt -- tie step.
        """

        self.DtoU()
        self.UtoD()

#---------------------------------------------------------------------------------------------------

    def PrintInfo(self):
        """
        Print information.
        """

        print('grid : %d cells, %d faces'
              % (len(self.Cells), len(self.Faces)))

#---------------------------------------------------------------------------------------------------

    def Draw(self, cell_draw_value_fun):
        """
        Draw grid.

        Arguments:
            cell_draw_value_fun -- function of draw value.
        """

        d = draw.Drawer(draw_area = (0.0, 0.0, self.SizeX, self.SizeY),
                        pic_size = (500, 500))

        # Draw cells.
        for cell in self.Cells:
            dv = int(cell_draw_value_fun(cell) * 255)
            color = (dv, dv, dv)
            pen = aggdraw.Pen(color, 1.0)
            brush = aggdraw.Brush(color)
            d.Rect((cell.Left, cell.Down), (cell.Right, cell.Up),
                   pen = pen, brush = brush)

        d.FSS()

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def cell_draw_value(cell):
    """
    Cell draw value from 0.0 to 1.0.

    Return:
        Value for drawing.
    """

    v = cell.D.r
    limit = 2.0

    if v < 0.0:
        return 1.0
    elif v > limit:
        return 0.0
    else:
        return 1.0 - v / limit

#---------------------------------------------------------------------------------------------------
# Main.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    print('HYDRO_2D:')

    # Init.
    grid = Grid(100.0, 100.0, 1.0)
    grid.CreateUniformGrid(100, 100, 1)
    grid.InitD()

    # Print.
    grid.PrintInfo()

    # Calc.
    grid.Step(0.001)

    # Draw.
    grid.Draw(cell_draw_value)

#---------------------------------------------------------------------------------------------------
