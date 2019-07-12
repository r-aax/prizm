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

# Count of rectangle faces.
RectangleFacesCount = 6

#---------------------------------------------------------------------------------------------------
# Cell class.
#---------------------------------------------------------------------------------------------------

class Cell:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor:
        """

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

    def __init__(self):
        """
        Constructor.
        """

        self.Cells = []

#---------------------------------------------------------------------------------------------------
# Boundary condition class.
#---------------------------------------------------------------------------------------------------

class BoundaryCondition:

    # Hard border.
    TypeHard = 0

    # Free border.
    TypуFree = 1

#---------------------------------------------------------------------------------------------------

    def __init__(self, t):
        """
        Constructor.

        Arguments:
            type -- border type.
        """

        self.Type = t

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
        self.BoundaryConditions = []

#---------------------------------------------------------------------------------------------------

    def SetBoundaryCondition(self, cell, cell_face_index, t):
        """
        Set boundary condition.

        Arguments:
            cell -- cell,
            cell_face_index -- face index,
            t -- type.
        """

        bc = BoundaryCondition(t)
        self.BoundaryConditions.append(bc)
        cell.Faces[cell_face_index] = bc

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionL(self, cell):
        """
        Set left boundary condition.

        Arguments:
            cell -- cell.
        """

        self.SetBoundaryCondition(cell, 0, BoundaryCondition.TypуFree)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionR(self, cell):
        """
        Set right boundary condition.

        Arguments:
            cell -- cell.
        """

        self.SetBoundaryCondition(cell, 1, BoundaryCondition.TypуFree)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionD(self, cell):
        """
        Set down boundary condition.

        Arguments:
            cell -- cell.
        """

        self.SetBoundaryCondition(cell, 2, BoundaryCondition.TypуFree)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionU(self, cell):
        """
        Set up boundary condition.

        Arguments:
            cell -- cell.
        """

        self.SetBoundaryCondition(cell, 3, BoundaryCondition.TypуFree)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionB(self, cell):
        """
        Set back boundary condition.

        Arguments:
            cell -- cell.
        """

        self.SetBoundaryCondition(cell, 4, BoundaryCondition.TypуFree)

#---------------------------------------------------------------------------------------------------

    def SetBoundaryConditionF(self, cell):
        """
        Set front boundary condition.

        Arguments:
            cell -- cell.
        """

        self.SetBoundaryCondition(cell, 5, BoundaryCondition.TypуFree)

#---------------------------------------------------------------------------------------------------

    def Link(self,
             cell1, cell1_face_index,
             cell2, cell2_face_index):
        """
        Link two cells with face.

        Arguments:
            cell1 -- first cell,
            cell1_face_index -- index of face in cell1 faces list,
            cell2 -- second face,
            cell2_face_index -- index of face in cell2 faces list.
        """

        face = Face()
        face.Cells = [cell1, cell2]
        self.Faces.append(face)
        cell1.Faces[cell1_face_index] = face
        cell2.Faces[cell2_face_index] = face

#---------------------------------------------------------------------------------------------------

    def LinkLR(self, cell_l, cell_r):
        """
        Link L-R cells.

        Arguments:
            cell_l -- left cell,
            cell_r -- right cell.
        """

        self.Link(cell_l, 1, cell_r, 0)

#---------------------------------------------------------------------------------------------------

    def LinkDU(self, cell_d, cell_u):
        """
        Link D-U cells.

        Arguments:
            cell_d -- down cell,
            cell_u -- upper cell.
        """

        self.Link(cell_d, 3, cell_u, 2)

#---------------------------------------------------------------------------------------------------

    def LinkBF(self, cell_b, cell_f):
        """
        Link B-F cells.

        Arguments:
            cell_b -- back cell,
            cell_f -- front cell.
        """

        self.Link(cell_b, 5, cell_f, 4)

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

    def SetPressure(self, pressure_fun):
        """
        Set pressure.

        Arguments:
            pressure_fun -- function of pressure.
        """

        for cell in self.Cells:
            cell.p = pressure_fun(cell.Center())

#---------------------------------------------------------------------------------------------------

    def PrintInfo(self):
        """
        Print information.
        """

        print('grid : %d cells, %d faces, %d boundary condidtions'
              % (len(self.Cells), len(self.Faces), len(self.BoundaryConditions)))

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
            print(dv)
            color = (dv, dv, dv)
            pen = aggdraw.Pen(color, 1.0)
            brush = aggdraw.Brush(color)
            d.Rect((cell.Left, cell.Down), (cell.Right, cell.Up),
                   pen = pen, brush = brush)

        d.FSS()

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def pressure_fun(point):
    """
    Set pressure function.

    Arguments:
        point -- point.

    Return:
        Pressure value.
    """

    (x, y, z) = point
    p = x * x + y * y

    return p

#---------------------------------------------------------------------------------------------------

def cell_draw_value(cell):
    """
    Cell draw value from 0.0 to 1.0.

    Return:
        Value for drawing.
    """

    p = cell.p
    limit = 10000.0

    if p < 0.0:
        return 0.0
    elif p > limit:
        return 1.0
    else:
        return p / limit

#---------------------------------------------------------------------------------------------------
# Main.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    print('HYDRO_2D:')
    grid = Grid(100.0, 100.0, 1.0)
    grid.CreateUniformGrid(100, 100, 1)
    grid.SetPressure(pressure_fun)
    grid.PrintInfo()
    grid.Draw(cell_draw_value)

#---------------------------------------------------------------------------------------------------
