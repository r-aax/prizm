# -*- coding: utf-8 -*-
"""
Grid realization (Grid Surface Unstructured).

Created on Wed Jun  5 13:13:50 2019

@author: Rybakov
"""

from geom import Vector
import lst

#---------------------------------------------------------------------------------------------------
# Class node.
#---------------------------------------------------------------------------------------------------

class Node:

#---------------------------------------------------------------------------------------------------

    def __init__(self, v = Vector()):
        """
        Constructor.

        Arguments:
            v -- vector.
        """

        self.Vector = v
        self.Index = -1

#---------------------------------------------------------------------------------------------------

    def ToString(self):
        """
        Convert to string.

        Result:
            String.
        """

        return self.Vector.ToString()

#---------------------------------------------------------------------------------------------------
# Class grid.
#---------------------------------------------------------------------------------------------------

class Edge:

#---------------------------------------------------------------------------------------------------

    def __init__(self, a, b):
        """
        Constructor.

        Arguments:
            a -- first node,
            b -- second node.
        """

        self.Nodes = [a, b]

#---------------------------------------------------------------------------------------------------

    def ToString(self):
        """
        Convert to string.

        Result:
            String.
        """

        return '[%s - %s]' % (self.Nodes[0].ToString(), self.Nodes[1].ToString())

#---------------------------------------------------------------------------------------------------
# Class grid.
#---------------------------------------------------------------------------------------------------

class Face:

#---------------------------------------------------------------------------------------------------

    def __init__(self, a, b, c):
        """
        Constructor.

        Arguments:
            a -- first node,
            b -- second node,
            c -- third node.
        """

        self.Nodes = [a, b, c]

#---------------------------------------------------------------------------------------------------

    def ToString(self):
        """
        Convert to string.

        Result:
            String.
        """

        return '<%s, %s, %s>' % (self.Nodes[0].ToString(),
                                 self.Nodes[1].ToString(),
                                 self.Nodes[2].ToString())

#---------------------------------------------------------------------------------------------------
# Class grid.
#---------------------------------------------------------------------------------------------------

class Grid:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.Nodes = []
        self.Edges = []
        self.Faces = []

#---------------------------------------------------------------------------------------------------

    def NodesCount(self):
        """
        Count of nodes.

        Result:
            Nodes count.
        """

        return len(self.Nodes)

#---------------------------------------------------------------------------------------------------

    def EdgesCount(self):
        """
        Count of edges.

        Result:
            Edges count.
        """

        return len(self.Edges)

#---------------------------------------------------------------------------------------------------

    def FacesCount(self):
        """
        Count of faces.

        Result:
            Faces count.
        """

        return len(self.Faces)

#---------------------------------------------------------------------------------------------------

    def Print(self):
        """
        Print information.
        """

        print('Grid:')

        print('Nodes (%d):' % self.NodesCount())
        for n in self.Nodes:
            print('    %s' % n.ToString())

        print('Edges (%d):' % self.EdgesCount())
        for e in self.Edges:
            print('    %s' % e.ToString())

        print('Face (%d):' % self.FacesCount())
        for f in self.Faces:
            print('    %s' % f.ToString())

#---------------------------------------------------------------------------------------------------

    def AddNode(self, v):
        """
        Add new node.

        Arguments:
            v -- vector.
        """

        self.Nodes.append(Node(v))

#---------------------------------------------------------------------------------------------------

    def AddEdgeFromNodesIndices(self, ai, bi):
        """
        Connect nodes with edge.

        Arguments:
            ai -- first node index,
            bi -- second node index.
        """

        self.Edges.append(Edge(self.Nodes[ai],
                               self.Nodes[bi]))

#---------------------------------------------------------------------------------------------------

    def AddFaceFromNodesIndices(self, ai, bi, ci):
        """
        Group nodes with face.

        Arguments:
            ai -- first node index,
            bi -- second node index,
            ci -- third node index.
        """

        self.Faces.append(Face(self.Nodes[ai],
                               self.Nodes[bi],
                               self.Nodes[ci]))

#---------------------------------------------------------------------------------------------------

    def ConstructFromVectorsMatrix(self, m):
        """
        Construct grid from vectors matrix.
        Add nodes, edges, faces.

        Arguments:
            m -- matrix.
        """

        h = len(m)
        w = len(m[0])
        off = self.NodesCount()

        # Add nodes.
        for row in m:
            for v in row:
                self.AddNode(Vector(v))

        # Lambda function for index detection.
        idx = lambda r, c: r * w + c + off

        # Add horizontal edges.
        for r in range(h):
            for c in range(w - 1):
                self.AddEdgeFromNodesIndices(idx(r, c), idx(r, c + 1))
        # Add vertical edges.
        for r in range(h - 1):
            for c in range(w):
                self.AddEdgeFromNodesIndices(idx(r, c), idx(r + 1, c))

        # Diagonal edges.
        for r in range(h - 1):
            for c in range(w - 1):
                self.AddEdgeFromNodesIndices(idx(r, c + 1), idx(r + 1, c))

        # Add faces.
        for r in range(h - 1):
            for c in range(w - 1):
                self.AddFaceFromNodesIndices(idx(r, c), idx(r, c + 1), idx(r + 1, c))
                self.AddFaceFromNodesIndices(idx(r, c + 1), idx(r + 1, c), idx(r + 1, c + 1))

#---------------------------------------------------------------------------------------------------

    def ConstructFromVectorsFlatMatrix(self, h, w, m):
        """
        Construct grid from vectors matrix, written as list.

        Arguments:
            h -- height,
            w -- width,
            m -- flat matrix.
        """

        assert w * h == len(m), 'wrong flat matrix sizes'
        self.ConstructFromVectorsMatrix(lst.slice_rows(m, w))

#---------------------------------------------------------------------------------------------------

    def SetNodesIndices(self):
        """
        Set nodes indices.
        """

        for i in range(self.NodesCount()):
            self.Nodes[i].Index = i

#---------------------------------------------------------------------------------------------------

    def ExportToTecplot(self, filename):
        """
        Export to tecplot.

        Arguments:
            filename -- name of file.
        """

        f = open(filename, 'w')

        f.write('TITLE = "FE Surface Data ASCII"\n')
        f.write('VARIABLES = "X", "Y", "Z"\n')
        f.write('ZONE T="TRIANGLES", NODES="%d", ELEMENTS="%d", DATAPACKING="BLOCK", ZONETYPE="FETRIANGLE"\n'
                % (self.NodesCount(), self.FacesCount()))

        # Print coordinates.
        for n in self.Nodes:
            f.write('%f ' % n.Vector.X)
        f.write('\n')
        for n in self.Nodes:
            f.write('%f ' % n.Vector.Y)
        f.write('\n')
        for n in self.Nodes:
            f.write('%f ' % n.Vector.Z)
        f.write('\n')

        # Print faces.
        self.SetNodesIndices()
        for fc in self.Faces:
            f.write('%d %d %d\n' % (fc.Nodes[0].Index + 1,
                                    fc.Nodes[1].Index + 1,
                                    fc.Nodes[2].Index + 1))

        f.close()

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    g = Grid()
    g.ConstructFromVectorsFlatMatrix(3, 3, [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 2.0, 0.0),
                                            (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (1.0, 2.0, 0.0),
                                            (2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (2.0, 2.0, 2.0)])
    g.ExportToTecplot('air_inlet_2.dat')

#---------------------------------------------------------------------------------------------------
