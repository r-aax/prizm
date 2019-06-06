# -*- coding: utf-8 -*-
"""
Grid realization (Grid Surface Unstructured).

Created on Wed Jun  5 13:13:50 2019

@author: Rybakov
"""

from geom import Vector
import lst
import fun

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
        self.Ref = None

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return str(self.Vector)

#---------------------------------------------------------------------------------------------------

    def __lt__(self, n):
        """
        Overload '<'

        Arguments:
            n -- node.

        Result:
            Boolean value.
        """

        return self.Vector < n.Vector

#---------------------------------------------------------------------------------------------------

    def GRef(self):
        """
        Global recurent reference.

        Result:
            Reference.
        """

        if self.Ref == None:
            return self
        else:
            return self.Ref.GRef()

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

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return '[%s - %s]' % (str(self.Nodes[0]), str(self.Nodes[1]))

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

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return '<%s, %s, %s>' % (str(self.Nodes[0]), str(self.Nodes[1]), str(self.Nodes[2]))

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
            print('    %s' % str(n))

        print('Edges (%d):' % self.EdgesCount())
        for e in self.Edges:
            print('    %s' % str(e))

        print('Faces (%d):' % self.FacesCount())
        for f in self.Faces:
            print('    %s' % str(f))

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
                nul = idx(r, c)
                nur = idx(r, c + 1)
                self.AddEdgeFromNodesIndices(nul, nur)
        # Add vertical edges.
        for r in range(h - 1):
            for c in range(w):
                nul = idx(r, c)
                ndl = idx(r + 1, c)
                self.AddEdgeFromNodesIndices(nul, ndl)

        # Diagonal edges.
        for r in range(h - 1):
            for c in range(w - 1):
                nur = idx(r, c + 1)
                ndl = idx(r + 1, c)
                self.AddEdgeFromNodesIndices(nur, ndl)

        # Add faces.
        for r in range(h - 1):
            for c in range(w - 1):
                nul = idx(r, c)
                nur = idx(r, c + 1)
                ndl = idx(r + 1, c)
                ndr = idx(r + 1, c + 1)
                self.AddFaceFromNodesIndices(nul, ndl, nur)
                self.AddFaceFromNodesIndices(ndr, nur, ndl)

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

        f.write('# EXPORT MODE: CHECK_POINT\n')
        f.write('TITLE = "FE Surface Data ASCII"\n')
        f.write('VARIABLES = "X", "Y", "Z", '
                '"Node_HTC", "Node_Beta", '
                '"Node_TauX", "Node_TauY", "Node_TauZ", '
                '"T", "Hw", "Hi"\n')
        f.write('ZONE T="TRIANGLES", '
                'NODES=%d, '
                'ELEMENTS=%d, '
                'DATAPACKING="BLOCK", '
                'ZONETYPE="FETRIANGLE" '
                'VARLOCATION=([9-11]=CELLCENTERED)\n'
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

        # 8 zeroed values (5 in nodes, 3 in faces).
        for iteration in range(5):
            for n in self.Nodes:
                f.write('0.0 ')
            f.write('\n')
        for iteration in range(3):
            for fc in self.Faces:
                f.write('0.0 ')
            f.write('\n')

        # Print faces.
        self.SetNodesIndices()
        for fc in self.Faces:
            f.write('%d %d %d\n' % (fc.Nodes[0].Index + 1,
                                    fc.Nodes[1].Index + 1,
                                    fc.Nodes[2].Index + 1))

        f.close()

#---------------------------------------------------------------------------------------------------

    def FindPairsOfNearPoints(self, eps):
        """
        Find pairs of points near to each other.

        Arguments:
            eps -- epsilon (if distance between points is less than eps, then points are the same).
        """

        # Sort nodes.
        ns = self.Nodes.copy()
        ns.sort()

        # Check each pair.
        count = 0
        for i in range(len(ns) - 1):
            n1 = ns[i]
            n2 = ns[i + 1]
            if n1.Vector.DistTo(n2.Vector) < eps:
                n2.Ref = n1
                count = count + 1

        # Report.
        print('FindPairsOfNearPoints : %d pairs are found.' % count)

#---------------------------------------------------------------------------------------------------

    def MergeNodes(self):
        """
        Merge nodes.
        """

        # Redirect links from b node to a node.
        for e in self.Edges:
            for i in range(2):
                if e.Nodes[i].Ref != None:
                    e.Nodes[i] = e.Nodes[i].GRef()
        for f in self.Faces:
            for i in range(3):
                if f.Nodes[i].Ref != None:
                    f.Nodes[i] = f.Nodes[i].GRef()

        # Refresh nodes.
        self.Nodes = [n for n in self.Nodes if n.Ref == None]

        # Restore nodes indices again.
        self.SetNodesIndices();

#---------------------------------------------------------------------------------------------------

    def SewFaces(self, eps = 1.0e-6):
        """
        Sew faces.

        Arguments:
            eps -- epsilon for similarity detect.
        """

        # Prepare.
        self.SetNodesIndices()
        for n in self.Nodes:
            n.Flag = False

        # Find nodes pairs for merge.
        self.FindPairsOfNearPoints(eps)

        # And merge them.
        self.MergeNodes()

        # TODO:
        # Delete extra edges.

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    g = Grid()
    g.ConstructFromVectorsFlatMatrix(3, 3, [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 2.0, 0.0),
                                            (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (1.0, 2.0, 0.0),
                                            (2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (2.0, 2.0, 0.0)])
    g.ConstructFromVectorsFlatMatrix(2, 2, [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0),
                                            (-1.0, 0.0, 0.0), (-1.0, 1.0, 0.0)])
    g.Print()
    g.SewFaces()
    g.Print()

    g.ExportToTecplot('test.dat')

#---------------------------------------------------------------------------------------------------
