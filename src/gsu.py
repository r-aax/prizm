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

    def __init__(self, v):
        """
        Constructor.

        Arguments:
            v -- vector.
        """

        self.Vector = v
        self.HTC = 0.0
        self.Beta = 0.0
        self.Tau = Vector()
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
        self.T = 0.0
        self.Hw = 0.0
        self.Hi = 0.0

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return '<%s, %s, %s>' % (str(self.Nodes[0]),
                                 str(self.Nodes[1]),
                                 str(self.Nodes[2]))

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
        self.ConstructFromVectorsMatrix(lst.chop(m, w))

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
                '"T", "Hw", "Hi", '
                '"Node_HTC", "Node_Beta", '
                '"Node_TauX", "Node_TauY", "Node_TauZ"\n')
        f.write('ZONE T="TRIANGLES", '
                'NODES=%d, '
                'ELEMENTS=%d, '
                'DATAPACKING=BLOCK, '
                'ZONETYPE=FETRIANGLE '
                'VARLOCATION=([4-6]=CELLCENTERED)\n'
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

        # Values (3 in faces, 5 in nodes).
        for fc in self.Faces:
            f.write('%f ' % fc.T)
        f.write('\n')
        for fc in self.Faces:
            f.write('%f ' % fc.Hw)
        f.write('\n')
        for fc in self.Faces:
            f.write('%f ' % fc.Hi)
        f.write('\n')
        for n in self.Nodes:
            f.write('%f ' % n.HTC)
        f.write('\n')
        for n in self.Nodes:
            f.write('%f ' % n.Beta)
        f.write('\n')
        for n in self.Nodes:
            f.write('%f ' % n.Tau.X)
        f.write('\n')
        for n in self.Nodes:
            f.write('%f ' % n.Tau.Y)
        f.write('\n')
        for n in self.Nodes:
            f.write('%f ' % n.Tau.Z)
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

    def ImportFromTecplot(self, filename):
        """
        Import from tecplot.

        Arguments:
            filename -- name of file.
        """

        is_zone_str = lambda s: 'ZONE' in s

        with open(filename, 'r') as f:

            # Comment.
            mode_l = f.readline()
            if mode_l != '# EXPORT MODE: CHECK_POINT\n':
                raise Exception('wrong mode line : %s' % mode_l)

            # Type.
            title_l = f.readline()
            if title_l != 'TITLE = "FE Surface Data ASCII"\n':
                raise Exception('wrong title line : %s' % title_l)

            # Variables.
            variables_l = f.readline()
            if variables_l != 'VARIABLES = "X", "Y", "Z", "T", "Hw", "Hi", "Node_HTC", "Node_Beta", "Node_TauX", "Node_TauY", "Node_TauZ"\n':
                raise Exception('wrong variables line : %s' % variables_l)

            l = f.readline()
            while l:

                if is_zone_str(l):
                    zone_l_split = l.split()
                    nodes_count = int(zone_l_split[2][6:-1])
                    faces_count = int(zone_l_split[3][9:-1])
                    print('imported block information : nodes_count = %d, faces_count = %d'
                          % (nodes_count, faces_count))

                    # Read nodes and faces data.
                    x_s = f.readline().split()
                    y_s = f.readline().split()
                    z_s = f.readline().split()
                    t_s = f.readline().split()
                    hw_s = f.readline().split()
                    hi_s = f.readline().split()
                    node_htc_s = f.readline().split()
                    node_beta_s = f.readline().split()
                    node_tau_x_s = f.readline().split()
                    node_tau_y_s = f.readline().split()
                    node_tau_z_s = f.readline().split()

                    # Add nodes.
                    cur_nodes_count = self.NodesCount()
                    cur_faces_count = self.FacesCount()
                    for _ in range(nodes_count):
                        self.AddNode(Vector())
                    for i in range(nodes_count):
                        node = self.Nodes[cur_nodes_count + i]
                        node.Vector.X = float(x_s[i])
                        node.Vector.Y = float(y_s[i])
                        node.Vector.Z = float(z_s[i])
                        node.HTC = float(node_htc_s[i])
                        node.Beta = float(node_beta_s[i])
                        node.Tau.X = float(node_tau_x_s[i])
                        node.Tau.Y = float(node_tau_y_s[i])
                        node.Tau.Z = float(node_tau_z_s[i])

                    # Add faces.
                    for _ in range(faces_count):
                        l = f.readline()
                        s = l.split()
                        self.AddFaceFromNodesIndices(cur_nodes_count + int(s[0]) - 1,
                                                     cur_nodes_count + int(s[1]) - 1,
                                                     cur_nodes_count + int(s[2]) - 1)
                    for i in range(faces_count):
                        face = self.Faces[cur_faces_count + i]
                        face.T = float(t_s[i])
                        face.Hw = float(hw_s[i])
                        face.Hi = float(hi_s[i])
                
                else:
                    raise Exception('unexpected line : %s' % l)

                l = f.readline()

            f.close()

#---------------------------------------------------------------------------------------------------

    def NearestNodeInSortedNodes(self, x, y, z):
        """
        Find nearest node.

        Arguments:
            x -- x coordinate,
            y -- y coordinate,
            z -- z coordinate.

        Result:
            Nearest node.
        """

        eps = 0.001

        # Search in x.
        li = 0
        hi = len(self.Nodes) - 1
        while self.Nodes[hi].Vector.X - self.Nodes[li].Vector.X > eps:
            mi = (li + hi) // 2
            if x <= self.Nodes[mi].Vector.X:
                if hi == mi:
                    break
                hi = mi
            else:
                if li == mi:
                    break
                li = mi

        v = Vector((x, y, z))

        # Total search.
        i = li
        res = self.Nodes[i]
        dist = res.Vector.DistTo(v)
        i += 1
        while i <= hi:
            node = self.Nodes[i]
            new_dist = node.Vector.DistTo(v)
            if new_dist < dist:
                res = node
                dist = new_dist
            i += 1

        return (res, dist)

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
