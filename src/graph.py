"""
Graph realization.
"""

import itertools
import geom

# ==================================================================================================

EPS = 0.001

# ==================================================================================================


class Vertex:
    """Vertex class.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, p, color=0):
        """Constructor.

        Parameters
        ----------
        p : geom.Vector
            Coordinates.
        color : int
            Color.
        """

        self.Id = 0
        self.P = p
        self.Color = color
        self.Edges = []

    # ----------------------------------------------------------------------------------------------

    def __repr__(self):
        """Representation.

        Returns
        -------
        str
            Representation.
        """

        return f'V {self.P}, {self.Color}'

    # ----------------------------------------------------------------------------------------------

    def neighbour(self, e):
        """Get neighbour.

        Parameters
        ----------
        e : Edge
            Edge to find neighbour.

        Returns
        -------
        Vertex
            Neighbour.
        """

        v0, v1 = e.Vertices[0], e.Vertices[1]

        if self == v0:
            return v1
        elif self == v1:
            return v0
        else:
            raise Exception('Vertex.neighbour : internal error')

    # ----------------------------------------------------------------------------------------------

    def neighbours_colors(self):
        """Get list of neighbours colors.

        Returns
        -------
        list(int)
            List of neighbours colors.
        """

        return [self.neighbour(e).Color for e in self.Edges]

# ==================================================================================================


class Edge:
    """Edge class.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, v0, v1):
        """Constructor.

        Parameters
        ----------
        v0 : Vertex
            First vertex.
        v1 : Vertex
            Second vertex.
        """

        self.Vertices = [v0, v1]

    # ----------------------------------------------------------------------------------------------

    def __repr__(self):
        """Representation.

        Returns
        -------
        str
            Representation.
        """

        return f'{self.Vertices[0]} - {self.Vertices[1]}'

# ==================================================================================================


class Graph:
    """Graph class.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self):
        """Constructor.
        """

        self.Vertices = []
        self.Edges = []

    # ----------------------------------------------------------------------------------------------

    def new_vertex(self, p, color=0):
        """Create new vertex.

        Parameters
        ----------
        p : geom.Vector
            Coordinates.
        color : int
            Color.

        Returns
        -------
        Vertex
            New vertex.
        """

        v = Vertex(p=p, color=color)
        self.Vertices.append(v)

        return v

    # ----------------------------------------------------------------------------------------------

    def find_vertex(self, p):
        """Find vertex.

        Parameters
        ----------
        p : geom.Vector
            Point to find.

        Returns
        -------
            Found vertex or None.
        """

        for v in self.Vertices:
            if p.dist_to(v.P) < EPS:
                return v

        return None

    # ----------------------------------------------------------------------------------------------

    def find_or_new_vertex(self, p):
        """Find or new vertex.

        Parameters
        ----------
        p : geom.Vector
            Point.

        Returns
        -------
            Found or new vertex.
        """

        v = self.find_vertex(p)

        if v is None:
            v = self.new_vertex(p=p)

        return v

    # ----------------------------------------------------------------------------------------------

    def new_edge(self, v0, v1):
        """Create new edge.

        Parameters
        ----------
        v0 : Vertex
            First vertex.
        v1 : Vertex
            Second vertex.

        Returns
        -------
        Edge
            New edge.
        """

        e = Edge(v0, v1)
        v0.Edges.append(e)
        v1.Edges.append(e)
        self.Edges.append(e)

        return e

    # ----------------------------------------------------------------------------------------------

    def vertices_count(self):
        """Vertices count.

        Returns
        -------
        int
            Vertices count.
        """

        return len(self.Vertices)

    # ----------------------------------------------------------------------------------------------

    def edges_count(self):
        """Edges count.

        Returns
        -------
        int
            Edges count.
        """

        return len(self.Edges)

    # ----------------------------------------------------------------------------------------------

    def clear(self):
        """Clear graph.
        """

        self.Vertices.clear()
        self.Edges.clear()

    # ----------------------------------------------------------------------------------------------

    def set_vertices_ids(self):
        """Set vertices identifiers.
        """

        for (i, v) in enumerate(self.Vertices):
            v.Id = i

    # ----------------------------------------------------------------------------------------------

    def print(self):
        """Print graph.
        """

        print(f'Graph : {self.vertices_count()} vertices, {self.edges_count()} edges')
        for v in self.Vertices:
            print(f'  {v}')
        for e in self.Edges:
            print(f'  {e}')

    # ----------------------------------------------------------------------------------------------

    def init_ecg_for_2d_rect_mesh(self, cells_x, cells_y):
        """Init edges conflict graph for 2d rectangular mesh.

        Parameters
        ----------
        cells_x : int
            Cells count along OX axis.
        cells_y : int
            Cells count along OY axis.
        """

        self.clear()

        for xi in range(cells_x):
            for yi in range(cells_y):

                # Process one cell (xi, xi + 1) * (yi, yi + 1)
                # First create vertices in centers.
                l = self.find_or_new_vertex(geom.Vector(xi, yi + 0.5, 0.0))
                r = self.find_or_new_vertex(geom.Vector(xi + 1.0, yi + 0.5, 0.0))
                d = self.find_or_new_vertex(geom.Vector(xi + 0.5, yi, 0.0))
                u = self.find_or_new_vertex(geom.Vector(xi + 0.5, yi + 1.0, 0.0))

                # Now add conflicts.
                for (v0, v1) in itertools.combinations([l, r, d, u], 2):
                    self.new_edge(v0, v1)

    # ----------------------------------------------------------------------------------------------

    def init_ecg_for_3d_rect_mesh(self, cells_x, cells_y, cells_z):
        """Init edges conflict graph for 3d rectangular mesh.

        Parameters
        ----------
        cells_x : int
            Cells count along OX axis.
        cells_y : int
            Cells count along OY axis.
        cells_z : int
            Cells count along OZ axis.
        """

        self.clear()

        for xi in range(cells_x):
            for yi in range(cells_y):
                for zi in range(cells_z):

                    # Process one cell (xi, xi + 1) * (yi, yi + 1) * (zi, zi + 1)
                    # First create vertices in centers.
                    l = self.find_or_new_vertex(geom.Vector(xi, yi + 0.5, zi + 0.5))
                    r = self.find_or_new_vertex(geom.Vector(xi + 1.0, yi + 0.5, zi + 0.5))
                    d = self.find_or_new_vertex(geom.Vector(xi + 0.5, yi, zi + 0.5))
                    u = self.find_or_new_vertex(geom.Vector(xi + 0.5, yi + 1.0, zi + 0.5))
                    f = self.find_or_new_vertex(geom.Vector(xi + 0.5, yi + 0.5, zi))
                    b = self.find_or_new_vertex(geom.Vector(xi + 0.5, yi + 0.5, zi + 1.0))

                    # Now add conflicts.
                    for (v0, v1) in itertools.combinations([l, r, d, u, f, b], 2):
                        self.new_edge(v0, v1)

    # ----------------------------------------------------------------------------------------------

    def load(self, filename):
        """Load graph.

        Parameters
        ----------
        filename : str
            Name of file.
        """

        self.clear()

        with open(filename, 'r') as f:

            # Read head.
            f.readline()  # comment
            f.readline()  # title
            f.readline()  # variables
            f.readline()  # zone
            ns = int(f.readline().split('=')[-1])
            es = int(f.readline().split('=')[-1])
            f.readline()  # datapacking
            f.readline()  # zonetype

            # Create vertices.
            for _ in range(ns):
                self.new_vertex(geom.Vector())

            # Set coordinates and colors.
            xs = [float(xi) for xi in f.readline().split()]
            ys = [float(xi) for xi in f.readline().split()]
            zs = [float(xi) for xi in f.readline().split()]
            cs = [int(xi) for xi in f.readline().split()]
            for i in range(ns):
                v = self.Vertices[i]
                v.P.X, v.P.Y, v.P.Z, v.C = xs[i], ys[i], zs[i], cs[i]

            # Load data for edges.
            for _ in range(es):
                ds = [int(i) for i in f.readline().split()]
                self.new_edge(self.Vertices[ds[0] - 1], self.Vertices[ds[1] - 1])

            f.close()

        self.set_vertices_ids()

    # ----------------------------------------------------------------------------------------------

    def save(self, filename):
        """Save graph.

        Parameters
        ----------
        filename : str
            Name of file.
        """

        self.set_vertices_ids()

        with open(filename, 'w') as f:

            # Head.
            f.write('# Graph.\n')
            f.write('TITLE="Graph coloring"\n')
            f.write('VARIABLES="X", "Y", "Z", "C"\n')
            f.write('ZONE T="Graph coloring"\n')
            f.write(f'NODES={self.vertices_count()}\n')
            f.write(f'ELEMENTS={self.edges_count()}\n')
            f.write('DATAPACKING=BLOCK\n')
            f.write('ZONETYPE=FETRIANGLE\n')

            # Vertices.
            f.write(' '.join([str(v.P.X) for v in self.Vertices]) + '\n')
            f.write(' '.join([str(v.P.Y) for v in self.Vertices]) + '\n')
            f.write(' '.join([str(v.P.Z) for v in self.Vertices]) + '\n')
            f.write(' '.join([str(v.Color) for v in self.Vertices]) + '\n')

            # Edges.
            for e in self.Edges:
                id0, id1 = e.Vertices[0].Id + 1, e.Vertices[1].Id + 1
                f.write(f'{id0} {id1} {id1}\n')

            f.close()

    # ----------------------------------------------------------------------------------------------

    def decolor(self):
        """Reset colors.
        """

        for v in self.Vertices:
            v.Color = 0

    # ----------------------------------------------------------------------------------------------

    def coloring_greedy(self):
        """Greedy coloring algorithm.

        Trying to color vertex by vertex in greedy manner.

        Returns
        -------
        int
            Colors count.
        """

        self.decolor()
        max_c = 0

        for v in self.Vertices:
            cs = v.neighbours_colors()
            c = 1
            while c in cs:
                c += 1
            max_c = max(max_c, c)
            v.Color = c

        print(f'Graph.coloring_greedy : {max_c} colors')

        return max_c

# ----------------------------------------------------------------------------------------------

    def coloring_greedy2(self):
        """Greedy coloring algorithm.

        First try to color maxinum number of vertices with color 1, 2, and so on..

        Returns
        -------
        int
            Colors count.
        """

        self.decolor()
        max_c = 0

        # Put all vertices to new array.
        a = [v for v in self.Vertices]

        # Try to color all possible vertices with new color.
        while len(a) > 0:
            max_c += 1
            for v in a:
                cs = v.neighbours_colors()
                if not (max_c in cs):
                    v.Color = max_c
            a = [v for v in a if v.Color == 0]

        print(f'Graph.coloring_greedy2 : {max_c} colors')

        return max_c

# ==================================================================================================


if __name__ == '__main__':

    g = Graph()

    # 2d regular mesh.
    # g.init_ecg_for_2d_rect_mesh(4, 3)

    # 3d regular mesh.
    # g.init_ecg_for_3d_rect_mesh(4, 3, 2)

    # Make colorings for bunny.
    grid_file = 'bunny_ecg.dat'
    g.load(grid_file)

    g.coloring_greedy()
    g.save('greedy.dat')
    g.coloring_greedy2()
    g.save('greedy2.dat')

# ==================================================================================================
