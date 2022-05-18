"""
Graph realization.
"""

import itertools
import random
import plotly.graph_objects as go
import networkx as nx
import geom

# ==================================================================================================

# Small value.
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

        self.P = p
        self.Color = color
        self.Edges = []
        self.Faces = []

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

    def neighbours(self):
        """Get all neighbours.

        Returns
        -------
        list(Vertex)
            List of neighbours.
        """

        return [self.neighbour(e) for e in self.Edges]

    # ----------------------------------------------------------------------------------------------

    def neighbours_colors(self):
        """Get list of neighbours colors.

        Returns
        -------
        list(int)
            List of neighbours colors.
        """

        return [n.Color for n in self.neighbours()]

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
        self.Faces = []

    # ----------------------------------------------------------------------------------------------

    def __repr__(self):
        """Representation.

        Returns
        -------
        str
            Representation.
        """

        return f'{self.Vertices[0]} - {self.Vertices[1]}'

    # ----------------------------------------------------------------------------------------------

    @property
    def A(self):
        """A vertex.

        Returns
        -------
        Vertex
            First vertex.
        """

        return self.Vertices[0]

    # ----------------------------------------------------------------------------------------------

    @property
    def B(self):
        """B vertex.

        Returns
        -------
        Vertex
            Second vertex.
        """

        return self.Vertices[1]

# ==================================================================================================


class Face:
    """Face class.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, vs):
        """Create new face.

        Parameters
        ----------
        vs : list(Vertex)
            List of vertices.
        """

        self.Vertices = vs
        self.Edges = []

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
        self.Faces = []

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
        Vertex | None
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
        Vertex
            Found or new vertex.
        """

        v = self.find_vertex(p)

        if v is None:
            v = self.new_vertex(p=p)

        return v

    # ----------------------------------------------------------------------------------------------

    def new_edge(self, a, b):
        """Create new edge.

        Parameters
        ----------
        a : Vertex
            First vertex.
        b : Vertex
            Second vertex.

        Returns
        -------
        Edge
            New edge.
        """

        e = Edge(a, b)
        a.Edges.append(e)
        b.Edges.append(e)
        self.Edges.append(e)

        return e

    # ----------------------------------------------------------------------------------------------

    def find_edge(self, a, b):
        """Find vertex.

        Parameters
        ----------
        a : Vertex
            First vertex.
        b : Vertex
            Second vertex.

        Returns
        -------
        Edge | None
            Found edge or None.
        """

        for e in a.Edges:
            if a.neighbour(e) == b:
                return e

        return None

    # ----------------------------------------------------------------------------------------------

    def find_or_new_edge(self, a, b):
        """Find or new edge.

        Parameters
        ----------
        a : Vertex
            First vertex.
        b : Vertex
            Second vertex.

        Returns
        -------
        Edge
            Found or new edge.
        """

        e = self.find_edge(a, b)

        if e is None:
            e = self.new_edge(a, b)

        return e

    # ----------------------------------------------------------------------------------------------

    def new_face(self, vs):
        """Add new face based on vertices list.

        Parameters
        ----------
        vs : list(Vertex)
            Vertices list.

        Returns
        -------
        Face
            New face.
        """

        f = Face(vs)

        # Add Vertex-Face links.
        for v in vs:
            v.Faces.append(f)

        # Add edges.
        for i in range(len(vs)):
            a, b = vs[i], vs[(i + 1) % len(vs)]
            e = self.find_or_new_edge(a, b)
            e.Faces.append(f)
            f.Edges.append(e)

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

    def load_dat_ecg(self, filename):
        """Load ECG graph from dat format.

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

    # ----------------------------------------------------------------------------------------------

    def save_dat_ecg(self, filename):
        """Save ECG graph in dat format.

        Parameters
        ----------
        filename : str
            Name of file.
        """

        # Aux index.
        for (i, v) in enumerate(self.Vertices):
            v.Id = i + 1

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
                id0, id1 = e.Vertices[0].Id, e.Vertices[1].Id
                f.write(f'{id0} {id1} {id1}\n')

            f.close()

    # ----------------------------------------------------------------------------------------------

    def load_dat_mesh(self, filename):
        """Load surface mesh from dat file.

        Parameters
        ----------
        filename : str
            Name of file.
        """

        self.clear()

        with open(filename, 'r') as f:

            # Variables info.
            vc = len(f.readline().split('=')[1].split())

            # Zone. Ingore it - we work with one zone.
            f.readline()

            # Nodes.
            nc = int(f.readline().split('=')[1])

            # Elements.
            ec = int(f.readline().split('=')[1])

            # Datapacking, zonetype, varlocation. Ignore it.
            f.readline()
            f.readline()
            f.readline()

            # Read coordinates of nodes.
            xs = [float(x) for x in f.readline().split()]
            ys = [float(y) for y in f.readline().split()]
            zs = [float(z) for z in f.readline().split()]
            assert (len(xs) == nc) and (len(ys) == nc) and (len(zs) == nc)

            # Add nodes.
            for i in range(nc):
                self.new_vertex(geom.Vector(xs[i], ys[i], zs[i]))

            # Ignore extra variables.
            for _ in range(vc - 3):
                f.readline()

            # Read links.
            for _ in range(ec):
                link = f.readline().split()
                idxs = [int(link[i]) for i in range(3)]
                self.new_face([self.Vertices[idx - 1] for idx in idxs])

            f.close()

    # ----------------------------------------------------------------------------------------------

    def construct_and_show_networkx_graph(self):
        """Construct NetworkX graph and show it.
        """

        # Aux index.
        for (i, v) in enumerate(self.Vertices):
            v.Id = i

        # Create graph.
        g = nx.Graph()

        # We need to collect list of positions.
        positions = []

        # Add nodes.
        for (i, v) in enumerate(self.Vertices):
            pos = (v.P.X, v.P.Y)
            positions.append(pos)
            g.add_node(i, pos=pos)

        # Add edges.
        for e in self.Edges:
            g.add_edge(e.A.Id, e.B.Id)

        nx.draw(g, pos=positions)

    # ----------------------------------------------------------------------------------------------

    def decolor_vertices(self):
        """Reset colors.
        """

        for v in self.Vertices:
            v.Color = 0

    # ----------------------------------------------------------------------------------------------

    def vertices_of_color(self, color):
        """Get list of vertices of a given color.

        Parameters
        ----------
        color : int
            Color.

        Returns
        -------
        list
            List of vertices.
        """

        return filter(lambda v: v.Color == color, self.Vertices)

    # ----------------------------------------------------------------------------------------------

    def calculate_max_vertex_color(self):
        """Calculate max color.

        Returns
        -------
        int
            Max color.

        """

        return max([v.Color for v in self.Vertices])

    # ----------------------------------------------------------------------------------------------

    def is_vertex_coloring_correct(self):
        """Check coloring for correctness.

        Returns
        -------
        bool
            Is coloring is correct.
        """

        if len(self.vertices_of_color(0)) > 0:
            return False

        if len(filter(lambda e: e.A.Color == e.B.Color, self.Edges)) > 0:
            return False

        return True

    # ----------------------------------------------------------------------------------------------

    def vertex_coloring_greedy(self, verbose=True):
        """Greedy vertex coloring algorithm.

        Trying to color vertex by vertex in greedy manner.

        Parameters
        ----------
        verbose : bool
            Verbose mode.

        Returns
        -------
        int
            Colors count.
        """

        self.decolor_vertices()
        max_c = 0

        for v in self.Vertices:
            cs = v.neighbours_colors()
            c = 1
            while c in cs:
                c += 1
            max_c = max(max_c, c)
            v.Color = c

        if verbose:
            print(f'Graph.coloring_greedy : {max_c} colors')

        return max_c

# ----------------------------------------------------------------------------------------------

    def vertex_coloring_greedy2(self, verbose=True):
        """Greedy vertex coloring algorithm.

        First try to color maxinum number of vertices with color 1, 2, and so on..

        Parameters
        ----------
        verbose : bool
            Verbose mode.

        Returns
        -------
        int
            Colors count.
        """

        self.decolor_vertices()
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

        if verbose:
            print(f'Graph.coloring_greedy2 : {max_c} colors')

        return max_c

# ----------------------------------------------------------------------------------------------

    def vertex_coloring_recolor5(self, verbose=True):
        """Vertex coloring with trying to recolor vertices of color 5.

        Parameters
        ----------
        verbose : bool
            Verbose mode.

        Returns
        -------
        int
            Colors count.
        """

        c = self.vertex_coloring_greedy(verbose=False)

        # This coloring doesn't work with more than 5 colors.
        assert c <= 5

        if c == 5:
            vs = self.vertices_of_color(5)

            # Try to recolor all color 5 vertices.
            while len(vs) > 0:
                v = vs[0]
                neigh_cs = v.neighbours_colors()
                is_recolored = False

                # Try to recolor this vertex.
                while not is_recolored:

                    # Try to recolor in 1, 2, 3, 4.
                    for new_c in range(1, 5):
                        if new_c not in neigh_cs:
                            v.Color = new_c
                            # print(f'vertex {v.Id} is recolored')
                            is_recolored = True
                            break

                    # Move to first direction.
                    if not is_recolored:
                        goto_v = v.neighbour(v.Edges[random.randint(0, 3)])
                        v.Color, goto_v.Color = goto_v.Color, v.Color
                        v = goto_v
                        neigh_cs = v.neighbours_colors()

                # Get color 5 vertices again.
                vs = self.vertices_of_color(5)

            # Recalc max color.
            assert self.is_vertex_coloring_correct()
            c = self.calculate_max_color()

        if verbose:
            print(f'Graph.coloring_recolor5 : {c} colors')

        return c

# ==================================================================================================


if __name__ == '__main__':

    # Example for vertex coloring.
    # g = Graph()
    # g.load_dat_ecg('../data/dat/ecg/wing_1_ecg.dat')
    # g.vertex_coloring_recolor5()
    # g.save_dat_ecg('t.dat')

    # Load graph from dat file.
    g = Graph()
    g.load_dat_mesh('../data/dat/meshes/wing_1.dat')
    g.construct_and_show_networkx_graph()

# ==================================================================================================
