"""
Graph realization.
"""

import itertools
import random
import plotly.graph_objects as go
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
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
        self.Parent = None

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

        a, b = e.A, e.B

        if self == a:
            return b
        elif self == b:
            return a
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

        self.Color = 0
        self.Vertices = [v0, v1]
        self.Faces = []
        self.Parent = None

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

    # ----------------------------------------------------------------------------------------------

    def is_adjacent(self, e):
        """Check if edge adjacent to another edge.

        Parameters
        ----------
        a : Edge
            Edge.

        Returns
        -------
        bool
            True - if adjacent, False - otherwise.
        """

        a, b = self.A, self.B

        return (a == e.A) or (a == e.B) or (b == e.A) or (b == e.B)

    # ----------------------------------------------------------------------------------------------

    def split(self):
        """Split edge into 2 edges.

        Returns
        -------
        Vertex
            New vertex.
        """

        old_b = self.B
        new_v = self.Parent.new_vertex(p=0.5 * (self.A.P + self.B.P))

        # Replace second vertex of current edge with new_vertex.
        self.Vertices[1] = new_v
        new_v.Edges.append(self)

        # Remove link from old_b to this edges.
        old_b.Edges.remove(self)

        # Add new edge.
        new_e = self.Parent.new_edge(new_v, old_b)

        # Link faces with new objects.
        for f in self.Faces:
            f.append_vertex_between(new_v, self.A, old_b)
            f.append_edge(new_e)
            new_v.Faces.append(f)
            new_e.Faces.append(f)

        return new_v

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
        self.Parent = None

    # ----------------------------------------------------------------------------------------------

    def is_even(self):
        """Check if face even.

        Returns
        -------
        bool
            True - if face is even, False - otherwise.
        """

        return len(self.Edges) % 2 == 0

    # ----------------------------------------------------------------------------------------------

    def is_odd(self):
        """Check if face odd.

        Returns
        -------
            True - if face is odd, False - otherwise.
        """

        return not self.is_even()

    # ----------------------------------------------------------------------------------------------

    def neighbour(self, e):
        """Get neighbour.

        Parameters
        ----------
        e : Edge
            Edge to find neighbour.

        Returns
        -------
        Face
            Neighbour or None.
        """

        if len(e.Faces) == 1:
            return None

        f1, f2 = e.Faces[0], e.Faces[1]

        if self == f1:
            return f2
        elif self == f2:
            return f1
        else:
            raise Exception('Face.neighbour : internal error')

    # ----------------------------------------------------------------------------------------------

    def append_vertex_between(self, new_v, a, b):
        """Append vertex between two that already in face.

        Parameters
        ----------
        new_v : Vertex
            New vertex.
        a : Vertex
            First vertex.
        b : Vertex
            Second vertex.
        """

        ai = self.Vertices.index(a)
        bi = self.Vertices.index(b)
        min_i = min(ai, bi)
        max_i = max(ai, bi)
        if (min_i == 0) and (max_i == len(self.Vertices) - 1):
            self.Vertices.append(new_v)
        else:
            self.Vertices.insert(max_i, new_v)

    # ----------------------------------------------------------------------------------------------

    def append_edge(self, e):
        """Append edge while keeping order of traversal edges.

        Parameters
        ----------
        e : Edge
            Edge to append.
        """

        if e.is_adjacent(self.Edges[0]) and e.is_adjacent(self.Edges[-1]):
            # Our edge is between first and last, so we can just append it.
            self.Edges.append(e)
            return
        else:
            for i in range(1, len(self.Edges)):
                if e.is_adjacent(self.Edges[i - 1]) and e.is_adjacent(self.Edges[i]):
                    self.Edges.insert(i, e)
                    return

        raise Exception('internal error')

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
        v.Parent = self
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
        e.Parent = self
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

    def remove_edge(self, e):
        """Remove edge.

        Parameters
        ----------
        e : Edge
            Edge.
        """

        # Removing edge from graph with faces requires special corrections,
        # so don't do it in a common manner.
        assert self.Faces == []

        self.Edges.remove(e)

        for v in self.Vertices:
            if e in v.Edges:
                v.Edges.remove(e)

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
        f.Parent = self
        self.Faces.append(f)

        # Add Vertex-Face links.
        for v in vs:
            v.Faces.append(f)

        # Add edges.
        for i in range(len(vs)):
            a, b = vs[i], vs[(i + 1) % len(vs)]
            e = self.find_or_new_edge(a, b)
            e.Faces.append(f)
            f.Edges.append(e)

        return f

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

    def init_2d_rect_mesh(self, xn, yn):
        """Init 2 dimensional mesh graph.

        Parameters
        ----------
        xn : int
            Number of cells along x direction.
        yn : int
            Number of cells along y direction.
        """

        self.clear()

        for xi in range(xn):
            for yi in range(yn):
                ul = self.find_or_new_vertex(geom.Vector(xi, yi, 0.0))
                ur = self.find_or_new_vertex(geom.Vector(xi + 1.0, yi, 0.0))
                dl = self.find_or_new_vertex(geom.Vector(xi, yi + 1.0, 0.0))
                dr = self.find_or_new_vertex(geom.Vector(xi + 1.0, yi + 1.0, 0.0))
                self.find_or_new_edge(ul, ur)
                self.find_or_new_edge(ul, dl)
                self.find_or_new_edge(ur, dr)
                self.find_or_new_edge(dl, dr)

    # ----------------------------------------------------------------------------------------------

    def init_3d_rect_mesh(self, xn, yn, zn):
        """Init 3 dimensional mesh graph.

        Parameters
        ----------
        xn : int
            Number of cells along x direction.
        yn : int
            Number of cells along y direction.
        zn : int
            Number of cells along z direction.
        """

        self.clear()

        for xi in range(xn):
            for yi in range(yn):
                for zi in range(zn):
                    ulf = self.find_or_new_vertex(geom.Vector(xi + 0.1 * zi,
                                                              yi + 0.1 * zi, 0.0))
                    urf = self.find_or_new_vertex(geom.Vector(xi + 1.0 + 0.1 * zi,
                                                              yi + 0.1 * zi, 0.0))
                    dlf = self.find_or_new_vertex(geom.Vector(xi + 0.1 * zi,
                                                              yi + 1.0 + 0.1 * zi, 0.0))
                    drf = self.find_or_new_vertex(geom.Vector(xi + 1.0 + 0.1 * zi,
                                                              yi + 1.0 + 0.1 * zi, 0.0))
                    ulb = self.find_or_new_vertex(geom.Vector(xi + 0.1 * (zi + 1),
                                                              yi + 0.1 * (zi + 1), 0.0))
                    urb = self.find_or_new_vertex(geom.Vector(xi + 1.0 + 0.1 * (zi + 1),
                                                              yi + 0.1 * (zi + 1), 0.0))
                    dlb = self.find_or_new_vertex(geom.Vector(xi + 0.1 * (zi + 1),
                                                              yi + 1.0 + 0.1 * (zi + 1), 0.0))
                    drb = self.find_or_new_vertex(geom.Vector(xi + 1.0 + 0.1 * (zi + 1),
                                                              yi + 1.0 + 0.1 * (zi + 1), 0.0))
                    self.find_or_new_edge(ulf, urf)
                    self.find_or_new_edge(ulf, dlf)
                    self.find_or_new_edge(urf, drf)
                    self.find_or_new_edge(dlf, drf)
                    self.find_or_new_edge(ulb, urb)
                    self.find_or_new_edge(ulb, dlb)
                    self.find_or_new_edge(urb, drb)
                    self.find_or_new_edge(dlb, drb)
                    self.find_or_new_edge(ulf, ulb)
                    self.find_or_new_edge(urf, urb)
                    self.find_or_new_edge(dlf, dlb)
                    self.find_or_new_edge(drf, drb)

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

    def construct_and_show_networkx_graph(self, is_draw_labels=False):
        """Construct NetworkX graph and show it.

        Parameters
        ----------
        is_draw_labels : bool
            Is needed to draw labels.
        """

        figure(figsize=(8, 4))
        plt.axes([0, 0, 1, 1])

        # Aux index.
        for (i, v) in enumerate(self.Vertices):
            v.Id = i

        # Sort edges.
        self.Edges.sort(key=lambda e: 1000 * min(e.A.Id, e.B.Id) + max(e.A.Id, e.B.Id))

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

        vertex_colors = ['#000000', '#ffffff', '#555555', '#dddddd', '#999999']
        vertex_shapes = ['o', 'h', '^', 'v', 's']

        # Draw nodes.
        # If no vertex coloring found then color it in neutral manner.
        if self.min_vertex_color() == 0:
            if is_draw_labels:
                node_size = 200
            else:
                node_size = 20
            nx.draw_networkx_nodes(g, pos=positions, node_color='black', node_size=node_size)
        else:
            for c in range(self.max_vertex_color()):
                c = c + 1
                vertices = filter(lambda v: v.Color == c, self.Vertices)
                idxs = [v.Id for v in vertices]
                node_color = vertex_colors[(c - 1) % len(vertex_colors)]
                node_shape = vertex_shapes[(c - 1) % len(vertex_shapes)]
                nx.draw_networkx_nodes(g,
                                       pos=positions,
                                       nodelist=idxs,
                                       node_color=node_color,
                                       node_shape=node_shape,
                                       node_size=120, edgecolors='black')

        # Draw edges.
        # If no coloring color it in neutral color.
        if self.min_edge_color() == 0:
            nx.draw_networkx_edges(g, pos=positions, edge_color='silver', style=':', width=1)
        else:
            edges_colors = ['black']
            edges_colors_map = [edges_colors[(e.Color - 1) % len(edges_colors)] for e in self.Edges]
            edges_styles = ['-', '--', '-.', ':', '-', '--']
            edges_styles_map = [edges_styles[(e.Color - 1) % len(edges_styles)] for e in self.Edges]
            edges_widths = [3, 2, 1, 3, 1, 1]
            edges_widths_map = [edges_widths[(e.Color - 1) % len(edges_widths)] for e in self.Edges]
            nx.draw_networkx_edges(g,
                                   pos=positions,
                                   edge_color=edges_colors_map,
                                   style=edges_styles_map,
                                   width=edges_widths_map)

        # Draw vertices labels.
        if is_draw_labels:
            nx.draw_networkx_labels(g, pos=positions,
                                    font_color='white', font_size=11, font_weight='bold')

        plt.axis('off')

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

        return [v for v in self.Vertices if v.Color == color]

    # ----------------------------------------------------------------------------------------------

    def max_vertex_color(self):
        """Calculate max color.

        Returns
        -------
        int
            Max color.
        """

        return max([v.Color for v in self.Vertices])

    # ----------------------------------------------------------------------------------------------

    def min_vertex_color(self):
        """Calculate min color.

        Returns
        -------
        int
            Min color.
        """

        return min([v.Color for v in self.Vertices])

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

        if len([e for e in self.Edges if e.A.Color == e.B.Color]) > 0:
            return False

        return True

    # ----------------------------------------------------------------------------------------------

    def vertex_coloring_random(self, n, verbose=True):
        """Random coloring of vertices.

        Parameters
        ----------
        n : int
            Count of colors.
        verbose : bool
            Verbose mode.

        Returns
        -------
        int
            Colors count.
        """

        for v in self.Vertices:
            v.Color = random.randint(1, n)

        max_c = self.max_vertex_color()

        if verbose:
            print(f'Graph.vertex_coloring_random : {max_c} colors')

        return max_c

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
            print(f'Graph.vertex_coloring_greedy : {max_c} colors')

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
            print(f'Graph.vertex_coloring_greedy2 : {max_c} colors')

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
            c = self.max_vertex_color()

        if verbose:
            print(f'Graph.vertex_coloring_recolor5 : {c} colors')

        return c

    # ----------------------------------------------------------------------------------------------

    def decolor_edges(self):
        """Reset colors for edges.
        """

        for e in self.Edges:
            e.Color = 0

    # ----------------------------------------------------------------------------------------------

    def edges_of_color(self, color):
        """Get list of edges of a given color.

        Parameters
        ----------
        color : int
            Color.

        Returns
        -------
        list
            List of edges.
        """

        return [e for e in self.Edges if e.Color == color]

    # ----------------------------------------------------------------------------------------------

    def max_edge_color(self):
        """Calculate max edge color.

        Returns
        -------
        int
            Max edge color.
        """

        return max([e.Color for e in self.Edges])

    # ----------------------------------------------------------------------------------------------

    def min_edge_color(self):
        """Calculate min edge color.

        Returns
        -------
        int
            Min edge color.
        """

        return min([e.Color for e in self.Edges])

    # ----------------------------------------------------------------------------------------------

    def is_edge_coloring_correct(self):
        """Check coloring for correctness.

        Returns
        -------
        bool
            Is coloring is correct.
        """

        if len(self.edges_of_color(0)) > 0:
            return False

        for v in self.Vertices:
            for e1 in v.Edges:
                for e2 in v.Edges:
                    if e1 != e2:
                        if e1.Color == e2.Color:
                            return False

        return True

    # ----------------------------------------------------------------------------------------------

    def edge_coloring_random(self, n, verbose=True):
        """Random coloring.

        Parameters
        ----------
        n : int
            Count of colors.
        verbose : bool
            Verbose mode.

        Returns
        -------
        int
            Colors count.
        """

        for e in self.Edges:
            e.Color = random.randint(1, n)

        max_c = self.max_edge_color()

        if verbose:
            print(f'Graph.edge_coloring_random : {max_c} colors')

        return max_c

    # ----------------------------------------------------------------------------------------------

    def edge_coloring_greedy(self, verbose=True):
        """Greedy edge coloring.

        Returns
        -------
        int
            Colors count.
        """

        self.decolor_edges()
        max_c = 0

        for e in self.Edges:
            if e.Color == 0:
                a, b = e.A, e.B
                cs = [e.Color for e in a.Edges + b.Edges]
                print(cs)
                c = 1
                while c in cs:
                    c += 1
                e.Color = c

        max_c = self.max_edge_color()

        if verbose:
            print(f'Graph.edge_coloring_greedy : {max_c} colors')

        return max_c

    # ----------------------------------------------------------------------------------------------

    def edge_coloring_even_loops_3colors(self, verbose=True):
        """Coloring for graph with only even loops

        Parameters
        ----------
        verbose : bool
            Verbose mode.

        Returns
        -------
        int
            Colors count.
        """

        self.decolor_edges()

        for f in self.Faces:
            assert f.is_even()

        # Set colors for any loop.
        for (i, e) in enumerate(self.Faces[0].Edges):
            e.Color = (i % 2) + 1

        def neighbours_without_color(v):
            colors0 = [e.Color for e in v.Edges if e.Color == 0]
            return len(colors0)

        # Color while all edges are colored.
        while self.min_edge_color() == 0:
            start_vertices = [v for v in self.Vertices if neighbours_without_color(v) == 1]
            start_vertex = start_vertices[0]
            edge0 = [e for e in start_vertex.Edges if e.Color == 0][0]
            edgen = [e for e in start_vertex.Edges if e.Color != 0][0]
            new_color = 6 - sum([e.Color for e in start_vertex.Edges])
            f_common = None
            for f0 in edge0.Faces:
                for fn in edgen.Faces:
                    if f0 == fn:
                        f_common = f0
            idx = f_common.Edges.index(edgen)
            for (i, e) in enumerate(f_common.Edges):
                if (i - idx) % 2 == 0:
                    e.Color = edgen.Color
                else:
                    e.Color = new_color

        max_c = self.max_edge_color()

        if verbose:
            print(f'Graph.edge_coloring_even_loops_3colors : {max_c} colors')

        return max_c

    # ----------------------------------------------------------------------------------------------

    def new_aux_edge(self, a, b):
        """New aux edge.
        """

        # Find common face.
        f_common = None
        for fa in a.Faces:
            for fb in b.Faces:
                if fa == fb:
                    f_common = fa

        # Add new face.
        f_new = self.new_face([])

        # Find all vertices between a and b in common_face and put them into new face.
        ai = f_common.Vertices.index(a)
        bi = f_common.Vertices.index(b)
        min_i = min(ai, bi)
        max_i = max(ai, bi)
        v_cnt_all = len(f_common.Vertices)
        v_cnt_new = max_i - min_i + 1
        v_cnt_com = v_cnt_all - v_cnt_new + 2
        v_new = []
        for i in range(v_cnt_new):
            ind_of_new = min_i + i
            v_new.append(f_common.Vertices[ind_of_new])
        v_com = []
        for i in range(v_cnt_com):
            ind_of_com = (max_i + i) % v_cnt_all
            v_com.append(f_common.Vertices[ind_of_com])

        # a, b get links to f_new
        a.Faces.append(f_new)
        b.Faces.append(f_new)
        # links of new vertices to new face
        for v in v_new[1:-1]:
            v.Faces.remove(f_common)
            v.Faces.append(f_new)
        f_common.Vertices = v_com
        f_new.Vertices = v_new

        # Correct links of edges.
        for i in range(len(v_new) - 1):
            e = self.find_edge(v_new[i], v_new[i + 1])
            if e in f_common.Edges:
                f_common.Edges.remove(e)
            if f_common in e.Faces:
                e.Faces.remove(f_common)
            f_new.Edges.append(e)
            e.Faces.append(f_new)

        # Add new edge.
        new_e = self.new_edge(a, b)
        f_common.append_edge(new_e)
        f_new.append_edge(new_e)
        new_e.Faces = [f_common, f_new]

    # ----------------------------------------------------------------------------------------------

    def add_edges_for_odd_faces_elimination(self):
        """Add new edges to eliminate add faces
        """

        # Hardcode.
        a = [[(2, 8), (3, 16)], [(2, 3), (20, 21)], [(4, 5), (0, 1)], [(0, 5), (24, 25)]]
        for i in range(len(a)):
            a0 = a[i][0]
            a1 = a[i][1]
            e1 = self.find_edge(self.Vertices[a0[0]], self.Vertices[a0[1]])
            e2 = self.find_edge(self.Vertices[a1[0]], self.Vertices[a1[1]])
            new_a = e1.split()
            new_b = e2.split()
            self.new_aux_edge(new_a, new_b)

    # ----------------------------------------------------------------------------------------------

    def delete_aux_edge(self, new_color_a, new_color_b):
        """Delete aux edges.
        """

        # Delete faces - we don't need them.
        self.Faces = []
        for v in self.Vertices:
            v.Faces = []
        for e in self.Edges:
            e.Faces = []

        # Hardcode.
        a, b = self.Vertices[-1], self.Vertices[-2]
        ab_e = self.find_edge(a, b)
        a_edges = [ae for ae in a.Edges if ae != ab_e]
        b_edges = [be for be in b.Edges if be != ab_e]
        a_nbs = [a.neighbour(ae) for ae in a_edges]
        b_nbs = [b.neighbour(be) for be in b_edges]
        # remove edges.
        for e in [ab_e] + a_edges + b_edges:
            self.remove_edge(e)
        # remove vertices.
        self.Vertices.remove(a)
        self.Vertices.remove(b)
        # add needed edges back.
        a_new = self.new_edge(a_nbs[0], a_nbs[1])
        a_new.Color = new_color_a
        b_new = self.new_edge(b_nbs[0], b_nbs[1])
        b_new.Color = new_color_b

    # ----------------------------------------------------------------------------------------------

    def rotate_colors(self, el):
        """Rotate colors of edges.

        Parameters
        ----------
        el : list(Edge)
            List of edges.
        """

        edges = [g.find_edge(self.Vertices[ei[0]], self.Vertices[ei[1]]) for ei in el]
        c0, c1 = edges[0].Color, edges[1].Color

        for (i, e) in enumerate(edges):
            if i % 2 == 1:
                e.Color = c0
            else:
                e.Color = c1

# ==================================================================================================


if __name__ == '__main__':

    test = 1

    if test == 1:

        # Example for vertex coloring.
        g = Graph()
        g.load_dat_ecg('../data/dat/ecg/bunny_ecg.dat')
        g.vertex_coloring_recolor5()
        g.save_dat_ecg('t.dat')

    elif test == 2:

        # Load graph from dat file.
        g = Graph()
        cs = \
            [
                [0.0, 10.0, 0.0],
                [5.0, 5.0, 0.0],
                [5.0, -5.0, 0.0],
                [0.0, -10.0, 0.0],
                [-5.0, -5.0, 0.0],
                [-5.0, 5.0, 0.0],
                [10.0, 10.0, 0.0],
                [15.0, 0.0, 0.0],
                [10.0, -10.0, 0.0],
                [-10.0, 10.0, 0.0],
                [-15.0, 0.0, 0.0],
                [-10.0, -10.0, 0.0],
                [0.0, 20.0, 0.0],
                [15.0, 15.0, 0.0],
                [20.0, 0.0, 0.0],
                [15.0, -15.0, 0.0],
                [0.0, -20.0, 0.0],
                [-15.0, -15.0, 0.0],
                [-20.0, 0.0, 0.0],
                [-15.0, 15.0, 0.0]
            ]
        for c in cs:
            g.new_vertex(geom.Vector(c[0], c[1], c[2]))
        fs = \
            [
                [0, 1, 2, 3, 4, 5],
                [1, 6, 7, 8, 2],
                [5, 9, 10, 11, 4],
                [12, 13, 6, 1, 0],
                [13, 14, 7, 6],
                [7, 14, 15, 8],
                [15, 8, 2, 3, 16],
                [16, 3, 4, 11, 17],
                [17, 11, 10, 18],
                [18, 10, 9, 19],
                [19, 12, 0, 5, 9],
                [12, 13, 14, 15, 16, 17, 18, 19]
            ]
        for f in fs:
            g.new_face([g.Vertices[i] for i in f])
        #
        g.add_edges_for_odd_faces_elimination()
        g.edge_coloring_even_loops_3colors()
        assert g.is_edge_coloring_correct()
        #
        g.rotate_colors([(5, 26), (26, 27), (24, 27), (5, 24)])
        g.delete_aux_edge(3, 3)
        g.rotate_colors([(5, 9), (9, 10), (10, 11), (11, 4), (4, 24)])
        assert g.is_edge_coloring_correct()
        #
        g.rotate_colors([(0, 25), (25, 24), (24, 4), (4, 11), (11, 10), (10, 9), (9, 5), (5, 0)])
        g.delete_aux_edge(1, 3)
        g.rotate_colors([(0, 12), (12, 19), (19, 9), (9, 5)])
        assert g.is_edge_coloring_correct()
        #
        g.rotate_colors([(2, 22), (22, 23), (23, 20), (20, 2)])
        g.delete_aux_edge(3, 3)
        g.rotate_colors([(2, 1), (1, 6), (6, 7), (7, 8), (8, 20)])
        assert g.is_edge_coloring_correct()
        #
        g.rotate_colors([(20, 21), (21, 3), (3, 2), (2, 1), (1, 6), (6, 7), (7, 8), (8, 20)])
        g.delete_aux_edge(3, 1)
        g.rotate_colors([(16, 15), (15, 8)])
        assert g.is_edge_coloring_correct()
        #
        g.construct_and_show_networkx_graph()

    elif test == 3:

        g = Graph()
        g.init_2d_rect_mesh(4, 3)
        g.edge_coloring_greedy()
        g.construct_and_show_networkx_graph()

    elif test == 4:

        g = Graph()
        g.init_3d_rect_mesh(4, 3, 3)
        g.edge_coloring_greedy()
        g.construct_and_show_networkx_graph()

# ==================================================================================================
