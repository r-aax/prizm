"""
Graph realization.
"""

# ==================================================================================================


class Vertex:
    """Vertex class.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, p=(0.0, 0.0, 0.0), color=0):
        """Constructor.

        Parameters
        ----------
        p : (float, float, float)
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
        basestring
            Representation.
        """

        return f'V {self.P}, {self.Color}'

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
        basestring
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

    def new_vertex(self, p=(0.0, 0.0, 0.0), color=0):
        """Create new vertex.

        Parameters
        ----------
        p : (float, float, float)
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

    def save(self, filename):
        """Save graph.

        Parameters
        ----------
        filename : basestring
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
            f.write(' '.join([str(v.P[0]) for v in self.Vertices]) + '\n')
            f.write(' '.join([str(v.P[1]) for v in self.Vertices]) + '\n')
            f.write(' '.join([str(v.P[2]) for v in self.Vertices]) + '\n')
            f.write(' '.join([str(v.Color) for v in self.Vertices]) + '\n')

            # Edges.
            for e in self.Edges:
                id0, id1 = e.Vertices[0].Id + 1, e.Vertices[1].Id + 1
                f.write(f'{id0} {id1} {id1}\n')

            f.close()

# ==================================================================================================


if __name__ == '__main__':

    g = Graph()
    v0 = g.new_vertex(p=(0.0, 0.0, 0.0), color=0)
    v1 = g.new_vertex(p=(1.0, 0.0, 0.0), color=1)
    v2 = g.new_vertex(p=(0.0, 1.0, 0.0), color=2)
    e0 = g.new_edge(v0, v1)
    e1 = g.new_edge(v0, v2)
    e2 = g.new_edge(v1, v2)
    g.print()
    g.save('test_graph.dat')

# ==================================================================================================
