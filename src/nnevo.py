# -*- coding: utf-8 -*-
"""
Neuroevolution system (NeuroNet EVOlution).

Created on Fri Jun 21 11:35:26 2019

@author: Rybakov
"""

import lst
import mth
import fun
import time

#---------------------------------------------------------------------------------------------------
# Class settings.
#---------------------------------------------------------------------------------------------------

class Settings:

    # Default weight.
    DefaultEdgeWeight = 0.01

    # Default bias.
    DefaultNodeBias = 0.0

#---------------------------------------------------------------------------------------------------
# Class Node (neuron).
#---------------------------------------------------------------------------------------------------

class Node:

    def __init__(self):
        """
        Constructor.
        """

        self.Id = None
        self.IEdges = []
        self.OEdges = []
        self.Signals = []
        self.SavedSignals = []
        self.Errors = []
        self.Bias = Settings.DefaultNodeBias
        self.Z = None
        self.A = None
        self.E = None
        self.Mark = False

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return 'Node %s : Sg/Er/B/Z/A/E = %s/%s/%s/%s/%s/%s' % (self.Id,
                                                                self.Signals,
                                                                self.Errors,
                                                                self.Bias,
                                                                self.Z,
                                                                self.A,
                                                                self.E)

#---------------------------------------------------------------------------------------------------

    def IsSignalsReady(self):
        """
        Check if signals are ready.

        Result:
            True - if signals are ready,
            Faslse - otherwise.
        """

        return fun.is_all(self.Signals, lambda x: x != None)

#---------------------------------------------------------------------------------------------------

    def IsErrorsReady(self):
        """
        Check if error are ready.

        Result:
            True - if errors  are ready,
            False - otherwise.
        """

        return fun.is_all(self.Errors, lambda x: x != None)

#---------------------------------------------------------------------------------------------------

    def ForwardPropagation(self):
        """
        Forward propagation of the signal.
        """

        # Calculate out signal.
        self.Z = sum(fun.zipwith(self.Signals,
                                 self.IEdges,
                                 lambda s, e: s * e.Weight))
        self.A = mth.sigmoid(self.Z)

        # Propagate.
        for oe in self.OEdges:
            oe.Dst.Signals[oe.IIndex] = self.A

        # Save signals.
        self.SavedSignals = self.Signals
        self.Signals = [None] * len(self.SavedSignals)

#---------------------------------------------------------------------------------------------------

    def BackPropagation(self):
        """
        Back propagation of the error.
        """

        # Calculate total error.
        self.E = sum(self.Errors) * self.A * (1.0 - self.A)

        # Propagate.
        for ie in self.IEdges:
            ie.Src.Errors[ie.OIndex] = self.E * ie.Weight

#---------------------------------------------------------------------------------------------------
# Class Edge.
#---------------------------------------------------------------------------------------------------

class Edge:

    def __init__(self):
        """
        Constructor.
        """

        self.Id = None
        self.Src = None
        self.Dst = None
        self.OIndex = None
        self.IIndex = None
        self.Weight = Settings.DefaultEdgeWeight

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return 'Edge %s : [n%s/o%s -> n%s/i%s] W = %s' % (self.Id,
                                                          self.Src.Id,
                                                          self.OIndex,
                                                          self.Dst.Id,
                                                          self.IIndex,
                                                          self.Weight)

#---------------------------------------------------------------------------------------------------
# Class Net.
#---------------------------------------------------------------------------------------------------

class Net:

    def __init__(self):
        """
        Constructor.
        """

        self.Nodes = []
        self.FirstLayer = []
        self.LastLayer = []
        self.Edges = []

#---------------------------------------------------------------------------------------------------

    def NodesCount(self):
        """
        Count of nodes.

        Result:
            Count of nodes.
        """

        return len(self.Nodes)

#---------------------------------------------------------------------------------------------------

    def EdgesCount(self):
        """
        Count of edges.

        Result:
            Count of edges.
        """

        return len(self.Edges)

#---------------------------------------------------------------------------------------------------

    def CreateMultilayer(self, LayersSizes):
        """
        Create multilayer net.

        Arguments:
            LayersSizes -- layers sizes.
        """

        # Check the layers count.
        if len(LayersSizes) < 2:
            raise Exception('not enough layers (must be >= 2)')

        # Create all nodes.
        nodes = [[Node() for i in range(LayerSize)] for LayerSize in LayersSizes]

        # Set links to nodes.
        self.Nodes = lst.flatten(nodes)
        self.FirstLayer = nodes[0]
        self.LastLayer = lst.last(nodes)

        # First layer has signals.
        for node in self.FirstLayer:
            node.Signals = [None]

        # Last layer has errors.
        for node in self.LastLayer:
            node.Errors = [None]

        # Nodes ids.
        for i in range(len(self.Nodes)):
            self.Nodes[i].Id = i

        # Create edges.
        eid = 0
        for i in range(len(nodes) - 1):

            src_layer = nodes[i]
            dst_layer = nodes[i + 1]

            # Edges between i-th and (i + 1)-th layers.
            for src in src_layer:
                for dst in dst_layer:

                    e = Edge()
                    e.Id = eid
                    eid = eid + 1

                    # Links.
                    e.Src = src
                    e.Dst = dst
                    src.OEdges.append(e);
                    src.Errors.append(None);
                    dst.IEdges.append(e);
                    dst.Signals.append(None);
                    self.Edges.append(e);

        # Correct indices.
        for node in self.Nodes:
            for i in range(len(node.IEdges)):
                node.IEdges[i].IIndex = i
            for i in range(len(node.OEdges)):
                node.OEdges[i].OIndex = i

#---------------------------------------------------------------------------------------------------

    def Print(self):
        """
        Print neuronet.
        """

        print('Nodes (%d):' % self.NodesCount())
        for node in self.Nodes:
            print(str(node))

        print('Edges (%d)' % self.EdgesCount())
        for edge in self.Edges:
            print(str(edge))

#---------------------------------------------------------------------------------------------------

    def SenseForward(self, x):
        """
        Sense neuronet forward.

        Arguments:
            x -- signal.

        Result:
            Response.
        """

        if len(x) != len(self.FirstLayer):
            raise Exception('wrong input signal size')

        # Clean old data.
        for node in self.Nodes:
            node.Z = None
            node.A = None
            node.E = None
            node.Signals = [None] * len(node.Signals)
            node.Errors = [None] * len(node.Errors)
            node.Mark = False

        # Propagate forward.
        for i in range(len(x)):
            node = self.FirstLayer[i]
            node.Signals[0] = x[i]
            node.Mark = True
            
        traversal = self.FirstLayer.copy()

        # Propagate signal for all neuronet.
        while traversal != []:

            # Split list.
            h, t = traversal[0], traversal[1:]

            # Check for signals ready.
            if not h.IsSignalsReady():
                raise Exception('no signal on neuron')

            # Forward propagation for head.
            h.ForwardPropagation()

            # Add all successors of head to the end of tail.
            for oe in h.OEdges:
                succ = oe.Dst
                if not succ.Mark:
                    t.append(succ)
                    succ.Mark = True

            # Shift traversal.
            traversal = t

        return [node.A for node in self.LastLayer]

#---------------------------------------------------------------------------------------------------

    def A(self):
        """
        Result value.

        Result:
            Result.
        """

        return [node.A for node in self.LastLayer]

#---------------------------------------------------------------------------------------------------

    def Cost(self, y):
        """
        Calculate cost function.

        Arguments:
            y -- true result.

        Result:
            Cost function value.
        """

        return 0.5 * sum(fun.zipwith(y, self.A(), lambda a, b: (a - b) * (a - b)))

#---------------------------------------------------------------------------------------------------

    def SenseBack(self, y):
        """
        Sense neuronet back.

        Arguments:
            y -- right answer.
        """

        # Clean old data.
        for node in self.Nodes:
            node.Z = None
            #node.A = None
            node.E = None
            node.Signals = [None] * len(node.Signals)
            node.Errors = [None] * len(node.Errors)
            node.Mark = False

        # Propagate back.
        for i in range(len(self.LastLayer)):
            node = self.LastLayer[i]
            node.Errors[0] = y[i] - node.A
            node.Mark = True
            
        traversal = self.LastLayer.copy()

        # Propagate signal for all neuronet in back direction.
        while traversal != []:

            # Split list.
            h, t = traversal[0], traversal[1:]

            # Check for errors ready.
            if not h.IsErrorsReady():
                raise Exception('no signal on neuron')

            # Forward propagation for head.
            h.BackPropagation()

            # Add all successors of head to the end of tail.
            for ie in h.IEdges:
                pred = ie.Src
                if not pred.Mark:
                    t.append(pred)
                    pred.Mark = True

            # Shift traversal.
            traversal = t

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    net = Net()
    net.CreateMultilayer([800, 15, 10])
    t0 = time.clock()
    res = net.SenseForward([0.1] * 800)
    net.SenseBack([0.2] * 10)
    t1 = time.clock()
    print('res = %s, time = %s' % (str(res), t1 - t0))

#---------------------------------------------------------------------------------------------------
