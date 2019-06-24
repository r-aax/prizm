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
import random

#---------------------------------------------------------------------------------------------------
# Class settings.
#---------------------------------------------------------------------------------------------------

class Settings:

    # Default weight.
    DefaultEdgeWeight = 0.01

    # Default bias.
    DefaultNodeBias = 0.0

    # Default learning rate.
    DefaultLearningRate = 0.1

#---------------------------------------------------------------------------------------------------
# Class Node (neuron).
#---------------------------------------------------------------------------------------------------

class Node:

#---------------------------------------------------------------------------------------------------

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
        self.DBias = 0.0
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

        return 'Node %s : B = %s, dB = %s' % (self.Id, self.Bias, self.DBias)

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

#---------------------------------------------------------------------------------------------------

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
        self.DWeight = 0.0

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return 'Edge %s : [n%s/o%s -> n%s/i%s] : W = %s, dW = %s' % (self.Id,
                                                                     self.Src.Id,
                                                                     self.OIndex,
                                                                     self.Dst.Id,
                                                                     self.IIndex,
                                                                     self.Weight,
                                                                     self.DWeight)

#---------------------------------------------------------------------------------------------------
# Class Net.
#---------------------------------------------------------------------------------------------------

class Net:

#---------------------------------------------------------------------------------------------------

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

    def SetNodesTraversalOrder(self):
        """
        Set nodes to right traversal order.
        """

        # First remove all signals.
        for node in self.Nodes:
            node.Signals = [None] * len(node.Signals)
            node.Mark = False

        # Set 0.0 signal to first layer.
        for node in self.FirstLayer:
            node.Signals = [0.0]
            node.Mark = True

        # First layer to new order list.
        order = self.FirstLayer.copy()

        # Begin traversal.
        i = 0
        while i < len(order):
            node = order[i]
            for oe in node.OEdges:
                succ = oe.Dst
                succ.Signals[oe.IIndex] = 0.0
                if succ.IsSignalsReady():
                    order.append(succ)
                    succ.Mark = True
            i += 1

        # Set new nodes order.
        self.Nodes = order

        # Clean.
        for node in self.Nodes:
            node.Signals = [None] * len(node.Signals)
            node.Mark = False

#---------------------------------------------------------------------------------------------------

    def Print(self, is_print_nodes = True, is_print_edges = False):
        """
        Print neuronet.

        Arguments:
            is_print_nodes -- flag for nodes print,
            is_print_edges -- flag for edges print.
        """

        if is_print_nodes:
            print('Nodes (%d):' % self.NodesCount())
            for node in self.Nodes:
                print(str(node))

        if is_print_edges:
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
            node.Signals = [None] * len(node.Signals)

        # Propagate forward.
        for i in range(len(self.FirstLayer)):
            self.FirstLayer[i].Signals[0] = x[i]
        for node in self.Nodes:
            node.ForwardPropagation()

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
            node.E = None
            node.Errors = [None] * len(node.Errors)

        # Propagate back.
        for i in range(len(self.LastLayer)):
            self.LastLayer[i].Errors[0] = y[i] - node.A
        for node in self.Nodes.__reversed__():
            node.BackPropagation()

#---------------------------------------------------------------------------------------------------

    def ZeroDWeightsAndDBiases(self):
        """
        Set all DWeights and DBiases to zero.
        """

        for node in self.Nodes:
            node.DBias = 0.0
            for edge in node.IEdges:
                edge.DWeight = 0.0

#---------------------------------------------------------------------------------------------------

    def StoreDWeightsAndDBiases(self):
        """
        Store all DWeigths and DBiases.
        """

        for node in self.Nodes:
            node.DBias += node.E
            for i in range(len(node.IEdges)):
                node.IEdges[i].DWeight += (node.SavedSignals[i] * node.E)

#---------------------------------------------------------------------------------------------------

    def CorrectWeightsAndBiases(self, eta):
        """
        Correct weights and biases.

        Arguments:
            Learning speed.
        """

        for node in self.Nodes:
            node.Bias += eta * node.DBias
            for i in range(len(node.IEdges)):
                node.IEdges[i].Weight += eta * node.IEdges[i].DWeight

#---------------------------------------------------------------------------------------------------

    def SingleLearn(self, x, y):
        """
        Learn on single case.

        Arguments:
            x -- input,
            y -- right output.
        """

        while True:

            t0 = time.clock()
            self.SenseForward(x)
            c = self.Cost(y)

            if c < 0.001:
                print('single learn : learning is finished')
                return
            else:
                self.SenseBack(y)
                self.CorrectWeightsAndBiases()
                t1 = time.clock()
                print('cost = %s, iter time = %s' % (c, t1 - t0))

#---------------------------------------------------------------------------------------------------
# Class MNIST tests.
#---------------------------------------------------------------------------------------------------

class MNISTTests:

#---------------------------------------------------------------------------------------------------

    def __init__(self,
                 imgs_file = '../data/mnist/t10k-images.idx3-ubyte',
                 labs_file = '../data/mnist/t10k-labels.idx1-ubyte'):
        """
        Parser initialization.

        Arguments:
            imgs_file -- binary file of images,
            labs_file -- binary file of labels.
        """

        with open(imgs_file, 'rb') as imgs_bf:
            with open(labs_file, 'rb') as labs_bf:
                imgs_b, labs_b = imgs_bf.read(), labs_bf.read()

                # Cut first 4 bytes.
                imgs_b_803, labs_b_801 = imgs_b[0:4], labs_b[0:4]
                imgs_v_803 = int.from_bytes(imgs_b_803, byteorder = 'big')
                labs_v_801 = int.from_bytes(labs_b_801, byteorder = 'big')
                if imgs_v_803 != 0x803:
                    raise Exception('img file is corrupted')
                if labs_v_801 != 0x801:
                    raise Exception('lab file is corrupted')

                # Count.
                imgs_b_count, labs_b_count = imgs_b[4:8], labs_b[4:8]
                self.Count = int.from_bytes(imgs_b_count, byteorder = 'big')
                labs_count = int.from_bytes(labs_b_count, byteorder = 'big')
                if self.Count != labs_count:
                    raise Exception('data is corrupted')

                # Correct binaries.
                imgs_b, labs_b = imgs_b[8:], labs_b[8:]

                # Cut binaries.
                s = 784
                self.Cases = []
                for i in range(self.Count):
                    img_v = [v / 255.0 for v in list(imgs_b[i * s : (i + 1) * s])]
                    lab_v = labs_b[i]
                    if lab_v > 9:
                        raise Exception('wrong label value')
                    lab_a = [0.0] * 10
                    lab_a[lab_v] = 1.0
                    self.Cases.append((img_v, lab_a))

#---------------------------------------------------------------------------------------------------

    def Cut(self, n):
        """
        Cur tests.

        Arguments:
            n -- target tests count.
        """

        if (self.Count > n):
            self.Count = n
            self.Cases = self.Cases[:n]

#---------------------------------------------------------------------------------------------------

    def RandomCase(self):
        """
        Get random case.

        Result:
            Random case.
        """

        return self.Cases[random.randrange(self.Count)]

#---------------------------------------------------------------------------------------------------

    def Batch(self, n):
        """
        Get test batch.

        Arguments:
            n -- batch size.
        """

        return [self.RandomCase() for i in range(n)]

#---------------------------------------------------------------------------------------------------
# Odd-even tests.
#---------------------------------------------------------------------------------------------------

class XorTests:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Initialize tests.
        """

        self.Count = 4
        self.Cases = \
        [
            ([0.0, 0.0], [0.0]),
            ([0.0, 1.0], [1.0]),
            ([1.0, 0.0], [1.0]),
            ([1.0, 1.0], [0.0])
        ]

#---------------------------------------------------------------------------------------------------

    def Cut(self, n):
        """
        Cur tests.

        Arguments:
            n -- target tests count.
        """

        if (self.Count > n):
            self.Count = n
            self.Cases = self.Cases[:n]

#---------------------------------------------------------------------------------------------------

    def RandomCase(self):
        """
        Get random case.

        Result:
            Random case.
        """

        return self.Cases[random.randrange(self.Count)]

#---------------------------------------------------------------------------------------------------

    def Batch(self, n):
        """
        Get test batch.

        Arguments:
            n -- batch size.
        """

        return [self.RandomCase() for i in range(n)]

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    # Create net.
    net = Net()
    net.CreateMultilayer([2, 1])
    net.SetNodesTraversalOrder()

    # Create parser.
    par = XorTests()
    par.Cut(100)

    t0 = time.clock()

    stages = 1
    tests = 1
    iterations = 5

    for stage in range(stages):

        # Get batch.
        batch = par.Cases[1:2] #par.Batch(tests)

        for iteration in range(iterations):
            net.ZeroDWeightsAndDBiases()
            for (x, y) in batch:
                net.SenseForward(x)
                net.SenseBack(y)
                net.StoreDWeightsAndDBiases()
            net.CorrectWeightsAndBiases(Settings.DefaultLearningRate / tests)
        
        t1 = time.clock()
        dt = t1 - t0
        print('time = %s' %dt)

        # Check.
        net.Print(is_print_edges = True)
        for (x, y) in batch:
            a = net.SenseForward(x)
            print('check   : ', x, y, [round(ai, 2) for ai in a])
        (x, y) = ([0.0, 0.0], [0.0])
        az = net.SenseForward(x)
        print('check z : ', x, y, [round(ai, 2) for ai in az])


#---------------------------------------------------------------------------------------------------
