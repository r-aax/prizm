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

    # Default learning rate.
    DefaultLearningRate = 3.0

#---------------------------------------------------------------------------------------------------
# Class Node (neuron).
#---------------------------------------------------------------------------------------------------

class Node:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.IEdges = []
        self.OEdges = []
        self.B = 0.0
        self.dB = 0.0
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

        return 'Node : B = %s' % self.B

#---------------------------------------------------------------------------------------------------

    def Signals(self):
        """
        Signals vector.

        Result:
            Signals.
        """

        return [e.S for e in self.IEdges]

#---------------------------------------------------------------------------------------------------

    def Errors(self):
        """
        Errors vector.

        Result:
            Errors.
        """

        return [e.E for e in self.OEdges]

#---------------------------------------------------------------------------------------------------


    def IsSignalsReady(self):
        """
        Check if signals are ready.

        Result:
            True - if signals are ready,
            Faslse - otherwise.
        """

        return fun.is_all(self.Signals(), lambda x: x != None)

#---------------------------------------------------------------------------------------------------

    def IsErrorsReady(self):
        """
        Check if error are ready.

        Result:
            True - if errors  are ready,
            False - otherwise.
        """

        return fun.is_all(self.Errors(), lambda x: x != None)

#---------------------------------------------------------------------------------------------------

    def ForwardPropagation(self):
        """
        Forward propagation of the signal.
        """

        # We need calculate self.A only for work nodes.
        if self.IEdges != []:
            self.A = mth.sigmoid(sum([e.S * e.W for e in self.IEdges]) + self.B)

        # Propagate.
        for oe in self.OEdges:
            oe.S = self.A

#---------------------------------------------------------------------------------------------------

    def BackPropagation(self):
        """
        Back propagation of the error.
        """

        # Calculate error for nodes not from the last layer.
        if self.OEdges != []:
            self.E = sum(self.Errors()) * self.A * (1.0 - self.A)

        # Propagate.
        for ie in self.IEdges:
            ie.E = self.E * ie.W

#---------------------------------------------------------------------------------------------------
# Class Edge.
#---------------------------------------------------------------------------------------------------

class Edge:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.Src = None
        self.Dst = None
        self.S = None
        self.W = 1.0
        self.dW = 0.0
        self.E = None

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return 'Edge : W = %s' % self.W

#---------------------------------------------------------------------------------------------------
# Class Net.
#---------------------------------------------------------------------------------------------------

class Net:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        # Nodes - all nodes of the net.
        #
        # | FirstLayer |           HiddenNodes          |  LastLayer |
        # + ---------- + ------------------------------ + ---------- +
        # |            |                  WorkNodes                  |

        # Nodes initialization.
        self.Nodes = []
        self.FirstLayer = []
        self.HiddenNodes = []
        self.LastLayer = []
        self.WorkNodes = []

        # Edges initialization.
        self.Edges = []

#---------------------------------------------------------------------------------------------------

    def AddEdge(self, src, dst):
        """
        Add edge between two nodes.

        Arguments:
            src -- source node,
            dst -- destination node.
        """

        e = Edge()

        # Links.
        e.Src = src
        e.Dst = dst
        src.OEdges.append(e);
        dst.IEdges.append(e);
        self.Edges.append(e);

#---------------------------------------------------------------------------------------------------

    def CreateMultilayer(self, layers_sizes):
        """
        Create multilayer net.

        Arguments:
            LayersSizes -- layers sizes.
        """

        # Check the layers count.
        if len(layers_sizes) < 2:
            raise Exception('not enough layers (must be >= 2)')

        # Create all nodes.
        layers = [[Node() for _ in range(layer_size)] for layer_size in layers_sizes]
        nodes = lst.flatten(layers)

        # First and last layers sizes.
        first_layer_size = layers_sizes[0]
        last_layer_size = layers_sizes[-1]

        # Set links to nodes.
        self.Nodes = nodes
        self.FirstLayer = nodes[:first_layer_size]
        self.HiddenNodes = nodes[first_layer_size:-last_layer_size]
        self.LastLayer = nodes[-last_layer_size:]
        self.WorkNodes = nodes[first_layer_size:]

        # Create edges.
        for i in range(len(layers) - 1):
            for src in layers[i]:
                for dst in layers[i + 1]:
                    self.AddEdge(src, dst)

        # Correct nodes biases and edges weights.
        for n in self.WorkNodes:
            n.B = random.gauss(0.0, 1.0)
        for e in self.Edges:
            e.W = random.gauss(0.0, 1.0)

#---------------------------------------------------------------------------------------------------

    def CleanSingals(self):
        """
        Clean signals.
        """

        for e in self.Edges:
            e.S = None

#---------------------------------------------------------------------------------------------------

    def SetNodesTraversalOrder(self):
        """
        Set nodes to right traversal order.
        """

        # Set 0.0 signal to first layer.
        for n in self.FirstLayer:
            n.Mark = True

        # First remove all signals.
        for n in self.WorkNodes:
            n.Mark = False

        self.CleanSingals()

        # First layer to new order list.
        order = self.FirstLayer.copy()

        # Begin traversal.
        i = 0
        while i < len(order):
            n = order[i]
            for oe in n.OEdges:
                oe.S = 0.0
                succ = oe.Dst
                if succ.Mark == False:
                    if succ.IsSignalsReady():
                        order.append(succ)
                        succ.Mark = True
            i += 1

        # Set new nodes order.
        self.Nodes = order

        self.CleanSingals()

        # Clean.
        for n in self.Nodes:
            n.Mark = False

#---------------------------------------------------------------------------------------------------

    def Print(self, is_print_nodes = True, is_print_edges = False):
        """
        Print neuronet.

        Arguments:
            is_print_nodes -- flag for nodes print,
            is_print_edges -- flag for edges print.
        """

        if is_print_nodes:
            print('Nodes (%d):' % len(self.Nodes))
            for i, n in enumerate(self.Nodes):
                print('[%d] : %s' % (i, str(n)))

        if is_print_edges:
            print('Edges (%d)' % len(self.Edges))
            for i, e in enumerate(self.Edges):
                print('[%d] [%d - %d] : %s' % (i,
                                               self.Nodes.index(e.Src),
                                               self.Nodes.index(e.Dst),
                                               str(e)))

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

        # Propagate forward.
        # First set signals to the first layer.
        for i, n in enumerate(self.FirstLayer):
            n.A = x[i]
        for n in self.Nodes:
            n.ForwardPropagation()

        return self.A()

#---------------------------------------------------------------------------------------------------

    def A(self):
        """
        Result value.

        Result:
            Result.
        """

        return [n.A for n in self.LastLayer]

#---------------------------------------------------------------------------------------------------

    def Cost(self, y):
        """
        Calculate cost function.

        Arguments:
            y -- true result.

        Result:
            Cost function value.
        """

        return 0.5 * sum([(yi - ai) * (yi - ai) for yi, ai in zip(y, self.A())])

#---------------------------------------------------------------------------------------------------

    def TotalCost(self, batch):
        """
        Calculate total cost.

        Arguments:
            batch -- batch of tests.

        Result:
            Total cost.
        """

        c = 0.0

        for (x, y) in batch:
            self.SenseForward(x)
            c += self.Cost(y)

        return c / len(batch)

#---------------------------------------------------------------------------------------------------

    def SenseBack(self, y):
        """
        Sense neuronet back.

        Arguments:
            y -- right answer.
        """

        # Propagate back.
        for i, n in enumerate(self.LastLayer):
            n.E = (n.A - y[i]) * n.A * (1.0 - n.A)
        for n in self.Nodes.__reversed__():
            n.BackPropagation()

#---------------------------------------------------------------------------------------------------

    def ZeroDWeightsAndDBiases(self):
        """
        Set all DWeights and DBiases to zero.
        """

        for n in self.Nodes:
            n.dB = 0.0
        for e in self.Edges:
            e.dW = 0.0

#---------------------------------------------------------------------------------------------------

    def StoreDWeightsAndDBiases(self):
        """
        Store all DWeigths and DBiases.
        """

        for n in self.Nodes:
            n.dB += n.E
            for i, e in enumerate(n.IEdges):
                e.dW += (e.S * n.E)

#---------------------------------------------------------------------------------------------------

    def CorrectWeightsAndBiases(self, eta):
        """
        Correct weights and biases.

        Arguments:
            Learning speed.
        """

        for n in self.WorkNodes:
            n.B -= eta * n.dB
        for e in self.Edges:
            e.W -= eta * e.dW

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

    def MiniBatches(self, size = 30):
        """
        Return test in mini batches.

        Arguments:
            size -- mini batch size.

        Return:
            Mini batches array.
        """

        return lst.chop(self.Cases, size)

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

    def MiniBatches(self, size = 4):
        """
        Return test in mini batches.

        Arguments:
            size -- mini batch size.

        Return:
            Mini batches array.
        """

        return [self.Cases]

#---------------------------------------------------------------------------------------------------
# Trainer.
#---------------------------------------------------------------------------------------------------

class Trainer:

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        pass

#---------------------------------------------------------------------------------------------------

    def Train(self,
              net,
              mini_batches,
              max_epochs_count = 1000000,
              min_total_cost = 0.0001,
              print_step = -1):
        """
        Train.

        Arguments:
            net -- neuronet to train,
            mini_batches -- mini batches of tests,
            max_epochs_count -- max epochs count for training,
            min_total_cost -- total cost for training finish,
            print_step -- step of epochs for print (if -1 - then print is forbidden).
        """

        epoch = 0

        # Calculate total cost.
        total_cost = sum([net.TotalCost(batch) for batch in mini_batches]) / len(mini_batches)

        while True:

            # Check max epochs count is reached.
            if epoch >= max_epochs_count:
                print('max epochs count is reached', max_epochs_count)
                return (epoch, total_cost);

            # The net if educated enough.
            if total_cost <= min_total_cost:
                print('min total cost is reached', min_total_cost)
                return (epoch, total_cost);

            t0 = time.clock()

            # Studying.
            for batch in mini_batches:
                net.ZeroDWeightsAndDBiases()
                for (x, y) in batch:
                    net.SenseForward(x)
                    net.SenseBack(y)
                    net.StoreDWeightsAndDBiases()
                net.CorrectWeightsAndBiases(Settings.DefaultLearningRate / len(batch))

            # Calculate total cost.
            total_cost = sum([net.TotalCost(batch) for batch in mini_batches]) / len(mini_batches)

            t1 = time.clock()
            dt = t1 - t0

            if epoch % print_step == 0:
                print('Epoch %d : time = %f, total_cost = %f' % (epoch, dt, total_cost))
            epoch += 1

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    # Create net.
    net = Net()
    net.CreateMultilayer([784, 15, 10])
    net.SetNodesTraversalOrder()

    # Create parser.
    par = MNISTTests()

    # Train.
    trainer = Trainer()
    res = trainer.Train(net, par.MiniBatches(20)[:1],
                        max_epochs_count = 1000, print_step = 1)
    print('Result : ', res)

    for (x, y) in par.MiniBatches(20)[0]:
        a = net.SenseForward(x)
        print(y, 'vs', [round(ai, 2) for ai in a])

#---------------------------------------------------------------------------------------------------
