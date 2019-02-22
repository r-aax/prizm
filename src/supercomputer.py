# -*- coding: utf-8 -*-
"""
Some supercomputer functions.

Created on Thu Feb 21 15:06:46 2019

@author: Rybakov
"""

from ftree import FTree

#---------------------------------------------------------------------------------------------------
# CPUs.
#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5450():
    """
    Harpertown microprocessor (Intel Xeon E5450).

    Result:
        Harpertown microprocessor.
    """

    t = FTree("cpu", "ht", "Intel Xeon E5450 Harpertown")
    t.Set("cores_count", 4)
    t.Set("freq", 3.0)
    t.Set("tfs", 0.048)

    return t;

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_X5675():
    """
    Westmere microprocessor (Intel Xeon X5675).

    Result:
        Westmere microprocessor.
    """

    t = FTree("cpu", "wm", "Intel Xeon X5675 Westmere")
    t.Set("cores_count", 6)
    t.Set("freq", 3.06)
    t.Set("tfs", 0.14488)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2667():
    """
    Ivy Bridge microprocessor (Intel Xeon E5-2667).

    Result:
        Ivy Bridge microprocessor.
    """

    t = FTree("cpu", "ib", "Inte Xeon E5-2667 Ivy Bridge")
    t.Set("cores_count", 8)
    t.Set("freq", 3.3)
    t.Set("tfs", 0.2112)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2690():
    """
    Sandy Bridge microprocessor (Intel Xeon E5-2690).

    Result:
        Sandy Bridge microprocessor.
    """

    t = FTree("cpu", "sb", "Intel Xeon E5-2690 Sandy Bridge")
    t.Set("cores_count", 8)
    t.Set("freq", 2.9)
    t.Set("tfs", 0.1856)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2697v3():
    """
    Haswell microprocessor (Intel Xeon E5-2697v3).

    Result:
        Haswell microprocessor.
    """

    t = FTree("cpu", "hw", "Intel Xeon E5-2697v3 Haswell")
    t.Set("cores_count", 14)
    t.Set("freq", 2.6)
    t.Set("tfs", 0.5824)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2697Av4():
    """
    Broadwell microprocessor (Intel Xeon E5-2697Av4).

    Result:
        Broadwell microprocessor.
    """

    t = FTree("cpu", "bw", "Intel Xeon E5-2697Av4 Broadwell")
    t.Set("cores_count", 16)
    t.Set("freq", 2.6)
    t.Set("tfs", 0.6656)

    return t;

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_Phi_7120D():
    """
    Knights Corner microprocessor in Petastream (Intel Xeon Phi 7120D).

    Result:
        Knights Corner microprocessor.
    """

    t = FTree("cpu", "knc", "Intel Xeon Phi 7120D KNC")
    t.Set("cores_count", 61)
    t.Set("freq", 1.238)
    t.Set("tfs", 1.208)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_Phi_7110X():
    """
    Knights Corner microprocessor in Tornado (Intel Xeon Phi 7110X).

    Result:
        Knights Corner microprocessor.
    """

    t = FTree("cpu", "knc", "Intel Xeon Phi 7110X KNC")
    t.Set("cores_count", 61)
    t.Set("freq", 1.1)
    t.Set("tfs", 1.0736)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_Phi_7290():
    """
    Knights Landing microprocessor (Intel Xeon Phi 7290).

    Result:
        Knights Landing microprocessor.
    """

    t = FTree("cpu", "knl", "Intel Xeon Phi 7290 KNL")
    t.Set("cores_count", 72)
    t.Set("freq", 1.5)
    t.Set("tfs", 3.456)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_NVIDIA_Tesla_M2090():
    """
    Tesla microprocessor (NVIDIA Tesla M2090).

    Result:
        Tesla microprocessor.
    """

    t = FTree("cpu", "nv", "NVIDIA Tesla M2090")
    t.Set("cores_count", 512)
    t.Set("freq", 1.3)
    t.Set("tfs", 0.665)

    return t

#---------------------------------------------------------------------------------------------------
# Supercomputer center.
#---------------------------------------------------------------------------------------------------

def center_jscc():
    """
    JSCC center.

    Result:
        JSCC center.
    """

    t = FTree("center", "jscc", "Joint Supercomputer Center");

    # Init.
    s = t.AddChildTN("segment", "100k", "MVS-100K")
    n = s.AddChildTN("node", "100k")
    n.AddChildTree(cpu_Intel_Xeon_E5450())
    #
    s = t.AddChildTN("segment", "ps", "Petastream")
    n = s.AddChildTN("node", "ps")
    n.AddChildTree(cpu_Intel_Xeon_E5_2667())
    n.AddChildTree(cpu_Intel_Xeon_Phi_7120D())
    #
    s = t.AddChildTN("segment", "tr", "Tornado")
    n = s.AddChildTN("node", "tr")
    n.AddChildTree(cpu_Intel_Xeon_E5_2690())
    n.AddChildTree(cpu_Intel_Xeon_Phi_7110X())
    #
    s = t.AddChildTN("segment", "hw", "Haswell")
    n = s.AddChildTN("node", "hw")
    n.AddChildTree(cpu_Intel_Xeon_E5_2697v3())
    #
    s = t.AddChildTN("segment", "bw", "Broadwell")
    n = s.AddChildTN("node", "bw")
    n.AddChildTree(cpu_Intel_Xeon_E5_2697Av4())
    #
    s = t.AddChildTN("segment", "knl", "Knights Landing")
    n = s.AddChildTN("node", "knl")
    n.AddChildTree(cpu_Intel_Xeon_Phi_7290())
    #
    s = t.AddChildTN("segment", "nv", "NVIDIA")
    n = s.AddChildTN("node", "nv")
    n.AddChildTree(cpu_Intel_Xeon_X5675())
    n.AddChildTree(cpu_NVIDIA_Tesla_M2090())

    # Add special print function for cpus.
    t.Apply(lambda t: t.Set("to_string_fun",
                            lambda x: x.BaseStr()
                                      + " (c" + str(x.Get("cores_count"))
                                      + "/f" + str(x.Get("freq"))
                                      + "/t" + str(x.Get("tfs")) + ")"),
            lambda t: t.IsType("cpu"))

    # Add special print function for node.
    t.Apply(lambda t: t.Set("to_string_fun",
                            lambda x: x.BaseStr()),
            lambda t: t.IsType("node"))

    # Add special print function for segments.
    t.Apply(lambda t: t.Set("to_string_fun",
                            lambda x: x.BaseStr()),
            lambda t: t.IsType("segment"))

    return t;

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    jscc = center_jscc()
    jscc.PrintTree()

#---------------------------------------------------------------------------------------------------
