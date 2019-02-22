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
    s.Set("watt", 36.0)
    s.Set("pue", 2.0)
    s.Set("interconnect", 56.0)
    n = s.AddChildTN("node", "100k")
    s.SetOuter(n, "count", 110)
    c = n.AddChildTree(cpu_Intel_Xeon_E5450())
    n.SetOuter(c, "count", 2)
    n.SetOuter(c, "ram", 8)
    #
    s = t.AddChildTN("segment", "ps", "Petastream")
    s.Set("watt", 15.0)
    s.Set("pue", 1.25)
    s.Set("interconnect", 56.0)
    n = s.AddChildTN("node", "ps")
    s.SetOuter(n, "count", 8)
    c = n.AddChildTree(cpu_Intel_Xeon_E5_2667())
    n.SetOuter(c, "count", 1)
    n.SetOuter(c, "ram", 8)
    c = n.AddChildTree(cpu_Intel_Xeon_Phi_7120D())
    n.SetOuter(c, "count", 8)
    n.SetOuter(c, "ram", 16)
    #
    s = t.AddChildTN("segment", "tr", "Tornado")
    s.Set("watt", 223.0)
    s.Set("pue", 1.25)
    s.Set("interconnect", 56.0)
    n = s.AddChildTN("node", "tr")
    s.SetOuter(n, "count", 207)
    c = n.AddChildTree(cpu_Intel_Xeon_E5_2690())
    n.SetOuter(c, "count", 2)
    n.SetOuter(c, "ram", 64)
    c = n.AddChildTree(cpu_Intel_Xeon_Phi_7110X())
    n.SetOuter(c, "count", 2)
    n.SetOuter(c, "ram", 16)
    #
    s = t.AddChildTN("segment", "hw", "Haswell")
    s.Set("watt", 28.0)
    s.Set("pue", 1.06)
    s.Set("interconnect", 100.0)
    n = s.AddChildTN("node", "hw")
    s.SetOuter(n, "count", 42)
    c = n.AddChildTree(cpu_Intel_Xeon_E5_2697v3())
    n.SetOuter(c, "count", 2)
    n.SetOuter(c, "ram", 128)
    #
    s = t.AddChildTN("segment", "bw", "Broadwell")
    s.Set("watt", 91.0)
    s.Set("pue", 1.06)
    s.Set("interconnect", 100.0)
    n = s.AddChildTN("node", "bw")
    s.SetOuter(n, "count", 136)
    c = n.AddChildTree(cpu_Intel_Xeon_E5_2697Av4())
    n.SetOuter(c, "count", 2)
    n.SetOuter(c, "ram", 128)
    #
    s = t.AddChildTN("segment", "knl", "Knights Landing")
    s.Set("watt", 29.0)
    s.Set("pue", 1.06)
    s.Set("interconnect", 100.0)
    n = s.AddChildTN("node", "knl")
    s.SetOuter(n, "count", 38)
    c = n.AddChildTree(cpu_Intel_Xeon_Phi_7290())
    n.SetOuter(c, "count", 1)
    n.SetOuter(c, "ram", 96)
    #
    s = t.AddChildTN("segment", "nv", "NVIDIA")
    s.Set("watt", 19.0)
    s.Set("pue", 2.0)
    s.Set("interconnect", 56.0)
    n = s.AddChildTN("node", "nv")
    s.SetOuter(n, "count", 6)
    c = n.AddChildTree(cpu_Intel_Xeon_X5675())
    n.SetOuter(c, "count", 2)
    n.SetOuter(c, "ram", 192)
    c = n.AddChildTree(cpu_NVIDIA_Tesla_M2090())
    n.SetOuter(c, "count", 8)
    n.SetOuter(c, "ram", 48)

    return t;

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    jscc = center_jscc()
    jscc.ApplyUpward(lambda t: t.GatherTacticSumWithCount("tfs"))
    jscc.Apply(lambda t: t.GatherTacticSumOuterProperties("ram"))
    jscc.ApplyUpward(lambda t: t.GatherTacticSumWithCount("ram"))
    jscc.PrintTree()

#---------------------------------------------------------------------------------------------------
