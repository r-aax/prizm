# -*- coding: utf-8 -*-
"""
Some supercomputer functions.

Created on Thu Feb 21 15:06:46 2019

@author: Rybakov
"""

from ftree.tree import Tree

#---------------------------------------------------------------------------------------------------
# CPUs.
#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5450():
    """
    Harpertown microprocessor (Intel Xeon E5450).

    Result:
        Harpertown microprocessor.
    """

    t = Tree('cpu', 'ht', 'Intel Xeon E5450 Harpertown')
    t.Set('cores_count', 4)
    t.Set('freq', 3.0)
    t.Set('tfs', 0.048)

    return t;

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_X5675():
    """
    Westmere microprocessor (Intel Xeon X5675).

    Result:
        Westmere microprocessor.
    """

    t = Tree('cpu', 'wm', 'Intel Xeon X5675 Westmere')
    t.Set('cores_count', 6)
    t.Set('freq', 3.06)
    t.Set('tfs', 0.14488)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2667():
    """
    Ivy Bridge microprocessor (Intel Xeon E5-2667).

    Result:
        Ivy Bridge microprocessor.
    """

    t = Tree('cpu', 'ib', 'Intel Xeon E5-2667 Ivy Bridge')
    t.Set('cores_count', 8)
    t.Set('freq', 3.3)
    t.Set('tfs', 0.2112)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2690():
    """
    Sandy Bridge microprocessor (Intel Xeon E5-2690).

    Result:
        Sandy Bridge microprocessor.
    """

    t = Tree('cpu', 'sb', 'Intel Xeon E5-2690 Sandy Bridge')
    t.Set('cores_count', 8)
    t.Set('freq', 2.9)
    t.Set('tfs', 0.1856)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2697v3():
    """
    Haswell microprocessor (Intel Xeon E5-2697v3).

    Result:
        Haswell microprocessor.
    """

    t = Tree('cpu', 'hw', 'Intel Xeon E5-2697v3 Haswell')
    t.Set('cores_count', 14)
    t.Set('freq', 2.6)
    t.Set('tfs', 0.5824)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_E5_2697Av4():
    """
    Broadwell microprocessor (Intel Xeon E5-2697Av4).

    Result:
        Broadwell microprocessor.
    """

    t = Tree('cpu', 'bw', 'Intel Xeon E5-2697Av4 Broadwell')
    t.Set('cores_count', 16)
    t.Set('freq', 2.6)
    t.Set('tfs', 0.6656)

    return t;

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_Phi_7120D():
    """
    Knights Corner microprocessor in Petastream (Intel Xeon Phi 7120D).

    Result:
        Knights Corner microprocessor.
    """

    t = Tree('cpu', 'knc', 'Intel Xeon Phi 7120D KNC')
    t.Set('cores_count', 61)
    t.Set('freq', 1.238)
    t.Set('tfs', 1.208)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_Phi_7110X():
    """
    Knights Corner microprocessor in Tornado (Intel Xeon Phi 7110X).

    Result:
        Knights Corner microprocessor.
    """

    t = Tree('cpu', 'knc', 'Intel Xeon Phi 7110X KNC')
    t.Set('cores_count', 61)
    t.Set('freq', 1.1)
    t.Set('tfs', 1.0736)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_Intel_Xeon_Phi_7290():
    """
    Knights Landing microprocessor (Intel Xeon Phi 7290).

    Result:
        Knights Landing microprocessor.
    """

    t = Tree('cpu', 'knl', 'Intel Xeon Phi 7290 KNL')
    t.Set('cores_count', 72)
    t.Set('freq', 1.5)
    t.Set('tfs', 3.456)

    return t

#---------------------------------------------------------------------------------------------------

def cpu_NVIDIA_Tesla_M2090():
    """
    Tesla microprocessor (NVIDIA Tesla M2090).

    Result:
        Tesla microprocessor.
    """

    t = Tree('cpu', 'nv', 'NVIDIA Tesla M2090')
    t.Set('cores_count', 512)
    t.Set('freq', 1.3)
    t.Set('tfs', 0.665)

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

    t = Tree('center', 'jscc', 'Joint Supercomputer Center');

    # Init.
    with t.AddChildTN('segment', '100k', 'MVS-100K') as s:
        s.Set('watt', 36.0)
        s.Set('pue', 2.0)
        s.Set('interconnect', 56.0)
        with s.AddChildTN('node', '100k') as n:
            n.Pred.Set('count', 110)
            with n.AddChildTN('cpus', 'ht') as cs:
                cs.Set('ram', 8.0)
                with cs.AddChild(cpu_Intel_Xeon_E5450()) as c:
                    c.Pred.Set('count', 2)
    with t.AddChildTN('segment', 'ps', 'Petastream') as s:
        s.Set('watt', 15.0)
        s.Set('pue', 1.25)
        s.Set('interconnect', 56.0)
        with s.AddChildTN('node', 'ps') as n:
            n.Pred.Set('count', 8)
            with n.AddChildTN('cpus', 'ib') as cs:
                cs.Set('ram', 8.0)
                with cs.AddChild(cpu_Intel_Xeon_E5_2667()) as c:
                    c.Pred.Set('count', 1)
            with n.AddChildTN('cpus', 'knc') as cs:
                cs.Set('ram', 16.0)
                with cs.AddChild(cpu_Intel_Xeon_Phi_7120D()) as c:
                    c.Pred.Set('count', 8)
    with t.AddChildTN('segment', 'tr', 'Tornado') as s:
        s.Set('watt', 223.0)
        s.Set('pue', 1.25)
        s.Set('interconnect', 56.0)
        with s.AddChildTN('node', 'tr') as n:
            n.Pred.Set('count', 207)
            with n.AddChildTN('cpus', 'sb') as cs:
                cs.Set('ram', 64.0)
                with cs.AddChild(cpu_Intel_Xeon_E5_2690()) as c:
                    c.Pred.Set('count', 2)
            with n.AddChildTN('cpus', 'knc') as cs:
                cs.Set('ram', 16.0)
                with cs.AddChild(cpu_Intel_Xeon_Phi_7110X()) as c:
                    c.Pred.Set('count', 2)
    with t.AddChildTN('segment', 'hw', 'Haswell') as s:
        s.Set('watt', 28.0)
        s.Set('pue', 1.06)
        s.Set('interconnect', 100.0)
        with s.AddChildTN('node', 'hw') as n:
            n.Pred.Set('count', 42)
            with n.AddChildTN('cpus', 'hw') as cs:
                cs.Set('ram', 128.0)
                with cs.AddChild(cpu_Intel_Xeon_E5_2697v3()) as c:
                    c.Pred.Set('count', 2)
    with t.AddChildTN('segment', 'bw', 'Broadwell') as s:
        s.Set('watt', 91.0)
        s.Set('pue', 1.06)
        s.Set('interconnect', 100.0)
        with s.AddChildTN('node', 'bw') as n:
            n.Pred.Set('count', 136)
            with n.AddChildTN('cpus', 'bw') as cs:
                cs.Set('ram', 128.0)
                with cs.AddChild(cpu_Intel_Xeon_E5_2697Av4()) as c:
                    c.Pred.Set('count', 2)
    with t.AddChildTN('segment', 'knl', 'Knights Landing') as s:
        s.Set('watt', 29.0)
        s.Set('pue', 1.06)
        s.Set('interconnect', 100.0)
        with s.AddChildTN('node', 'knl') as n:
            n.Pred.Set('count', 38)
            with n.AddChildTN('cpus', 'knl') as cs:
                cs.Set('ram', 96.0)
                with cs.AddChild(cpu_Intel_Xeon_Phi_7290()) as c:
                    c.Pred.Set('count', 1)
    with t.AddChildTN('segment', 'nv', 'NVIDIA') as s:
        s.Set('watt', 19.0)
        s.Set('pue', 2.0)
        s.Set('interconnect', 56.0)
        with s.AddChildTN('node', 'nv') as n:
            n.Pred.Set('count', 6)
            with n.AddChildTN('cpus', 'wm') as cs:
                cs.Set('ram', 192.0)
                with cs.AddChild(cpu_Intel_Xeon_X5675()) as c:
                    c.Pred.Set('count', 2)
            with n.AddChildTN('cpus', 'nv') as cs:
                cs.Set('ram', 48.0)
                with cs.AddChild(cpu_NVIDIA_Tesla_M2090()) as c:
                    c.Pred.Set('count', 8)

    return t;

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    jscc = center_jscc()
    jscc.ApplyUpward(lambda t: t.GatherTacticSumWithCount('tfs'))
    jscc.PrintTree()

#---------------------------------------------------------------------------------------------------
