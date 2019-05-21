# -*- coding: utf-8 -*-
"""
Hierarchical tree.

Created on Tue May 21 16:40:37 2019

@author: Rybakov
"""

class HTree:
    """
    Hierarchical tree.
    """

#---------------------------------------------------------------------------------------------------
# Constructor.
#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        # Links to parent and children.
        self.Parent = None
        self.Children= []

#---------------------------------------------------------------------------------------------------
