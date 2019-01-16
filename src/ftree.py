# -*- coding: utf-8 -*-
"""
Functional tree realization.
Tree consists of nodes linked by edges.
Each node is dictionary.

Predefined node keys are:
    - type
    - name
    - children

Created on Tue Jan 15 16:05:25 2019

@author: Rybakov
"""

#---------------------------------------------------------------------------------------------------

def get_type(el):
    """
    Get special "type" field value of the element.
    
    Arguments:
        el -- element.
        
    Result:
        Value of the field "type".
    """
    
    return el["type"]

#---------------------------------------------------------------------------------------------------

def set_type(el, v):
    """
    Set special "type" field of the element.
    
    Arguments:
        el -- element,
        v -- value.
    """
    
    el["type"] = v
    
#---------------------------------------------------------------------------------------------------

def is_type(el, v):
    """
    Check if the element has given type.
    
    Arguments:
        el -- element,
        v -- value of the field "type".
        
    Result:
        True -- if the element has given type,
        False -- otherwise.
    """
    
    return get_type(el) == v
    
#---------------------------------------------------------------------------------------------------

def get_name(el):
    """
    Get special "name" field value of the element.
    
    Arguments:
        el -- element.
        
    Result:
        Value of the field "name".
    """
    
    return el["name"]

#---------------------------------------------------------------------------------------------------

def set_name(el, v):
    """
    Set special "name" field of the element.
    
    Arguments:
        el -- element,
        v -- value.
    """
    
    el["name"] = v
    
#---------------------------------------------------------------------------------------------------

def is_name(el, v):
    """
    Check if the element has given name.
    
    Arguments:
        el -- element,
        v -- value of the field "name".
        
    Result:
        True -- if the element has given name,
        False -- otherwise.
    """
    
    return get_name(el) == v

#---------------------------------------------------------------------------------------------------
    
def get_children(el):
    """
    Get special "children" field value of the element.
    
    Arguments:
        el -- element.
        
    Result:
        Value of the field "children".
    """
    
    return el["children"]

#---------------------------------------------------------------------------------------------------

def empty_children(el):
    """
    Empty children list.
    
    Arguments:
        el -- element.
    """
    
    el["children"] = []

#---------------------------------------------------------------------------------------------------
    
def add_child(el, c):
    """
    Add new child.
    
    Arguments:
        el -- element,
        c -- new child.
    """
    
    get_children(el).append(c)
    
#---------------------------------------------------------------------------------------------------

def new_element(tp, nm):
    """
    Create new element from type and name.
    
    Argumnets:
        tp -- type,
        nm -- name.
    
    Result:
        New element.
    """
    
    el = {}
    set_type(el, tp)
    set_name(el, nm)
    empty_children(el)
    
    return el

#---------------------------------------------------------------------------------------------------
    
def print_tree(tree, level=0):
    """
    Recursive print tree.
    
    Arguments:
        t -- tree,
        sh - shift of print.
    """
    
    if level == 0:
        sh_str = "# " # root
    else:
        sh= 4 * (level - 1)
        sh_str = (" " * sh) + "|--->"
    
    # Print current element.
    print(sh_str + get_type(tree) + " : " + get_name(tree))
    
    # Print all children.
    for c in get_children(tree):
        print_tree(c, level + 1)
    
#---------------------------------------------------------------------------------------------------
        
def find_element(tree, fun):
    """
    Find element in the tree.
    
    Arguments:
        tree -- tree,
        fun -- search function.
        
    Result:
        Found element, if it is in th tree,
        None, if element is not found.
    """
    
    if fun(tree):
        return tree
    else:
        
        # Check all children.
        for c in get_children(tree):
            if find_element(c, fun) != None:
                return c
        
        # Nothing is found.
        return None
        
#---------------------------------------------------------------------------------------------------
    
def find_element_by_name(tree, nm):
    """
    Find element by its name.

    Arguments:
        tree -- tree,
        nm -- name.
    
    Result:
        Found element, if it is in the tree,
        None, if element is not found.
    """
    
    return find_element(tree, lambda x: is_name(x, nm))
    
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    print("ftree tests:")
    jscc = new_element("resources", "JSCC")
    k100 = new_element("supercomputer", "100K")
    p10 = new_element("supercomputer", "10P")
    add_child(jscc, k100)
    add_child(jscc, p10)
    print_tree(jscc)
    
    print(find_element_by_name(jscc, "100K"))
    
    
    
    
    
    
    
    
    
    