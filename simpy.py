# -*- coding: utf-8 -*-
"""
Created on Thu May 18 21:34:37 2023

Here we give some functions and classes for computing the Alexander
polynomial of a (1,1) knot in a lens space, particularly simple knots and
Hedden knots

@author: jacob
"""
from collections import defaultdict

def add_laurent(poly1, poly2):
    
    """
    Function adds two Laurent polynomials represented by dictionaries whose
    keys are integers and poly[i] is the coefficient of t^i in the polynomial.

    Parameters
    ----------
    poly1 : dict
    poly2 : dict

    Returns
    -------
    poly1 + poly2 as a dict

    """
    
    poly = defaultdict(int)
    for power in poly1.keys():
        poly[power] += poly1[power]
    for power in poly2.keys():
        poly[power] += poly2[power]
    return {key:poly[key] for key in poly.keys() if poly[key] != 0}

def multiply_laurent(poly1, poly2):
    
    """
    Parameters
    ----------
    poly1 : dict
    poly2 : dict

    Returns
    -------
    poly1*poly2 as a dict

    """
    poly = defaultdict(int)
    for power1 in poly1.keys():
        for power2 in poly2.keys():
            poly[power1 + power2] += poly1[power1]*poly2[power2]
    return {power:poly[power] for power in poly.keys() if poly[power] != 0}

def simplify_laurent(poly):
    
    """
    Removes keys in poly.keys() whose values are 0

    Parameters
    ----------
    poly : dict

    Returns
    -------
    poly : dict

    """
    
    for k in poly.keys():
        if poly[k] == 0:
            poly.pop(k)
    return poly

def symmetrize_laurent(poly):
    
    """
    Normalize a symmetric Laurent polynomial to be symmetric about degree 0
    
    Parameters
    ----------
    poly : dict

    Returns
    -------
    dict
    
    """
    
    k = (max(poly) + min(poly))//2
    return multiply_laurent({-k:1}, poly)
 
def positive_laurent(poly):
    """
    
    Normalize a Laurent polynomial to be an honest polynomial in Z[t]
    
    Parameters
    ----------
    poly : dict

    Returns
    -------
    dict

    """
    k = min(poly)
    return multiply_laurent({-k:1},poly)

def divide_laurent(poly1, poly2, quotient = defaultdict(int)):
    """
    Returns the Laurent polynomial quotient poly1/poly2

    Parameters
    ----------
    poly1 : dict
    poly2 : dict
    quotient : dict

    Returns
    -------
    dict

    """
    k = max(poly1)-max(poly2)
    if k == 0:
        quotient[0] = 1
        return quotient
    else:
        quotient[k] = int(poly1[max(poly1.keys())]/poly2[max(poly2.keys())])
        return divide_laurent(add_laurent(poly1, multiply_laurent({k:-quotient[k]}, poly2)), poly2, quotient)
        
def allones(p):
    """
    Returns the polynomial 1 + t + ... + t^{p-1} as a dict

    Parameters
    ----------
    p : int

    Returns
    -------
    dict

    """
    return {power:1 for power in range(p)}
 

def laurent_expression(poly):
    """
    Writes a Laurent polynomial as a string

    Parameters
    ----------
    poly : dict

    Returns
    -------
    string

    """
    keys = sorted(poly.keys())
    expression = ''
    for key in keys:
        expression += str(poly[key]) + "t^(" + str(key) + ") + "
    return expression[:-3]

def simple_relator(p, q, k):
    """
    Returns the relator associated to the simple knot in the lens space
    L(p,q) representing the homology class k in Z/pZ

    Parameters
    ----------
    p : int > 0
    q : 0 < int < p co-prime with p
    k : 0 < int < p

    Returns
    -------
    relator : string
    
    """
    relator = ""
    for i in range(p):
        j = (i*q) % p
        if j < k:
            relator += "ab"
        else:
            relator += "a"
    return relator

def hedden_relator(p,q):
    """
    Return the relator associated to the Hedden knot T_L in the lens space
    L(p,q)

    Parameters
    ----------
    p : 0 < int
    q : 0 < int < p co-prime with p

    Returns
    -------
    relator : string

    """
    relator = "abAba"
    for i in range(p-1):
        if (i*q)%p + q >= p:
            relator += "ba"
        else:
            relator += "a"
    return relator

def abelianization(relator):
    """
    Given a relator (a string in 'a', 'b', 'A', 'B'), return a dictionary
    containing the images of the generators 'a', 'b', 'A', 'B' of the 
    fundamental group in the abelianization Z

    Parameters
    ----------
    relator : string in 'a', 'b', 'A', 'B'

    Returns
    -------
    dict

    """
    j,k = 0,0
    for char in relator:
        if char == "a":
            j += 1
        else:
            k += 1
    return {"a":k, "b":-j, "A":-k, "B":j}


class hedden_knot:
    def __init__(self, p, q):
        self.p = p
        self.q = q
        self.k = q + 1
        self.relator = hedden_relator(p,q)
        self.abelianization = abelianization(self.relator)
        self.alexander_polynomial = symmetrize_laurent(positive_laurent(simplify_laurent(self.free_differential_a(self.relator))))


    def free_differential_a(self, relator):
        """
        Implements the Fox calculus algorithm for producing the Alexander
        polynomial of a knot 

        Parameters
        ----------
        relator : string in 'a', 'b', 'A', 'B'

        Returns
        -------
        dict (Alexander (Laurent) polynomial)

        """
        if len(relator) == 0:
            return
        if len(relator) == 1:
            if relator == "a":
                return {0:1}
            elif relator == "b":
                return {0:0}
            else:
                return {self.abelianization["a"]:-1}
        else:
            return add_laurent(self.free_differential_a(relator[0]), multiply_laurent({self.abelianization[relator[0]]:1}, self.free_differential_a(relator[1:])))
        
# def relator_pattern(relator):
#     pattern = []
#     i = 0
#     while i < len(relator):
#         j = 0
#         while j < len(relator)-i and relator[i] == relator[i+j]:
#             j+=1
#         pattern.append(j)
#         i += j
#     return pattern

# def find_1(p,q):
#     for i in range(p):
#         if (i*q) % p == 1:
#             return i
        
        
class simple_knot:
    def __init__(self, p, q, k):
        self.p = p
        self.q = q
        self.k = k
        self.relator = simple_relator(p,q,k)
        self.abelianization = abelianization(self.relator)
        self.alexander_polynomial = symmetrize_laurent(divide_laurent(positive_laurent(simplify_laurent(self.free_differential_a(self.relator))), allones(self.p), {0:0}))
        
    def free_differential_a(self, relator):
        if len(relator) == 0:
            return
        if len(relator) == 1:
            return {0:1} if relator == "a" else {0:0}
        else:
            return add_laurent(self.free_differential_a(relator[0]), multiply_laurent({self.abelianization[relator[0]]:1}, self.free_differential_a(relator[1:])))
        
        
        