# -*- coding: utf-8 -*-
"""
Geom 3D point.

Created on Wed Jun  5 13:55:30 2019

@author: Rybakov
"""

import mth
import math

#---------------------------------------------------------------------------------------------------
# Class point.
#---------------------------------------------------------------------------------------------------

class Vector:

#---------------------------------------------------------------------------------------------------

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        """
        Constructor.

        Arguments:
            x -- x coordinate,
            y -- y coordinate,
            z -- z coordinate.
        """

        self.X = x
        self.Y = y
        self.Z = z

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Convert to string.

        Result:
            String.
        """

        return '(%f, %f, %f)' % (self.X, self.Y, self.Z)

#---------------------------------------------------------------------------------------------------

    def __neg__(self):
        """
        Negated vector.

        Result:
            Negated vector.
        """

        return Vector(-self.X, -self.Y)

#---------------------------------------------------------------------------------------------------

    def __lt__(self, v):
        """
        Overload '<'.

        Arguments:
            v -- vector.

        Result:
            Boolean value.
        """

        return (self.X, self.Y, self.Z) < (v.X, v.Y, v.Z)

#---------------------------------------------------------------------------------------------------

    def __add__(self, v):
        """
        Overload '+'.

        Arguments:
            v -- vector.

        Result:
            Sum of vectors.
        """

        return Vector(self.X + v.X, self.Y + v.Y, self.Z + v.Z)

#---------------------------------------------------------------------------------------------------

    def __sub__(self, v):
        """
        Overload '-'.

        Arguments:
            v -- vector.

        Result:
            Sub of two vectors.
        """

        return Vector(self.X - v.X, self.Y - v.Y, self.Z - v.Z)

#---------------------------------------------------------------------------------------------------

    def __mul__(self, k):
        """
        Overload '*'.

        Arguments:
            k -- coefficient.

        Result:
            New vector.
        """

        return Vector(k * self.X, k * self.Y, k * self.Z)

#---------------------------------------------------------------------------------------------------

    def __rmul__(self, k):
        """
        Overload '*'.

        Arguments:
            k -- coefficient.

        Result:
            New vector.
        """

        return self * k

#---------------------------------------------------------------------------------------------------

    def Copy(self):
        """
        Get copy.

        Result:
            Copy.
        """

        return Vector(self.X, self.Y, self.Z)

#---------------------------------------------------------------------------------------------------

    def Tuple(self, n = 3):
        """
        Get tuple of coordinates.

        Arguments:
            n -- count of coordinates.

        Result:
            Tuple.
        """

        if n == 1:
            return (self.X)
        elif n == 2:
            return (self.X, self.Y)
        elif n == 3:
            return (self.X, self.Y, self.Z)
        else:
            raise Exception('wrong dimensiality')

#---------------------------------------------------------------------------------------------------

    def Mod2(self):
        """
        Square of module.

        Result:
            Square of module.
        """

        return self.X * self.X + self.Y * self.Y + self.Z * self.Z

#---------------------------------------------------------------------------------------------------

    def Mod(self):
        """
        Module.

        Result:
            Module.
        """

        return math.sqrt(self.Mod2())

#---------------------------------------------------------------------------------------------------

    def Orthogonal(self):
        """
        Orthogonal vector.

        Result:
            Orthogonal vector.
        """

        return Vector(-self.Y, self.X)

#---------------------------------------------------------------------------------------------------

    def Scale(self, k):
        """
        Scale vector.

        Arguments:
            k -- coefficient.
        """

        self.X *= k
        self.Y *= k
        self.Z *= k

#---------------------------------------------------------------------------------------------------

    def Scaled(self, k):
        """
        Get scaled vector.

        Arguments:
            k -- coefficient.

        Result:
            Scaled vector.
        """

        r = self.Copy()
        r.Scale(k)

        return r

#---------------------------------------------------------------------------------------------------

    def Normalize(self):
        """
        Normalize to size.

        Arguments:
            size -- vector size.
        """

        self.Scale(1.0 / self.Mod())

#---------------------------------------------------------------------------------------------------

    def Normalized(self):
        """
        Get normalized vector.

        Result:
            Normalized vector.
        """

        r = self.Copy()
        r.Normalize()

        return r

#---------------------------------------------------------------------------------------------------

    def Negate(self):
        """
        Negate.
        """

        self.X = -self.X
        self.Y = -self.Y
        self.Z = -self.Z

#---------------------------------------------------------------------------------------------------

    def dist_to(self, v):
        """
        Distance to.

        Arguments:
            v -- vector.

        Result:
            Distance.
        """

        return mth.hypot(self.X - v.X, self.Y - v.Y, self.Z - v.Z)

#---------------------------------------------------------------------------------------------------

    def DotProduct(self, v):
        """
        Dot product with another vector.

        Arguments:
            v -- vector.

        Result:
            Dot product.
        """

        return self.X * v.X + self.Y * v.Y + self.Z * v.Z

#---------------------------------------------------------------------------------------------------

    def CosAngle(self, v):
        """
        Cosine of angle with another vector.

        Arguments:
            v -- vector.

        Result:
            Angle cosine.
        """

        r = self.DotProduct(v) / (self.Mod() * v.Mod())

        if r > 1.0:
            r = 1.0

        return r

#---------------------------------------------------------------------------------------------------

    def VectorProduct(self, v):
        """
        Vector product.

        Arguments:
            v -- vector.

        Result:
            Vector product.
        """

        return Vector(self.Y * v.Z - self.Z * v.Y,
                      self.Z * v.X - self.X * v.Z,
                      self.X * v.Y - self.Y * v.X)

#---------------------------------------------------------------------------------------------------

    def NearestPoint(self, target):
        """
        Find nearest point to target.

        Argumnets:
            target -- target.

        Result:
            Point.
        """

        if isinstance(target, Sphere):
            return target.C + (self - target.C).Normalized().Scaled(target.R)
        else:
            raise Exception('unknown target type')

#---------------------------------------------------------------------------------------------------
# Class triangle.
#---------------------------------------------------------------------------------------------------

class Triangle:

#---------------------------------------------------------------------------------------------------

    def __init__(self, a, b, c):
        """
        Constructor.

        Arguments:
            a -- a point,
            b -- b point,
            c -- c point.
        """

        self.Points = [a, b, c]

#---------------------------------------------------------------------------------------------------

    def Area(self):
        """
        Area.

        Result:
            Area.
        """

        v1 = self.Points[1] - self.Points[0]
        v2 = self.Points[2] - self.Points[0]
        v = v1.VectorProduct(v2)

        return abs(0.5 * v.Mod())

#---------------------------------------------------------------------------------------------------

    def R(self):
        """
        Outer circle radius.

        Result:
            Outer circle radius.
        """

        s1 = self.Points[0].DistTo(self.Points[1])
        s2 = self.Points[0].DistTo(self.Points[2])
        s3 = self.Points[1].DistTo(self.Points[2])

        return (s1 * s2 * s3) / (4.0 * self.Area())

#---------------------------------------------------------------------------------------------------
# Sphere.
#---------------------------------------------------------------------------------------------------

class Sphere:

#---------------------------------------------------------------------------------------------------

    def __init__(self, c, r):
        """
        Constructor.

        Arguments:
            c -- center,
            r -- radius.
        """

        self.C = c
        self.R = r

#---------------------------------------------------------------------------------------------------

    def LoPoint(self):
        """
        Get lo point.

        Result:
            Lo point.
        """

        return self.C - Vector(self.R, self.R, self.R)

#---------------------------------------------------------------------------------------------------

    def HiPoint(self):
        """
        Get hi point.

        Result:
            Hi point.
        """

        return self.C + Vector(self.R, self.R, self.R)

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    pass

#---------------------------------------------------------------------------------------------------
