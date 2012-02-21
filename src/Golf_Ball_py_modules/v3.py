"""
This module provides a reasonably complete vector class for 3d calculations.
"""

__author__="Igor Volobouev (i.volobouev@ttu.edu)"
__version__="0.2"
__date__ ="Jan 29 2008"

import math

class V3(object):
    "Basic 3d vector class."
    def __init__(self, x=0.0, y=0.0, z=0.0):
        # Use float for internal representation.
        # This avoids various problems with integer division.
        self.x = x*1.0
        self.y = y*1.0
        self.z = z*1.0

    def lengthSquared(self):
        "Spatial length of the vector, squared."
        return self.x*self.x + self.y*self.y + self.z*self.z

    def length(self):
        "Spatial length of the vector."
        return math.sqrt(self.lengthSquared())

    def direction(self):
        "Direction of the vector as a 3d vector with unit length."
        if (self.__nonzero__()):
            return self / self.length()
        else:
            return V3(1.0, 0.0, 0.0)

    def phi(self):
        "Azimuthal angle."
        return math.atan2(self.y, self.x)

    def cosTheta(self):
        "Cosine of the polar angle."
        if (self.__nonzero__()):
            return self.z/self.length()
        else:
            return 0.0

    def theta(self):
        "Polar angle."
        cosTheta = self.cosTheta()
        if (math.fabs(cosTheta) < 0.99):
            return math.acos(cosTheta)
        else:
            # acos would loose too much numerical precision
            th = math.asin(math.sqrt((self.x*self.x + \
                                      self.y*self.y)/self.lengthSquared()))
            if (cosTheta > 0.0):
                return th
            else:
                return math.pi - th

    def dot(self, other):
        "Scalar product with another vector, unit metric."
        return self.x*other.x + self.y*other.y + self.z*other.z

    def cross(self, other):
        "Cross product with another vector, unit metric."
        return V3(self.y*other.z - self.z*other.y,
                  self.z*other.x - self.x*other.z,
                  self.x*other.y - self.y*other.x)

    def project(self, other):
        "Projection onto the direction of another vector."
        othermag2 = other.lengthSquared()
        assert othermag2 > 0.0
        return other * (self.dot(other)/othermag2)

    def angle(self, other):
        "Angle between two vectors in radians."
        u = self.direction()
        v = other.direction()
        cosa = u.dot(v)
        if (math.fabs(cosa) < 0.99):
            return math.acos(cosa)
        else:
            # acos would loose too much numerical precision
            if (cosa > 0.0):
                return 2.0*math.asin(abs(v - u)/2.0)
            else:
                return math.pi - 2.0*math.asin(abs(-v - u)/2.0)

    def __repr__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + \
               ', ' + str(self.z) + ')'

    def __nonzero__(self):
        return int(self.x != 0.0 or self.y != 0.0 or self.z != 0.0)

    def __neg__(self):
        return V3(-self.x, -self.y, -self.z)

    def __pos__(self):
        # This is useful if one wants to get a copy of the object
        return V3(self.x, self.y, self.z)

    def __abs__(self):
        return self.length()

    def __eq__(self, other):
        return self.x == other.x and \
               self.y == other.y and \
               self.z == other.z

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        return V3(self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        return V3(self.x-other.x, self.y-other.y, self.z-other.z)

    def __mul__(self, other):
        return V3(self.x*other, self.y*other, self.z*other)

    def __rmul__(self, other):
        return self*other

    def __div__(self, other):
        if (other == 0.0):
            raise ZeroDivisionError, "3d vector divided by zero"
        return V3(self.x/other, self.y/other, self.z/other)

    # Mutators
    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        self.z *= other
        return self

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __idiv__(self, other):
        if (other == 0.0):
            raise ZeroDivisionError, "3d vector divided by zero"
        self.x /= other
        self.y /= other
        self.z /= other
        return self
