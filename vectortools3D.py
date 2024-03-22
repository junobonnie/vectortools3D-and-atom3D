# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 21:58:32 2024

@author: replica
"""

import math as m
import numpy as np

class Tensor:
    def __init__(self, xx, xy, xz,
                        yx, yy, yz,
                        zx, zy, zz):
        self.xx = xx
        self.xy = xy
        self.xz = xz
        
        self.yx = yx
        self.yy = yy
        self.yz = yz
        
        self.zx = zx
        self.zy = zy
        self.zz = zz
        
    def __abs__(self):
        return self.xx*self.yy*self.zz + self.xy*self.yz*self.zx + self.xz*self.yx*self.zy - self.xx*self.yz*self.zy - self.xy*self.yx*self.zz - self.xz*self.yy*self.zx
    
    def __add__(self, other_Tensor):
        return Tensor(self.xx + other_Tensor.xx, self.xy + other_Tensor.xy, self.xz + other_Tensor.xz,
                      self.yx + other_Tensor.yx, self.yy + other_Tensor.yy, self.yz + other_Tensor.yz,
                      self.zx + other_Tensor.zx, self.zy + other_Tensor.zy, self.zz + other_Tensor.zz)
    def __sub__(self, other_Tensor):
        return Tensor(self.xx - other_Tensor.xx, self.xy - other_Tensor.xy, self.xz - other_Tensor.xz,
                      self.yx - other_Tensor.yx, self.yy - other_Tensor.yy, self.yz - other_Tensor.yz,
                      self.zx - other_Tensor.zx, self.zy - other_Tensor.zy, self.zz - other_Tensor.zz)
    def __mul__(self, num):
        return Tensor(self.xx*num, self.xy*num, self.xz*num,
                      self.yx*num, self.yy*num, self.yz*num,
                      self.zx*num, self.zy*num, self.zz*num)
    
    def __rmul__(self, num):
        return Tensor(self.xx*num, self.xy*num, self.xz*num,
                      self.yx*num, self.yy*num, self.yz*num,
                      self.zx*num, self.zy*num, self.zz*num)
    
    def __truediv__(self, num):
        return Tensor(self.xx/num, self.xy/num, self.xz/num,
                      self.yx/num, self.yy/num, self.yz/num,
                      self.zx/num, self.zy/num, self.zz/num)
    
    def __pow__(self, num):
        result = Tensor(1, 0, 0,
                        0, 1, 0,
                        0, 0, 1)
        if num < 0:
            pow_tensor = self.inverse()
        else:
            pow_tensor = self
        for i in range(abs(num)):
            result = result.dot(pow_tensor)
        return result
    
    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False
        
    def __ne__(self, other):
        if str(self) != str(other):
            return True
        else:
            return False
        
    def __str__(self):
        return 'Tensor(' + str(self.xx) + ', ' + str(self.xy) + ', ' + str(self.xz) + ', ' + str(self.yx) + ', ' + str(self.yy) + ', ' + str(self.yz) + ', ' + str(self.zx) + ', ' + str(self.zy) + ', ' + str(self.zz) + ')'
    
    def inverse(self):
        if abs(self) == 0:
            return None
        else:
            R = np.linalg.inv(self.array())
            return Tensor(R[0][0], R[0][1], R[0][2], 
                          R[1][0], R[1][1], R[1][2],
                          R[2][0], R[2][1], R[2][2])
        
    def T(self):
        return Tensor(self.xx, self.yx, self.zx,
                      self.xy, self.yy, self.zy,
                      self.xz, self.yz, self.zz)
        
    def dot(self, other_object):
        if isinstance(other_object, Tensor):
            R = np.dot(self.array(), other_object.array())
            return Tensor(R[0][0], R[0][1], R[0][2], 
                          R[1][0], R[1][1], R[1][2],
                          R[2][0], R[2][1], R[2][2])
            
        elif isinstance(other_object, Vector):
            return Vector(self.xx*other_object.x + self.xy*other_object.y + self.xz*other_object.z,
                          self.yx*other_object.x + self.yy*other_object.y + self.yz*other_object.z,
                          self.zx*other_object.x + self.zy*other_object.y + self.zz*other_object.z)
        
    def tuple(self):
        return ((self.xx, self.xy, self.xz),
                (self.yx, self.yy, self.yz),
                (self.zx, self.zy, self.zz))
    
    def list(self):
        return [[self.xx, self.xy, self.xz],
                [self.yx, self.yy, self.yz],
                [self.zx, self.zy, self.zz]]
    
    def array(self):
        return np.array([[self.xx, self.xy, self.xz],
                        [self.yx, self.yy, self.yz],
                        [self.zx, self.zy, self.zz]])

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        
    def __abs__(self):
        return m.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def __add__(self, other_Vector):
        return Vector(self.x + other_Vector.x, self.y + other_Vector.y, self.z + other_Vector.z)
    
    def __sub__(self, other_Vector):
        return Vector(self.x - other_Vector.x, self.y - other_Vector.y, self.z - other_Vector.z)
    
    def __mul__(self, num):
        return Vector(self.x*num, self.y*num, self.z*num)
    
    def __rmul__(self, num):
        return Vector(self.x*num, self.y*num, self.z*num)
    
    def __truediv__(self, num):
        return Vector(self.x/num, self.y/num, self.z/num)
    
    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False
        
    def __ne__(self, other):
        if str(self) != str(other):
            return True
        else:
            return False
    
    def __str__(self):
        return 'Vector(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'
    
    def dot(self, other_Vector):
        return self.x*other_Vector.x + self.y*other_Vector.y + self.z*other_Vector.z
    
    def cross(self, other_Vector):
        return Vector(self.y*other_Vector.z - self.z*other_Vector.y, self.z*other_Vector.x - self.x*other_Vector.z, self.x*other_Vector.y - self.y*other_Vector.x)
    
    def tuple(self):
        return (self.x, self.y, self.z)
    
    def list(self):
        return [self.x, self.y, self.z]
    
    def array(self):
        return np.array([self.x, self.y, self.z])
    
def SO3_z(theta):
    return Tensor(m.cos(theta), -m.sin(theta), 0,
                  m.sin(theta), m.cos(theta), 0,
                  0, 0, 1)

def SO3_y(theta):
    return Tensor(m.cos(theta), 0, m.sin(theta),
                  0, 1, 0,
                  -m.sin(theta), 0, m.cos(theta))

def SO3_x(theta):
    return Tensor(1, 0, 0,
                  0, m.cos(theta), -m.sin(theta),
                  0, m.sin(theta), m.cos(theta))
if __name__ == '__main__':
    t1 = Tensor(1, 2, 3,
                4, 5, 6,
                7, 8, 10)
    t2 = SO3_z(m.pi/4)
    t3 = Tensor(2, 3, 4, 
                5, 6, 7,
                8, 9, 10)
    print(abs(t1))
    print(t1 + t2)
    print(t1 - t2)
    print(t1*2)
    print(2*t1)
    print(t1/2)
    print(t1**3)
    print(t1**(-1))
    print(t1 == t3)
    print(t1 != t3)
    print(t1)
    print(t1.inverse())
    print(t1.T())
    print(t1.dot(t2))
    print(t1.tuple())
    print(t1.list())
    print(t1.array())
    
    v1 = Vector(3, 4, 5)
    v2 = Vector(5, 4, 3)
    print(t1.dot(v1))
    print(abs(v1))
    print(v1 + v2)
    print(v1 - v2)
    print(2*v1)
    print(v1*2)
    print(v1/2)
    print(v1 == v2)
    print(v1 != v2)
    print(v1)
    print(v1.dot(v2))
    print(v1.cross(v2))
    print(v1.tuple())
    print(v1.list())
    print(v1.array())
