import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Sphere:
    def __init__(self, paramsDict):
        self.h = paramsDict["h"]
        self.k = paramsDict["k"]
        self.l = paramsDict["l"]
        self.r = paramsDict["r"]
        self.genMesh()
        self.toCarts()
    def genMesh(self):
        self.u, self.v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    def toCarts(self):
        self.X = self.r * np.cos(self.u) * np.sin(self.v) + self.h
        self.Y = self.r * np.sin(self.u) * np.sin(self.v) + self.k
        self.Z = self.r * np.cos(self.v) + self.l
    def getXYZ(self):
        return{"x":self.X, "y":self.Y, "z":self.Z}
