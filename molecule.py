import atom
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import atomicradii as arad
import periodic_table as pt
import sphere as sph

# rgp12dec2024 this class getting too big
#  maybe pull all nontrivial algorithms out of this class into an external library #
#  apply each algorithm to molecule via strategy pattern  #

class Molecule (object):

    def __init__(self, fileN):
        self.atoms = {}
        self.filename = fileN
        self.connections = set() # was []
        self.maxBondLength = 1.6
        self.formulator = {}
        with open(fileN, mode = 'r') as file:
            next(file)  #don't read the first two lines of the .xyz files
            next(file)  #into the hash table. These are just headers
            for index, line in enumerate(file, 0):
                tmpList = line.split()
                tmpList.append(index)
                #print(tmpList)
                if (len(tmpList) == 5):
                    self.atoms[index] = atom.Atom(tmpList) #tmpList
        
        self.npdata = np.loadtxt(fileN, skiprows=2, usecols = (1,2,3))
        self.setConnections()
        self.thetax = 0
        self.thetay = 0
        self.thetaz = 0

    def origCoords():
        return np.loadtxt(self.filename, skiprows=2, usecols = (1,2,3))

    def plot(self):
        self.fig = plt.figure()
        self.canvas = FigureCanvasTkAgg(self.molecules[self.filenames[-1]].fig, master = self.frame1)
        self.ax  = self.fig.add_subplot(231, projection = '3d')
        plt.axis('off')
        self.ax.scatter(self.pData.x,self.pData.y,self.pData.z, marker = "o")
        for (i, j) in self.connections:
            self.ax.plot([self.atoms[i].getx(),self.atoms[j].getx()],[self.atoms[i].gety(),self.atoms[j].gety()],[self.atoms[i].getz(),self.atoms[j].getz()], color="black")
        plt.show()

    def projPlot(self):
        self.projAx = []
        self.projAx.append(self.fig.add_subplot(234))
        self.projAx.append(self.fig.add_subplot(235))
        self.projAx.append(self.fig.add_subplot(236))
        datxy = self.projXY()
        datxz = self.projXZ()
        datyz = self.projYZ()        
        self.projAx[0].scatter(datxy[0], datxy[1])
        self.projAx[0].title.set_text("xy projection")
        for (i, j) in self.connections:
            self.projAx[0].plot([datxy[0][i], datxy[0][j]], [datxy[1][i], datxy[1][j]], color="black", alpha = 0.25)
        self.projAx[1].scatter(datxz[0], datxz[2])
        self.projAx[1].title.set_text("xz projection")
        for (i, j) in self.connections:
            self.projAx[1].plot([datxz[0][i], datxz[0][j]], [datxz[2][i], datxz[2][j]], color="black", alpha = 0.25)
        self.projAx[2].scatter(datyz[1], datyz[2])
        self.projAx[2].title.set_text("yz projection")
        for (i, j) in self.connections:
            self.projAx[2].plot([datyz[1][i], datyz[1][j]], [datyz[2][i], datyz[2][j]], color="black", alpha = 0.25)
        plt.show()
        
    def projXY(self):
        Pxy = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
        return (Pxy*(self.npdata.T)).tolist()
        
    def projXZ(self):
        Pxz = np.matrix([[1, 0, 0], [0, 0, 0], [0, 0, 1]])
        return (Pxz*(self.npdata.T)).tolist()
        
    def projYZ(self):
        Pyz = np.matrix([[0, 0, 0], [0, 1, 0], [0, 0, 1]])
        return (Pyz*(self.npdata.T)).tolist()
        
    def findClosestAtom(self, xyzVec):
        d = 10.0
        for atm in self.atoms:
            tmpD = np.sqrt((self.atoms[atm].getx() - xyzVec[0])**2 + (self.atoms[atm].gety() - xyzVec[1])**2 + (self.atoms[atm].getz() - xyzVec[2])**2)
            if (tmpD < d):
                d = tmpD
                atmNo = atm
        return atmNo

    def findClosestBond(self, xyzVec):
        d = 10.0
        for i, j in self.connections:
            tmpD1 = np.sqrt((self.atoms[i].getx() - xyzVec[0])**2 + (self.atoms[i].gety() - xyzVec[1])**2 + (self.atoms[i].getz() - xyzVec[2])**2)
            tmpD2 = np.sqrt((self.atoms[j].getx() - xyzVec[0])**2 + (self.atoms[j].gety() - xyzVec[1])**2 + (self.atoms[j].getz() - xyzVec[2])**2)
            tmpD = (tmpD1 + tmpD2)/2
            if (tmpD < d):
                d = tmpD
                bond = (i, j)
        return bond

    def optInternalRadii(self):
        for key in self.formulator.keys():
            for q in np.linspace(0, 3, 31):#atom type being optimized
                for r in np.linspace(0, 3, 31):#
                    pass

    def calculateCenterOfMass(self):#, atomnumbers):
        self.centerofmass = {}
        x = 0.0
        y = 0.0
        z = 0.0
        for index, atm in enumerate(self.atoms):  
            x = x + float(self.atoms[atm].x)#[1]
            y = y + float(self.atoms[atm].y)#[2]
            z = z + float(self.atoms[atm].z)#[3]
        self.centerofmass["x"] = x/(index + 1)
        self.centerofmass["y"] = y/(index + 1)
        self.centerofmass["z"] = z/(index + 1)     


    def centerMoleculeAtOrigin(self):
        self.calculateCenterOfMass()
        for i in range(np.shape(self.npdata)[0]):
            self.npdata[i, 0] = self.npdata[i, 0] - self.centerofmass["x"]
            self.npdata[i, 1] = self.npdata[i, 1] - self.centerofmass["y"]
            self.npdata[i, 2] = self.npdata[i, 2] - self.centerofmass["z"]
        self.updateMoleculeFromNpdata()
        self.setConnections()
        return [[self.npdata[:,0].tolist(), self.npdata[:,1].tolist(), self.npdata[:,2].tolist()], self.projXY(), self.projXZ(), self.projYZ()]

    def alignZaxis(self, aPoint):
        apoint = np.array(aPoint[2:])

        pt = np.array([float(aPoint[2]), float(aPoint[3]), float(aPoint[4])])#.reshape(-1,1)

        unitZ = np.array([0.0000, 0.0000, 1.00]) # this could just be (0, 0 , 1) if only aligning here

        rotMat = self.rotation_matrix_from_vectors(apoint, unitZ)

        self.npdata = np.matmul(self.npdata, rotMat)

        self.updateMoleculeFromNpdata()

        self.setConnections()

        return [[self.npdata[:,0].tolist(), self.npdata[:,1].tolist(), self.npdata[:,2].tolist()], self.projXY(), self.projXZ(), self.projYZ()]

    

    def alignZaxis2(self, aPoint):
        apoint = np.array(aPoint[2:])
        pt = np.array([float(aPoint[2]), float(aPoint[3]), float(aPoint[4])])#.reshape(-1,1)
        Pxy = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
        Pxz = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1]])
        Pyz = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]])
        r = np.sqrt((pt[0])**2 + (pt[1])**2 + (pt[2])**2)
        newpoint = np.array([0.0000, 0.0000, 1.00]) # this could just be (0, 0 , 1) if only aligning here
        npxy = np.dot(Pxy, newpoint)
        apxy = np.dot(Pxy, pt)
        angz = np.arccos(np.clip(np.dot(npxy, apxy), -1.0, 1.0))*(180/np.pi)
        npxz = np.dot(Pxz, newpoint)
        apxz = np.dot(Pxz, pt)
        angy = np.arccos(np.clip(np.dot(npxz, apxz), -1.0, 1.0))*(180/np.pi)
        npyz = np.dot(Pyz, newpoint)
        apyz = np.dot(Pyz, pt)
        angx = np.arccos(np.clip(np.dot(npyz, apyz), -1.0, 1.0))*(180/np.pi)
        plottingData = self.genRotXYZ([angx, angy, 0])
        return plottingData

    def totalMarcusDifference(self, dicz):
        print("totalMarcusDifference() from molecule module")
        Vi = 0
        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                d = ((((self.atoms[i].getx())-(self.atoms[j].getx()))**2)+(((self.atoms[i].gety())-(self.atoms[j].gety()))**2)+(((self.atoms[i].getz())-(self.atoms[j].getz()))**2))**(0.50)
                ri = arad.crystalRadii[self.atoms[i].getSymbol()]
                rj = arad.crystalRadii[self.atoms[i].getSymbol()]
                if (ri >= rj):
                    R = ri
                    r = rj
                else:
                    R = rj
                    r = ri
                h = (ri + R - d)/2
                if h < 0.0:
                    print("found a negative h value")
                    continue
                if R <= (d+ri):
                    vi = ((np.pi)*((R +r -d)**2)*((d**2)+(2*d*r)-(3*(r**2))+(2*d*R)+(6*r*R)-(3*(R**2))))/(12*(d))
                else:
                    vi = (4/3)*np.pi*(r**3)
                Vi = Vi + vi
        dicz[str(h)+ "_" + str(k)] = Vi

    def CHinternalRadiiOpt(self, dicz):
        print("ran molecule version")
        f = open('/home/robert/Documents/Literature/manuscripts/histericall__1/StericVolumeFittingExperiments/rotFitLR8/volCalc.dat', 'w')
        rads = {}
        rads["C"] = np.linspace(0.6, 1.0, 5)
        rads["H"] = np.linspace(0.6, 1.0, 5)
        for C in rads["C"]:
            for H in rads["H"]:
                f.write("C radius: {}, H radius: {}\n".format(C, H))
                Vi = 0
                for i in range(len(self.atoms)):
                    for j in range(i + 1, len(self.atoms)):
                        d = ((((self.atoms[i].getx())-(self.atoms[j].getx()))**2)+(((self.atoms[i].gety())-(self.atoms[j].gety()))**2)+(((self.atoms[i].getz())-(self.atoms[j].getz()))**2))**(0.50)
                        if self.atoms[i].getSymbol() == "C":
                            ri = C
                        else:
                            ri = H
                        if self.atoms[j].getSymbol() == "C":
                            rj = C
                        else:
                            rj = H
                        if (ri >= rj):
                            R = ri
                            r = rj
                        else:
                            R = rj
                            r = ri
                        h = (r + R - d)/2
                        if h < 0.0:
                            vi = 0
                        elif R <= (d+r):
                            vi = ((np.pi)*((R + r -d)**2)*((d**2)+(2*d*r)-(3*(r**2))+(2*d*R)+(6*r*R)-(3*(R**2))))/(12*(d))
                        else:
                            vi = (4/3)*np.pi*(r**3)
                        f.write("\tAtom: {}\tvi: {}\n".format(i, vi))
                        Vi = Vi + vi
                f.write("\tTotal Atomic Overlap Volume: {}\n".format(Vi))
                dicz[str(np.round(C, 2))+ "_" + str(np.round(H, 2))] = Vi #{"v" : Vi, "h" : h}
        f.close()
    def printInternalScan(self, dicz):#ethane connectivity     0 -->  1, (3,4,7),    1 -->  0, (2,5,6)
        print("InternalScan() from molecule")
        sfig, sax = plt.subplots(subplot_kw={"projection": "3d"})
        spheres = {}
        for i in range(len(self.atoms)):
            spheres[i] = True
        n = 0
        rads = {}
        rads["C"] =  [0.2]    #np.linspace(0, 2, 11)
        rads["H"] =  [3.0]    #np.linspace(0, 2, 11)
        for C in rads["C"]:
            for H in rads["H"]:
                Vi = 0
                for i in range(len(self.atoms)):
                    print("i: {}".format(i))
                    for j in range(i + 1, len(self.atoms)):
                        if j != i:
                            print("j: {}".format(j))
                            d = ((((self.atoms[i].getx())-(self.atoms[j].getx()))**2)+(((self.atoms[i].gety())-(self.atoms[j].gety()))**2)+(((self.atoms[i].getz())-(self.atoms[j].getz()))**2))**(0.50)
                            if self.atoms[i].getSymbol() == "C":
                                colori = "blue"
                                ri = C
                            else:
                                ri = H
                                colori = "red"
                            if self.atoms[j].getSymbol() == "C":
                                rj = C
                                colorj = "blue"
                            else:
                                rj = H
                                colorj = "red"
                            # next two lines insure R is always the biggest    
                            if (ri >= rj):
                                R = ri
                                r = rj
                            else:
                                R = rj
                                r = ri
                            h = (r + R - d)/2
                            if h < 0.0:
                                continue                            
                            elif R <= (d+r):
                                vi = ((np.pi)*((R + r -d)**2)*((d**2)+(2*d*r)-(3*(r**2))+(2*d*R)+(6*r*R)-(3*(R**2))))/(12*(d))
                            else:
                                vi = (4/3)*np.pi*(r**3)
                            Vi = Vi + vi                                                        
                            if spheres[i]:
                                sp = sph.Sphere({"h" : self.atoms[i].getx(), "k" : self.atoms[i].gety(), "l" : self.atoms[i].getz(), "r" : ri})
                                spcoords = sp.getXYZ()
                                sax.plot_wireframe(spcoords["x"], spcoords["y"], spcoords["z"], color = colori, linewidth = 0.5)
                                spheres[i] = False
                            
                            if spheres[j]:
                                sphe = sph.Sphere({"h" : self.atoms[j].getx(), "k" : self.atoms[j].gety(), "l" : self.atoms[j].getz(), "r" : rj})                            
                                sphecoords = sphe.getXYZ()                            
                                sax.plot_wireframe(sphecoords["x"], sphecoords["y"], sphecoords["z"], color = colorj, linewidth = 0.5)
                                spheres[j] = False
                            n = n + 1

        tmpBonds = self.getBonds()
        for bond in tmpBonds:
            sax.plot(bond[0],bond[1],bond[2], color="b")

        sax.set_xlim(xmin = -2.0, xmax = 2.0)
        sax.set_ylim(ymin = -2.0, ymax = 2.0)
        sax.set_zlim(zmin = -2.0, zmax = 2.0)
        sax.annotate('VDW overlap volume: {}'.format(round(Vi, 2)), xy=(1, 0), xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')
        plt.axis("off")
        plt.show()
                                                
                
    def genDF(self):
        tmpDict = {}
        tmpDict["atom"] = []
        tmpDict["x"] = []
        tmpDict["y"] = []
        tmpDict["z"] = []
        for atm in self.atoms:
            tmpDict["atom"].append(self.atoms[atm].getSymbol())
            tmpDict["x"].append(self.atoms[atm].getx())
            tmpDict["y"].append(self.atoms[atm].gety())
            tmpDict["z"].append(self.atoms[atm].getz())
        return pd.DataFrame(tmpDict)

    def genX(self):
        tmpX = []
        for atm in self.atoms:
            tmpX.append(self.atoms[atm].getx())
        return tmpX

    def genY(self):
        tmpY = []
        for atm in self.atoms:
            tmpY.append(self.atoms[atm].gety())
        return tmpY

    def genZ(self):
        tmpZ = []
        for atm in self.atoms:
            tmpZ.append(self.atoms[atm].getz())
        return tmpZ

    def rotation_matrix_from_vectors(self, vec1, vec2):

        """ Find the rotation matrix that aligns vec1 to vec2

        :param vec1: A 3d "source" vector

        :param vec2: A 3d "destination" vector

        :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.

        """

        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)

        v = np.cross(a, b)

        c = np.dot(a, b)

        s = np.linalg.norm(v)

        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

        return rotation_matrix

    def genRotXYZ(self, angles):  
        self.thetax = (np.pi/180) * angles[0]
        self.thetay = (np.pi/180) * angles[1]
        self.thetaz = (np.pi/180) * angles[2]
        R = np.matmul(np.matmul(self.rotMatX(), self.rotMatY()), self.rotMatZ())
        self.npdata = np.matmul(self.npdata, R)
        self.updateMoleculeFromNpdata()
        return [[self.npdata[:,0].tolist(), self.npdata[:,1].tolist(), self.npdata[:,2].tolist()], self.projXY(), self.projXZ(), self.projYZ()]

    def updateMoleculeFromNpdata(self):
        for index, atm in enumerate(self.atoms):
            self.atoms[atm].x = self.npdata[index][0]
            self.atoms[atm].y = self.npdata[index][1]
            self.atoms[atm].z = self.npdata[index][2]
        self.setConnections()
        
    def genAtoms(self):
        tmpAt = []
        for atm in self.atoms:
            tmpAt.append(self.atoms[atm].getSymbol())
        return tmpAt

    def genColors(self):
        return np.array([arad.colors[element] for element in self.genAtoms()])

    def setConnections(self):
        for i in range(len(self.atoms)):
            if self.atoms[i].getSymbol() not in self.formulator.keys():
                self.formulator[self.atoms[i].getSymbol()] = 1
            else:
                self.formulator[self.atoms[i].getSymbol()] += 1
            for j in range(i+1, len(self.atoms)):
                d = ((((self.atoms[i].getx())-(self.atoms[j].getx()))**2)+(((self.atoms[i].gety())-(self.atoms[j].gety()))**2)+(((self.atoms[i].getz())-(self.atoms[j].getz()))**2))**(0.50)
                if d < self.maxBondLength:
                    self.connections.add((i, j))
        
    def getBonds(self):
        L = []
        for i, j in self.connections:            
            L.append([[self.atoms[i].getx(),self.atoms[j].getx()],[self.atoms[i].gety(),self.atoms[j].gety()],[self.atoms[i].getz(),self.atoms[j].getz()]])
        return L
        

    def rotMatX(self):
        return np.array([[1, 0, 0],
            [0, np.cos(self.thetax), -np.sin(self.thetax)],
            [0, np.sin(self.thetax), np.cos(self.thetax)]])

    def rotMatY(self):
        return np.array([[np.cos(self.thetay), 0, np.sin(self.thetay)],
            [0, 1, 0],
            [-np.sin(self.thetay), 0, np.cos(self.thetay)]])

    def rotMatZ(self):
        return np.array([[np.cos(self.thetaz), -np.sin(self.thetaz), 0],
            [np.sin(self.thetaz), np.cos(self.thetaz), 0],
            [0, 0, 1]])

    def getnumpyData(self):
        return npdata

    def printRotdata(self):
        for i in range(self.rotdata[:,0].size):
            print("{}   {}   {}".format(self.rotdata[i][0], self.rotdata[i][1], self.rotdata[i][2]))

    def calcMomentsOfInertia(self):
        self.momentsOfInertia = {}
        self.momentsOfInertia["Ix"] = self.getMolWt() * (np.sqrt((atm.gety())**2 + (atm.getz())**2))
        self.momentsOfInertia["Iy"] = self.getMolWt() * (np.sqrt((atm.getx())**2 + (atm.getz())**2))
        self.momentsOfInertia["Iz"] = self.getMolWt() * (np.sqrt((atm.getx())**2 + (atm.gety())**2))
        
    def getMolWt(self):
        molWt = 0
        for atm in self.atoms:
            molWt = molWt + pt.elements[atm.getSymbol()]["AtomicMass"]
        return molWt

    def dirCos(self):
        axes = ["x", "y", "z"]
        for atm in rotor:
            for i in axes:
                lmda = atm[i]/(np.sqrt(atm["x"]**2 + atm["y"]**2 + atm["z"]))
                tmp1 = (Am*lmda**2)/i
            


