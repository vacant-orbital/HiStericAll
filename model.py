import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
import molecule
import atom
import atomicradii as arad
from sklearn.preprocessing import StandardScaler  

class Main_Model:
    def __init__(self):
        self.data = {}

    def runLinReg(self):
        pass
        
        
class System_Model:
    def __init__(self, f):
        self.molecule = molecule.Molecule(f)
        self.defaultTmpAtom = atom.Atom()
        self.tmpAtom = self.defaultTmpAtom
        self.centroidCounter = 0

        self.centroidCalcList = []
        self.centroidCalcResult = {}
        self.rotorCalcList = []
        
        self.dummycount = 0
        self.dummies = {}
        self.excludedAtomCounter = 0
        self.excludedAtoms = []
        self.rm = np.linspace(0, 5, 51)


    def send3Dminmax(self):
        x ={}
        dif = 0.05 * np.abs(self.molecule.npdata[:,0].max() - self.molecule.npdata[:,0].min())
        x["xmin"] = self.molecule.npdata[:,0].min() - dif
        x["xmax"] = self.molecule.npdata[:,0].max() + dif
        dif = 0.05 * np.abs(self.molecule.npdata[:,1].max() - self.molecule.npdata[:,1].min())
        x["ymin"] = self.molecule.npdata[:,1].min() - dif
        x["ymax"] = self.molecule.npdata[:,1].max() + dif
        dif = 0.05 * np.abs(self.molecule.npdata[:,2].max() - self.molecule.npdata[:,2].min())
        x["zmin"] = self.molecule.npdata[:,2].min() - dif
        x["zmax"] = self.molecule.npdata[:,2].max() + dif
        return x
       
    def sendxyz(self):
        '''sends the xyz coords for plotting the molecule as a 3D projection'''
        return [self.molecule.genX(), self.molecule.genY(), self.molecule.genZ(), self.molecule.genColors()]

    def sendDF(self):
        return self.molecule.genDF()

    def centerMolecule(self):
        x = self.molecule.centerMoleculeAtOrigin()
        x.append(self.molecule.genColors())
        return x

    def alignZaxis(self, coords):
        x = self.molecule.alignZaxis(coords)
        x.append(self.molecule.genColors())
        return x

    def sendBonds(self):
        return self.molecule.getBonds()

    def sendxy(self):
        return self.molecule.projXY()

    def sendxz(self):
        return self.molecule.projXZ()

    def sendyz(self):
        return self.molecule.projYZ()

    def sendColors(self):
        return self.molecule.genColors()

    def defineTempAtom(self, atomNumber):
        self.tmpAtom = self.molecule.atoms[atomNumber]

    def sendClosestAtomInfo(self, inputCoords):
        atmNo = self.molecule.findClosestAtom(inputCoords)
        self.defineTempAtom(atmNo)
        return [self.tmpAtom.getIndex(), self.tmpAtom.getSymbol(), self.tmpAtom.getx(), self.tmpAtom.gety(), self.tmpAtom.getz()]

    def sendClosestBondInfo(self, inputCoords):
        atmNo = self.molecule.findClosestBond(inputCoords)
        self.defineTempAtom(atmNo)
        return [self.tmpAtom.getIndex(), self.tmpAtom.getSymbol(), self.tmpAtom.getx(), self.tmpAtom.gety(), self.tmpAtom.getz()]

    def sendResetTmpAtom(self, ar):
        self.centroidCalcList.append(ar)# this appends current atom, x, y, z from boxes frrom  def appendPointToCentroidCalculator(self):
        return [self.tmpAtom.getx(), self.tmpAtom.gety(), self.tmpAtom.getz(), self.tmpAtom.getSymbol()]

    def sendResetTmpAtomUpdateRotor(self, ar):
        self.rotorCalcList.append(ar)# this appends current atom, x, y, z from boxes from  def updateRotorCenterCalculator(ar):
        return [self.tmpAtom.getx(), self.tmpAtom.gety(), self.tmpAtom.getz(), self.tmpAtom.getSymbol()]

    def sendResetTmpAtomUpdateRotorCenter(self, ar):
        self.rotorCalcList.insert(0, ar)
        return [self.tmpAtom.getx(), self.tmpAtom.gety(), self.tmpAtom.getz(), self.tmpAtom.getSymbol()]
    
    def updateCentroidCalc(self, ar):
        self.centroidCalcList.append(ar)

    def updateRotorCalc(self, ar):
        self.centroidCalcList.append(ar)
    
    def calculateCentroid(self):#, atomnumbers):
        self.centroidCounter += 1
        self.centroidCalcResult[self.centroidCounter] = {}
        x = 0.0
        y = 0.0
        z = 0.0
        for coords in self.centroidCalcList:  #make this a dict of atoms!!!
            #print("during calc x: {}, y: {}, z: {}".format(x,y,z))
            x = x + float(coords[2])#[1]
            y = y + float(coords[3])#[2]
            z = z + float(coords[4])#[3]
            print("type of coords[4]: {}".format(type(coords[0])))
            print("coords[4]: {}".format(coords[0]))
            #self.excludedAtoms.append(map(float, coords[4][4]))
            self.excludedAtoms.append(coords[0])
        #print("after calc x: {}, y: {}, z: {}".format(x,y,z))
        self.centroidCalcResult[self.centroidCounter]["x"] = str(x/len(self.centroidCalcList))#atomnumbers)
        self.centroidCalcResult[self.centroidCounter]["y"] = str(y/len(self.centroidCalcList))
        self.centroidCalcResult[self.centroidCounter]["z"] = str(z/len(self.centroidCalcList))
        self.centroidCalcResult[self.centroidCounter]["id"] = "centroid_" + str(self.centroidCounter)
        #self.appendPointToCentroidResults()       
        for atm in self.centroidCalcList:
            print("atm: {}".format(atm))
        self.centroidCalcList.clear()
        print("Excluded Atoms: {}".format(self.excludedAtoms))
            
        return [self.centroidCalcResult[self.centroidCounter]["id"], self.centroidCalcResult[self.centroidCounter]["x"], self.centroidCalcResult[self.centroidCounter]["y"], self.centroidCalcResult[self.centroidCounter]["z"]]

    def calculateInternalRotor(self):
        pass
        
    def updateModelWithSelectedCentroids(self, kys):        
        for key in kys:
            self.dummycount += 1
            #self.dummies[self.dummycount] = {"symbol" : "centroid" + str(self.dummycount) , "x" : self.centroidCalcResult[-1]["x"], "y" : self.centroidCalcResult["y"], "z" : self.centroidCalcResult["z"]}
            self.dummies[self.dummycount] = self.centroidCalcResult.pop(int(key))

    def rotateMolecule(self, ar):
        x = self.molecule.genRotXYZ(ar)
        x.append(self.molecule.genColors())
        return x
    
    def rotation_matrix_from_vectors(vec1, vec2):
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

    def vanderwaalsDifSpheres(self, ): # CHECK BEFORE USED  before was only set to use the first dummy atom
        # put all algorithms in a function library to be incorporated in Strategy pattern WRT molecule
        print("vanderwaalsDifSpheres() from model module")
        for dum in self.dummies:
            x1 = float(self.dummies[dum]["x"])
            y1 = float(self.dummies[dum]["y"])
            z1 = float(self.dummies[dum]["z"])
            for R in self.rm:
                V = 0.0
                Vi = 0.0
                for atm in self.molecule.atoms:
                    if (str(self.molecule.atoms[atm].getIndex()) not in self.excludedAtoms):#[:,4]
                        #ri = arad.crystalRadii[self.molecule.atoms[atm].getSymbol()]
                        ri = arad.vanderwaals[self.molecule.atoms[atm].getSymbol()]
                        x2 = float(self.molecule.atoms[atm].getx())
                        y2 = float(self.molecule.atoms[atm].gety())
                        z2 = float(self.molecule.atoms[atm].getz())
                        d = ((((x2)-(x1))**2)+(((y2)-(y1))**2)+(((z2)-(z1))**2))**(0.50)
                        h = (ri + R - d)/2
                        if h < 0.0:
                            continue
                        if R <= (d+ri):
                            vi = ((np.pi)*((R +ri -d)**2)*((d**2)+(2*d*ri)-(3*(ri**2))+(2*d*R)+(6*ri*R)-(3*(R**2))))/(12*(d))
                        else:
                            vi = (4/3)*np.pi*(ri**3)
                        Vi = Vi + vi
                dicy[R]= Vi

    def vanderwaalsDifScan(self, dicy):
        #not used yet put all algorithms in a function library to be incorporated in Strategy pattern WRT molecule
        print("vanderwaalsDifScan() from model module")
        for dum in self.dummies:
            x1 = float(self.dummies[dum]["x"])
            y1 = float(self.dummies[dum]["y"])
            z1 = float(self.dummies[dum]["z"])
            for R in self.rm:
                dicz = dicy[R]
                V = 0.0
                Vi = 0.0
                rads["C"] = np.linspace(0.5, 3, 31)
                rads["H"] = np.linspace(0.5, 3, 31)
                #vdwCube = np.array(len(self.atoms), len(rads["C"]), len(rads["H"]))
                #grid = np.meshgrid(varR, varr)
                for C in rads["C"]:
                    for H in rads["H"]:
                        Vi = 0
                        for i in range(len(self.atoms)):
                            for j in range(1, len(self.atoms)):
                                if (str(self.molecule.atoms[atm].getIndex()) not in self.excludedAtoms):#[:,4]
                                    if j != i:
                                        d = ((((self.atoms[i].getx())-(self.atoms[j].getx()))**2)+(((self.atoms[i].gety())-(self.atoms[j].gety()))**2)+(((self.atoms[i].getz())-(self.atoms[j].getz()))**2))**(0.50)
                                        if self.atoms[i].getSymbol() == "C":
                                            ri = C
                                        elif self.atoms[i].getSymbol() == "O":
                                            ri = arad.vanderwaals[self.molecule.atoms[atm].getSymbol()]
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
                                            #print("found a negative h value")
                                            continue
                                        if R <= (d+r):
                                            vi = ((np.pi)*((R +r -d)**2)*((d**2)+(2*d*r)-(3*(r**2))+(2*d*R)+(6*r*R)-(3*(R**2))))/(12*(d))
                                        else:
                                            vi = (4/3)*np.pi*(r**3)
                                        Vi = Vi + vi
                        dicz[str(np.round(C, 2))+ "_" + str(np.round(H, 2))] = Vi #{"v" : Vi, "h" : h}





        
        for dum in self.dummies:
            x1 = float(self.dummies[dum]["x"])
            y1 = float(self.dummies[dum]["y"])
            z1 = float(self.dummies[dum]["z"])
            for R in self.rm:
                V = 0.0
                Vi = 0.0
                for atm in self.molecule.atoms:
                    if (str(self.molecule.atoms[atm].getIndex()) not in self.excludedAtoms):#[:,4]
                        #ri = arad.crystalRadii[self.molecule.atoms[atm].getSymbol()]
                        ri = arad.vanderwaals[self.molecule.atoms[atm].getSymbol()]
                        x2 = float(self.molecule.atoms[atm].getx())
                        y2 = float(self.molecule.atoms[atm].gety())
                        z2 = float(self.molecule.atoms[atm].getz())
                        d = ((((x2)-(x1))**2)+(((y2)-(y1))**2)+(((z2)-(z1))**2))**(0.50)
                        h = (ri + R - d)/2
                        if h < 0.0:
                            continue
                        if R <= (d+ri):
                            vi = ((np.pi)*((R +ri -d)**2)*((d**2)+(2*d*ri)-(3*(ri**2))+(2*d*R)+(6*ri*R)-(3*(R**2))))/(12*(d))
                        else:
                            vi = (4/3)*np.pi*(ri**3)
                        Vi = Vi + vi
                dicy[R]= Vi

    def getinternalVDWcalc(self, molvoldict): # This is triggered from the 'internal scan' option in the analysis menu
        self.molecule.CHinternalRadiiOpt(molvoldict)
        #self.molecule.printInternalScan(molvoldict)

                  

