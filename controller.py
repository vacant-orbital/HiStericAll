import model


class MainController:
    def __init__(self):
        self.model = model.Main_Model()
    def getlinearRegressionCalc(self):
        self.model.runLinReg()

class SystemController:
    
    def __init__(self):
        ''' No constructer used.  model object is created when this.open file is called'''
        pass

    def updatePlot(self):
        pass
       
    def openFile(self, fl):
        self.model = model.System_Model(fl)

    def get3Dminmax(self):
        return self.model.send3Dminmax()

    def getCoords(self):
        return self.model.sendxyz()

    def getDF(self):
        return self.model.sendDF()

    def getBonds(self):
        return self.model.sendBonds()

    def getxy(self):
        return self.model.sendxy()

    def getxz(self):
        return self.model.sendxz()

    def getyz(self):
        return self.model.sendyz()

    def getColors(self):
        return self.model.sendColors()

    def getRotation(self, angles):
        #print(angles)
        return self.model.rotateMolecule(angles)

    def centerAtOrigin(self):
        return self.model.centerMolecule()

    def alignWz(self, coords):
        return self.model.alignZaxis(coords)

    def getClosestAtomCoords(self, ar):
        return self.model.sendClosestAtomInfo(ar)

    def getClosestBond(self, ar):
        return self.model.sendClosestBondInfo(ar)

    def updateCentroidCalculator(self, ar):
        #self.model.excluded
        return self.model.sendResetTmpAtom(ar)

    def updateRotorCalculator(self, ar):
        #self.model.excluded
        return self.model.sendResetTmpAtomUpdateRotor(ar)
    
    def updateRotorCenterCalculator(self, ar):
        #self.model.excluded
        return self.model.sendResetTmpAtomUpdateRotorCenter(ar)
    
    def refreshCentroidCalculator(self, ar):
        return self.model.updateCentroidCalc(ar)

    def getCentroidCalc(self):
        return self.model.calculateCentroid()

    def requestUpdateDummyAtoms(self, ky):
        self.model.updateModelWithSelectedCentroids(ky)

##    def getVDWsphereCalc(self):
##        self.model.vanderwaalsDifSphere()

    def getVDWsphereCalcs(self, vdict):
        self.model.vanderwaalsDifSpheres(vdict)

    def getVDWscanCalcs(self, vdict):
        self.model.vanderwaalsDifScan(vdict)
        
    def getinternalVDWcalc(self, movodict):
        self.model.getinternalVDWcalc(movodict)

