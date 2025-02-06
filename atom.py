class Atom (object):
    #self.number = atomNumber ##not necesary now that Atoms are held in dict
    def __init__(self, atm = ["VT", 0.00, 0.00, 0.00, 609618]):
        #self.index = atm[0]
        self.symbol = atm[0]
        self.x = atm[1]
        self.y = atm[2]
        self.z = atm[3]
        self.index = atm[4]
    def getIndex(self):
        return self.index
    def getSymbol(self):
        return self.symbol
    def getx(self):
        return float(self.x)
    def gety(self):
        return float(self.y)
    def getz(self):
        return float(self.z)
