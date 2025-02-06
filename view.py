import tkinter as tk
import controller
import sphere as sph
from tkinter import messagebox
from tkinter.filedialog import askopenfilename
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import pandas as pd
import atomicradii as arad

class MainView(tk.Frame):
    def __init__(self, parent, controller): 
        self.parent = parent
        self.parent.geometry("1200x700")
        self.controller = controller
        super().__init__(self.parent)
        self.data = {}
        self.vdwDict = {}
        self.molvol = {}
        self.menu = tk.Menu(self.parent) 
        self.parent.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.parent.config(menu=self.menu)

        self.filemenu = tk.Menu(self.menu, tearoff=False)
        self.menu.add_cascade(label='File', menu=self.filemenu)
        self.filemenu.add_command(label='New')
        self.filemenu.add_command(label='Open', command = self.onOpenFile)
        self.filemenu.add_separator()
        self.filemenu.add_command(label='Exit', command=lambda:[self.parent.quit(), self.parent.destroy()])

        self.analysismenu = tk.Menu(self.menu, tearoff=False)
        self.vdwOptionsMenu = tk.Menu(self.analysismenu, tearoff = False)
        self.optimizeMenu = tk.Menu(self.analysismenu, tearoff = False)
        self.optimizationMenu = tk.Menu(self.analysismenu, tearoff = False)
        self.vdwOptionsMenu.add_command(label = "spherical model", command = self.calcAllVDW)
        self.vdwOptionsMenu.add_command(label = "spherical model scan", command = self.calcAllVDWscan)
        self.vdwOptionsMenu.add_command(label = "internal scan", command = self.requestInternalScan)

        
        self.menu.add_cascade(label='Analysis', menu=self.analysismenu)
        self.analysismenu.add_cascade(label='vdw models', menu=self.vdwOptionsMenu)
        self.analysismenu.add_cascade(label='optimization models', menu=self.optimizationMenu)
        self.optimizationMenu.add_command(label = "linear regression", command = self.requestlinearRegression)

        self.windowmenu = tk.Menu(self.menu)
        self.menu.add_cascade(label='Window', menu=self.windowmenu)
        self.helpmenu = tk.Menu(self.menu)
        self.menu.add_cascade(label='Help', menu=self.helpmenu)
        self.helpmenu.add_command(label='About', command = self.launchHelpWindow)

      
        self.grid(row = 0, column = 0, sticky = 'NSEW')
        self.columnconfigure(0, weight = 1)
        self.rowconfigure(0, weight = 1)
        
    def on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.parent.quit()
            self.parent.destroy()

    def onHelpClosing(self):
        self.help.t3.destroy()

    def launchHelpWindow(self):
        self.help.t3 = tk.Toplevel(self.parent)
        self.help.t3.protocol("WM_DELETE_WINDOW", self.onHelpClosing)
        self.help.helpText = tk.Text(self.t3, height = 10, width = 30, wrap = 'word')
        self.help.helpText.pack(expand = True, fill = tk.BOTH)
        about = 'This is a small program for calculating total volume of the intersection of'\
        'a set of atom centered spheres generated from a molecular xyz file.  The resulting '\
        'volume can be used to estimate the magnitude of steric interations related to the '\
        'atomic set.  Often this set will represent a moecule and could then be referred to'\
        'a molecular steric volume.  The program is also capable of identifying subsets of '\
        'atoms to measure the overlap of the subset with the remaining atoms of the original'\
        'set.  Radii of the sphere of a particular element can scanned for optimization'
        self.help.helpText.insert(tk.END, about)


    def onOpenFile(self):
        f = askopenfilename(parent=self.parent)
        currentFile = (f.split("/")[-1]).strip(".xyz")
        if currentFile in self.data:
            self.data[currentFile][len(self.data[currentFile].keys())] = {} # explain why not matching
            currentSys = self.data[currentFile][len(self.data[currentFile].keys()) - 1]# these indices
            currentFile = currentFile + "_" + str(len(self.data[currentFile].keys()) - 1)
        else:
            self.data[currentFile] = {}
            self.data[currentFile][0] = {}
            currentSys = self.data[currentFile][0]
            
        ctlr = controller.SystemController()
        currentSys["filename"] = currentFile # just here for a test
        currentSys["view"] = SystemView(self, ctlr)
        currentSys["view"].controller.openFile(f)
        currentSys["view"].make3Dfigure()
        currentSys["view"].make1DFigure()
        currentSys["view"].updateAtomiclistbox()
        self.windowmenu.add_command(label=currentFile, command = lambda: self.changeFile(currentSys))
        
        currentSys["vdwData"] = {}


    def requestlinearRegression(self):
        pass

    def calcAllVDW(self):  # CALLED WHEN USER CHOOSES SPHERICAL MODEL--runs model.vanderwaalsDifSpheres(self, dicy)
        for filename in self.data:
            print("filename: {} data type: {} ".format(filename, self.data[filename]))
            for syst in self.data[filename].keys():
                self.vdwDict[filename + "_" + str(syst)] = {}
                print("syst: {} data type: {} ".format(syst, self.data[filename][syst]))
                # the app vdwDict gets passed to requestCalcMultVDW method of every SystemView object open
                self.data[filename][syst]["view"].requestCalcMultVDW(self.vdwDict[filename + "_" + str(syst)])
                # should build a self.vdwDict
        allVDW = pd.DataFrame.from_dict(self.vdwDict, orient = 'index')
        for key in self.vdwDict.keys():
            print("{}".format(key))
            for ky in self.vdwDict[key].keys():
                print("\t {}".format(ky))

        allVDW.to_csv("/home/robert/Documents/python/dev/HiStericAll/stericEstimates.csv")

    def calcAllVDWscan(self):  # CALLED WHEN USER CHOOSES SPHERICAL MODEL SCAN--runs model.vanderwaalsDifScan(self, dicy)
        for filename in self.data:
            print("filename: {} data type: {} ".format(filename, self.data[filename]))
            for syst in self.data[filename].keys():
                self.vdwDict[filename + "_" + str(syst)] = {}
                print("syst: {} data type: {} ".format(syst, self.data[filename][syst]))
                # the app vdwDict gets passed to requestCalcMultVDW method of every SystemView object open
                self.data[filename][syst]["view"].requestCalcMultVDWscan(self.vdwDict[filename + "_" + str(syst)])
                # should build a self.vdwDict
        allVDW = pd.DataFrame.from_dict(self.vdwDict, orient = 'index')
        for key in self.vdwDict.keys():
            print("{}".format(key))
            for ky in self.vdwDict[key].keys():
                print("\t {}".format(ky))

        allVDW.to_csv("/home/robert/Documents/python/dev/HiStericAll/stericScanEstimates.csv")
        
    def requestInternalScan(self):
        for filename in self.data:
            #self.vdwDict[filename] = {}
            print("filename: {} data type: {} ".format(filename, self.data[filename]))
            for syst in self.data[filename].keys():
                self.molvol[filename + "_" + str(syst)] = {}
                self.data[filename][syst]["view"].requestinternalVDW(self.molvol[filename + "_" + str(syst)])
        self.allVDW = pd.DataFrame.from_dict(self.molvol, orient = 'index')
        #allVDW = pd.DataFrame.from_dict(self.molvol)
        for key in self.molvol.keys():
            print("{}".format(key))
            #for ky in self.molvol[key].keys():
             #   print("\t {}".format(ky))
        self.allVDW.to_csv("/home/robert/Documents/python/dev/HiStericAll/vdwScan.csv")

    def FitInternalScan(self):
        pass
                
    def changeFile(self, cs):
        #print(cs["filename"])
        cs["view"].mainFrame.tkraise()


class SystemView(): #tk.Frame

    def __init__(self, parent, controller): 
        self.parent = parent
        self.controller = controller
        self.frames = []
        self.mainFrame = tk.Frame(self.parent)            
        for i in range(10):
            self.frames.append(tk.Frame(self.mainFrame, highlightbackground="black", highlightthickness=1))
        self.mainFrame.grid(row = 0, column = 0, sticky = 'NSEW')
        for index, frame in enumerate(self.frames, 0):
            frame.grid(row = index//2, column = index % 2, sticky = "NSEW")

        self.centroidResults = {}
        self.dummycount = 0
        self.toggles = {}
        self.toggles["xyVar"] = tk.IntVar()
        self.toggles["xzVar"] = tk.IntVar()
        self.toggles["yzVar"] = tk.IntVar()
        self.toggles["toggleFrame"] = tk.Frame(self.frames[0])        
        self.toggles["xy"] = tk.Checkbutton(self.toggles["toggleFrame"], text="XY projection Active", variable=self.toggles["xyVar"], command = self.xyActivate)
        self.toggles["xz"] = tk.Checkbutton(self.toggles["toggleFrame"], text="XZ projection Active", variable=self.toggles["xzVar"], command = self.xzActivate)
        self.toggles["yz"] = tk.Checkbutton(self.toggles["toggleFrame"], text="YZ projection Active", variable=self.toggles["yzVar"], command = self.yzActivate)        
        self.atomSelector = {}        
        self.atomSelector["importButtons"] = tk.Frame(self.frames[4])
        self.atomSelector["oneDrotFrame"] = tk.Frame(self.frames[4])        
        self.atomSelector["rotLabelFrame"] = tk.Frame(self.atomSelector["oneDrotFrame"])
        self.atomSelector["rotFrame"] = tk.Frame(self.atomSelector["oneDrotFrame"])
        self.atomSelector["2DfigureToolbarFrame"] = tk.Frame(self.frames[4])# not used. toolbar caused flickering
        self.atomSelector["symmetry"] = {}
        self.atomSelector["symmetryFrame"] = tk.Frame(self.frames[4])
        self.atomSelector["symmetry"]["centerAtOrigin"] = tk.Button(self.atomSelector["symmetryFrame"], text = "center at origin", command = self.requestCenterAtOrigin)
        self.atomSelector["symmetry"]["allignNearestAtomZaxis"] = tk.Button(self.atomSelector["symmetryFrame"], text = "align nearest atom with z", command = self.requestAlignWithZ)
        self.atomSelector["rotXvar"] = tk.IntVar()
        self.atomSelector["rotYvar"] = tk.IntVar()
        self.atomSelector["rotZvar"] = tk.IntVar()
        self.atomSelector["rotXlabel"] = tk.Label(self.atomSelector["rotLabelFrame"], text = "    X", width = 5)
        self.atomSelector["rotYlabel"] = tk.Label(self.atomSelector["rotLabelFrame"], text = "      Y", width = 5)
        self.atomSelector["rotZlabel"] = tk.Label(self.atomSelector["rotLabelFrame"], text = "        Z", width = 5)
        self.atomSelector["rotX"] = tk.Scale(self.atomSelector["rotFrame"], variable = self.atomSelector["rotXvar"], from_ = 0, to = 360, orient = tk.VERTICAL)
        self.atomSelector["rotY"] = tk.Scale(self.atomSelector["rotFrame"], variable = self.atomSelector["rotYvar"], from_ = 0, to = 360, orient = tk.VERTICAL)
        self.atomSelector["rotZ"] = tk.Scale(self.atomSelector["rotFrame"], variable = self.atomSelector["rotZvar"], from_ = 0, to = 360, orient = tk.VERTICAL)
        self.atomSelector["rotButton"] = tk.Button(self.atomSelector["oneDrotFrame"], text = "rotate projections", command = self.requestRotation)
        self.atomSelector["statusbar"] = tk.Frame(self.atomSelector["symmetryFrame"])
        self.atomSelector["atomSelectorbuttons"] = tk.Frame(self.atomSelector["importButtons"])
        self.atomSelector["tkVars"] = {}
        self.atomSelector["tkVars"]["xvar"] = tk.StringVar(self.atomSelector["statusbar"], "0.0000") #self.frames[9]
        self.atomSelector["tkVars"]["yvar"] = tk.StringVar(self.atomSelector["statusbar"], "0.0000") #self.frames[9]
        self.atomSelector["tkVars"]["zvar"] = tk.StringVar(self.atomSelector["statusbar"], "0.0000") #self.frames[9]
        self.atomSelector["tkVars"]["symbolVar"] = tk.StringVar(self.atomSelector["statusbar"], "VT") #self.frames[9]
        self.atomSelector["atomLabel"] = tk.Label(self.atomSelector["statusbar"], text = "Atom: " + self.atomSelector["tkVars"]["symbolVar"].get(), width = 10) # frames[9]
        self.atomSelector["xlabel"] = tk.Label(self.atomSelector["statusbar"], text = "X: ", width = 3) #frames[9]
        self.atomSelector["xbox"] = tk.Entry(self.atomSelector["statusbar"], width = 12, textvariable = self.atomSelector["tkVars"]["xvar"])#frames[9]
        self.atomSelector["ylabel"] = tk.Label(self.atomSelector["statusbar"], text = "Y: ", width = 3) #frames[9]
        self.atomSelector["ybox"] = tk.Entry(self.atomSelector["statusbar"], width = 12, textvariable = self.atomSelector["tkVars"]["yvar"])#frames[9]
        self.atomSelector["zlabel"] = tk.Label(self.atomSelector["statusbar"], text = "Z: ", width = 3) #frames[9]
        self.atomSelector["zbox"] = tk.Entry(self.atomSelector["statusbar"], width = 12, textvariable = self.atomSelector["tkVars"]["zvar"]) #frames[9]
        self.atomSelector["nearestAtomButton"] = tk.Button(self.atomSelector["atomSelectorbuttons"], text = "nearest atom", command = self.requestClosestAtomInfo) #frames[9]
        self.atomSelector["AddToCentroidCalcButton"] = tk.Button(self.atomSelector["atomSelectorbuttons"], text = "import Nearest Atom to Centroid Calcuator", command = self.appendPointToCentroidCalculator)
        self.atomSelector["AddToRotorCalcButton"] = tk.Button(self.atomSelector["importButtons"], text = "import Nearest Atom to Internal Rotor Calcuator", command = self.appendPointToRotorCalculator)
        self.atomSelector["AddRotationalCenterButton"] = tk.Button(self.atomSelector["importButtons"], text = "import Nearest Atom as rotational center", command = self.appendCenterToRotorCalculator)
        self.centroidCalculator = {}
        self.centroidCalculator["xyzViewerFrame"] = tk.Frame(self.frames[8])
        self.centroidCalculator["xyzLabel"] = tk.Label(self.centroidCalculator["xyzViewerFrame"], text = "Atom Selector")
        self.centroidCalculator["xyzViewer"] = tk.Listbox(self.centroidCalculator["xyzViewerFrame"], height = 10, width = 32, selectmode = "single")
        self.centroidCalculator["xyzInstruction"] = tk.Label(self.centroidCalculator["xyzViewerFrame"], text = "Select an atom from above")
        self.centroidCalculator["addToCalculatorButton"] = tk.Button(self.centroidCalculator["xyzViewerFrame"], text = "add to calculator", command = self.moveAtomToCalculator)
        self.centroidCalculator["xyzViewer"].bind("<<ListboxSelect>>", self.onSelect)
        self.centroidCalculator["xyzViewerscrollbar"] = tk.Scrollbar(self.centroidCalculator["xyzViewerFrame"])
        self.centroidCalculator["xyzViewer"].config(yscrollcommand = self.centroidCalculator["xyzViewerscrollbar"].set)
        self.centroidCalculator["xyzViewerscrollbar"].config(command = self.centroidCalculator["xyzViewer"].yview)
        self.centroidCalculator["centroidCalcResultVar"] = tk.StringVar(self.frames[8], "select desired points and press calculate")
        self.centroidCalculator["centroidListLabel"] = tk.Label(self.frames[8], text = "Centroid Calculator")#frame6
        self.centroidCalculator["centroidListBox"] = tk.Listbox(self.frames[8], height = 10, width = 32)#frame6
        self.centroidCalculator["lastCalcLabel"] = tk.Label(self.frames[8], text = self.centroidCalculator["centroidCalcResultVar"].get())
        self.centroidCalculator["calculateCentroidButton"] = tk.Button(self.frames[8], text = "calculate centroid", command = self.requestCentroidCalculation)#command = self.calculateCentroid)#frame6
        self.centroidCalculator["removePointFromCentroidCalcButton"] = tk.Button(self.frames[8], text = "remove selected point from calculator")#frame6
        self.rotors = {}
        self.rotors["rotorsFrame"] = tk.Frame(self.frames[8])
        self.rotors["rotorWindowLabel"] = tk.Label(self.rotors["rotorsFrame"], text = "Selected Internal Rotors")#frame8
        self.rotors["rotorListBox"] = tk.Listbox(self.rotors["rotorsFrame"], height = 10, width = 65, selectmode='extended')#frame8
        self.rotors["calculateRotorParameters"] = tk.Button(self.rotors["rotorsFrame"], text = "add rotor to model", command = self.requestAppendPointToModel)##frame8
        self.rotors["removeRotorButton"] = tk.Button(self.rotors["rotorsFrame"], text = "remove selected rotor", command = self.removeRotor)#frame9
        self.centroids = {}
        self.centroids["CentroidResultsWindowLabel"] = tk.Label(self.frames[9], text = "Calculated Centroids")#frame8
        self.centroids["CentroidResultsWindow"] = tk.Listbox(self.frames[9], height = 10, width = 65, selectmode='extended')#frame8
        self.centroids["addCentroidsToModelButton"] = tk.Button(self.frames[9], text = "add centroids to model", command = self.requestAppendPointToModel)##frame8
        self.centroids["removeCentroidButton"] = tk.Button(self.frames[9], text = "remove selected centroid", command = self.removeCentroid)#frame9
        self.centroids["centroidCounter"] = 1       
        self.toggles["toggleFrame"].grid(row=0, column = 0)#, sticky = 'NSEW'
        self.toggles["xy"].grid(row=0, column=0)
        self.toggles["xz"].grid(row=0, column=1)
        self.toggles["yz"].grid(row=0, column=2)
        self.atomSelector["statusbar"].grid(row = 0,column = 0, columnspan = 1)
        self.atomSelector["symmetry"]["centerAtOrigin"].grid(row = 1, column = 0)
        self.atomSelector["symmetry"]["allignNearestAtomZaxis"].grid(row = 2, column = 0)                
        self.atomSelector["atomLabel"].grid(row = 0, column = 0, sticky = 'w')
        self.atomSelector["xlabel"].grid(row = 0, column = 1, sticky = 'w')#, padx = 10) #sticky = 'e'
        self.atomSelector["xbox"].grid(row = 0, column = 2)#, sticky = 'w')#, padx = 10)
        self.atomSelector["ylabel"].grid(row = 0, column = 3)#, sticky = 'w')#, padx = 10)
        self.atomSelector["ybox"].grid(row = 0, column = 4)#, sticky = 'w')#, padx = 10)
        self.atomSelector["zlabel"].grid(row = 0, column = 5)#, sticky = 'w')#, padx = 10)
        self.atomSelector["zbox"].grid(row = 0, column = 6)#, sticky = 'w')#, padx = 10)
        self.atomSelector["atomSelectorbuttons"].grid(row = 2, column = 0, columnspan = 1)
        self.atomSelector["symmetryFrame"].grid(row = 0, column = 0)
        self.atomSelector["oneDrotFrame"].grid(row = 0, column = 2) # frame that holds rotate projections button,
        self.atomSelector["importButtons"].grid(row = 0, column = 1)
        self.atomSelector["rotLabelFrame"].grid(row = 0, column = 0, columnspan = 1)  # 1dRotFrame
        self.atomSelector["rotXlabel"].grid(row = 0, column = 0, sticky = 'NSEW') 
        self.atomSelector["rotYlabel"].grid(row = 0, column = 1, sticky = 'NSEW')
        self.atomSelector["rotZlabel"].grid(row = 0, column = 2, sticky = 'NSEW')
        self.atomSelector["rotFrame"].grid(row = 1, column = 0, rowspan = 1)   #  # 1dRotFrame
        self.atomSelector["rotButton"].grid(row = 2, column = 0, rowspan = 1) # # 1dRotFrame
        self.atomSelector["rotX"].grid(row = 0, column = 0)
        self.atomSelector["rotY"].grid(row = 0, column = 1)
        self.atomSelector["rotZ"].grid(row = 0, column = 2)        
        self.atomSelector["nearestAtomButton"].grid(row = 0, column = 3) #, columnspan = 3
        self.atomSelector["AddToCentroidCalcButton"].grid(row = 0, column = 4)        
        self.atomSelector["AddRotationalCenterButton"].grid(row = 0, column = 0)
        self.atomSelector["AddToRotorCalcButton"].grid(row = 1, column = 0)
        self.centroidCalculator["xyzViewerFrame"].grid(row = 0, column = 1, rowspan = 4)
        self.centroidCalculator["xyzLabel"].grid(row = 0, column = 0)
        self.centroidCalculator["xyzViewer"].grid(row = 1, column = 0)
        self.centroidCalculator["xyzViewerscrollbar"].grid(row = 1, column = 1, columnspan = 2, sticky = 'w')
        self.centroidCalculator["xyzInstruction"].grid(row = 2, column = 0)
        self.centroidCalculator["addToCalculatorButton"].grid(row = 3, column = 0)
        self.centroidCalculator["centroidListLabel"].grid(row = 0, column = 2)
        self.centroidCalculator["centroidListBox"].grid(row = 1, column = 2)
        self.centroidCalculator["lastCalcLabel"].grid(row = 2, column = 2)
        self.centroidCalculator["calculateCentroidButton"].grid(row = 3, column = 2)
        self.centroidCalculator["removePointFromCentroidCalcButton"].grid(row = 4, column = 2)
        
        self.centroids["CentroidResultsWindowLabel"].grid()
        self.centroids["CentroidResultsWindow"].grid()
        self.centroids["addCentroidsToModelButton"].grid()
        self.centroids["removeCentroidButton"].grid()

        self.rotors["rotorsFrame"].grid(row = 0, column = 0, rowspan = 4)
        self.rotors["rotorWindowLabel"].grid()
        self.rotors["rotorListBox"].grid()
        self.rotors["calculateRotorParameters"].grid()
        self.rotors["removeRotorButton"].grid()
        
        for col in range(2):
            self.mainFrame.columnconfigure(col, weight = 1)
        for row in [1]: #range(1, 5):
            self.mainFrame.rowconfigure(row, weight = 1)

        expandableFrames = range(2,10)
        for frame in self.frames:    #expandableFrames:
            frame.columnconfigure(0, weight = 1)
            frame.rowconfigure(0, weight = 1)

    def make3Dfigure(self):
        self.figure3D = plt.figure(figsize = [3.2, 2.4])#default is [6.4, 4.8] figsize = [5.4, 2.2]
        self.rotLock = tk.Checkbutton(self.frames[1], text="show cartesian axes", variable=self.toggles["xyVar"], command = self.xyActivate)
        self.canvas3D = FigureCanvasTkAgg(self.figure3D, master = self.frames[3]) # master = self.rt ##frame1
        self.axs3D = self.figure3D.add_subplot(projection = '3d') # removed gspec1[0,:],
        coords = self.requestCoords()
        bonds = self.requestBonds()
        self.figure3D.subplots_adjust(left=0, right=1, bottom=0, top=1)
       
        self.axs3D.scatter(coords[0],coords[1],coords[2], marker = "o", color = coords[3], edgecolor="darkgrey")
        for bond in bonds:
            self.axs3D.plot(bond[0],bond[1],bond[2], color="b")
        #self.set_size(1, 1, self.axs3D) not yet working
        self.updatePlotLimits()
        self.axs3D.set_xlabel("x")
        self.axs3D.set_ylabel("y")
        self.axs3D.set_zlabel("z")
        self.rotLock.grid(row = 0, column = 0)#, sticky = 'NSEW'
        self.canvas3D.get_tk_widget().grid(row = 0, column = 0, sticky = 'NSEW')

    def make1DFigure(self):
        self.figure1D = plt.figure(figsize = [5.4, 2.4])
        self.canvas1D = FigureCanvasTkAgg(self.figure1D, master = self.frames[2]) #frame4
        gspec = self.figure1D.add_gridspec(nrows = 1, ncols = 3, wspace = 0.5)
        self.axs1D = {}
        bonds = self.requestBonds()
        
        xy = self.requestxy()
        clr = self.requestColors()
        self.axs1D["xy"] = self.figure1D.add_subplot(gspec[0,0])
        self.axs1D["xy"].scatter(xy[0], xy[1], color = clr, edgecolor="darkgrey")
        for bond in bonds:
            self.axs1D["xy"].plot(bond[0],bond[1], color="b")
        xz = self.requestxz()
        self.axs1D["xz"] = self.figure1D.add_subplot(gspec[0,1])
        self.axs1D["xz"].scatter(xz[0], xz[2], color = clr, edgecolor="darkgrey")
        for bond in bonds:
            self.axs1D["xz"].plot(bond[0],bond[2], color="b")
        yz = self.requestyz()
        self.axs1D["yz"] = self.figure1D.add_subplot(gspec[0,2])
        self.axs1D["yz"].scatter(yz[1], yz[2], color = clr, edgecolor="darkgrey")
        for bond in bonds:
            self.axs1D["yz"].plot(bond[1],bond[2], color="b")
        self.update2dPlotLimits()
        self.canvas1D.draw()
        self.canvas1D.get_tk_widget().grid(row = 0, column = 0, sticky = 'NSEW')


    def update1dPlots(self, dCoords):
        self.figure1D.clear()
        bonds = self.requestBonds()
        gspec = self.figure1D.add_gridspec(nrows = 1, ncols = 3, wspace = 0.5)
        xy = dCoords[1] #self.requestxy()
        clr = dCoords[4]        
        self.axs1D["xy"] = self.figure1D.add_subplot(gspec[0,0])
        self.axs1D["xy"].scatter(xy[0], xy[1], color = clr, edgecolor="darkgrey")
        for bond in bonds:
            self.axs1D["xy"].plot(bond[0],bond[1], color="b")
        xz = dCoords[2]  #self.requestxz()
        self.axs1D["xz"] = self.figure1D.add_subplot(gspec[0,1])
        self.axs1D["xz"].scatter(xz[0], xz[2], color = clr, edgecolor="darkgrey")
        for bond in bonds:
            self.axs1D["xz"].plot(bond[0],bond[2], color="b")
        yz = dCoords[3]      #self.requestyz()
        self.axs1D["yz"] = self.figure1D.add_subplot(gspec[0,2])
        self.axs1D["yz"].scatter(yz[1], yz[2], color = clr, edgecolor="darkgrey")
        for bond in bonds:
            self.axs1D["yz"].plot(bond[1],bond[2], color="b")
        self.canvas1D.draw()

    def update3dPlots(self, dCoords):
        self.figure3D.clear()
        self.axs3D = self.figure3D.add_subplot(projection = '3d')
        atoms = dCoords[0]
        bonds = self.requestBonds()
        self.axs3D.scatter(atoms[0],atoms[1],atoms[2], marker = "o", color = dCoords[4], edgecolor="darkgrey")
        for bond in bonds:
            self.axs3D.plot(bond[0],bond[1],bond[2], color="b")
        self.updatePlotLimits()
        self.axs3D.set_xlabel("x")
        self.axs3D.set_ylabel("y")
        self.axs3D.set_zlabel("z")
        self.canvas3D.draw()

    def updatePlotLimits(self):
        limits = self.request3Dminmax()
        low = min(limits["xmin"], limits["ymin"], limits["zmin"])
        high = max(limits["xmax"], limits["ymax"], limits["zmax"])
        self.axs3D.set_xlim(xmin = low, xmax = high)
        self.axs3D.set_ylim(ymin = low, ymax = high)
        self.axs3D.set_zlim(zmin = low, zmax = high)

    def update2dPlotLimits(self):
        limits = self.request3Dminmax()
        self.axs1D["xy"].set_xlim(xmin = min(limits["xmin"], limits["ymin"]), xmax = max(limits["xmax"], limits["ymax"]))
        self.axs1D["xy"].set_ylim(ymin = min(limits["xmin"], limits["ymin"]), ymax = max(limits["xmax"], limits["ymax"]))
        self.axs1D["xz"].set_xlim(xmin = min(limits["xmin"], limits["zmin"]), xmax = max(limits["xmax"], limits["zmax"]))
        self.axs1D["xz"].set_ylim(ymin = min(limits["xmin"], limits["zmin"]), ymax = max(limits["xmax"], limits["zmax"]))
        self.axs1D["yz"].set_xlim(xmin = min(limits["ymin"], limits["zmin"]), xmax = max(limits["ymax"], limits["zmax"]))
        self.axs1D["yz"].set_ylim(ymin = min(limits["ymin"], limits["zmin"]), ymax = max(limits["ymax"], limits["zmax"]))

    def updateAtomiclistbox(self):
        self.centroidCalculator["xyzViewer"].delete(0,tk.END)
        for indx, (atm, x, y, z) in  enumerate(self.requestDF().values.tolist()):
            self.centroidCalculator["xyzViewer"].insert(tk.END, str(indx) + "    " + atm + "    " + str(round(x, 4)) + "    " + str(round(y, 4)) + "    " + str(round(z, 4)))

    def updateRotorListBox(self):
        pass

    def updateCentroidCalcListBox(self):
        pass

    def updateCentroidListBox(self):
        pass
        
    def requestDF(self):
        return self.controller.getDF()

    def requestCenterAtOrigin(self):
        centeredCoords = self.controller.centerAtOrigin()
        self.update1dPlots(centeredCoords)
        self.update3dPlots(centeredCoords)
        self.canvas1D.draw()
        self.updateAtomiclistbox()

    def requestAlignWithZ(self):
        closestAtom = self.controller.getClosestAtomCoords(self.readFromXYZboxes())
        alignedCoords = self.controller.alignWz(closestAtom)
        self.update1dPlots(alignedCoords)
        self.update3dPlots(alignedCoords)
        self.canvas1D.draw()
        self.updateAtomiclistbox()
    
    def requestRotation(self):  # getRotation returns the array below with colors appended to end
        rotCoords = self.controller.getRotation([self.atomSelector["rotXvar"].get(), self.atomSelector["rotYvar"].get(), self.atomSelector["rotZvar"].get()])
        self.update1dPlots(rotCoords)
        self.update3dPlots(rotCoords)
        self.canvas1D.draw()
        self.updateAtomiclistbox()

    def requestSymmetryAdjustment(self):
        pass

    def requestClosestAtomInfo(self):
        self.closestAtom = self.controller.getClosestAtomCoords(self.readFromXYZboxes())
        self.atomSelector["tkVars"]["xvar"].set(self.closestAtom[2])
        self.atomSelector["tkVars"]["yvar"].set(self.closestAtom[3])
        self.atomSelector["tkVars"]["zvar"].set(self.closestAtom[4])
        self.atomSelector["tkVars"]["symbolVar"].set(self.closestAtom[1])
        self.atomSelector["atomLabel"].configure(text = "Atom: " + self.atomSelector["tkVars"]["symbolVar"].get())#self.atomSelector["atomLabel"].insert(0, closestAtom[3]))

        
    def request3Dminmax(self):
        return self.controller.get3Dminmax()

    def requestTestCoords(self):
        return self.controller.getTestCoords()

    def requestCoords(self):
        return self.controller.getCoords()
    
    def requestBonds(self):
        return self.controller.getBonds()

    def requestColors(self):
        return self.controller.getColors()
    
    def requestxy(self):
        return self.controller.getxy()
     
    def requestxz(self):
        return self.controller.getxz()

    def requestyz(self):
        return self.controller.getyz()

    def onclickXY(self, event):
        ix = event.xdata
        iy = event.ydata
        if (self.atomSelector["xbox"].get()):
            self.atomSelector["xbox"].delete(0, 'end')
        self.atomSelector["xbox"].insert(0,ix)
        if (self.atomSelector["ybox"].get()):
            self.atomSelector["ybox"].delete(0, 'end')
        self.atomSelector["ybox"].insert(0,iy)

    def xyActivate(self):
        if self.toggles["xyVar"].get() == 1:
            self.cidXY = self.canvas1D.mpl_connect('button_press_event', self.onclickXY)
        else:
            self.canvas1D.mpl_disconnect(self.cidXY)

    def xzActivate(self):
        if self.toggles["xzVar"].get() == 1:
            self.cidXZ = self.canvas1D.mpl_connect('button_press_event', self.onclickXZ)
        else:
            self.canvas1D.mpl_disconnect(self.cidXZ)

    def onclickXZ(self, event):
        ix = event.xdata
        iz = event.ydata
        if (self.atomSelector["xbox"].get()):
            self.atomSelector["xbox"].delete(0, 'end')
        self.atomSelector["xbox"].insert(0,ix)
        if (self.atomSelector["zbox"].get()):
            self.atomSelector["zbox"].delete(0, 'end')
        self.atomSelector["zbox"].insert(0,iz)

    def yzActivate(self):
        if self.toggles["yzVar"].get() == 1:
            self.cidYZ = self.canvas1D.mpl_connect('button_press_event', self.onclickYZ)
        else:
            self.canvas1D.mpl_disconnect(self.cidYZ)

    def onclickYZ(self, event):
        iy = event.xdata
        iz = event.ydata
        if (self.atomSelector["ybox"].get()):
            self.atomSelector["ybox"].delete(0, 'end')
        self.atomSelector["ybox"].insert(0,iy)
        if (self.atomSelector["zbox"].get()):
            self.atomSelector["zbox"].delete(0, 'end')
        self.atomSelector["zbox"].insert(0,iz)

    def onSelect(self, event):
        selection = event.widget.curselection()
        if selection:            
            index = selection[0]
            self.centroidCalculator["tmpTxt"] = event.widget.get(index).split()
            try:
                self.tmpSelection.remove()
            except AttributeError:
                pass
            self.tmpSelection = self.axs3D.scatter(float(self.centroidCalculator["tmpTxt"][2]), float(self.centroidCalculator["tmpTxt"][3]), float(self.centroidCalculator["tmpTxt"][4]), marker = "o", s = 42, color = "chartreuse")
            self.canvas3D.draw()
        else:
            self.tmpSelection.remove()
            self.canvas3D.draw()
            

    def readFromXYZboxes(self):
        return [float(self.atomSelector["xbox"].get()), float(self.atomSelector["ybox"].get()), float(self.atomSelector["zbox"].get())]
  
    def moveAtomToCalculator(self):
        try:
            print("moving to calculator: {}".format(self.centroidCalculator["tmpTxt"]))
            self.controller.refreshCentroidCalculator(self.centroidCalculator["tmpTxt"])
            self.centroidCalculator["centroidListBox"].insert(tk.END, self.centroidCalculator["tmpTxt"][1] + "    " + self.centroidCalculator["tmpTxt"][2] + "    " + self.centroidCalculator["tmpTxt"][3] + "    " + self.centroidCalculator["tmpTxt"][4])
        except AttributeError:
            tk.messagebox.showerror('Python Error', 'Error: Select an atom before requsting transfer!')

    def appendPointToCentroidCalculator(self):
        xyz = self.controller.getClosestAtomCoords(self.readFromXYZboxes())
        label = str(xyz[1]) + "    " + str(xyz[2]) + "    " + str(xyz[3]) + "    " + str(xyz[4])
        self.centroidCalculator["centroidListBox"].insert(tk.END, label)
        newxyz = self.controller.updateCentroidCalculator(xyz)
        self.atomSelector["tkVars"]["xvar"].set(newxyz[0])
        self.atomSelector["tkVars"]["yvar"].set(newxyz[1])
        self.atomSelector["tkVars"]["zvar"].set(newxyz[2])
        self.atomSelector["tkVars"]["symbolVar"].set("Atom: " + newxyz[3])
        self.atomSelector["atomLabel"].configure(text = self.atomSelector["tkVars"]["symbolVar"].get())        

    def appendPointToRotorCalculator(self):
        xyz = self.controller.getClosestAtomCoords(self.readFromXYZboxes())
        print("appending to calculator: {}".format(xyz))
        label = str(xyz[1]) + "    " + str(xyz[2]) + "    " + str(xyz[3]) + "    " + str(xyz[4])
        self.rotors["rotorListBox"].insert(tk.END, label)
        newxyz = self.controller.updateRotorCalculator(xyz)
        self.atomSelector["tkVars"]["xvar"].set(newxyz[0])
        self.atomSelector["tkVars"]["yvar"].set(newxyz[1])
        self.atomSelector["tkVars"]["zvar"].set(newxyz[2])
        self.atomSelector["tkVars"]["symbolVar"].set("Atom: " + newxyz[3])
        self.atomSelector["atomLabel"].configure(text = self.atomSelector["tkVars"]["symbolVar"].get())

    def appendCenterToRotorCalculator(self):#this one needs to be the axix...dot product of H's...through the center
        xyz = self.controller.getClosestAtomCoords(self.readFromXYZboxes()) # xyz has: [index, symbol, x, y, z]
        label = "center" + "    " + str(xyz[1]) + "    " + str(xyz[2]) + "    " + str(xyz[3]) + "    " + str(xyz[4])
        self.rotors["rotorListBox"].insert(tk.END, label)
        newxyz = self.controller.updateRotorCenterCalculator(xyz.append("center")) #whythe append???
        self.atomSelector["tkVars"]["xvar"].set(newxyz[0])
        self.atomSelector["tkVars"]["yvar"].set(newxyz[1])
        self.atomSelector["tkVars"]["zvar"].set(newxyz[2])
        self.atomSelector["tkVars"]["symbolVar"].set("Atom: " + newxyz[3])
        self.atomSelector["atomLabel"].configure(text = self.atomSelector["tkVars"]["symbolVar"].get())
        
    def removeCentroid(self):
        pass

    def removeRotor(self):
        pass

    def requestCentroidCalculation(self):
        x = self.controller.getCentroidCalc()
        self.centroidResults[x[0]] = [x[1], x[2], x[3]]
        self.centroidCalculator["centroidCalcResultVar"].set("last centroid calculation result: x = " + x[1][:8] + "  y = " + x[2][:8] + "  z = " + x[3][:8])
        self.centroidCalculator["lastCalcLabel"].configure(text = self.centroidCalculator["centroidCalcResultVar"].get())
        self.centroids["CentroidResultsWindow"].insert(tk.END, x[0] + "  x = " + x[1] + "  y = " + x[2] + "  z = " + x[3])
        self.centroidCalculator["centroidListBox"].delete(0, tk.END)

    def requestAppendPointToModel(self):
        keys = []
        if (len(self.centroids["CentroidResultsWindow"].curselection()) > 0):
            for centroid in self.centroids["CentroidResultsWindow"].curselection():
                self.dummycount += 1
                xyzs = self.centroids["CentroidResultsWindow"].get(centroid).split()
                keys.append(xyzs[0].split("_", 1)[1])
                x = xyzs[3]
                y = xyzs[6]
                z = xyzs[9]
                self.axs3D.scatter(float(x), float(y), float(z), marker = 'o', color = "chartreuse")
                self.genSphere([x, y, z])
            self.canvas3D.draw()
            self.controller.requestUpdateDummyAtoms(keys)
            self.centroids["CentroidResultsWindow"].delete(self.centroids["CentroidResultsWindow"].curselection())

    def requestAppendRotor(self):
        keys = []
        if (len(self.rotors["rotorWindowLabel"].curselection()) > 0):
            for centroid in self.rotors["rotorWindowLabel"].curselection():
                self.dummycount += 1
                xyzs = self.rotors["rotorWindowLabel"].get(centroid).split()
                keys.append(xyzs[0].split("_", 1)[1])
                x = xyzs[3]
                y = xyzs[6]
                z = xyzs[9]
                self.axs3D.scatter(float(x), float(y), float(z), marker = 'o', color = "chartreuse")
                self.genSphere([x, y, z])
            self.canvas3D.draw()
            self.controller.requestUpdateDummyAtoms(keys)
            self.centroids["CentroidResultsWindow"].delete(self.centroids["CentroidResultsWindow"].curselection())


    def set_size(self, w, h, ax=None):
        """ w, h: width, height in inches """
        if not ax: ax=plt.gca()
        l = ax.figure.subplotpars.left
        r = ax.figure.subplotpars.right
        t = ax.figure.subplotpars.top
        b = ax.figure.subplotpars.bottom
        figw = float(w)/(r-l)
        figh = float(h)/(t-b)
        ax.figure.set_size_inches(figw, figh)

    def df_to_list(self,df):
        self.tree["columns"] = df.columns.values.tolist()
        for x in range(len(df.columns.values)):
            self.tree.column(df.columns.values[x], width=100)
            self.tree.heading(df.columns.values[x], text=df.columns.values[x], command=self.populate_selection)

        for index, row in df.iterrows():
            self.tree.insert("",0,text=index,values=list(row))

        self.tree.grid(row=50,column=0,rowspan=1,columnspan=12,sticky=N+E+W+S)
        self.tree.bind("<<TreeviewSelect>>", self.populate_selection)


    def genSphere(self, ar):  #NeedsWorK!  put these in a model framework like scikitlearn!!!!
        sp = sph.Sphere({"h" : float(ar[0]), "k" : float(ar[1]), "l" : float(ar[2]), "r" : 1.0})
        self.spcoords = sp.getXYZ()
        self.axs3D.plot_wireframe(self.spcoords["x"], self.spcoords["y"], self.spcoords["z"], linewidth = 0.5)
        self.canvas3D.draw()

    def requestVDWsphereCalc(self):
        if self.dummycount == 0:
            pass
        else:
            self.controller.getVDWsphereCalc()

    def requestCalcMultVDW(self, vdwdict):
        if self.dummycount == 0:
            pass
        else:
            self.controller.getVDWsphereCalcs(vdwdict)

    def requestCalcMultVDWscan(self, vdwdict):
        if self.dummycount == 0:
            pass
        else:
            self.controller.getVDWscanCalcs(vdwdict)
            
    def requestinternalVDW(self, mvdict):
        self.controller.getinternalVDWcalc(mvdict)
        
    def on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.parent.quit()
            self.parent.destroy()


def main(controller):
    root = tk.Tk()
    root.columnconfigure(index=0, weight=1)
    root.rowconfigure(index=0, weight=1)
    root.title("MVC")
    v1 = MainView(root, controller)# removed controller
    root.mainloop()
    
if __name__ == '__main__':
    cMain = controller.MainController()
    main(cMain)#

