# 
#  -de-Convolution algorithm for data scanned with
#  miniature confocal microscope setup by Rob Heath
#                    rmh9@hw.ac.uk
#
#  Required inputs: filename, spotsize, rangeX, rangeY
#
#  !!!!                                           !!!!
#  !!!!  getting the rangeX or rangeY wrong will  !!!!
#  !!!!       invalidate the deconvolution        !!!!
#
#  !!!! the scaling for the spot size requires it !!!!
#  
#  Relevant documentation can be found at:
#   www.lmd.ens.fr/jori/pdfs/FFTzeropad.pdf
#   hebb.mit.edu/courses/9.29/2002/readings/c13-1.pdf
#

print 'Importing SciPy, NumPy, wxPython and PIL...',
from scipy import fftpack, array, signal
import csv
import os
from math import *
import numpy as np
try:
    import Image, ImageDraw
except ImportError:
    print 'Python Image Library required!'
    raise ImportError
try:
    import wx
except ImportError:
    print 'wxPython required!'
    raise ImportError
print 'Done!'

progname = "Convolver"
version = "1.0" #increment by one every time you update something
author = "Rob Heath"

##setup
filename = r'input.bmp' #the image to be convolved
filename_export = r'convolved.bmp'
spotsize = 2.0 #FWHM, in microns; assumed equal in X and Y
rangeX = 24.0 #in microns, the image's range in X
rangeY = 24.0 #in microns, the image's range in Y
##end setup

lastConvolved = ''

def importBitmap(filename):
    picture = Image.open(filename)
    pixels = list(picture.getdata())
    width, height = picture.size

    data_array = []
    pixel_maximum = 0
    for p, pixel in enumerate(pixels):
        if p%width==0:
            data_array.append([])
        data_array[-1].append(float(sum(pixel)))
        if sum(pixel) > pixel_maximum:
            pixel_maximum = sum(pixel)

    data_array = np.array(data_array)
    data_array /= float(pixel_maximum)
    return data_array, (width,height)

def exportBitmap(filename, image_data, shape, filetype):
    image_data = list(image_data)
    picture = Image.new('L',shape) #change L to RGB for colour (only pain)
    draw = ImageDraw.Draw(picture)
    maximum = 0
    for row in image_data:
        for cell in row:
            if cell > maximum:
                maximum = cell
    image_data = list(255*(np.array(image_data)/float(maximum)))
    for y, row in enumerate(image_data):
        for x, value in enumerate(row):
            image_data[y][x]=int(round(value))
            draw.point((x,y),fill=int(round(value)))
    
    picture.save(filename,filetype)
    print 'Wrote '+str(filetype)+' data to: '+str(filename)
    

def writeCsv(data, filename_csv):
    numfile = open('out-input.csv','wb')
    numcsv = csv.writer(numfile, delimiter=',')
    for row in data:
        numcsv.writerow(row)
    numfile.close()
    print 'Wrote ' + str(filename_csv)

def gaussian2(x, y, FWHM, xOffset, yOffset):
    x,y,FWHM,xOffset,yOffset = float(x),float(y),float(FWHM),float(xOffset),float(yOffset)
    max = 1.0
    sigma = FWHM / 2.3548
    output = max * exp(-1*((x - xOffset)**2 + ((y - yOffset)**2)) / (2*(sigma**2)))
    return output

def zeropad(fftImage):
    #split fftImage into quarters
    fftImageSplit=[]
    for fftHalf in np.array_split(fftImage, 2, 0):
        for fftQuarter in np.array_split(fftHalf, 2, 1):
            fftImageSplit.append(fftQuarter)
    #create zeroes
    fftZeros=[]
    for fftQuarter in fftImageSplit:
        fftZeros.append(np.zeros(fftQuarter.shape))
    #concatenate zeros and data "zero padding" the matrix
    #for why see:
    # www.lmd.ens.fr/jori/pdfs/FFTzeropad.pdf
    fftConcat=[]
    #spread to edges ##doesn't change anything
    fftConcat.append(np.concatenate((fftImageSplit[0],fftZeros[0],fftZeros[1],fftImageSplit[1]),1))
    fftConcat.append(np.concatenate((fftZeros[0],fftZeros[0],fftZeros[1],fftZeros[1]),1))
    fftConcat.append(np.concatenate((fftZeros[2],fftZeros[2],fftZeros[3],fftZeros[3]),1))
    fftConcat.append(np.concatenate((fftImageSplit[2],fftZeros[2],fftZeros[3],fftImageSplit[3]),1))
    result = np.concatenate((fftConcat[0],fftConcat[1],fftConcat[2],fftConcat[3]))
    return result

def convolve(filename,filename_export,spotsize,rangeX,rangeY,filetype,outputcsv=False,outputinterim=True):

    numericalData, shape = importBitmap(filename)
    print 'Read input.bmp'

    if outputcsv==True:
        writeCsv(numericalData, 'out-input.csv')
    if outputinterim==True:
        exportBitmap('out-input.bmp',numericalData,shape,filetype)

    ##create gaussian spot profile of correct size
    stepsX,stepsY = shape
    centerX = float(stepsX-1)/2.0
    centerY = float(stepsY-1)/2.0
    stepLengthX = float(rangeX)/float(stepsX)
    stepLengthY = float(rangeY)/float(stepsY)

    #spot zeropadding######
    spot=[]
    twoSigma = 2.0*spotsize / 2.3548

    for y in range(stepsY):
        spot.append([])
        for x in range(stepsX):
            r = (((x-centerX)*stepLengthX)**2+((y-centerY)*stepLengthY)**2)**0.5
            spot[y].append(gaussian2(x*stepLengthX,y*stepLengthY,spotsize,centerX*stepLengthX,centerY*stepLengthY))
    spot = array(spot)
    #####################

    print 'Generated Gaussian'

    if outputcsv==True:
        writeCsv(spot,'spot.csv')
    if outputinterim==True:
        exportBitmap('spot.bmp',spot,spot.shape,filetype)
    
    
        
    ##perform deconvolution
    #fft & multiply
    #zeropad inserts padding; not sure why any more :(

    #the magic happens here:
    fftImage = np.multiply(fftpack.fft2(numericalData),fftpack.fft2(spot))
    convolvedImage = fftpack.ifft2(fftImage).real #returns only real components
    #done

    ##reorganise data (weird indexes, look up fftpack module to explain)
    imageSplit=[]
    for half in np.array_split(convolvedImage, 2, 0):
        for quarter in np.array_split(half, 2, 1):
            imageSplit.append(quarter)
    concat=[]
    concat.append(np.concatenate((imageSplit[3],imageSplit[2]),1))
    concat.append(np.concatenate((imageSplit[1],imageSplit[0]),1))
    convolvedImage=np.concatenate((concat[0],concat[1]))

    print 'Convolved spot and input'

    ##output result
    
    if outputcsv==True:
        writeCsv(convolvedImage,'convolved.csv')
    exportBitmap(filename_export, convolvedImage, shape,filetype)

    print 'Done!'



class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self,parent,title=title, size=(600,480))
        #self.control = wx.TextCtrl(self)
        self.panel = wx.Panel(self,-1)
        self.CreateStatusBar()
        
        self.plt = False

        #menubar
        filemenu = wx.Menu()
        menuAbout = filemenu.Append(wx.ID_ABOUT, "&About","Information about "+str(progname))
        menuExit = filemenu.Append(wx.ID_EXIT, "E&xit","Terminate "+str(progname))
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File")
        self.SetMenuBar(menuBar)

        #sizer parts
        self.sizer1 = wx.GridBagSizer(2,5)
        self.sizer2 = wx.GridBagSizer(2,5)
        self.sizer3 = wx.GridBagSizer(2,5)
        self.sizer4 = wx.GridBagSizer(2,5)
        self.sizer5 = wx.GridBagSizer(2,5)

        #loading file
        self.button1 = wx.Button(self.panel, -1, "Load from...")
        self.filepath = wx.TextCtrl(self.panel, size=wx.Size(256,-1))
        
        self.sizer1.Add(self.button1,(0,0), wx.DefaultSpan, wx.ALL, 3)
        self.sizer1.Add(self.filepath, (0,1), (1,1), wx.ALIGN_CENTER)

        #X-Y bits
        self.textX = wx.StaticText(self.panel, label="X:")
        self.editX = wx.TextCtrl(self.panel, size=wx.Size(48,-1))
        self.editX.SetEditable(False)
        self.editX.SetBackgroundColour((212,212,212))
        self.textY = wx.StaticText(self.panel, label="Y:")
        self.editY = wx.TextCtrl(self.panel, size=wx.Size(48,-1))
        self.editY.SetEditable(False)
        self.editY.SetBackgroundColour((212,212,212))
        self.textXR = wx.StaticText(self.panel, label="X:")
        self.editXR = wx.TextCtrl(self.panel, size=wx.Size(48,-1), style=wx.TE_PROCESS_ENTER)
        self.textYR = wx.StaticText(self.panel, label="Y:")
        self.editYR = wx.TextCtrl(self.panel, size=wx.Size(48,-1), style=wx.TE_PROCESS_ENTER)
        self.textXRR = wx.StaticText(self.panel, label="X:")
        self.editXRR = wx.TextCtrl(self.panel, size=wx.Size(48,-1), style=wx.TE_PROCESS_ENTER)
        self.textYRR = wx.StaticText(self.panel, label="Y:")
        self.editYRR = wx.TextCtrl(self.panel, size=wx.Size(48,-1), style=wx.TE_PROCESS_ENTER)

        self.sizer2.Add(wx.StaticText(self.panel, label="Data points (pixels):"),(0,0),(1,4),wx.EXPAND)
        self.sizer2.Add(self.textX, (1,1), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.editX, (1,2), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.textY, (1,3), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.editY, (1,4), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(wx.StaticText(self.panel, label="Range (micron):"),(2,0),(1,4),wx.EXPAND)
        self.sizer2.Add(self.textXR, (3,1), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.editXR, (3,2), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.textYR, (3,3), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.editYR, (3,4), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(wx.StaticText(self.panel, label="Resolution (micron per pixel):"),(2,5),(1,4),wx.EXPAND)
        self.sizer2.Add(self.textXRR, (3,6), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.editXRR, (3,7), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.textYRR, (3,8), wx.DefaultSpan, wx.ALL, 5)
        self.sizer2.Add(self.editYRR, (3,9), wx.DefaultSpan, wx.ALL, 5)

        #sizer 3: the sizer strikes back
        self.textSpot = wx.StaticText(self.panel, label="Spot size (micron):")
        self.editSpot = wx.TextCtrl(self.panel, size=wx.Size(48,-1))
        self.outputCsv = wx.CheckBox(self.panel, -1, 'Output CSVs')
        self.outputOther = wx.CheckBox(self.panel, -1, 'Output input & spot')
        
        self.sizer3.Add(self.textSpot, (0,1), wx.DefaultSpan, wx.ALL, 5)
        self.sizer3.Add(self.editSpot, (0,2), wx.DefaultSpan, wx.ALL, 5)
        self.sizer3.Add(self.outputCsv, (0,3), (1,2), wx.ALIGN_CENTER)
        self.sizer3.Add(self.outputOther, (0,5), (1,2), wx.ALIGN_CENTER)

        #sizer4: it's go time
        self.saveButton = wx.Button(self.panel, -1, 'Save to...')
        self.saveLocation = wx.TextCtrl(self.panel, size=wx.Size(256,-1))
        
        self.sizer4.Add(self.saveButton,(0,0), wx.DefaultSpan, wx.ALL, 5)
        self.sizer4.Add(self.saveLocation, (0,1), wx.DefaultSpan, wx.ALL, 5)

        #sizer5: bigg buttonz
        self.button2 = wx.Button(self.panel, -1, "Convolve")
        self.button3 = wx.Button(self.panel, -1, "Show")
        self.sizer5.Add(self.button2,(0,0), wx.DefaultSpan, wx.ALL, 5)
        self.sizer5.Add(self.button3,(0,1), wx.DefaultSpan, wx.ALL, 5)
        
        #events
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        self.Bind(wx.EVT_BUTTON, self.OnOpen, self.button1)
        self.Bind(wx.EVT_BUTTON, self.OnSaveLocation, self.saveButton)
        self.Bind(wx.EVT_BUTTON, self.Convolute, self.button2)
        self.Bind(wx.EVT_BUTTON, self.ShowFile, self.button3)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEditRange, self.editXR)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEditRange, self.editYR)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEditReso, self.editXRR)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEditReso, self.editYRR)

        #set defaults
        self.editSpot.SetValue(str(spotsize))
        self.editXR.SetValue(str(rangeX))
        self.editYR.SetValue(str(rangeY))
        self.filepath.SetValue(os.getcwd()+'\\'+filename)
        self.saveLocation.SetValue(os.getcwd()+'\\'+filename_export)

        #wrapping it up into one
        self.sizer = wx.GridBagSizer(2,5)
        self.sizer.Add(self.sizer1, (0,0), (1,1), wx.ALIGN_CENTER)
        self.sizer.Add(self.sizer2, (1,0), (1,1), wx.ALIGN_CENTER)
        self.sizer.Add(self.sizer3, (2,0), (1,1), wx.ALIGN_CENTER)
        self.sizer.Add(self.sizer4, (3,0), (1,1), wx.ALIGN_CENTER)
        self.sizer.Add(self.sizer5, (4,0), (1,1), wx.ALIGN_CENTER)

        #and go!
        self.SetSizer(self.sizer)
        self.SetSizerAndFit(self.sizer)
        self.sizer.Fit(self)
        
        self.Show(True)
        
    #definitions
    def ShowFile(self,e):
        if lastConvolved!='':
            print 'Showing '+str(lastConvolved)
            picture = Image.open(lastConvolved)
            picture.show()
        else:
            print 'Nothing has been convolved yet!'
        
    def Convolute(self,e):
        filename = self.filepath.GetValue()
        filename_export = self.saveLocation.GetValue()
        spotsize = float(self.editSpot.GetValue())
        rangeX = float(self.editXR.GetValue())
        rangeY = float(self.editYR.GetValue())
        filetype = os.path.splitext(filename_export)[1][1:]
        outputcsv=True if self.outputCsv.GetValue() == True else False  
        outputinterim=True if self.outputOther.GetValue() == True else False

        picture = Image.open(filename)
        width, height = picture.size
        self.editX.SetValue(str(width))
        self.editY.SetValue(str(height))
        self.editXRR.SetValue(str(rangeX/width))
        self.editYRR.SetValue(str(rangeY/height))
        
        convolve(filename,filename_export,spotsize,rangeX,rangeY,filetype,outputcsv,outputinterim)
        global lastConvolved
        lastConvolved = filename_export
        
        
    def OnSaveLocation(self,e):
        self.dirname=''
        filters = 'BMP (*.bmp)|*.bmp|PNG (*.png)|*.png'
        dlg = wx.FileDialog(self,"Save destination", self.dirname, "", filters, wx.SAVE | wx.OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.saveLocation.SetValue(self.dirname+'\\'+self.filename)
        dlg.Destroy()
        
    def OnEditRange(self,e):
        X = int(self.editX.GetValue())
        Y = int(self.editY.GetValue())
        XR = float(self.editXR.GetValue())
        YR = float(self.editYR.GetValue())
        self.editXRR.SetValue(str(XR/X))
        self.editYRR.SetValue(str(YR/Y))

    def OnEditReso(self,e):
        X = int(self.editX.GetValue())
        Y = int(self.editY.GetValue())
        XRR = float(self.editXRR.GetValue())
        YRR = float(self.editYRR.GetValue())
        self.editXR.SetValue(str(XRR*X))
        self.editYR.SetValue(str(YRR*Y))
        
    def OnAbout(self,e):
        dlg = wx.MessageDialog(self, "Written by "+str(author), "About "+str(progname)+" "+str(version), wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def OnExit(self,e):
        self.Close(True)

    def OnOpen(self,e):
        self.dirname=''
        dlg = wx.FileDialog(self,"Pick source file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.filepath.SetValue(self.dirname+'\\'+self.filename)
            picture = Image.open(filename)
            width, height = picture.size
            self.editX.SetValue(str(width))
            self.editY.SetValue(str(height))
            XR = float(self.editXR.GetValue())
            YR = float(self.editYR.GetValue())
            self.editXRR.SetValue(str(XR/width))
            self.editYRR.SetValue(str(YR/height))
        dlg.Destroy()
        
        
app = wx.App(False)
frame = MainWindow(None,str(progname)+" "+str(version))
app.MainLoop()
