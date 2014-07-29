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

from scipy import fftpack, array, signal
import csv
from math import *
import numpy as np
try:
    import Image, ImageDraw
except ImportError:
    print 'Python Image Library required!'
    raise ImportError

##setup
filename = r'input.bmp' #the image to be convolved
filename_export = r'convolved.bmp'
spotsize = 2.0 #FWHM, in microns; assumed equal in X and Y
rangeX = 24.0 #in microns, the image's range in X
rangeY = 24.0 #in microns, the image's range in Y
##end setup

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

def exportBitmap(filename, image_data, shape):
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
    
    picture.save(filename,"BMP")

numericalData, shape = importBitmap(filename)

print 'Read input.bmp'

numfile = open('out-input.csv','wb')
numcsv = csv.writer(numfile, delimiter=',')
for row in numericalData:
    numcsv.writerow(row)
numfile.close()

exportBitmap('out-input.bmp',numericalData,shape)

print 'Wrote out-input.bmp'

##create gaussian spot profile of correct size
stepsX,stepsY = shape
centerX = float(stepsX-1)/2.0
centerY = float(stepsY-1)/2.0
stepLengthX = float(rangeX)/float(stepsX)
stepLengthY = float(rangeY)/float(stepsY)
def gaussian2(x, y, FWHM, xOffset, yOffset):
    x,y,FWHM,xOffset,yOffset = float(x),float(y),float(FWHM),float(xOffset),float(yOffset)
    max = 1.0
    sigma = FWHM / 2.3548
    output = max * exp(-1*((x - xOffset)**2 + ((y - yOffset)**2)) / (2*(sigma**2)))
    return output
#spot = []
#for y in range(stepsY):
#    spot.append([])
#    for x in range(stepsX):
#        spot[y].append(gaussian2(x*stepLengthX,y*stepLengthY,spotsize,1,centerX*stepLengthX,centerY*stepLengthY))
#spot = array(spot)

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

print 'Created Gaussian'

spotfile = open('spot.csv','wb')
outspot = csv.writer(spotfile,delimiter=',')
for row in spot:
    outspot.writerow(row)
spotfile.close()

exportBitmap('spot.bmp',spot,spot.shape)

print 'Wrote spot.bmp'
    
##perform deconvolution
#fft & multiply
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
outfile = open('convolved.csv','wb')
outcsv = csv.writer(outfile, delimiter=',')
for row in convolvedImage:
    outcsv.writerow(row)
outfile.close()

exportBitmap(filename_export, convolvedImage, shape)

print 'Wrote output, done!'
