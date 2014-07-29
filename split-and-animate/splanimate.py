import numpy as np
import matplotlib.pyplot as plt
import csv

fileName = r"08 Scanner_Counts_stepZ.txt"
gridSize = (42,42)
gridAxisX = (0,40)
gridAxisY = (0,40)
normalization = 0.7 #can use this factor to convert to seconds (normally 0.7) or to QE if you're a crafty bugger

fileHandle = open(fileName,'r')
csvHandle = csv.reader(fileHandle,delimiter=" ")

processed_data = []
maximum = 0

for r, row in enumerate(csvHandle):
	if r%(gridSize[0]-1) == 0:
		processed_data.append([])
	if gridSize[0] == len(row):
		processed_data[-1].append([])
		for cell in row[:-1]:
			processed_data[-1][-1].append(float(cell)/float(normalization))
			if float(cell)/float(normalization)>maximum:
				maximum = float(cell)/float(normalization)
	else:
		print 'gridSize does not accurately describe the width of this file!'
		print row
		print 'This is a typical row.'
		print 'There is a row length of',len(row),'while a length of',gridSize[0],'is specified.'
		raise IOError

x = np.linspace(gridAxisX[0],gridAxisX[1],gridSize[0]-1)
y = np.linspace(gridAxisY[0],gridAxisY[1],gridSize[1]-1)


for p, processed_datum in enumerate(processed_data):
	image = np.array(processed_datum)
	plt.pcolor(x,y,image,vmin=0,vmax=maximum)
	plt.xlabel('X displacement (nm)')
	plt.ylabel('Y displacement (nm)')
	plt.title('Z step:'+str(p))

	plt.axis([gridAxisX[0],gridAxisX[1],gridAxisY[0],gridAxisY[1]])

	plt.jet()
	plt.colorbar(format='%d')
	fileStart = str(p)
	while len(fileStart)<len(str(len(processed_data))):
		fileStart = '0'+fileStart
	plt.savefig(fileStart+'.png')
	plt.clf()
