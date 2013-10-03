import copy

#P = [['(1-a)','a','0'],['0','(1-a/2)','a/2'],['0','0','1']]
P = [['b','a','0'],['0','c','d'],['0','0','1']]
#P = [['(1-a)','a'],['0','1']]

def square_matrix(P,n):
	n = n-1
	origP = copy.deepcopy(P)
	for a in range(n):
		newP = copy.deepcopy(P)
		for i, row in enumerate(P):
			for j, cell in enumerate(row):
				temp = ''
				for k in range(len(P)):
					if origP[i][k]!='0' and P[k][j]!='0':
						if origP[i][k] != '1' and P[k][j] != '1':
							temp += '('+origP[i][k]+')*('+P[k][j]+')'
						elif origP[i][k] == '1' and P[k][j] != '1':
							temp += '('+P[k][j]+')'
						elif origP[i][k] != '1' and P[k][j] == '1':
							temp += '('+origP[i][k]+')'
						elif origP[i][k] == '1' and P[k][k] == '1':
							temp += '1'
						else:
							print origP[i][k],P[k][j],': both are 1'
						if k<len(P)-1:
							temp += ' + '
				if temp == '':
					temp = '0'
				elif temp[-3:] == ' + ':
					temp = temp[:-3]
				elif temp[-2:] == '*1':
					temp = temp[:-2]
				newP[i][j]=temp
		P = copy.deepcopy(newP)
	return P


#for row in sqMat:
#	print row

for a in range(1,6):
	sqMat = square_matrix(P,a)
	#for row in sqMat:
	#	print row
	print 'a:',a,':',sqMat[0][-1]