import csv
import random

random.seed() #herp derp

number_of_errors 	= 10 			#the exact number of errors to distribute across the film
film_size 			= 10000			#the size of the film; arbitrary, 1D, and used only to determine the number of steps
film_size_step 		= 100 			#on balance not required, but film_size/film_step_size = number of data points to take
trials 				= 10000 		#number of trials per film_size; higher = higher quality output, lower = faster
filename 			= 'output.txt'  #output filename


def createErrors(number_of_errors, film_size):
	errors = []
	for error in range(number_of_errors):
		errors.append(random.random()*film_size)
	return errors

def isItPerfect(value,errors):
	result = 0
	for error in errors:
		if error < value:
			return 0
	return 1


data = []
for size in range(0,film_size,film_size_step)+[film_size]:
	print 'length',size,'of',film_size
	result = 0
	for trial in range(trials):
		errors = createErrors(number_of_errors,film_size)
		result += isItPerfect(float(size),errors)
	perfection = float(result)/float(trials)
	data.append([size,perfection])



file_handle = open(filename,'w')
file_writer = csv.writer(file_handle,delimiter=',')
for row in data:
	file_writer.writerow(row)
file_handle.close()

print '-fin-'