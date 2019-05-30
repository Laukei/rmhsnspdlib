import csv
from math import exp

#datapoints are 10 per ns

OUTPUT = 'snspd_signal.csv'
I_BIAS = 20e-6 #bias current, amps (typically tens of microamps)
L_K = 1000e-9 #kinetic inductance, Henrys (typically several hundred nanohenrys)
R_N = 1000 #detector resistance (typically few kohm)
R_SHUNT = 50 #shunt resistance (typically 50 ohm)
T_START = 0 #start time for pulse (seconds)
V_PEAK = 100e-6 # pre-amplification peak voltage (100s of microvolts volts)
ITERATIONS = 70000 #number of datapoints to generate for signal
TIME_FACTOR = 10 #number of points per nanosecond
NORMALIZE = True #whether to set max output to 1

class PulseGenerator:
	def __init__(self,*args,**kwargs):
		self._I_bias = kwargs.get('I_bias',I_BIAS)
		self._L_k = kwargs.get('L_k',L_K)
		self._R_n = kwargs.get('R_n',R_N)
		self._R_shunt = kwargs.get('R_shunt',R_SHUNT)
		self._T_start = kwargs.get('T_start',T_START)
		self._V_peak = kwargs.get('V_peak',V_PEAK)
		self._T_stop = None
		self._time_factor = kwargs.get('t_factor',TIME_FACTOR)
		self.stage = 0
		self.tick = 0
		self.voltage = 0
		self.normalization_peak = None
		self.data = []

	def _Tau_fall(self):
		return self._L_k / (self._R_shunt + self._R_n)

	def _Tau_rise(self):
		return self._L_k / self._R_shunt

	def step(self):
		self._check_stage()
		data = self._evolve()
		if not self.normalization_peak or data > self.normalization_peak:
			self.normalization_peak = data
		self._increment()
		return data

	def _increment(self):
		self.tick += 1

	def _return_superconducting(self):
		'''
		device returns to superconductivity when I_snspd**2 * R_n < P_cooling
		and I_snspd depends on I_bias redistribution through shunt which
		in turn depends on L_k
		'''
		if self.voltage >= self._V_peak:
			return True
		else:
			return False

	def _check_stage(self):
		if self.stage == 0 and self._current_time() >= self._T_start:
			self.stage = 1
		elif self.stage == 1 and self._return_superconducting():
			self.stage = 2
			self._T_stop = self._current_time()

	def _current_time(self):
		return self.tick * (1e-9 / self._time_factor)

	def _evolve(self):
		if self.stage == 1:
			self.voltage = exp((self._current_time() - self._T_start)/self._Tau_fall()) - 1
		elif self.stage == 2:
			self.voltage = self.normalization_peak * exp(-(self._current_time() - self._T_stop)/self._Tau_rise())
		if self.voltage > self._V_peak:
			self.voltage = self._V_peak
		return self.voltage

	def normalization_factor(self):
		return self.normalization_peak

def create_dataset(pulse_generator,iterations,updown,normalize=NORMALIZE):
	dataset = []
	for i in range(iterations):
		dataset.append([pulse_generator.step(),int(i//updown / 2 == int (i//updown / 2)),0])
	if normalize == True:
		for r, row in enumerate(dataset):
			dataset[r][0] = row[0]/pulse_generator.normalization_factor()
	return dataset

def write_dataset(dataset,filename):
	with open(filename,'w',newline='\n') as f:
		writer = csv.writer(f,delimiter=',')
		writer.writerows(dataset)
	print('done')


if __name__ == "__main__":
	p = PulseGenerator()
	dataset = create_dataset(p,ITERATIONS,125)
	write_dataset(dataset,OUTPUT)