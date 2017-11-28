import random
import math
import numpy as np


class dynamicKidney():

	def __init__(self):
		self.c = 1
		self.rho = 0.05
		self.max_time = 1. 
		self.lam = 100
		self.elements = ['A', 'AB', 'B', 'O']
		self.freqs = np.array([224806., 26036., 98893., 329133.])
		self.probs = probs = self.freqs / sum(self.freqs)
		self.od_pairs = [['AB','B'],['A','O'],['B','O'],['AB','O'],['AB','A']]
		self.p_c = 0.11 # tissue incompatibility probability
		self.mult_thresh = 2

	def init_params(self):
		self.t = 0
		self.state_BA = 0
		self.matches = 0
		self.matches_dict = {0.: 0}
	
	def calc_surplus(self):
		sorted_times = sorted(self.matches_dict)
		surplus = 0.
		for i in range(len(sorted_times)):
			L = sorted_times[i]
			if i == (len(sorted_times) - 1):
				U = self.max_time
			else:
				U = sorted_times[i+1]
			surplus += (self.c * (self.matches_dict[L] / self.rho) * (np.exp(-1 * self.rho * L) - np.exp(-1 * self.rho * U)))
		return surplus
	
	def two_way_sim(self):
		self.init_params()
		while (self.t < self.max_time):
			i_t = random.expovariate(self.lam)
			self.t += i_t
			types = list(np.random.choice(self.elements, 2, replace=False, p=self.probs))
			patient_type = types[0]
			donor_type = types[1]
			if (types == ['B', 'A']):
				if self.state_BA < 0:
					self.matches += 2
					self.matches_dict[self.t] = self.matches
				self.state_BA += 1
			elif (types == ['A', 'B']):
				if self.state_BA > 0:
					self.matches += 2
					self.matches_dict[self.t] = self.matches
				self.state_BA -= 1
			elif types in self.od_pairs:
				r_unif = np.random.uniform(0., 1.)
				if r_unif < self.p_c:
					self.matches += 2
					self.matches_dict[self.t] = self.matches
		print self.matches,self.matches_dict
		return (self.matches, self.matches_dict)
	
	def maximal_matching(self, pair):
		## TODO 
		return 2
	## paper assumes BA is the relevant state variable
	def multi_way_sim(self):
		self.init_params()
		while (self.t < self.max_time):
			i_t = random.expovariate(self.lam)
			self.t += i_t
			types = list(np.random.choice(self.elements, 2, replace=False, p=self.probs))
			patient_type = types[0]
			donor_type = types[1]
			if (types == ['B', 'A']):
				## match with AB if possible
				if self.state_BA < 0:
					self.matches += 2
					self.matches_dict[self.t] = self.matches
				self.state_BA += 1
			elif (types == ['A', 'B']):
				## match with BA if possible
				if self.state_BA > 0:
					self.matches += 2
					self.matches_dict[self.t] = self.matches
				self.state_BA -= 1
			elif (types in self.od_pairs):
				r_unif = np.random.uniform(0., 1.)
				if r_unif < self.p_c:
					# matching AB-O, O-B (O-A), (B-A) (A-B), A-AB (B-AB)
					if types == ['AB', 'O']:
						if self.state_BA < 0:
							self.matches += 4
							self.state_BA += 1
						elif self.state_BA > self.mult_thresh:
							self.matches +=4
							self.state_BA -=1
						else:
							self.matches += 3
						self.matches_dict[self.t] = self.matches
					# matching AB-A, (A-B), A-AB (B-AB)
					if types == ['AB', 'A']:
						if self.state_BA < 0:
							self.matches += 3
							self.state_BA +=1
						else:
							self.matches += 2
						self.matches_dict[self.t] = self.matches
					# matching AB-B, (B-A), B-AB (A-AB)
					if types == ['AB', 'B']:
						if self.state_BA > self.mult_thresh:
							self.matches += 3
							self.state_BA -= 1
						else: 
							self.matches += 2
						self.matches_dict[self.t] = self.matches
					## compare with threshold to decide whether to use AB pair (always use BA)


mySim = dynamicKidney()
(matches_2, matches_dict_2) = mySim.two_way_sim()
surplus_2 = mySim.calc_surplus()
print(surplus_2)

			

