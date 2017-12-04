from __future__ import print_function

import random
import math
import numpy as np
import sys


class dynamicKidney():

	def __init__(self, thresh = 2, perturb = 0.1):
		self.c = 1
		self.rho = 0.05
		self.max_time = 1. 
		self.lam = 10000
		self.elements = ['A', 'AB', 'B', 'O']
		self.freqs = np.array([224806., 26036., 98893., 329133.])
		self.probs = probs = self.freqs / sum(self.freqs)
		self.od_pairs = [['AB','B'],['A','O'],['B','O'],['AB','O'],['AB','A']]
		self.ud_pairs = [['AB','O'],['O','A'],['O','B'],['A','AB'],['B','AB']]
		self.p_c = 0.11 # tissue incompatibility probability
		self.mult_thresh = thresh
		self.perturb = perturb

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
	
	def flip_types(self, types):
		return types
				# sometimes switch to AB if we get BA
		if (types == ['B', 'A']):	
			r_unif = np.random.uniform(0., 1.)
			if r_unif < self.perturb:
				types = ['A', 'B']
		return types

	def two_way_sim(self):
		self.init_params()
		while (self.t < self.max_time):
			i_t = random.expovariate(self.lam)
			self.t += i_t
			types = list(np.random.choice(self.elements, 2, replace=True, p=self.probs))
			types = self.flip_types(types)

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
		# print self.matches
		return (self.matches, self.matches_dict)
	
	def maximal_matching(self, pair):
		## TODO 
		return 2
	## paper assumes BA is the relevant state variable
	def mult_way_sim(self):
		self.init_params()
		while (self.t < self.max_time):
			i_t = random.expovariate(self.lam)
			self.t += i_t
			types = list(np.random.choice(self.elements, 2, replace=True, p=self.probs))
			types = self.flip_types(types)
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
					# matching B-O, O-B (O-A), A-B
					if types == ['B', 'O']:
						if self.state_BA < 0:
							self.matches += 3
							self.state_BA += 1
						else:
							self.matches += 2
						self.matches_dict[self.t] = self.matches
					# matching A-O, O-A (O-B), B-A
					if types == ['A', 'O']:
						if self.state_BA > self.mult_thresh:
							self.matches += 3
							self.state_BA -= 1
						else:
							self.matches += 2
						self.matches_dict[self.t] = self.matches
		# print self.matches
		return (self.matches, self.matches_dict)		

	# run two-way sim starting from empty world to find lower bound on steady state
	def empty_world_sim(self):
		self.init_params()
		current_surplus = 0
		max_time_limit = 2
		ud_pairs = {('O','AB'): 0, ('O','A'): 0, ('O','B'): 0, ('A','AB'): 0, ('B','AB'): 0}
		od_pairs = {('AB','O'): 0, ('A','O'): 0, ('B','O'): 0, ('AB','A'): 0, ('AB','B'): 0}
		# while current not within some epsilon of steady state surplus
		while (self.t < max_time_limit):
			i_t = random.expovariate(self.lam)
			self.t += i_t
			types = list(np.random.choice(self.elements, 2, replace=True, p=self.probs))
			types = self.flip_types(types)

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
			# check if ud reciprocal available
			elif types in self.od_pairs:
				r_unif = np.random.uniform(0., 1.)
				if r_unif < self.p_c:
					# look for reciprocal
					if (types == ['AB','O']):
						if ud_pairs[('O','AB')] > 0:
							ud_pairs[('O','AB')] -= 1
							self.matches += 2
							self.matches_dict[self.t] = self.matches
						else:
							od_pairs[('AB','O')] += 1
					elif (types == ['A','O']):
						if ud_pairs[('O','A')] > 0:
							ud_pairs[('O','A')] -= 1
							self.matches += 2
							self.matches_dict[self.t] = self.matches
						else:
							od_pairs[('A','O')] += 1
					elif (types == ['B','O']):
						if ud_pairs[('O','B')] > 0:
							ud_pairs[('O','B')] -= 1
							self.matches += 2
							self.matches_dict[self.t] = self.matches
						else:
							od_pairs[('B','O')] += 1
					elif (types == ['AB','A']):
						if ud_pairs[('A','AB')] > 0:
							ud_pairs[('A','AB')] -= 1
							self.matches += 2
							self.matches_dict[self.t] = self.matches
						else:
							od_pairs[('AB','A')] += 1
					elif (types == ['AB','B']):
						if ud_pairs[('B','AB')] > 0:
							ud_pairs[('B','AB')] -= 1
							self.matches += 2
							self.matches_dict[self.t] = self.matches
						else:
							od_pairs[('AB','B')] += 1
			elif types in self.ud_pairs:
				# look for reciprocal
				if (types == ['O','AB']):
					if od_pairs[('AB','O')] > 0:
						od_pairs[('AB','O')] -= 1
						self.matches += 2
						self.matches_dict[self.t] = self.matches
					else:
						ud_pairs[('O','AB')] += 1					
				elif (types == ['O','A']):
					if od_pairs[('A','O')] > 0:
						od_pairs[('A','O')] -= 1
						self.matches += 2
						self.matches_dict[self.t] = self.matches
					else:
						ud_pairs[('O','A')] += 1
				elif (types == ['O','B']):
					if od_pairs[('B','O')] > 0:
						od_pairs[('B','O')] -= 1
						self.matches += 2
						self.matches_dict[self.t] = self.matches
					else:
						ud_pairs[('O','B')] += 1
				elif (types == ['A','AB']):
					if od_pairs[('AB','A')] > 0:
						od_pairs[('AB','A')] -= 1
						self.matches += 2
						self.matches_dict[self.t] = self.matches
					else:
						ud_pairs[('A','AB')] += 1
				elif (types == ['B','AB']):
					if od_pairs[('AB','B')] > 0:
						od_pairs[('AB','B')] -= 1
						self.matches += 2
						self.matches_dict[self.t] = self.matches
					else:
						ud_pairs[('B','AB')] += 1
		return (self.matches, self.matches_dict)

	def find_steady_state(self):
		sorted_times = sorted(self.matches_dict)
		surplus = 0.
		start_time = 0
		end_time = 0
		end_time_index = 0 # index of last time before 1
		# find the last time before 1
		for i in range(len(sorted_times)):
			if(sorted_times[i] > 1):
				end_time = sorted_times[i - 1]
				end_time_index = i - 1
				break
		# calculate first surplus (up to 1)
		for i in range(end_time_index):
			L = sorted_times[i]
			if i == end_time_index:
				U = 1
			else:
				U = sorted_times[i + 1]
			surplus += (self.c * (self.matches_dict[L] / self.rho) * (np.exp(-1 * self.rho * L) - np.exp(-1 * self.rho * U)))
		# subtract 1, add 1 until within steady state surplus
		current_start_index = 1
		current_end_time = 1 + sorted_times[current_start_index]
		end_time_index += 1
		# end_time_index is currently next index after 1.0
		while ((surplus > 765 or surplus < 755) and current_end_time < 2):
			# print ("Current End", current_end_time)
			# print ("Previous End", sorted_times[end_time_index])
			# print ("Current Start", sorted_times[current_start_index])
			# subtract 1
			L = sorted_times[current_start_index - 1]
			U = sorted_times[current_start_index]
			surplus -= (self.c * (self.matches_dict[L] / self.rho) * (np.exp(-1 * self.rho * L) - np.exp(-1 * self.rho * U)))
			# maybe add some
			while (end_time_index < len(sorted_times) and sorted_times[end_time_index] < current_end_time):
				L = sorted_times[end_time_index - 1]
				U = sorted_times[end_time_index]
				surplus += (self.c * (self.matches_dict[L] / self.rho) * (np.exp(-1 * self.rho * L) - np.exp(-1 * self.rho * U)))
				end_time_index += 1
			# move interval
			current_start_index += 1
			current_end_time = 1 + sorted_times[current_start_index]
			# print ("Surplus", surplus)
		return (sorted_times[current_start_index], current_start_index)

num_trials = 50
trials_surplus = []
# steady_state_times = []
# num_arrivals = []

mySim = dynamicKidney()
thresholds = [i for i in xrange(0, 4)]
for threshold in thresholds:
	mySim.mult_thresh = threshold
	print ("Threshold: ", threshold)
	for i in xrange(num_trials):
		(matches_2, matches_dict_2) = mySim.mult_way_sim()
		trials_surplus.append(mySim.calc_surplus())
	print (np.mean(trials_surplus))
	print (np.var(trials_surplus))	

# for i in xrange(num_trials):
# 	mySim.empty_world_sim()
# 	result = mySim.find_steady_state()
# 	steady_state_times.append(result[0])
# 	num_arrivals.append(result[1])
# print("Times Mean:", np.mean(steady_state_times))
# print("Times Variance:", np.var(steady_state_times))
# print("Arrivals Mean:", np.mean(num_arrivals))
# print("Arrivals Variance:", np.var(num_arrivals))

# trials_surplus = []

# for i in xrange(num_trials):
# 	(matches_2, matches_dict_2) = mySim.two_way_sim()
# 	trials_surplus.append(mySim.calc_surplus())
# print np.mean(trials_surplus)
# print np.var(trials_surplus)

# mySim = dynamicKidney()
# n = 50
# n = 10
# for _ in range(n):
# 	avg_surplus = 0
# 	(matches_mult, matches_dict_mult) = mySim.mult_way_sim()
# 	surplus_mult = mySim.calc_surplus()
# 	avg_surplus += surplus_mult
# print(avg_surplus/n)
# for perturb in [0.0, 0.05, 0.1, 0.15]:
# 	mySim.perturb = perturb
# 	print("Perturb", perturb)
# 	for thresh in [(perturb * 100 + 1) * x for x in range(15)]:
# 		mySim.mult_thresh = thresh
# 		avg_surplus = 0
# 		for _ in range(n):
# 			# print(".", end = "")
# 			# ys.stdout.flush()
# 			(matches_mult, matches_dict_mult) = mySim.mult_way_sim()
# 			surplus_mult = mySim.calc_surplus()
# 			avg_surplus += surplus_mult
# 		print("Threshold: ", thresh, ", Surplus:", avg_surplus/n)

# mySim = dynamicKidney(thresh = 2, perturb = 0.05)
# (matches_2, matches_dict_2) = mySim.two_way_sim()
# surplus_2 = mySim.calc_surplus()

# for i in xrange(num_trials):
# 	(matches_2, matches_dict_2) = mySim.mult_way_sim()
# 	trials_surplus.append(mySim.calc_surplus())
# print (np.mean(trials_surplus))
# print (np.var(trials_surplus))		
