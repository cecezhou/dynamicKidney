import random
import math
import numpy as np

max_time = 1. 
lam = 10000
t = 0.
state_BA = 0
matches = 0
matches_dict = {}
od_pairs = [['AB','B'],['A','O'],['B','O'],['AB','O'],['AB','A']]
while (t < max_time):
	i_t = random.expovariate(lam)
	t = t + i_t
	elements = ['A', 'AB', 'B', 'O']
	freqs = np.array([224806., 26036., 98893., 329133.])
	probs = freqs / sum(freqs)
	types = list(np.random.choice(elements, 2, replace=False, p=probs))
	patient_type = types[0]
	donor_type = types[1]
	if (types == ['B', 'A']):
		if state_BA < 0:
			matches += 2
			matches_dict[t] = matches
		state_BA += 1
	elif (types == ['A', 'B']):
		if state_BA > 0:
			matches += 2
			matches_dict[t] = matches
		state_BA -= 1
	elif types in od_pairs:
		matches += 2
		matches_dict[t] = matches

print matches
# print matches_dict




			

