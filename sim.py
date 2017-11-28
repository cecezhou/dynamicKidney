import random
import math
import numpy as np

c = 1
rho = 0.05
max_time = 1. 
lam = 5
init_t = 0.
state_BA = 0
matches = 0
matches_dict = {0.: 0}
od_pairs = [['AB','B'],['A','O'],['B','O'],['AB','O'],['AB','A']]
p_c = 0.11 # tissue incompatibility probability

def two_way_sim():
	t = init_t
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
			r_unif = np.random.uniform(0., 1.)
			if r_unif < p_c:
				matches += 2
				matches_dict[t] = matches

	print matches
	print matches_dict
	return (matches, matches_dict)

def calc_surplus(matches_dict):
	sorted_times = sorted(matches_dict)
	surplus = 0.
	for i in range(len(sorted_times)):
		L = sorted_times[i]
		if i == (len(sorted_times) - 1):
			U = max_time
		else:
			U = sorted_times[i+1]
		surplus += (c * (matches_dict[L] / rho) * (np.exp(-1 * rho * L) - np.exp(-1 * rho * U)))
	return surplus




(matches_2, matches_dict_2) = two_way_sim()
surplus_2 = calc_surplus(matches_dict_2)


			

