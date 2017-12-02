import random
import math
import numpy as np


class dynamicKidney():
    
    def __init__(self, lam=10000, rho=0.05, mu=1., sigma2=1.):
        self.lam = lam
        self.mu = mu # mean of marginal costs
        self.sigma2 = sigma2 # variance of marginal costs
        self.rho = rho
        self.max_time = 1.
        self.elements = ['A', 'AB', 'B', 'O']
        self.freqs = np.array([224806., 26036., 98893., 329133.])
        self.probs = probs = self.freqs / sum(self.freqs)
        self.od_pairs = [['AB','B'],['A','O'],['B','O'],['AB','O'],['AB','A']]
        self.p_c = 0.11 # tissue incompatibility probability
        self.mult_thresh = 2 # threshold for multi-way mechanism
    
    def init_params(self):
        self.t = 0
        self.state_BA = 0
        self.matches = 0
        self.matches_dict = {}
        self.cum_matches_dict = {0.: 0}
    
    def calc_surplus_random_mc(self):
        sorted_times = sorted(self.matches_dict)
        surplus = 0.
        for tau in sorted_times:
            nmatches = self.matches_dict[tau]
            m_costs = np.random.gamma(shape=(self.mu ** 2)/(self.sigma2), scale=(self.sigma2 / self.mu), size=nmatches) 
            surplus += (sum(m_costs) / self.rho * (np.exp(-1 * self.rho * tau) - np.exp(-1 * self.rho * self.max_time)))
        return surplus
    
    def calc_surplus(self):
        sorted_times = sorted(self.cum_matches_dict)
        surplus = 0.
        for i in range(len(sorted_times)):
            L = sorted_times[i]
            if i == (len(sorted_times) - 1):
                U = self.max_time
            else:
                U = sorted_times[i+1]
            surplus += (self.mu * (self.cum_matches_dict[L] / self.rho) * (np.exp(-1 * self.rho * L) - np.exp(-1 * self.rho * U)))
        return surplus
    
    def two_way_sim(self):
        self.init_params()
        while (self.t < self.max_time):
            i_t = random.expovariate(self.lam)
            self.t += i_t
            types = list(np.random.choice(self.elements, 2, replace=True, p=self.probs))
            patient_type = types[0]
            donor_type = types[1]
            if (types == ['B', 'A']):
                if self.state_BA < 0:
                    self.update_matches(2)
                self.state_BA += 1
            elif (types == ['A', 'B']):
                if self.state_BA > 0:
                    self.update_matches(2)
                self.state_BA -= 1
            elif types in self.od_pairs:
                r_unif = np.random.uniform(0., 1.)
                if r_unif < self.p_c:
                    self.update_matches(2)
        return (self.matches, self.matches_dict, self.cum_matches_dict)
    
    def maximal_matching(self, pair):
        ## TODO
        return 2

    def update_matches(self, nmatches):
        self.matches += nmatches
        self.matches_dict[self.t] = nmatches
        self.cum_matches_dict[self.t] = self.matches

    ## paper assumes BA is the relevant state variable
    def mult_way_sim(self):
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
                    self.update_matches(2)
                self.state_BA += 1
            elif (types == ['A', 'B']):
                ## match with BA if possible
                if self.state_BA > 0:
                    self.update_matches(2)
                self.state_BA -= 1
            elif (types in self.od_pairs):
                r_unif = np.random.uniform(0., 1.)
                if r_unif < self.p_c:
                    # matching AB-O, O-B (O-A), (B-A) (A-B), A-AB (B-AB)
                    if types == ['AB', 'O']:
                        if self.state_BA < 0:
                            self.update_matches(4)
                            self.state_BA += 1
                        elif self.state_BA > self.mult_thresh:
                            self.update_matches(4)
                            self.state_BA -=1
                        else:
                            self.update_matches(3)
                    # matching AB-A, (A-B), A-AB (B-AB)
                    if types == ['AB', 'A']:
                        if self.state_BA < 0:
                            self.update_matches(3)
                            self.state_BA +=1
                        else:
                            self.update_matches(2)
                    # matching AB-B, (B-A), B-AB (A-AB)
                    if types == ['AB', 'B']:
                        if self.state_BA > self.mult_thresh:
                            self.update_matches(3)
                            self.state_BA -= 1
                        else:
                            self.update_matches(2)
                    # matching B-O, O-B (O-A), A-B
                    if types == ['B', 'O']:
                        if self.state_BA < 0:
                            self.update_matches(3)
                            self.state_BA += 1
                        else:
                            self.update_matches(2)
                            # matching A-O, O-A (O-B), B-A
                    if types == ['A', 'O']:
                        if self.state_BA > self.mult_thresh:
                            self.update_matches(3)
                            self.state_BA -= 1
                        else:
                            self.update_matches(2)
        return (self.matches, self.matches_dict, self.cum_matches_dict)

for var in [0.000001]:
    n_sims = 1
    surpluses_2 = []
    surpluses_2_random_mc = []
    surpluses_mult = []
    surpluses_mult_random_mc = []
    for sim in range(n_sims):
        if (sim % 10 == 0):
            print 'Simulation', sim
        mySim = dynamicKidney(lam=10000, sigma2=var)
        (matches_2, matches_dict_2, cum_matches_dict_2) = mySim.two_way_sim()
        surpluses_2.append(mySim.calc_surplus())
        surpluses_2_random_mc.append(mySim.calc_surplus_random_mc())
        (matches_mult, matches_dict_mult, cum_matches_dict_mult) = mySim.mult_way_sim()
        surpluses_mult.append(mySim.calc_surplus())
        surpluses_mult_random_mc.append(mySim.calc_surplus_random_mc())
    print 'Var =', var
    print 'Two way regular mean, stdev:', np.mean(surpluses_2), np.std(surpluses_2)
    print 'Two way random MC mean, stdev:', np.mean(surpluses_2_random_mc), np.std(surpluses_2_random_mc)
    print 'Multi-way regular mean, stdev:', np.mean(surpluses_mult), np.std(surpluses_mult)
    print 'Multi-way random MC mean, stdev:', np.mean(surpluses_mult_random_mc), np.std(surpluses_mult_random_mc)





