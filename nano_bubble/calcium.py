from math import exp

R = 8.314   #ideal gas constant
T = 311.650 #average human brain temperature in Kelvin
z = 2.000       #calcium ion valence
F = 96485.332   #Faraday constant
nernst_constant = (R*T)/(z*F)


class CalciumCluster:

    def __init__(self, Vrest, Vmax, init_extra_conc):

        #the SER's membrane potential determine's the
        # nernst ratio which is what controls the change 
        # of extra/intracellular concentration over time
        self.extra_conc = init_extra_conc
        self.intra_conc = init_extra_conc/exp(Vrest/nernst_constant)
        self.total_conc = self.extra_conc + self.intra_conc
        self.ratio = self.extra_conc/self.intra_conc

        #calculate the total amount of required aequoerin molecules () 
        ratio_max = exp(Vmax/nernst_constant)
        self.Pout_max = self.total_conc - (self.total_conc/(ratio_max + 1))


    def updateConcentrations(self,Vt):

        ratio = exp(Vt/nernst_constant)
        self.intra_conc = self.total_conc/(ratio + 1)
        self.extra_conc = self.total_conc - self.intra_conc

        return (ratio, self.intra_conc, self.extra_conc)


    def getPoutMax(self):
        return self.Pout_max


