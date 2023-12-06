from math import pi, factorial, exp

 
def getProbabilityOfEmission(flux, t):

    C = (flux*t)
    q = pow(C,3)/factorial(3)

    return q*exp(-C)

def getProbabilityOfNoEmission(flux, t):

    C = (flux*t) 

    return (1+C)*exp(-C)

def getRate(prob, Ta):

   Tlem = 3*Ta

   return (1/Tlem)*(1-prob)  

def getMaxLightIntensity(concentration):
    
    Na_hc = 0.1196
    lamb = 4.8E-7

    return (concentration * Na_hc)/(3*lamb)

    
class PhotonEmission:
        
    def __init__(self, init_Pout, Pout_max, aequ_radius, diffusion_coef, dt):
        
        #setup init photon emission variables
        self.diffusion_coef = diffusion_coef
        self.flux = 4*pi*self.diffusion_coef*aequ_radius*(init_Pout*6.02214076E23)
        self.Ta = 1/self.flux

        C = (self.flux*dt)
        self.prob_no_emission = (1+C)*exp(-C)

        Tlem = 3*self.Ta
        self.radius = aequ_radius
        self.rate = (1/Tlem)*(1-self.prob_no_emission) 
       
        #total amount of aequorin is about 1/3 of max extracellular concentration
        self.total_aequorin = ((Pout_max)/3)*6.02214076E23
        
        
        self.total_emission_per_dt = self.rate * self.total_aequorin / 6.02214076E23

        self.output_intensity = self.total_emission_per_dt * 249222.0104166667  


    def getFlux(self, Cion):

        
        flux = 4*pi*self.diffusion_coef*self.radius*(Cion*6.02214076E23)

        return (flux,(1/flux))
 

    def getLightIntensity(self, concentration, dt):
        #calculate flux
        self.flux, self.Ta = self.getFlux(concentration)

        #probability of no photons being emited
        self.prob_no_emission = getProbabilityOfNoEmission(self.flux, dt)

        #rate of photons emited within dt time interval
        self.rate = getRate(self.prob_no_emission, self.Ta)

        #total photon molecules emited within dt time interval
        self.total_emission_per_dt = self.rate * self.total_aequorin / 6.02214076E23

        #equivalent light intensity of total emited photons
        self.output_intensity = self.total_emission_per_dt * 249222.0104166667  

        return self.output_intensity


    