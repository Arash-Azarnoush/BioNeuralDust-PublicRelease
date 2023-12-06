from math import exp

# {num_channels: (Rca [Ohms], tau [seconds])}
channel_resistance_tau = {1: (5.28E14,0.237),   
                          50: (2.11E12,9.5E-5),
                          100: (5.28E11,2.37E-5),
                          500: (2.11E10,9.5E-7),
                          1000: (5.28E9,2.37E-7)}



class SmoothEndoRet:
        
    def __init__(self, num_chan, Vrest, Vmax, dt):

        self.tau = channel_resistance_tau[num_chan][1]
        self.Vrest = Vrest
        self.Vmax = Vmax
        self.Vb = 0
        self.steady_state_V = Vrest
        self.dt = dt
        self.Vt = 0

    def updateVoltage(self, Vm13):
        self.Vb = Vm13
        self.steady_state_V = self.Vrest + Vm13

    def nextTimeStep(self):        
        if self.Vt >= self.Vmax:  
            self.Vt = self.Vrest
        else:
            self.Vt = self.steady_state_V + (self.Vt - self.steady_state_V)*exp(-self.dt/self.tau)
            
        return self.Vt
