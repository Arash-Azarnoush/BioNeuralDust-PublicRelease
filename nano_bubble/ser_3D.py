from ser import SmoothEndoRet
from vpython import dot, mag, cross, cylinder, vector, sphere, compound, norm

class SmoothEndoRet3D(SmoothEndoRet):

    def __init__(self, membrane_radius, channel_radius, length, num_chan, Vrest, Vmax, dt):

        SmoothEndoRet.__init__(self, num_chan, Vrest, Vmax, dt)

        self.membrane=cylinder(color=vector(0,1,0),
                               pos=vector(0,0,0),
                               axis=vector(1,0,0),
                               radius=membrane_radius, #25nm 
                               length=length, #880nm
                               opacity=0.6, 
                               shininess=5)

        num_channels = 20
        self.channels = []

        if((num_channels%2)==0):
            num_channels_hor=num_channels/2
            num_channels_vert=num_channels/2
        else:
            num_channels_hor=(num_channels+1)/2
            num_channels_vert=num_channels-num_channels_hor
        
        interval=length/(num_channels_vert+1)
        x=interval

        for _ in range(int(num_channels_vert)):
            for i in [-1,1]:
                x_ve=x
                y_ve=self.membrane.radius*i
                z_ve=0
                self.channels.append(sphere(color=vector(0,1,1), pos=vector(x_ve,y_ve,z_ve),radius=channel_radius,opacity=1,shininess=1))
            x = x + interval
            
        interval2=length/(num_channels_hor+1)
        y=interval2
        
        for _ in range(int(num_channels_hor)):
            for i in [-1,1]:
                x_ve=y
                y_ve=0
                z_ve=self.membrane.radius*i
                self.channels.append(sphere(color=vector(0,1,1), pos=vector(x_ve,y_ve,z_ve),radius=channel_radius,opacity=1,shininess=1))
            y = y + interval

        self.ser = compound((self.channels+[self.membrane]), origin=self.membrane.pos)
        self.end_point = self.ser.pos+(norm(self.ser.axis)*self.ser.length)


    #function to determine whether a molecule is inside the SER or not
    def isInside(self,coord):
        
        qp1=coord - self.ser.pos
        qp2=coord - self.end_point
        p1p2=self.end_point - self.ser.pos
        
        # math stuff (it works trust me I'm a pro)
        return ((dot(qp1,p1p2)>= 0) and (dot(qp2,p1p2)<= 0) and (mag(cross(qp1,p1p2))/mag(p1p2) <= (self.membrane.radius)))
        
        
    #these following class methods are not used by the main scripts but might be useful in the future
    def changePosition(self, new_pos):
        self.ser.pos = new_pos
        self.end_point = self.ser.pos+(norm(self.ser.axis)*self.ser.length)
    
    def changeAxis(self, new_axis):
        self.ser.axis = new_axis
        self.end_point = self.ser.pos+(norm(self.ser.axis)*self.ser.length)
        
    def getStartPosition(self):
        return self.ser.pos
        
    def getEndPosition(self):
        return self.end_point


    #these were initialy used to control the ingress/egress flows 
    #of the ions accross the membrane as they move randomly
    def toggleChannel(self, b):
        
        if b.text[:4]=="Open":
            self.open_channels=True
            b.text="Close Ionic Sesame"
        else:
            self.open_channels=False
            b.text="Open Ionic Sesame"
            
    def toggleIngress(self, b):

        if b.text[5:7]=="on":
           self.ingress_on=True
           b.text="Turn off absorbtion"
        else:
           self.ingress_on=False
           b.text="Turn on absorbtion" 
           
    def toggleEgress(self, b):

        if b.text[5:7]=="on":
           self.egress_on=True
           b.text="Turn off diffusion"
        else:
           self.egress_on=False
           b.text="Turn on diffusion"  