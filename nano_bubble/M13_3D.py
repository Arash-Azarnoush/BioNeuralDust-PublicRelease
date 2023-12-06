from M13 import M13Phage
from vpython import cylinder, vector


class M13Phage3D(M13Phage):

    def __init__(self, file, radius, length, position):

        M13Phage.__init__(self, file)

        self.bacteriophage=cylinder(color=vector(0,0,1),
                        pos=position,
                        axis=vector(1,-1,1),
                        radius=radius,   #6nm 
                        length=length, #880nm
                        opacity=1, 
                        shininess=5)