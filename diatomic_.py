from vpython import *

size,  k_bond = 5e-3/2,  18600.0 # These numbers are all made up
m_o, m_c,m_n,m_h =16E-3/6E23, 12E-3/6E23, 14E-3/6E23, 1E-3/6E23 
mass_dict = {"CO" : (m_o,m_c), "O2":(m_o,m_o), "N2" : (m_n,m_n),"H2":(m_h,m_h)}
d = 2.5*size
dt = 1E-16
class diatomic_molecule:
    def __init__(self, pos,size=size, axis = d*vec(1,0,0), type_keyword = "O2"):
        self.A1 = sphere(pos = pos, radius = size, color = vec(255,215,0)/255)
        self.A2 = sphere(pos = pos+axis, radius = size, color = vec(162,205,90)/255)
        self.bond = cylinder(pos = pos, axis = axis, radius = size/3.0, color = color.white)
        self.A1.m = mass_dict[type_keyword][0]
        self.A2.m = mass_dict[type_keyword][1]
        self.A1.v = vector(0, 0, 0)
        self.A2.v = vector(0, 0, 0)
        self.bond.k = k_bond
    def bond_force_on_A1(self): # return bond force acted on the O atom
        return self.bond.k*(mag(self.bond.axis)-d)*norm(self.bond.axis)

    def time_lapse(self, dt): # by bond's force, calculate a, v and pos of C and O, and bond's pos and axis after dt
        self.A2.a = - self.bond_force_on_A1() / self.A2.m
        self.A1.a = self.bond_force_on_A1() / self.A1.m
        self.A2.v += self.A2.a * dt
        self.A1.v += self.A1.a * dt
        self.A2.pos += self.A2.v * dt
        self.A1.pos += self.A1.v * dt
        self.bond.axis = self.A2.pos - self.A1.pos
        self.bond.pos = self.A1.pos
    def com(self): # return position of center of mass
        return (self.A1.pos*self.A1.m+self.A2.pos*self.A2.m)/(self.A1.m + self.A2.m)
    def com_v(self): # return velocity of center of mass
        return (self.A1.v*self.A1.m+self.A2.v*self.A2.m)/(self.A1.m + self.A2.m)
    def v_P(self): # return potential energy of the bond for the vibration motion
        return 1/2*self.bond.k*(mag(self.bond.axis)-d)**2
    def v_K(self): # return kinetic energy of the vibration motion
        O_re_v = self.A1.v - self.com_v()
        C_re_v = self.A2.v - self.com_v()
        O_re_v_N = proj(O_re_v,self.bond.axis.hat)
        C_re_v_N = proj(C_re_v,self.bond.axis.hat)
        return 1/2*self.A1.m*mag2(O_re_v_N)+1/2*self.A2.m*mag2(C_re_v_N)
    def r_K(self): # return kinetic energy of the rotational motion
        O_re_v = self.A1.v - self.com_v()
        C_re_v = self.A2.v - self.com_v()
        O_re_v_T = O_re_v-proj(O_re_v,self.bond.axis.hat)
        C_re_v_T = C_re_v-proj(C_re_v,self.bond.axis.hat)
        return 1/2*self.A1.m*mag2(O_re_v_T)+1/2*self.A2.m*mag2(C_re_v_T)
    def com_K(self): #return kinetic energy of the translational motion of the center of mass
        return 1/2*(self.A1.m+self.A2.m)*mag2(self.com_v())
    def check_collide(self,another):
        '''if collide, change two atoms' velocity '''
        if(mag(self.A1.pos-another.A1.pos)<another.A1.radius*2) and dot(self.A1.v-another.A1.v,self.A1.pos-another.A1.pos)<0:
                self.A1.v,another.A1.v = collision(self.A1,another.A1)
        if(mag(self.A1.pos-another.A2.pos)<(another.A1.radius+another.A2.radius)) and dot(self.A1.v-another.A2.v,self.A1.pos-another.A2.pos)<0:
            self.A1.v,another.A2.v = collision(self.A1,another.A2)
        if(mag(self.A2.pos-another.A1.pos)<(another.A1.radius+another.A2.radius)) and dot(self.A2.v-another.A1.v,self.A2.pos-another.A1.pos)<0:
            self.A2.v,another.A1.v = collision(self.A2,another.A1)
        if(mag(self.A2.pos-another.A2.pos)<another.A2.radius*2) and dot(self.A2.v-another.A2.v,self.A2.pos-another.A2.pos)<0:
            self.A2.v,another.A2.v = collision(self.A2,another.A2)
        
def collision(a1, a2):
    v1prime = a1.v - 2 * a2.m/(a1.m+a2.m) *(a1.pos-a2.pos) * dot (a1.v-a2.v, a1.pos-a2.pos) / mag(a1.pos-a2.pos)**2
    v2prime = a2.v - 2 * a1.m/(a1.m+a2.m) *(a2.pos-a1.pos) * dot (a2.v-a1.v, a2.pos-a1.pos) / mag(a2.pos-a1.pos)**2
    return v1prime, v2prime

