from vpython import *
from diatomic_ import *

N = 20 # 20 molecules
L = ((24.4E-3/(6E23))*N)**(1/3.0)/50 # 2L is the length of the cubic container box, the number is made up
m = 14E-3/6E23 # average mass of O and C
k, T = 1.38E-23, 298.0 # some constants to set up the initial speed
initial_v = (3*k*T/m)**0.5 # some constant
class group_diatom:

    def __init__(self,CM_pos,L,H,W,_type_keyword = "N2",_N = N):
        self.molecules = []
        self.num = _N
        CM_p = 0
        _CM_Pos = 0
        for i in range(_N):
            rand_pos = vec((random()-0.5)*L, (random()-0.5)*H, (random()-0.5)*W)
            atom = diatomic_molecule(pos = rand_pos,type_keyword=_type_keyword)
            atom.A1.v = initial_v*vec(random(), random(), random())
            atom.A2.v = initial_v*vec(random(), random(), random())
            CM_p += atom.A1.v*atom.A1.m + atom.A2.v*atom.A2.m
            _CM_Pos += atom.com()*(atom.A1.m + atom.A2.m)
            self.molecules.append(atom)
        total_M = _N*(self.molecules[0].A1.m+self.molecules[0].A2.m)
        CM_v = CM_p/total_M
        d_CM_pos = _CM_Pos - CM_pos
        for i in range(_N):
            self.molecules[i].A1.v -= CM_v
            self.molecules[i].A2.v -= CM_v
            self.molecules[i].A1.pos -= d_CM_pos
            self.molecules[i].A2.pos -= d_CM_pos
        self.L = L
        self.H = H
        self.W = W
    def check_collide(self):
        for i in range(self.num-1):
            for j in range(i,self.num):
                self.molecules[i].check_collide(self.molecules[j])
        