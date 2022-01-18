from vpython import *
from diatomic_ import *
from fix_ball import fix_balls
N = 2 # 20 molecules
L = ((24.4E-3)*20)**(1/3.0)/40 # 2L is the length of the cubic container box, the number is made up
H = 5*L/4
W = 3*L/4
nL,nH,nW = 4,5,3
unit = 2*L/(nL-1)
m = 1e-2/3 # average mass of O and C
atom_size = 0.5
initial_v = 100 # some constant
# scene = canvas(width = 400, height =400, align = 'left', background = vec(1, 1, 1))

class group_diatom:

    def __init__(self,CM_pos,L,H,W,_type_keyword = "N2",_N = N):
        self.molecules = []
        self.num = _N
        _CM_p = vec(0,0,0)
        _CM_Pos = vec(0,0,0)
        for i in range(_N):
            rand_pos = vec((random()-0.5)*L, (random()-0.5)*H, (random()-0.5)*W)
            rand_axis = d*vec(random(),random(),random()).hat
            atom = diatomic_molecule(pos = rand_pos,axis = rand_axis,type_keyword=_type_keyword)
            atom.A1.v = initial_v*vec(random(), random(), random())
            atom.A2.v = initial_v*vec(random(), random(), random())
            _CM_p += atom.A1.v*atom.A1.m + atom.A2.v*atom.A2.m
            _CM_Pos += atom.com()*(atom.A1.m + atom.A2.m)
            self.molecules.append(atom)
        total_M = _N*(self.molecules[0].A1.m+self.molecules[0].A2.m)
        CM_v = _CM_p/total_M
        d_CM_pos = _CM_Pos - CM_pos
        for i in range(_N):
            self.molecules[i].A1.v -= CM_v
            self.molecules[i].A2.v -= CM_v
            self.molecules[i].A1.pos -= d_CM_pos
            self.molecules[i].A2.pos -= d_CM_pos
            self.molecules[i].bond.pos -= d_CM_pos
        self.L = L
        self.H = H
        self.W = W
    def check_collide(self):
        for i in range(self.num-1):
            for j in range(i,self.num):
                self.molecules[i].check_collide(self.molecules[j])

    def time_elapse(self,dt):
        for i in range(self.num):
            self.molecules[i].time_elapse(dt)
    

container_ball = fix_balls(1,unit,nL,nH,nW,size=7e-3,hollow=True)
# container = box(pos=container_ball._CM,length = 2*L, height = 2*H, width = 2*W, opacity=0.4, color = color.yellow )
gas = group_diatom(container_ball._CM,L,H,W)

dt = 5e-5
t=0
prev_t = 0
times = 0  #polt graph per 1000 times 

# print("K =",K)
while t<50:
    rate(1000)
    t+=dt
    times += 1
    gas.check_collide()
    for _atom_ in gas.molecules:
        for (i,j,k) in container_ball.my_index[0]:
            if(mag(_atom_.A1.pos-container_ball.elements[0][i][j][k].pos)<(container_ball.elements[0][i][j][k].radius + _atom_.A1.radius)) and dot(_atom_.A1.v-container_ball._v[0][i][j][k],_atom_.A1.pos-container_ball.elements[0][i][j][k].pos)<0:
                v1prime = _atom_.A1.v - 2 * container_ball.elements[0][i][j][k].m/(_atom_.A1.m+container_ball.elements[0][i][j][k].m) *(_atom_.A1.pos-container_ball.elements[0][i][j][k].pos) * dot (_atom_.A1.v-container_ball._v[0][i][j][k], _atom_.A1.pos-container_ball.elements[0][i][j][k].pos) / mag(_atom_.A1.pos-container_ball.elements[0][i][j][k].pos)**2
                v2prime = container_ball._v[0][i][j][k] - 2 * _atom_.A1.m/(_atom_.A1.m+container_ball.elements[0][i][j][k].m) *(container_ball.elements[0][i][j][k].pos-_atom_.A1.pos) * dot (container_ball._v[0][i][j][k]-_atom_.A1.v, container_ball.elements[0][i][j][k].pos-_atom_.A1.pos) / mag(container_ball.elements[0][i][j][k].pos-_atom_.A1.pos)**2
                _atom_.A1.v , container_ball._v[0][i][j][k] = v1prime,v2prime
            if(mag(_atom_.A2.pos-container_ball.elements[0][i][j][k].pos)<(container_ball.elements[0][i][j][k].radius + _atom_.A2.radius)) and dot(_atom_.A2.v-container_ball._v[0][i][j][k],_atom_.A2.pos-container_ball.elements[0][i][j][k].pos)<0:
                v1prime = _atom_.A2.v - 2 * container_ball.elements[0][i][j][k].m/(_atom_.A1.m+container_ball.elements[0][i][j][k].m) *(_atom_.A1.pos-container_ball.elements[0][i][j][k].pos) * dot (_atom_.A2.v-container_ball._v[0][i][j][k], _atom_.A1.pos-container_ball.elements[0][i][j][k].pos) / mag(_atom_.A1.pos-container_ball.elements[0][i][j][k].pos)**2
                v2prime = container_ball._v[0][i][j][k] - 2 * _atom_.A1.m/(_atom_.A1.m+container_ball.elements[0][i][j][k].m) *(container_ball.elements[0][i][j][k].pos-_atom_.A1.pos) * dot (container_ball._v[0][i][j][k]-_atom_.A2.v, container_ball.elements[0][i][j][k].pos-_atom_.A1.pos) / mag(container_ball.elements[0][i][j][k].pos-_atom_.A1.pos)**2
                _atom_.A2.v , container_ball._v[0][i][j][k] = v1prime,v2prime
    gas.time_elapse(dt)
    container_ball.time_elapse(dt,t)
    
    # if(times%50==0):
    container_ball.scene_move()
    if(times%1000==0):
        container_ball.plot_graph(t)
