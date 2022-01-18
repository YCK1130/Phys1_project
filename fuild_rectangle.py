from vpython import *
from diatomic_ import *
from fix_ball import fix_balls

g_K = graph(title='kinetic energy', width=450, height=300,align = 'right', background=vec(0.5,0.5,0),
               xtitle="<i>t</i>(s)", ytitle="K")
p_K_gas = gcurve(color = color.red,width = 4,graph=g_K)
p_K_box = gcurve(color = color.blue,width = 4,graph=g_K)
N = 10 # 20 molecules
L = ((24.4E-3)*20)**(1/3.0)/40 *2 # 2L is the length of the cubic container box, the number is made up
H = 10*L/8
W = 6*L/8
nL,nH,nW = 8,10,6
unit = 2*L/(nL-1)
m = 1e-2/3 # average mass of O and C
atom_size = 0.5
initial_v = 10 # some constant
# scene = canvas(width = 400, height =400, align = 'left', background = vec(1, 1, 1))

class group_diatom:

    def __init__(self,CM_pos,L,H,W,_type_keyword = "N2",_N = N):
        self.molecules = []
        self.num = _N
        _CM_p = vec(0,0,0)
        _CM_Pos = vec(0,0,0)
        for i in range(_N):
            rand_pos = vec((random()-0.5)*L, (random()-0.5)*H, (random()-0.5)*W)
            rand_axis = d*vec.random().hat
            atom = diatomic_molecule(pos = rand_pos,axis = rand_axis,type_keyword=_type_keyword)
            atom.A1.v = initial_v*vec.random().hat
            atom.A2.v = initial_v*vec.random().hat
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
K_energy = 0
for _atom_ in gas.molecules:
    K_energy += _atom_.total_K()
print(K_energy,container_ball.total_K())
p_K_gas.plot(pos=(t,K_energy))
p_K_box.plot(pos=(t,container_ball.total_K()))
# print("K =",K)
while t<50:
    rate(1000)
    t+=dt
    times += 1
    gas.check_collide()

    # box_cm = container_ball._CM
    box_mass = container_ball.ball_num*container_ball.elements[0][0][0][0].m
    for _atom_ in gas.molecules:
        for (i,j,k) in container_ball.my_index[0]:
            if(mag(_atom_.A1.pos-container_ball._pos[0][i][j][k])<(container_ball.elements[0][i][j][k].radius + _atom_.A1.radius)) and dot(_atom_.A1.v-container_ball._v[0][i][j][k],_atom_.A1.pos-container_ball._pos[0][i][j][k])<0:
                v1prime = - 2 * box_mass/(_atom_.A1.m+box_mass) *(_atom_.A1.pos-container_ball._pos[0][i][j][k]) * dot (_atom_.A1.v-container_ball._v[0][i][j][k], _atom_.A1.pos-container_ball._pos[0][i][j][k]) / mag2(_atom_.A1.pos-container_ball._pos[0][i][j][k])
                # v2prime = - 2 * _atom_.A1.m/(_atom_.A1.m+box_mass) *(container_ball._pos[0][i][j][k]-_atom_.A1.pos) * dot (container_ball._v[0][i][j][k]-_atom_.A1.v, container_ball._pos[0][i][j][k]-_atom_.A1.pos) / mag2(container_ball._pos[0][i][j][k]-_atom_.A1.pos)
                # v2prime = -v1prime*_atom_.A1.m/box_mass
                p2prime = -v1prime*_atom_.A1.m
                _atom_.A1.v += v1prime
                # _atom_.A1.v , container_ball._v[0][i][j][k] = v1prime,v2prime
                
                r_CM = container_ball._pos[0][i][j][k]-container_ball._CM # position vec with respect to CM
                p_on_r_CM = proj(p2prime,r_CM)
                p_T_r_CM = p2prime - p_on_r_CM

                r_axis = cross(r_CM,p2prime.hat)
                I_axis = container_ball.cal_I(r_axis,container_ball._CM)
                
                w_axis = cross(r_CM,p2prime)/I_axis
                for (ni,nj,nk) in container_ball.my_index[0]:
                    container_ball._v[0][ni][nj][nk] += cross(w_axis,container_ball._pos[0][ni][nj][nk]-container_ball._CM) + p_on_r_CM/box_mass

            
            
            
            if(mag(_atom_.A2.pos-container_ball._pos[0][i][j][k])<(container_ball.elements[0][i][j][k].radius + _atom_.A2.radius)) and dot(_atom_.A2.v-container_ball._v[0][i][j][k],_atom_.A2.pos-container_ball._pos[0][i][j][k])<0:
                v1prime = - 2 * box_mass/(_atom_.A2.m+box_mass) *(_atom_.A2.pos-container_ball._pos[0][i][j][k]) * dot (_atom_.A2.v-container_ball._v[0][i][j][k], _atom_.A2.pos-container_ball._pos[0][i][j][k]) / mag(_atom_.A2.pos-container_ball._pos[0][i][j][k])**2
                # v2prime = - 2 * _atom_.A2.m/(_atom_.A2.m+container_ball.elements[0][i][j][k].m) *(container_ball._pos[0][i][j][k]-_atom_.A2.pos) * dot (container_ball._v[0][i][j][k]-_atom_.A2.v, container_ball._pos[0][i][j][k]-_atom_.A2.pos) / mag(container_ball._pos[0][i][j][k]-_atom_.A2.pos)**2
                p2prime = -v1prime*_atom_.A2.m
                _atom_.A2.v += v1prime
                # _atom_.A2.v , container_ball._v[0][i][j][k] = v1prime,v2prime
                
                r_CM = container_ball._pos[0][i][j][k]-container_ball._CM # position vec with respect to CM
                p_on_r_CM = proj(p2prime,r_CM)
                p_T_r_CM = p2prime - p_on_r_CM

                r_axis = cross(r_CM,p2prime.hat)
                I_axis = container_ball.cal_I(r_axis,container_ball._CM)
                
                w_axis = cross(r_CM,p2prime)/I_axis
                for (ni,nj,nk) in container_ball.my_index[0]:
                    container_ball._v[0][ni][nj][nk] += cross(w_axis,container_ball._pos[0][ni][nj][nk]-container_ball._CM) + p_on_r_CM/box_mass

    gas.time_elapse(dt)
    container_ball.time_elapse(dt,t)
    
    # if(times%50==0):
    container_ball.scene_move()

    
    
    if(times%100==0):
        container_ball.plot_graph(t)
        
        K_energy = 0
        for _atom_ in gas.molecules:
            K_energy += _atom_.total_K()
        p_K_gas.plot(pos=(t,K_energy))
        p_K_box.plot(pos=(t,container_ball.total_K()))
