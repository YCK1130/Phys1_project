from vpython import *
scene = canvas(width = 600, height =600, align = 'left', background = vec(102, 167, 223)/255)#,forward=vec(0,1,0))
g_AM = graph(title='Angular Momentum', width=450, height=300,align = 'right', background=vec(0.5,0.5,0),
               xtitle="<i>t</i>(s)", ytitle="I (kg m<sup>2</sup>)")
L_L = 4
L_W = 2
L_H = 8

S_L = 4
S_W = 2
S_H = 2

size = 1e-3
m = 1e-2 #1e-2
K = 1e6 #1e6
a = 1
b_coef = 1
neighbor = [] 
#sqrt(2)#上下左右
for i in (-1,0,1):
    for j in (-1,0,1):
        for k in (-1,0,1):
            if sqrt(i**2+j**2+k**2)>1.42 or sqrt(i**2+j**2+k**2)<0.99:continue 
            neighbor.append((i,j,k))

class fix_balls():
    '''1 for rectangle, 2 for T shape'''
    def __init__(self,num, unit = 1,_L_L = L_L,_L_H = L_H,_L_W = L_W,
     _S_L = S_L,_S_H=S_H,_S_W = S_W,size=size,hollow = False):
        self.elements_num = num
        self.elements = []
        self.my_index = []
        self.unit = unit
        self._a = []
        self._v = []
        self._pos = []

        self.L_L = _L_L
        self.L_H = _L_H
        self.L_W = _L_W

        self.S_L = _S_L
        self.S_H = _S_H
        self.S_W = _S_W

        LongSide = []
        LongSide_a = []
        LongSide_v = []
        LongSide_pos = []
        L_index = []
        #element 1, long side 
        for i in range(_L_L):
            side_y = []
            side_y_a = []
            side_y_v = []
            side_y_pos = []
            for j in range(_L_H):
                side_z = []
                side_z_a = []
                side_z_v = []
                side_z_pos = []

                for k in range(_L_W):
                    if hollow and (i not in (0,_L_L-1) and j not in (0,_L_H-1) and k not in (0,_L_W-1)):
                        side_z.append(None)
                        side_z_a.append(None)
                        side_z_v.append(None)
                        side_z_pos.append(None)
                        continue

                    side_z.append(sphere(pos = vec(i*self.unit,(j-_L_H/2)*self.unit,k*self.unit), radius = size, color = color.blue))
                    if hollow:
                        side_z[-1].opacity=0.1 
                    side_z[-1].opos=vec(i*self.unit,(j-_L_H/2)*self.unit,k*self.unit)
                    side_z[-1].m = m
                    side_z_a.append(vec(0,0,0))
                    side_z_v.append(vec(0,0,0))
                    side_z_pos.append(vec(i*self.unit,(j-_L_H/2)*self.unit,k*self.unit))
                    L_index.append((i,j,k))

                side_y.append(side_z)
                side_y_a.append(side_z_a)
                side_y_v.append(side_z_v)
                side_y_pos.append(side_z_pos)

            LongSide.append(side_y)
            LongSide_a.append(side_y_a)
            LongSide_v.append(side_y_v)
            LongSide_pos.append(side_y_pos)
        
        self.elements.append(LongSide)
        self._a.append(LongSide_a)
        self._v.append(LongSide_v)
        self._pos.append(LongSide_pos)
        self.my_index.append(L_index)
        
        #element 2 , short side
        if num >1:
            ShortSide = []
            S_index = []

            ShortSide_a = []
            ShortSide_v = []
            ShortSide_pos = []

            for i in range(_S_L):
                side_y = []
                side_y_a = []
                side_y_v = []
                side_y_pos = []
                for j in range(_S_H):
                    side_z = []
                    side_z_a = []
                    side_z_v = []
                    side_z_pos = []

                    for k in range(_S_W):
                        if hollow and (i not in (0,_S_L-1) and j not in (0,_S_H-1) and k not in (0,_S_W-1)):
                            side_z.append(None)
                            side_z_a.append(None)
                            side_z_v.append(None)
                            side_z_pos.append(None)
                            continue
                        side_z.append(sphere(pos = vec(-(i+1)*self.unit,(j-_S_H/2)*self.unit,k*self.unit), radius = size, color = color.blue))
                        if hollow:
                            side_z[-1].opacity=0.1    
                        side_z[-1].m = m
                        side_z[-1].opos=vec(-(i+1)*self.unit,(j-_S_H/2)*self.unit,k*self.unit)
                        side_z_a.append(vec(0,0,0))
                        side_z_v.append(vec(0,0,0))
                        side_z_pos.append(vec(-(i+1)*self.unit,(j-_S_H/2)*self.unit,k*self.unit))
                        S_index.append((i,j,k))
                    side_y.append(side_z)
                    side_y_a.append(side_z_a)
                    side_y_v.append(side_z_v)
                    side_y_pos.append(side_z_pos)
                ShortSide.append(side_y)
                ShortSide_a.append(side_y_a)
                ShortSide_v.append(side_y_v)
                ShortSide_pos.append(side_y_pos)

            self.elements.append(ShortSide)
            self._a.append(ShortSide_a)
            self._v.append(ShortSide_v)
            self._pos.append(ShortSide_pos)
            self.my_index.append(S_index)
        if hollow:
            self.ball_num = _L_L*_L_H*2+_L_H*_L_W*2+_L_L*_L_W*2 - _L_H*4 - _L_L*4 - _L_W*4 + 8
        else:
            self.ball_num = (_L_L*_L_H*_L_W)
            if self.elements_num>1:
                self.ball_num += (_S_L*_S_H*_S_W)
        self.set_CM()
        self.set_init_v()
        self.set_init_arr()
        self.curve_set(g_AM)
        self.plot_graph(0)

    def set_CM(self):
        C_M_pos = vec(0,0,0)
        
        for num in range(self.elements_num):
            for (i,j,k) in self.my_index[num]:
                C_M_pos += self._pos[num][i][j][k]/self.ball_num 
        self._CM = C_M_pos

    def set_init_v(self,velocity = 100):
        '''set init velocity with a small perturbation'''
        perturbation = vec((random()+1e-3)/100,0,0)
        for x in range(self.L_L):
            self._v[0][x][0][0] = velocity*cross(self._CM-self._pos[0][0][0][0],vec(-1,0,0)) + perturbation
            self._v[0][x][-1][-1] = velocity*cross(self._CM-self._pos[0][-1][-1][-1],vec(-1,0,0)) - perturbation
        print("perturbation =",perturbation)

    def set_init_arr(self):
        self.arr_list = []
        '''init rotating axis'''
        init_arr = arrow(pos = self._CM , axis = vec(-self.S_L-2,0,0)*self.unit, color=color.green,shaftwidth = 0.1*self.unit)

        '''three orthogonal axis'''
        self.myaxis = []
        if self.elements_num>1:
            R_axis3 =  (self._pos[0][0][-1][0]+self._pos[0][-1][-1][-1])/2 - (self._pos[0][0][0][0]+self._pos[0][-1][0][-1])/2
            R_axis2 =  (self._pos[1][-1][0][0]+self._pos[1][-1][-1][-1])/2 - (self._pos[1][0][0][0]+self._pos[1][0][-1][-1])/2
            R_axis1 = norm(cross(R_axis2,R_axis3))
        else:
            R_axis3 =  (self._pos[0][0][-1][0]+self._pos[0][-1][-1][-1])/2 - (self._pos[0][0][0][0]+self._pos[0][-1][0][-1])/2
            R_axis2 =  (self._pos[0][-1][0][0]+self._pos[0][-1][-1][-1])/2 - (self._pos[0][0][0][0]+self._pos[0][0][-1][-1])/2
            R_axis1 = norm(cross(R_axis2,R_axis3))
        self.myaxis.append(R_axis1.hat*self.L_W*self.unit/2)
        self.myaxis.append(R_axis2.hat*self.L_L*self.unit/2)
        self.myaxis.append(R_axis3.hat*self.L_H*self.unit/2)

        axis1_arr = arrow(pos = self._CM , axis = R_axis1.hat*self.L_W*self.unit/2, color=color.red,shaftwidth = 0.1*self.unit)
        axis2_arr = arrow(pos = self._CM , axis = R_axis2.hat*self.L_L*self.unit/2, color=color.blue,shaftwidth = 0.1*self.unit)
        axis3_arr = arrow(pos = self._CM , axis = R_axis3.hat*self.L_H*self.unit/2, color=color.yellow,shaftwidth = 0.1*self.unit)

        self.arr_list.append(init_arr)
        self.arr_list.append(axis1_arr)
        self.arr_list.append(axis2_arr)
        self.arr_list.append(axis3_arr)

        print("Moment of Inertia axis1 : {:.3f}".format(self.cal_I(R_axis1,self._CM)))
        print("Moment of Inertia axis2 : {:.3f}".format(self.cal_I(R_axis2,self._CM)))
        print("Moment of Inertia axis3 : {:.3f}".format(self.cal_I(R_axis3,self._CM)))

    def set_arr(self):
        
        if self.elements_num>1:
            R_axis3 =  (self._pos[0][0][-1][0]+self._pos[0][-1][-1][-1])/2 - (self._pos[0][0][0][0]+self._pos[0][-1][0][-1])/2
            R_axis2 =  (self._pos[1][-1][0][0]+self._pos[1][-1][-1][-1])/2 - (self._pos[1][0][0][0]+self._pos[1][0][-1][-1])/2
            R_axis1 = norm(cross(R_axis2,R_axis3))
        else:
            R_axis3 =  (self._pos[0][0][-1][0]+self._pos[0][-1][-1][-1])/2 - (self._pos[0][0][0][0]+self._pos[0][-1][0][-1])/2
            R_axis2 =  (self._pos[0][-1][0][0]+self._pos[0][-1][-1][-1])/2 - (self._pos[0][0][0][0]+self._pos[0][0][-1][-1])/2
            R_axis1 = norm(cross(R_axis2,R_axis3))
        self.myaxis[0] = R_axis1.hat*self.L_W*self.unit/2
        self.myaxis[1] = R_axis2.hat*self.L_L*self.unit/2
        self.myaxis[2] = R_axis3.hat*self.L_H*self.unit/2

    def curve_set(self,graph_):
        self._curve = []
        p_A1 = gcurve(color = color.red,width = 4,graph=graph_)
        p_A2 = gcurve(color = color.blue,width = 4,graph=graph_)
        p_A3 = gcurve(color = color.yellow,width = 4,graph=graph_)
        p_init = gcurve(color = color.green,width = 4,graph=graph_)
        self._curve.append(p_A1)
        self._curve.append(p_A2)
        self._curve.append(p_A3)
        self._curve.append(p_init)

    def plot_graph(self,t):
        for i in range(3):
            self._curve[i].plot(pos=(t,self.myaxis[i].hat.dot(self.cal_L(self.myaxis[i],self._CM))))
        self._curve[3].plot(pos=(t,self.arr_list[0].axis.hat.dot(self.cal_L(self.arr_list[0].axis.hat,self._CM))))
    
    def cal_I(self,rotation_axis,ref_point):
        Moment_of_Inertia = 0
        for num in range(self.elements_num):
            for (i,j,k) in self.my_index[num]:
                r = mag(cross(self._pos[num][i][j][k]-ref_point,rotation_axis.hat))
                Moment_of_Inertia += m*r**2
    
        return Moment_of_Inertia
    
    def cal_L(self,rotation_axis,ref_point):
        L_AM = vec(0,0,0)
        for num in range(self.elements_num):
            for (i,j,k) in self.my_index[num]:
                r = self._pos[num][i][j][k]-ref_point
                p = m*self._v[num][i][j][k]
                #print(i,j,k,r,p)
                L_AM += cross(r-r.dot(rotation_axis.hat)*rotation_axis.hat,p-p.dot(rotation_axis.hat)*rotation_axis.hat)
        return L_AM

    def time_elapse(self,dt,t):
        for num in range(self.elements_num):
            for (i, j, k) in self.my_index[num]:
                self._a[num][i][j][k] = vec(0,0,0)
                for delta in neighbor:
                    ni, nj, nk= i+delta[0],j+delta[1],k+delta[2]
                    if (ni,nj,nk) not in self.my_index[num]:continue
                    natural_L = mag(self.elements[num][i][j][k].opos-self.elements[num][ni][nj][nk].opos)
                    self._a[num][i][j][k] += -K/m*(mag(self._pos[num][i][j][k]-self._pos[num][ni][nj][nk])-natural_L)*norm(self._pos[num][i][j][k]-self._pos[num][ni][nj][nk])
        if self.elements_num>1:
            num1, num2= 0 ,1
            for (i1, j1, k1) in self.my_index[num1]:
                for (i2,j2,k2) in self.my_index[num2]:
                    o_distance = mag(self.elements[num2][i2][j2][k2].opos-self.elements[num1][i1][j1][k1].opos)
                    if  (o_distance>0.99*self.unit and o_distance <1.01*self.unit) or (o_distance>1.41*self.unit and o_distance <1.415*self.unit):
                        # print("hi",t)
                        self._a[num1][i1][j1][k1] += -K/m*(mag(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])-o_distance)*norm(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])
                        self._a[num2][i2][j2][k2] += K/m*(mag(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])-o_distance)*norm(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])
        for num in range(self.elements_num):
            for (i, j, k) in self.my_index[num]:
                # self._a[num][i][j][k] += -b_coef*self._v[num][i][j][k]
                self._v[num][i][j][k] += self._a[num][i][j][k]*dt
                self._pos[num][i][j][k] += self._v[num][i][j][k]*dt
        self.set_CM()
        self.set_arr()
        self.count_period(t)

    def count_period(self,_t):
        if(_t<0.01):
            self.now_I = 10
            self.prev_I = self.now_I
            self.prev_t = 0
        self.now_I = self.myaxis[1].hat.dot(self.cal_L(self.myaxis[1],self._CM))
        if self.now_I*self.prev_I<0 and abs(_t-self.prev_t) >= 1e-2 and _t>0.5:
            print("Period : {:.5f} sec @ {:.5f} , axis: {:}".format(2*(_t-self.prev_t),_t,self.myaxis[1].hat))
            self.prev_t = _t
        self.prev_I = self.now_I

    def scene_move(self):
        for num in range(self.elements_num):
            for (i, j, k) in self.my_index[num]:
                self.elements[num][i][j][k].pos = self._pos[num][i][j][k]       
        for i in range(3):
            self.arr_list[i+1].pos = self._CM
            self.arr_list[i+1].axis = self.myaxis[i]

    def total_K(self):
        _K_ = 0
        for num in range(self.elements_num):
            for (i,j,k) in self.my_index[num]:
                _K_ += 1/2*self.elements[0][0][0][0].m*mag2(self._v[num][i][j][k])
        return _K_


if __name__ == '__main__':
    '''1 for rectangle, 2 for T shape'''
    T_shape = fix_balls(1)#,_L_W=3,hollow=True) 

    dt = 5e-5
    t=0
    prev_t = 0
    times = 0  #plot graph per 1000 times 

    print("K =",K)
    while t<50:
        rate(5000)
        t+=dt
        times += 1
        T_shape.time_elapse(dt,t)
        
        T_shape.scene_move()
        if(times%100==0):
            T_shape.plot_graph(t)