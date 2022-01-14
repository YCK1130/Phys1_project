from vpython import *
scene = canvas(width = 600, height =600, align = 'left', background = vec(102, 167, 223)/255)#,forward=vec(0,1,0))
g_AM = graph(title='Angular Momentum', width=450, height=300,align = 'right', background=vec(0.5,0.5,0),
               xtitle="<i>t</i>(s)", ytitle="<i>I</i> (kg m<sup>2</sup>)")
L_L = 2 
L_W = 2
L_H = 8

S_L = 4
S_W = 2
S_H = 2

size = 2e-1
m = 1e-2
K = 1e6
a = 1
neighbor = [] 
#sqrt(2)#上下左右
for i in (-1,0,1):
    for j in (-1,0,1):
        for k in (-1,0,1):
            if sqrt(i**2+j**2+k**2)>1.42 or sqrt(i**2+j**2+k**2)<0.99:continue 
            neighbor.append((i,j,k))

class fix_balls():
    
    def __init__(self,num):
        self.elements_num = num
        self.elements = []
        self.my_index = []
        self._a = []
        self._v = []
        self._pos = []

        LongSide = []
        LongSide_a = []
        LongSide_v = []
        LongSide_pos = []
        L_index = []
        #element 1, long side 
        for i in range(L_L):
            side_y = []
            side_y_a = []
            side_y_v = []
            side_y_pos = []
            for j in range(L_H):
                side_z = []
                side_z_a = []
                side_z_v = []
                side_z_pos = []

                for k in range(L_W):
                    side_z.append(sphere(pos = vec(i*a,(j-L_H/2)*a,k*a), radius = size, color = color.blue))
                    side_z[-1].opos=vec(i*a,(j-L_H/2)*a,k*a)
                    side_z_a.append(vec(0,0,0))
                    side_z_v.append(vec(0,0,0))
                    side_z_pos.append(vec(i*a,(j-L_H/2)*a,k*a))

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

            for i in range(S_L):
                side_y = []
                side_y_a = []
                side_y_v = []
                side_y_pos = []
                for j in range(S_H):
                    side_z = []
                    side_z_a = []
                    side_z_v = []
                    side_z_pos = []

                    for k in range(S_W):
                        side_z.append(sphere(pos = vec(-(i+1)*a,(j-S_H/2)*a,k*a), radius = size, color = color.blue))
                        side_z[-1].opos=vec(-(i+1)*a,(j-S_H/2)*a,k*a)
                        side_z_a.append(vec(0,0,0))
                        side_z_v.append(vec(0,0,0))
                        side_z_pos.append(vec(-(i+1)*a,(j-S_H/2)*a,k*a))
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

        self.set_CM()
        self.set_init_v()
        self.set_init_arr()
        self.curve_set(g_AM)
        self.polt_graph(0)

    def set_CM(self):
        C_M_pos = vec(0,0,0)
        ball_num = (L_L*L_H*L_W)
        if self.elements_num>1:
            ball_num += (S_L*S_H*S_W)
        for num in range(self.elements_num):
            for (i,j,k) in self.my_index[num]:
                C_M_pos += self._pos[num][i][j][k]/ball_num 
        self._CM = C_M_pos

    def set_init_v(self,velocity = 100):
        
        perturbation = vec((random()+1e-3)/100,0,0)
        self._v[0][0][0][0] = velocity*cross(self._CM-self._pos[0][0][0][0],vec(-1,0,0)) + perturbation
        self._v[0][1][0][0] = velocity*cross(self._CM-self._pos[0][1][0][0],vec(-1,0,0))
        self._v[0][0][L_H-1][1] = velocity*cross(self._CM-self._pos[0][0][L_H-1][1],vec(-1,0,0))
        self._v[0][1][L_H-1][1] = velocity*cross(self._CM-self._pos[0][1][L_H-1][1],vec(-1,0,0)) - perturbation
        print("perturbation =",perturbation)

    def set_init_arr(self):
        self.arr_list = []
        '''init rotating axis'''
        init_arr = arrow(pos = self._CM-vec(1,0,0) , axis = vec(-S_L-2,0,0), color=color.green,shaftwidth = 0.1)

        '''counting C.M. , three orthogonal axis'''
        self.myaxis = []
        R_axis3 =  (self._pos[0][0][L_H-1][0]+self._pos[0][L_L-1][L_H-1][L_W-1])/2 - (self._pos[0][0][0][0]+self._pos[0][L_L-1][0][L_W-1])/2
        R_axis2 =  (self._pos[1][S_L-1][0][0]+self._pos[1][S_L-1][S_H-1][S_W-1])/2 - (self._pos[1][0][0][0]+self._pos[1][0][S_H-1][S_W-1])/2
        R_axis1 = norm(cross(R_axis2,R_axis3))

        self.myaxis.append(R_axis1)
        self.myaxis.append(R_axis2)
        self.myaxis.append(R_axis3)

        axis1_arr = arrow(pos = self._CM , axis = R_axis1, color=color.red,shaftwidth = 0.1)
        axis2_arr = arrow(pos = self._CM , axis = R_axis2, color=color.blue,shaftwidth = 0.1)
        axis3_arr = arrow(pos = self._CM , axis = R_axis3, color=color.yellow,shaftwidth = 0.1)

        self.arr_list.append(init_arr)
        self.arr_list.append(axis1_arr)
        self.arr_list.append(axis2_arr)
        self.arr_list.append(axis3_arr)

        print("Moment of Inertia axis1 : {:.3f}".format(self.cal_I(R_axis1,self._CM)))
        print("Moment of Inertia axis2 : {:.3f}".format(self.cal_I(R_axis2,self._CM)))
        print("Moment of Inertia axis3 : {:.3f}".format(self.cal_I(R_axis3,self._CM)))

    def set_arr(self):
        
        R_axis3 =  (self._pos[0][0][L_H-1][0]+self._pos[0][L_L-1][L_H-1][L_W-1])/2 - (self._pos[0][0][0][0]+self._pos[0][L_L-1][0][L_W-1])/2
        R_axis2 =  (self._pos[1][S_L-1][0][0]+self._pos[1][S_L-1][S_H-1][S_W-1])/2 - (self._pos[1][0][0][0]+self._pos[1][0][S_H-1][S_W-1])/2
        R_axis1 = norm(cross(R_axis2,R_axis3))

        self.myaxis[0] = R_axis1
        self.myaxis[1] = R_axis2
        self.myaxis[2] = R_axis3
        for i in range(3):
            self.arr_list[i+1].pos = self._CM
            self.arr_list[i+1].axis = self.myaxis[i]

    def curve_set(self,graph_):
        self._curve = []
        p_A1 = gcurve(color = color.red,width = 4,graph=graph_)
        p_A2 = gcurve(color = color.blue,width = 4,graph=graph_)
        p_A3 = gcurve(color = color.yellow,width = 4,graph=graph_)
        self._curve.append(p_A1)
        self._curve.append(p_A2)
        self._curve.append(p_A3)

    def polt_graph(self,t):
        for i in range(3):
            self._curve[i].plot(pos=(t,self.myaxis[i].hat.dot(self.cal_L(self.myaxis[i],self._CM))))

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
                L_AM += cross(r-r.dot(rotation_axis)*rotation_axis,p-p.dot(rotation_axis)*rotation_axis)
        return L_AM

    def time_elapse(self,dt):
        for num in range(self.elements_num):
            for (i, j, k) in self.my_index[num]:
                self._a[num][i][j][k] = vec(0,0,0)
                for delta in neighbor:
                    ni, nj, nk= i+delta[0],j+delta[1],k+delta[2]
                    if (ni,nj,nk) not in self.my_index[num]:continue
                    natural_L = mag(self.elements[num][i][j][k].opos-self.elements[num][ni][nj][nk].opos)
                    self._a[num][i][j][k] += -K/m*(mag(self._pos[num][i][j][k]-self._pos[num][ni][nj][nk])-natural_L*a)*norm(self._pos[num][i][j][k]-self._pos[num][ni][nj][nk])
        num1 = 0
        num2 = 1
        for (i1, j1, k1) in self.my_index[num1]:
            for (i2,j2,k2) in self.my_index[num2]:
                o_distance = mag(self.elements[num2][i2][j2][k2].opos-self.elements[num1][i1][j1][k1].opos)
                if  (o_distance>0.99*a and o_distance <1.01*a) or (o_distance>1.41*a and o_distance <1.415*a):
                    # print("hi",t)
                    self._a[num1][i1][j1][k1] += -K/m*(mag(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])-o_distance*a)*norm(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])
                    self._a[num2][i2][j2][k2] += K/m*(mag(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])-o_distance*a)*norm(self._pos[num1][i1][j1][k1]-self._pos[num2][i2][j2][k2])
        for num in range(self.elements_num):
            for (i, j, k) in self.my_index[num]:
                self._v[num][i][j][k] += self._a[num][i][j][k]*dt
                self._pos[num][i][j][k] += self._v[num][i][j][k]*dt
        self.set_CM()

    def count_period(self,t):
        if(t<0.01):
            self.now_I = 10
            self.prev_I = self.now_I
            self.prev_t = 0

        if self.now_I*self.prev_I<0 and abs(t-self.prev_t) >= 1e-2 and t>0.5:
            print("Period : {:.5f} sec".format(2*(t-self.prev_t)))
            self.prev_t = t

    def scene_move(self):
        for num in range(self.elements_num):
            for (i, j, k) in self.my_index[num]:
                self.elements[num][i][j][k].pos = self._pos[num][i][j][k] 
        self.set_arr()



if __name__ == '__main__':
    T_shape = fix_balls(2) 

    dt = 5e-5
    t=0
    prev_t = 0
    times = 0  #polt graph per 1000 times 

    print("K =",K)
    while t<50:
        rate(5000)
        t+=dt
        times += 1
        T_shape.time_elapse(dt)

        if(times%50==0):
            T_shape.scene_move()
        if(times%1000==0):
            T_shape.polt_graph(t)