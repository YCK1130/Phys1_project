from vpython import *
scene = canvas(width = 600, height =600, align = 'left', background = vec(102, 167, 223)/255)#,forward=vec(0,1,0))
g_AM = graph(title='Angular Momentum', width=450, height=300,align = 'right', background=vec(0.5,0.5,0),
               xtitle="<i>t</i>(s)", ytitle="<i>I</i> (kg m<sup>2</sup>)")
p_A1 = gcurve(color = color.red,width = 4,graph=g_AM)
p_A2 = gcurve(color = color.blue,width = 4,graph=g_AM)
p_A3 = gcurve(color = color.yellow,width = 4,graph=g_AM)


def cal_I(rotation_axis,ref_point):
    Moment_of_Inertia = 0
    for i in range(L_L):
        for j in range(L_H):
            for k in range(L_W):
                r = mag(cross(LongSide[i][j][k].pos-ref_point,rotation_axis.hat))
                Moment_of_Inertia += m*r**2
    for i in range(S_L):
        for j in range(S_H):
            for k in range(S_W):
                r = mag(cross(ShortSide[i][j][k].pos-ref_point,rotation_axis.hat))
                Moment_of_Inertia += m*r**2
    return Moment_of_Inertia



def cal_L(rotation_axis,ref_point):
    L_AM = vec(0,0,0)
    for i in range(L_L):
        for j in range(L_H):
            for k in range(L_W):
                r = LongSide[i][j][k].pos-ref_point
                p = m*LongSide[i][j][k].v
                L_AM += cross(r-r.dot(rotation_axis)*rotation_axis,p-p.dot(rotation_axis)*rotation_axis)
    for i in range(S_L):
        for j in range(S_H):
            for k in range(S_W):
                r = ShortSide[i][j][k].pos-ref_point
                p = m*ShortSide[i][j][k].v
                L_AM += cross(r-r.dot(rotation_axis)*rotation_axis,p-p.dot(rotation_axis)*rotation_axis)
    return L_AM

L_L = 2 
L_W = 2
L_H = 8

S_L = 4 
S_W = 2
S_H = 2

size = 2e-1
m = 1e-2
K = 1e6
a=1
natural_L = 1
neighbor = [] 
#sqrt(2)#上下左右
for i in (-1,0,1):
    for j in (-1,0,1):
        for k in (-1,0,1):
            if sqrt(i**2+j**2+k**2)>1.42 or sqrt(i**2+j**2+k**2)<0.99:continue 
            neighbor.append((i,j,k))



#print(neighbor)
LongSide = []
L_index = []
for i in range(L_L):
    side_y = []
    for j in range(L_H):
        side_z = []
        for k in range(L_W):
            side_z.append(sphere(pos = vec(i*a,(j-L_H/2)*a,k*a),v = vec(0,0,0), a =vec(0,0,0), radius = size, color = color.blue))
            side_z[-1].opos=vec(i*a,(j-L_H/2)*a,k*a)
            L_index.append((i,j,k))
        side_y.append(side_z)
    LongSide.append(side_y)
middle = (LongSide[0][0][0].pos + LongSide[1][L_H-1][1].pos)/2

velocity = 100
perturbation = vec((random()+1e-3)/100,0,0)
LongSide[0][0][0].v = velocity*cross(middle-LongSide[0][0][0].pos,vec(-1,0,0)) + perturbation
LongSide[1][0][0].v = velocity*cross(middle-LongSide[1][0][0].pos,vec(-1,0,0))
LongSide[0][L_H-1][1].v = velocity*cross(middle-LongSide[0][L_H-1][1].pos,vec(-1,0,0))
LongSide[1][L_H-1][1].v = velocity*cross(middle-LongSide[1][L_H-1][1].pos,vec(-1,0,0)) - perturbation

print("perturbation =",perturbation)


# print(LongSide[L_L-1][L_H-1][L_W-1].a)
#bigstone = sphere(pos=vec(L/2,2,W/2),v =vec(0,0,0), a=g,radius=size*3,color=color.red)
ShortSide = []
S_index = []
for i in range(S_L):
    side_y = []
    for j in range(S_H):
        side_z = []
        for k in range(S_W):
            side_z.append(sphere(pos = vec(-(i+1)*a,(j-S_H/2)*a,k*a),v = vec(0,0,0), a =vec(0,0,0), radius = size, color = color.blue))
            side_z[-1].opos=vec(-(i+1)*a,(j-S_H/2)*a,k*a)
            S_index.append((i,j,k))
        side_y.append(side_z)
    ShortSide.append(side_y)

init_arr = arrow(pos = middle-vec(1,0,0) , axis = vec(-S_L-2,0,0), color=color.green,shaftwidth = 0.1)
C_M_pos = vec(0,0,0)
for i in range(L_L):
    for j in range(L_H):
        for k in range(L_W):
            C_M_pos += LongSide[i][j][k].pos/((L_L*L_H*L_W)+(S_L*S_H*S_W))
for i in range(S_L):
    for j in range(S_H):
        for k in range(S_W):
            C_M_pos += ShortSide[i][j][k].pos/((L_L*L_H*L_W)+(S_L*S_H*S_W))
R_axis3 =  (LongSide[0][L_H-1][0].pos+LongSide[L_L-1][L_H-1][L_W-1].pos)/2 - (LongSide[0][0][0].pos+LongSide[L_L-1][0][L_W-1].pos)/2
R_axis2 =  (ShortSide[S_L-1][0][0].pos+ShortSide[S_L-1][S_H-1][S_W-1].pos)/2 - (ShortSide[0][0][0].pos+ShortSide[0][S_H-1][S_W-1].pos)/2
R_axis1 = norm(cross(R_axis2,R_axis3))

axis1_arr = arrow(pos = C_M_pos , axis = R_axis1, color=color.red,shaftwidth = 0.1)
axis2_arr = arrow(pos = C_M_pos , axis = R_axis2, color=color.blue,shaftwidth = 0.1)
axis3_arr = arrow(pos = C_M_pos , axis = R_axis3, color=color.yellow,shaftwidth = 0.1)
dt = 5e-5
t=0
times = 0
p_A1.plot(pos=(t,R_axis1.hat.dot(cal_L(R_axis1,C_M_pos))))
p_A2.plot(pos=(t,R_axis2.hat.dot(cal_L(R_axis2,C_M_pos))))
p_A3.plot(pos=(t,R_axis3.hat.dot(cal_L(R_axis3,C_M_pos))))

print("K =",K)
while t<100:
    rate(4000)
    t+=dt
    times += 1

    for i in range(L_L):
        for j in range(L_H):
            for k in range(L_W):
                LongSide[i][j][k].a = vec(0,0,0)
                
                for delta in neighbor:
                    ni, nj, nk= i+delta[0],j+delta[1],k+delta[2]
                    #print(delta[1])
                    #print(ni,nj,nk)
                    if (ni,nj,nk) not in L_index:continue
                    
                    natural_L = mag(LongSide[i][j][k].opos-LongSide[ni][nj][nk].opos)
                    # print(natural_L)
                    LongSide[i][j][k].a += -K/m*(mag(LongSide[i][j][k].pos-LongSide[ni][nj][nk].pos)-natural_L*a)*norm(LongSide[i][j][k].pos-LongSide[ni][nj][nk].pos)
                #print(LongSide[i][j][k].a)
    for i in range(S_L):
        for j in range(S_H):
            for k in range(S_W):
                ShortSide[i][j][k].a = vec(0,0,0)
                for delta in neighbor:
                    ni, nj, nk= i+delta[0],j+delta[1],k+delta[2] 
                    if (ni,nj,nk) not in S_index:continue
                    natural_L = mag(ShortSide[i][j][k].opos-ShortSide[ni][nj][nk].opos)
                    ShortSide[i][j][k].a += -K/m*(mag(ShortSide[i][j][k].pos-ShortSide[ni][nj][nk].pos)-natural_L*a)*norm(ShortSide[i][j][k].pos-ShortSide[ni][nj][nk].pos)


    #interface
    for Lj in range(L_H):
        for Lk in range(L_W):
            for Sj in range(S_H):
                for Sk in range(S_W):
                    for nat_L in (1,sqrt(2)):
                        if mag(ShortSide[0][Sj][Sk].opos-LongSide[0][Lj][Lk].opos) >0.99*a*nat_L and mag(ShortSide[0][Sj][Sk].opos-LongSide[0][Lj][Lk].opos) <1.01*a*nat_L:
                            # print("hi",t)
                            LongSide[0][Lj][Lk].a += -K/m*(mag(LongSide[0][Lj][Lk].pos-ShortSide[0][Sj][Sk].pos)-nat_L*a)*norm(LongSide[0][Lj][Lk].pos-ShortSide[0][Sj][Sk].pos)
                            ShortSide[0][Sj][Sk].a += K/m*(mag(LongSide[0][Lj][Lk].pos-ShortSide[0][Sj][Sk].pos)-nat_L*a)*norm(LongSide[0][Lj][Lk].pos-ShortSide[0][Sj][Sk].pos)

    #print(LongSide[0][0][0].a)
    C_M_pos = vec(0,0,0)
    for i in range(L_L):
        for j in range(L_H):
            for k in range(L_W):
                LongSide[i][j][k].v += LongSide[i][j][k].a*dt
                LongSide[i][j][k].pos += LongSide[i][j][k].v*dt
                C_M_pos += LongSide[i][j][k].pos/((L_L*L_H*L_W)+(S_L*S_H*S_W))
    for i in range(S_L):
        for j in range(S_H):
            for k in range(S_W):
                ShortSide[i][j][k].v += ShortSide[i][j][k].a*dt
                ShortSide[i][j][k].pos += ShortSide[i][j][k].v*dt
                C_M_pos += ShortSide[i][j][k].pos/((L_L*L_H*L_W)+(S_L*S_H*S_W))
    
    R_axis3 = (LongSide[0][L_H-1][0].pos+LongSide[L_L-1][L_H-1][L_W-1].pos)/2 - (LongSide[0][0][0].pos+LongSide[L_L-1][0][L_W-1].pos)/2
    R_axis2 =  (ShortSide[S_L-1][0][0].pos+ShortSide[S_L-1][S_H-1][S_W-1].pos)/2 - (ShortSide[0][0][0].pos+ShortSide[0][S_H-1][S_W-1].pos)/2
    R_axis1 = norm(cross(R_axis2,R_axis3))

    axis1_arr.pos = C_M_pos
    axis1_arr.axis = R_axis1
    axis2_arr.pos = C_M_pos
    axis2_arr.axis = R_axis2
    axis3_arr.pos = C_M_pos
    axis3_arr.axis = R_axis3
    if(times%1000==0):
        p_A1.plot(pos=(t,R_axis1.hat.dot(cal_L(R_axis1,C_M_pos))))
        p_A2.plot(pos=(t,R_axis2.hat.dot(cal_L(R_axis2,C_M_pos))))
        p_A3.plot(pos=(t,R_axis3.hat.dot(cal_L(R_axis3,C_M_pos))))