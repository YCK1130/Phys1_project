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
                r = mag(cross(LongSide_pos[i][j][k]-ref_point,rotation_axis.hat))
                Moment_of_Inertia += m*r**2
    for i in range(S_L):
        for j in range(S_H):
            for k in range(S_W):
                r = mag(cross(ShortSide_pos[i][j][k]-ref_point,rotation_axis.hat))
                Moment_of_Inertia += m*r**2
    return Moment_of_Inertia



def cal_L(rotation_axis,ref_point):
    L_AM = vec(0,0,0)
    for i in range(L_L):
        for j in range(L_H):
            for k in range(L_W):
                r = LongSide_pos[i][j][k]-ref_point
                p = m*LongSide_v[i][j][k]
                L_AM += cross(r-r.dot(rotation_axis)*rotation_axis,p-p.dot(rotation_axis)*rotation_axis)
    for i in range(S_L):
        for j in range(S_H):
            for k in range(S_W):
                r = ShortSide_pos[i][j][k]-ref_point
                p = m*ShortSide_v[i][j][k]
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
'''Long Side setting'''
LongSide = []

LongSide_a = []
LongSide_v = []
LongSide_pos = []

L_index = []

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
# print(LongSide_pos)
middle = (LongSide[0][0][0].pos + LongSide[1][L_H-1][1].pos)/2

velocity = 100
perturbation = vec((random()+1e-3)/100,0,0)
LongSide_v[0][0][0] = velocity*cross(middle-LongSide_pos[0][0][0],vec(-1,0,0)) + perturbation
LongSide_v[1][0][0] = velocity*cross(middle-LongSide_pos[1][0][0],vec(-1,0,0))
LongSide_v[0][L_H-1][1] = velocity*cross(middle-LongSide_pos[0][L_H-1][1],vec(-1,0,0))
LongSide_v[1][L_H-1][1] = velocity*cross(middle-LongSide_pos[1][L_H-1][1],vec(-1,0,0)) - perturbation

print("perturbation =",perturbation)


# print(LongSide[L_L-1][L_H-1][L_W-1].a)
'''Short Side setting'''
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

'''init rotating axis'''
init_arr = arrow(pos = middle-vec(1,0,0) , axis = vec(-S_L-2,0,0), color=color.green,shaftwidth = 0.1)

'''counting C.M. , three orthogonal axis'''
C_M_pos = vec(0,0,0)

for i,j,k in L_index:
    C_M_pos += LongSide_pos[i][j][k]/((L_L*L_H*L_W)+(S_L*S_H*S_W)) 

for (i, j, k) in S_index:
    C_M_pos +=ShortSide_pos[i][j][k]/((L_L*L_H*L_W)+(S_L*S_H*S_W))

R_axis3 =  (LongSide_pos[0][L_H-1][0]+LongSide_pos[L_L-1][L_H-1][L_W-1])/2 - (LongSide_pos[0][0][0]+LongSide_pos[L_L-1][0][L_W-1])/2
R_axis2 =  (ShortSide_pos[S_L-1][0][0]+ShortSide_pos[S_L-1][S_H-1][S_W-1])/2 - (ShortSide_pos[0][0][0]+ShortSide_pos[0][S_H-1][S_W-1])/2
R_axis1 = norm(cross(R_axis2,R_axis3))

axis1_arr = arrow(pos = C_M_pos , axis = R_axis1, color=color.red,shaftwidth = 0.1)
axis2_arr = arrow(pos = C_M_pos , axis = R_axis2, color=color.blue,shaftwidth = 0.1)
axis3_arr = arrow(pos = C_M_pos , axis = R_axis3, color=color.yellow,shaftwidth = 0.1)
print("Moment of Inertia axis1 :",cal_I(R_axis1,C_M_pos))
print("Moment of Inertia axis2 :",cal_I(R_axis2,C_M_pos))
print("Moment of Inertia axis3 :",cal_I(R_axis3,C_M_pos))


dt = 5e-5
t=0
prev_t = 0
times = 0  #polt graph per 1000 times 
p_A1.plot(pos=(t,R_axis1.hat.dot(cal_L(R_axis1,C_M_pos))))
p_A2.plot(pos=(t,R_axis2.hat.dot(cal_L(R_axis2,C_M_pos))))
p_A3.plot(pos=(t,R_axis3.hat.dot(cal_L(R_axis3,C_M_pos))))

print("K =",K)
while t<50:
    rate(5000)
    t+=dt
    times += 1
    

    for (i, j, k) in L_index:
        LongSide_a[i][j][k] = vec(0,0,0)
        for delta in neighbor:
            ni, nj, nk= i+delta[0],j+delta[1],k+delta[2]
            if (ni,nj,nk) not in L_index:continue
            natural_L = mag(LongSide[i][j][k].opos-LongSide[ni][nj][nk].opos)
            LongSide_a[i][j][k] += -K/m*(mag(LongSide_pos[i][j][k]-LongSide_pos[ni][nj][nk])-natural_L*a)*norm(LongSide_pos[i][j][k]-LongSide_pos[ni][nj][nk])
    
    for (i, j, k) in S_index:
        ShortSide_a[i][j][k] = vec(0,0,0)
        for delta in neighbor:
            ni, nj, nk= i+delta[0],j+delta[1],k+delta[2] 
            if (ni,nj,nk) not in S_index:continue
            natural_L = mag(ShortSide[i][j][k].opos-ShortSide[ni][nj][nk].opos)
            ShortSide_a[i][j][k] += -K/m*(mag(ShortSide_pos[i][j][k]-ShortSide_pos[ni][nj][nk])-natural_L*a)*norm(ShortSide_pos[i][j][k]-ShortSide_pos[ni][nj][nk])

    #interface
    for Lj in range(L_H):
        for Lk in range(L_W):
            for Sj in range(S_H):
                for Sk in range(S_W):
                    for nat_L in (1,sqrt(2)):
                        if mag(ShortSide[0][Sj][Sk].opos-LongSide[0][Lj][Lk].opos) >0.99*a*nat_L and mag(ShortSide[0][Sj][Sk].opos-LongSide[0][Lj][Lk].opos) <1.01*a*nat_L:
                            # print("hi",t)
                            LongSide_a[0][Lj][Lk] += -K/m*(mag(LongSide_pos[0][Lj][Lk]-ShortSide_pos[0][Sj][Sk])-nat_L*a)*norm(LongSide_pos[0][Lj][Lk]-ShortSide_pos[0][Sj][Sk])
                            ShortSide_a[0][Sj][Sk] += K/m*(mag(LongSide_pos[0][Lj][Lk]-ShortSide_pos[0][Sj][Sk])-nat_L*a)*norm(LongSide_pos[0][Lj][Lk]-ShortSide_pos[0][Sj][Sk])

    #print(LongSide[0][0][0].a)
    for (i, j, k) in L_index:
        LongSide_v[i][j][k] += LongSide_a[i][j][k]*dt
        LongSide_pos[i][j][k] += LongSide_v[i][j][k]*dt
    for (i, j, k) in S_index:
        ShortSide_v[i][j][k] += ShortSide_a[i][j][k]*dt
        ShortSide_pos[i][j][k] += ShortSide_v[i][j][k]*dt

    C_M_pos = vec(0,0,0)
    for (i, j, k) in L_index:
        C_M_pos += LongSide_pos[i][j][k]/((L_L*L_H*L_W)+(S_L*S_H*S_W)) 
    for(i, j, k) in S_index:
        C_M_pos += ShortSide_pos[i][j][k]/((L_L*L_H*L_W)+(S_L*S_H*S_W))
    
    R_axis3 =  (LongSide_pos[0][L_H-1][0]+LongSide_pos[L_L-1][L_H-1][L_W-1])/2 - (LongSide_pos[0][0][0]+LongSide_pos[L_L-1][0][L_W-1])/2
    R_axis2 =  (ShortSide_pos[S_L-1][0][0]+ShortSide_pos[S_L-1][S_H-1][S_W-1])/2 - (ShortSide_pos[0][0][0]+ShortSide_pos[0][S_H-1][S_W-1])/2
    R_axis1 = norm(cross(R_axis2,R_axis3))

    axis1_arr.pos = C_M_pos
    axis1_arr.axis = R_axis1
    axis2_arr.pos = C_M_pos
    axis2_arr.axis = R_axis2
    axis3_arr.pos = C_M_pos
    axis3_arr.axis = R_axis3
    
    now_I = R_axis2.hat.dot(cal_L(R_axis2,C_M_pos))
    if now_I*prev_I<0 and abs(t-prev_t) >= 1e-2 and t>0.5:
        print("Period : {:.5f} sec".format(2*(t-prev_t)))
        prev_t = t
    prev_I = now_I
    if(times%50==0):
        for (i, j, k) in L_index:
            LongSide[i][j][k].pos = LongSide_pos[i][j][k]
        for (i, j, k) in S_index:
            ShortSide[i][j][k].pos = ShortSide_pos[i][j][k]
    
    if(times%1000==0):
        p_A1.plot(pos=(t,R_axis1.hat.dot(cal_L(R_axis1,C_M_pos))))
        p_A2.plot(pos=(t,R_axis2.hat.dot(cal_L(R_axis2,C_M_pos))))
        p_A3.plot(pos=(t,R_axis3.hat.dot(cal_L(R_axis3,C_M_pos))))