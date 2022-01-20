from vpython import *
# oscillation = graph(width = 450, align = 'right') 
# omega1 = gcurve(graph = oscillation, color=color.blue, width=4)
# omega2 = gcurve(graph = oscillation, color=color.red, width=4)
# omega3 = gcurve(graph = oscillation, color=color.green, width=4)

M=0.1 
l=4
m=3
s=2

v=100
u=0.01

I1=1.360 #最小轉動慣量
I2=1.840 #中間轉動慣量
I3=2.960 #最大轉動慣量

#o2=m*M*sqrt((l-1)**2+(s-1)**2)*v/I2
#o1=s*M*sqrt((l-1)**2+(m-1)**2)*u/I1
#o3=0
#o1=8.350046525398392e-05/I1
#o2=-15.0/I2
#o3=-0.00025050139576234274/I3
L3,L2,L1 =  -1.2942961240442396e-06,14.14213562373095,0.0018305960314603615
o1=L1/I1 #放角動量/角動慣量
o2=L2/I2    
o3=L3/I3 

O1=0
O2=0
O3=0
t=0
dt=0.0001
T0=0
while True:
    rate(10000)
    t+=dt
    oo1=O1
    oo2=O2
    oo3=O3
    O1=o1
    O2=o2
    O3=o3

    temp=o2
    
    o1+=(I2-I3)/I1*O2*O3*dt
    o2+=(I3-I1)/I2*O3*O1*dt
    o3+=(I1-I2)/I3*O1*O2*dt

   ## omega1.plot(pos=(t,o1))
   ## omega2.plot(pos=(t,o2))
   ## omega3.plot(pos=(t,o3))
    if o2*temp<0:
        print(t-T0)
        T0=t
    if abs(O3)<abs(o3) and abs(O3)<abs(oo3):
        print("o3=",o3)
    if abs(O1)<abs(o1) and abs(O1)<abs(oo1):
        print("o1=",o1)