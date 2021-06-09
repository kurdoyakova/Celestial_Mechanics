from scipy.optimize import fsolve
import scipy.optimize as opt
from math import atan, tan, asin, cos, sin, atan2,sqrt, acos, pi
import numpy as np

delta=[0,6.88010/57.295779513, 8.88970/57.295779513, 10.61325/57.295779513]
alfa=[0,11.19621*15/57.295779513, 10.94933*15/57.295779513, 10.71958*15/57.295779513]

x=[0,0.82609444, 0.96646769, 0.99096251]
y=[0,-0.49685422, -0.20695689, 0.10783609]
z=[0,-0.21538952, -0.08972156, 0.04674478]
t=[0,0,20,40]
eps=23.436742/57.295779513
kappa=0.017202
kappa2=kappa**2
l=[0]
m=[0]
nu=[0]

tau1=20
tau2=20
tau=40

n01=tau2/tau
n02=tau1/tau
c1=kappa2*tau1*tau2*(1+n01)/6
c2=kappa2*tau1*tau2*(1+n02)/6
print('c1=',c1)
print('c2=',c2)

for i in range(1,4):
    l.append(cos(delta[i])*cos(alfa[i]))
    m.append(cos(delta[i])*sin(alfa[i]))
    nu.append(sin(delta[i]))
print('landa=',l)
print('m=',m)
print('nu=',nu)

l13=m[1]*nu[3]-nu[1]*m[3]
m13=nu[1]*l[3]-nu[3]*l[1]
nu13=l[1]*m[3]-l[3]*m[1]
print('lamda13=',l13)
print('m13=',m13)
print('nu13=',nu13)

u=[0]
for i in range(1,4):
    u.append(x[i]*l13+y[i]*m13+z[i]*nu13)

dd=[[l[2],l[1],l[3]],[m[2],m[1],m[3]],[nu[2],nu[1],nu[3]]]
d=np.linalg.det(dd)

p=(u[2]-n01*u[1]-n02*u[3])/d
q=(c1*u[1]+c2*u[3])/d
cc=-(l[2]*x[2]+m[2]*y[2]+nu[2]*z[2])

r2=(x[2])**2+(y[2])**2+(z[2])**2
print('D=',d)
print('P=',p)
print('Q=',q)
print('C=',cc)
print('R^2=',r2)

def f(variables) :
    (ro0,r) = variables

    first_eq = -ro0+p-q*r**(-3)
    second_eq = -(r**2)+ro0**2+2*cc*ro0+r2
    return [first_eq, second_eq]

solution = opt.fsolve(f, (0.1,10) )
r=solution[1]
ro0=solution[0]
print('r=',r,'ro=',ro0)

n1=n01+c1*r**(-3)
n2=n02+c2*r**(-3)

def f2(variables) :
    (ro1,ro2, ro3) = variables

    first_eq = ro1*n1*l[1]-ro2*l[2]+ro3*n2*l[3]-n1*x[1]+x[2]-n2*x[3]
    second_eq = ro1*n1*m[1]-ro2*m[2]+ro3*n2*m[3]-n1*y[1]+y[2]-n2*y[3]
    third_eq = ro1*n1*nu[1]-ro2*nu[2]+ro3*n2*nu[3]-n1*z[1]+z[2]-n2*z[3]
    return [first_eq, second_eq, third_eq]

solution2 = opt.fsolve(f2, (1,1,1) )

ro=[0]
xx=[0]
yy=[0]
zz=[0]
rrx=[0]
rry=[0]
rrz=[0]
for i in range(1,4):
    ro.append(solution2[i-1])
    xx.append(ro[i]*l[i]-x[i])
    yy.append(ro[i]*m[i]-y[i])
    zz.append(ro[i]*nu[i]-z[i])
print('ro=',ro)
print('x=',xx)
print('y=',yy)
print('x=',zz)

rrr1=[xx[1],yy[1],zz[1]]
rrr2=[xx[2],yy[2],zz[2]]
rrr3=[xx[3],yy[3],zz[3]]

def vec(x,y):
    return [x[1]*y[2]-x[2]*y[1], x[0]*y[2]-x[2]*y[0], x[0]*y[1]-x[1]*y[0]]

c13=(np.linalg.norm(vec(rrr1,rrr2)))/(np.linalg.norm(vec(rrr1,rrr3)))
c11=(np.linalg.norm(vec(rrr2,rrr3)))/(np.linalg.norm(vec(rrr1,rrr3)))
print('c1=',c13)
print('c3=',c11)

r1=np.linalg.norm(rrr1)
r22=np.linalg.norm(rrr2)
r3=np.linalg.norm(rrr3)

pp=(c11*r1+c13*r3-r22)/(c11+c13-1)
print('p=',pp)

f=atan2(np.linalg.norm(vec(rrr1,rrr3)),np.dot(rrr1,rrr3))
#f=f2*2
print('2f=',f)

e=sqrt((((pp/r1-1)*cos(f)-(pp/r3-1))/sin(f))**2+(pp/r1-1)**2)
a=pp/(1-e*e)

teta1=atan2(((pp/r1-1)*cos(f)-(pp/r3-1))/(sin(f)),(pp/r1-1))
teta3=atan2((-(pp/r3-1)*cos(f)+(pp/r1-1))/(sin(f)),(pp/r3-1))

xh=[0]
yh=[0]
zh=[0]

for i in range(1,4):
    xh.append(xx[i])
    yh.append(zz[i]*sin(eps)+yy[i]*cos(eps))
    zh.append(zz[i]*cos(eps)-yy[i]*sin(eps))

lx=yh[1]*zh[3]-yh[3]*zh[1]
ly=zh[1]*xh[3]-zh[3]*xh[1]
lz=xh[1]*yh[3]-yh[1]*xh[3]

om=atan2(lx,-ly)
i=atan2((lx/sin(om)),lz)
uu=atan2(zh[1]/sin(i),xh[1]*cos(om)+yh[1]*sin(om))
g=uu-teta1

print('a=',a,'e=',e,'i=',i*57.295779513,'OMEGA=',om*57.295779513,'teta1=',teta1*57.295779513,'teta3=',teta3*57.295779513,'u=',uu,'g=',(g+2*pi)*57.295779513 ,)
