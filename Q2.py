#importing libraries
import numpy as np
import matplotlib.pyplot as plt
import math as mt

#Fresh air
t1=43+273
t1w=27.5+273
#Room conditions
t2=25+273
phi2=.5
ti=t2

#using psychometric charts
w1=16.911*.001 
w2=9.921*.001
wi=w2
h1=84.532
hi=50.407

#mixing conditions. m1/m2=1/4
wm=(w1+4*w2)/5
tm=(t1*w2- t2*w1 -wm*(t1-t2))/(w2-w1) 
#using these, we get wm=11.31 and tm=28.6

#some constants
qs=20
ql=5
r=qs/ql
c1=.0204
c2=50
c=c1/c2
d=r/c
Ra=287.058
P=101325

#phi=1 line
def p_s(t):
    return(mt.exp(77.3450+0.0057*(t+273.15)-7235/(t+273.15))/pow(t+273.15,8.2))
P=101325
DBT= np.arange(start=5, stop=100, step=.1)
x=[]
y=[]
for j in range(0,len(DBT)):
    w=.622/( (P/(p_s(DBT[j]))) -1)
    x.append(DBT[j])
    y.append(w)


#function to find intersection point of phi=1 curve and line
def f(k1,k2):
    xx=[]
    yy=[]
    for j in range(0,len(DBT)):
        wca=k1*(DBT[j]+273)+k2
        xx.append(DBT[j])
        yy.append(wca)
        
    #plt.plot(xx,yy)
    #plt.xlim(5,100)
    yx=np.array(y)
    yy=np.array(yy)
    idx = np.argwhere(np.diff(np.sign(yx - yy))).flatten()
    if len(idx)!=0:
        return DBT[idx[0]]+273
    else:
        return -1

tca=[]
ts=[]
xaxi=[]
qtotal=[]
vflow=[]
#using by-pass factor of the coil
bb=np.arange(start=0, stop=1, step=.1)
for i in range(0,len(bb)):
    b=bb[i]
    k1=1/d
    k2=(wi-b*wm)/(1-b)+ (b*tm-ti)/(d*(1-b))
    #print(k2)
    temp=f(k1,k2)
    #print(temp)
    if temp!=-1:
        #print(temp)
        tca.append(temp-273)
        wca=k1*temp+k2
        ws=b*wm+(1-b)*wca
        ts1=b*tm+(1-b)*temp
        ts.append(ts1-273)
        vspec=Ra*ts1/P*(1+ws/.622)
        cmm=qs/(.0204*(ti-ts1))
        v=qs/(.0204*(ti-ts1))
        ms=cmm/(60*vspec)
        mo=.2*ms
        houtside=mo*(h1-hi)
        qtotal.append(houtside+qs+ql)
        vflow.append(v)
        xaxi.append(b)   

#plotting the graphs        
plt.plot(xaxi,tca, label='Tca')
plt.plot(xaxi,ts, label='Ts')
plt.legend(loc='upper center',bbox_to_anchor=(0.55,-0.2),shadow='True',ncol=2)
plt.xlabel('By-pass factor')
plt.ylabel('Temperature (degree C)')
plt.title('Variation of suply temperature and coil ADP vs bypass factor')
plt.show()

fig, ax = plt.subplots()
ax.plot(xaxi,qtotal)
ax.set_title('Total load(KW) vs bypass factor')

fig, ax = plt.subplots()
ax.plot(xaxi,vflow)
ax.set_title('Volume flow rate(CMM)  vs bypass factor')

