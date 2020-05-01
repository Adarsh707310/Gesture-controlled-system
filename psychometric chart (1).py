#ADARSH JAIN
#B17ME002

#libraries required 
import math as mh
import numpy as np 
import matplotlib.pyplot as plt

#saturation pressure(Pa) at a given temperature(degree C)
def p_s(t):
    return(mh.exp(77.3450+0.0057*(t+273.15)-7235/(t+273.15))/pow(t+273.15,8.2))

#Data given
specific_volume= [.8,.85,.9,.95]
relative_humidity= np.arange(start=0.1, stop=1.1, step=0.1)
DBT= np.arange(start=5, stop=55, step=5)

#Some constants
P=101325  #Atmospheric pressure in Pa 
Ra=287.058 # Gas constant for air in J/Kg*K
p=101.325 #Atmospheric pressure in KPa

x=[]
y=[]

#Plotting constant relative humidity lines
f=1
for i in range(0,len(relative_humidity)):
    phi=relative_humidity[i]
    for j in range(0,len(DBT)):
        w=.622/( (P/(p_s(DBT[j])*phi)) -1)
        x.append(DBT[j])
        y.append(w)
    if f==1:
        plt.plot(x,y,color='green', label='Constant relative humidity (%)')
        f=0
    else:
        plt.plot(x,y,color='green')
    x.clear()
    y.clear()
    
#plotting the line w=(t+5)/1000 which acts as a bound to constant enthalpy lines
#W in KJ/Kgda and t in degree C
DBT2=np.arange(start=5, stop=52, step=.1)
for i in range (0,len(DBT2)):
    w=(DBT2[i]+5)/1000
    x.append(DBT2[i])
    y.append(w)
plt.plot(x,y,color='red')

x.clear()
y.clear()
f=1

#plotting constant enthaply(KJ/Kg(da)) lines
h1=[]
temp=np.arange(start=5, stop=35, step=5)
for i in range(0,len(temp)):
    w=.622/( (P/p_s(temp[i])) -1)
    enth=(1.206*(DBT[i])+2500*w)
    h1.append(enth)
h1=np.array(h1)


x.clear()
y.clear()
for i in range(0,len(h1)):
    hh=h1[i]
    for j in range(0,len(DBT2)):
        jj=len(DBT2)-1-j
        w=(hh-1.0216*DBT2[jj])/2500
        if(w*1000<=5+DBT2[jj]):
            x.append(DBT2[jj])
            y.append(w)
        elif w<.025:
            plt.text(DBT2[jj]-1.8,w+.0001,'%.2f'%hh)
            break
    if f==1:
        plt.plot(x,y,color='red', linestyle='-', label='Constant enthalpy (KJ/Kg(da))')
        f=0
    else: 
        plt.plot(x,y,color='red', linestyle='-')
    x.clear()
    y.clear()

plt.text(38,.003,'10%')
plt.text(37,.0065,'20%')
f=1

#plotting constant specific volume(m^3/Kg(da))
DBT2=np.arange(start=5, stop=52, step=.001)
for i in range(0, len(specific_volume)):
    v=specific_volume[i]
    phi=1
    for j in range(0,len(DBT2)):
        jj=len(DBT2)-j-1
        t=DBT2[jj]+273
        w=.622 *((P*v)/(Ra*t)-1)
        if P-Ra*t/v <p_s(t-273):
            x.append(DBT2[jj])
            y.append(w)
        else:
            w1=min(.622/( (P/(p_s(DBT2[jj])*phi)) -1), .028)
            if v!=.95:
                plt.text(DBT2[jj],w1,'%.2f'%v)
            else:
                plt.text(44,w1,'%.2f'%v)
            break
    if f==1:
        plt.plot(x,y,color='blue', linestyle='--', label='Constant specific volume (m^3/Kg(da))')
        f=0
    else:
        plt.plot(x,y,color='blue', linestyle='--')
    x.clear()
    y.clear()

#Defining the plot parameters
plt.xlim(5, 50)
plt.ylim(0, .03)
plt.title('Psychometric chart')
plt.tick_params(axis='y',which='both',labelleft=1, labelright=1) 
plt.xlabel('Dry Bulb Temperature (degree Celcius)')
plt.ylabel('Specific Humidity (kg/kg(da))')
plt.legend(loc='upper center',bbox_to_anchor=(0.55,-0.2),shadow='True',ncol=2)
plt.show()
