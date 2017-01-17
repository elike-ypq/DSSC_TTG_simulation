import math
import matplotlib.pyplot as plt
#arfa=arfa_p-arfa_n
#i=arfa*Ig/Kg #is the dimentionless electric current
#Z=arfa**2/(Kg*Rg) #is the figure of merit of the thermoelectric materials
Z=1/300
T_L=300
#maximizing the power output of the TTG
def Output(n,x_opt,i,T_H=340) :
    #DSSC's operation temperature
    sita=T_H/T_L #is the temperature ratio between the DSSC and the low temperature reservoir.
    #x=m/n #is the ratio of thermoelectric element numbers between the top stag and bottom stage
    d1=sita*Z*T_L
    d2=Z*T_L*(sita+1)
    d3=Z*T_L*(2*sita**2*Z*T_L+sita+3)
    d4=Z*T_L*(sita-1)
    Kg=0.09 #Heat conductivity of a thermoelectric element,Kg(W*k^-1*m^-1)
    P_TTG_m=Kg*T_L*n*(x_opt+1)*((x_opt*sita-1)*i/(x_opt+1)-i**2/(2*Z*T_L)+(sita*x_opt+1)/(x_opt+1)+((x_opt+1)*i**2/(2*Z*T_L)+(x_opt*sita+1))/((x_opt-1)*i-(x_opt+1)))
    f1=(1-sita)*x_opt
    f2=-2*x_opt
    f3=x_opt-1-(1+x_opt)/(Z*T_L)
    f4=(x_opt-1)/(2*Z*T_L)
    f5=-2*sita*x_opt
    f6=x_opt*(x_opt-1)*sita+(x_opt*(1+x_opt))/(Z*T_L)
    f7= -x_opt*(x_opt-1)/(2*Z*T_L)
    yita_TTG_m=(1-(f1+f2*i+f3*i**2+f4*i**3)/(f1+f5*i+f6*i*i+f7*i**3))
    return P_TTG_m,yita_TTG_m