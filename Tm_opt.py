import math
import matplotlib.pyplot as plt
#arfa=arfa_p-arfa_n
#i=arfa*Ig/Kg #is the dimentionless electric current
#Z=arfa**2/(Kg*Rg) #is the figure of merit of the thermoelectric materials
Z=1/300
T_L=300
T_H=340 #DSSC's operation temperature
sita=T_H/T_L #is the temperature ratio between the DSSC and the low temperature reservoir.
#x=m/n #is the ratio of thermoelectric element numbers between the top stag and bottom stage
#maximizing the power output of the TTG
i_max=Z*(T_H-T_L)

d1=sita*Z*T_L
d2=Z*T_L*(sita+1)
d3=Z*T_L*(2*sita**2*Z*T_L+sita+3)
d4=Z*T_L*(sita-1)
eps=i_max/1000
j_max=eps
j=i_max
while math.fabs(j_max-j) > eps :
    j=j_max
    #print(j_max)
    x_opt=(j**3-(2*d1-1)*j*j-2*d1*j+2*math.sqrt(j*(-j**4+2*d1*j**3-d2*j*j+d3*j+2*d1*d4)))/(j*(j-1)*(j-2*d1))
    b=1/x_opt
    Tm=(j**2/(2*Z)*(1+b)+T_H+b*T_L)/(b*j-j+b+1)
    j_max=Z*(T_H-Tm)
print(Tm)
N=100
i=[(x+1)*j_max/N for x in range(N)]
x_opt=list(map(lambda i : (i**3-(2*d1-1)*i*i-2*d1*i+2*math.sqrt(i*(-i**4+2*d1*i**3-d2*i*i+d3*i+2*d1*d4)))/(i*(i-1)*(i-2*d1)),i))
b=list(map(lambda x_opt : 1/x_opt,x_opt))
Tm=list(map(lambda i,b : (i**2/(2*Z)*(1+b)+T_H+b*T_L)/(b*i-i+b+1),i,b))
#j_max=list(map(lambda Tm : Z*(T_H-Tm),Tm))
Kg=0.09 #Heat conductivity of a thermoelectric element,Kg(W*k^-1*m^-1)
n=10 #Number of TEGs in the second stage, n

P_TTG_m=list(map(lambda i,x_opt :Kg*T_L*n*(x_opt+1)*((x_opt*sita-1)*i/(x_opt+1)-i**2/(2*Z*T_L)+(sita*x_opt+1)/(x_opt+1)+((x_opt+1)*i**2/(2*Z*T_L)+(x_opt*sita+1))/((x_opt-1)*i-(x_opt+1))),i,x_opt))
f1=list(map(lambda x_opt : (1-sita)*x_opt,x_opt))
f2=list(map(lambda x_opt : -2*x_opt,x_opt))
f3=list(map(lambda x_opt : x_opt-1-(1+x_opt)/(Z*T_L),x_opt))
f4=list(map(lambda x_opt : (x_opt-1)/(2*Z*T_L),x_opt))
f5=list(map(lambda x_opt : -2*sita*x_opt,x_opt))
f6=list(map(lambda x_opt : x_opt*(x_opt-1)*sita+(x_opt*(1+x_opt))/(Z*T_L),x_opt))
f7=list(map(lambda x_opt : -x_opt*(x_opt-1)/(2*Z*T_L),x_opt))
yita_TTG_m=list(map(lambda f1,f2,f3,f4,f5,f6,f7,i : (1-(f1+f2*i+f3*i**2+f4*i**3)/(f1+f5*i+f6*i*i+f7*i**3)),f1,f2,f3,f4,f5,f6,f7,i))
select1=input()	
if select1=='y' :
    plt.figure(figsize=(8,4))
#    plt.plot(J, P_D_ref, 'b*')#,label="$sin(x)$",color="red",linewidth=2)
    plt.plot(i,Tm,'r')
#    plt.plot(J,P_D)
    plt.xlabel("current density (i)")
    plt.ylabel("middle temperature(K)")
    plt.ylim(min(Tm), max(Tm))
    plt.title('TTG')
    plt.legend()
    plt.show()
select2=input()	
if select2=='y' :
    plt.figure(figsize=(8,4))
#    plt.plot(J, P_D_ref, 'b*')#,label="$sin(x)$",color="red",linewidth=2)
    plt.plot(i,P_TTG_m,'r')
#    plt.plot(J,P_D)
    plt.xlabel("current density (i)")
    plt.ylabel("TTG's output power(W*m^-2)")
    plt.ylim(min(P_TTG_m), max(P_TTG_m))
    plt.title('TTG')
    plt.legend()
    plt.show()
select3=input()	
if select3=='y' :
    plt.figure(figsize=(8,4))
#    plt.plot(J, P_D_ref, 'b*')#,label="$sin(x)$",color="red",linewidth=2)
    plt.plot(i,yita_TTG_m,'r')
#    plt.plot(J,P_D)
    plt.xlabel("current density (i)")
    plt.ylabel("output effeciency")
    plt.ylim(min(yita_TTG_m), max(yita_TTG_m))
    plt.title('TTG')
    plt.legend()
    plt.show()
