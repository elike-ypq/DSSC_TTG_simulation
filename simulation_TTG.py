import math
import matplotlib.pyplot as plt
#arfa=arfa_p-arfa_n
#i=arfa*Ig/Kg #is the dimentionless electric current
#Z=arfa**2/(Kg*Rg) #is the figure of merit of the thermoelectric materials
Z=1/300
T_L=300
Q1=0
n=1
x_opt=1
Kg=0.09 #Heat conductivity of a thermoelectric element,Kg(W*k^-1*m^-1)
#maximizing the power output of the TTG
def Output(Q,T_H) :
    i_max=Z*(T_H-T_L)
    #DSSC's operation temperature
    sita=T_H/T_L #is the temperature ratio between the DSSC and the low temperature reservoir.
    #x=m/n #is the ratio of thermoelectric element numbers between the top stag and bottom stage
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
    #N=100
    #i=[(x+1)*j_max/N for x in range(N)]
    def PowerYita(i) :
        global x_opt
        x_opt=(i**3-(2*d1-1)*i*i-2*d1*i+2*math.sqrt(i*(-i**4+2*d1*i**3-d2*i*i+d3*i+2*d1*d4)))/(i*(i-1)*(i-2*d1))
        b=1/x_opt
        Tm=(i**2/(2*Z)*(1+b)+T_H+b*T_L)/(b*i-i+b+1)
        #j_max=list(map(lambda Tm : Z*(T_H-Tm),Tm))
        
        #n=10 #Number of TEGs in the second stage, 
        global n
        n=Q/(x_opt*Kg*(T_H-Tm))
        global Q1
        Q1=x_opt*n*(i*Kg*T_H-Kg/2*i**2/Z+Kg*(T_H-Tm))
        #print(n)
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
    return PowerYita
def WriteToFlie(J,P,yita) :
    f1 = open("OutputTTG.txt", "w+")
    f1.write("J\tP_TTG\tyita_TTG\n")
    for i in range(0, len(J)):
        f1.write(str(J[i]) + "\t" + str(P[i]) +"\t"+ str(yita[i]) + "\n")
    f1.close()
def DrawGraphics(x,y) :
    plt.figure(figsize=(8,4))
#    plt.plot(J, P_D_ref, 'b*')#,label="$sin(x)$",color="red",linewidth=2)
    plt.plot(x ,y , 'r')
    plt.xlabel("current density (A/cm^2)")
    plt.ylabel("y")
    plt.ylim(min(y)*0.99, max(y)*1.01)
    plt.title('the output charactoristics of TTG')
    plt.legend()
    plt.show()
