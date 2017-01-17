import math
import matplotlib.pyplot as plt
yita_opt=0.970
fai0=1.00e17
n0=1.00e16  #cm^-3
d=10.0e-4   #
fai_b=0.5 #eV
M=4.5 #ideality factor or m=2
#h=6.626e-34 #m^2*kg*s^-1
A_=6.71e2 #A*cm^-2*K^-2
#me=9.10938291e-31
#m_=5.6*me

tao=10e-3
D=5.00e-4
q=1.6e-19
k_B=1.38e-23
T_ref=T_a=300 #K
G=100e-3 #W/cm^2
Pc=0.76
P=0.4
a=4e-4 #cm^-2*s^-1
miu=0.82 
arfa=2568*(1-P)*(P+2.89)
D=a*math.pow(math.fabs(P-Pc),miu)
l=math.sqrt(D*tao)
#print (l)
#D is the electron diffusion parameter,tao is the electron lifetime,fai0 is the optical efficiency of the glass cover,and arfa is the light absorption coefficient of the porous electrode
fai=fai0*yita_opt
#kB is the Boltzmann constant,d is the thin film thickness, Tref is the operating temprature at the reference condition, J is the current density of the DSSC.
#G=100mW/cm^2;
#yita_opt is the optical efficient of the glass cover;
#because of the effeciency of tempreture

J_sc=q*fai*l*arfa/(1-l*l*arfa*arfa)*(-l*arfa+math.tanh(d/l)+l*arfa*math.exp(-d*arfa)/math.cosh(d/l))
#find the max current density J_sc
#N=100
#J=[x*J_sc/N for x in range(N+1)]
def OutputPower(J,T=340) : 
    V=k_B*T_ref*M/q*math.log1p(l*(J_sc-J)/(q*D*n0*math.tanh(d/l)))-k_B*T_ref/q*math.log1p(J/(A_*T_ref**2*math.exp(-q*fai_b/(k_B*T_ref))))
    P_D_ref=J*V  #W/cm^2
    yita_D_ref=P_D_ref/(G)
    ###################################################
    #the cm^-2 shuold be here, so it has been added to these two parameter.
    lambda0=0.00506e-7  #W*K^-1^cm^-2
    beta=0.0278e-4    #K^-1*cm^-2
    ###################################################
    P_D=P_D_ref-lambda0*(T-T_ref)
    yita_D=yita_D_ref*(1-beta*(T-T_ref))
    return P_D,yita_D
#print(J,P_D_ref,P_D)
def WriteToFlie(J,P,yita) :
    f1 = open("OutputDSSC.txt", "w+")
    f1.write("J\tP_DSSC\tyita_DSSC\n")
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
    plt.title('the output charactoristics of DSSC')
    plt.legend()
    plt.show()