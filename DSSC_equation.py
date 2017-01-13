import math
import matplotlib.pyplot as plt
yita_opt=0.970
fai0=1.00e17
arfa=5.00e3
n0=1.00e16
d=10.0e-4
M=4.5
tao=10e-3
D=5.00e-4
q=1.6e-19
k_B=1.38e-23
T_ref=T_a=300
G=100e-3
A_D=100
l=math.sqrt(D*tao)
print (l)
#D is the electron diffusion parameter,tao is the electron lifetime,fai0 is the optical efficiency of the glass cover,and arfa is the light absorption coefficient of the porous electrode
fai=fai0*yita_opt
#kB is the Boltzmann constant,d is the thin film thickness, Tref is the operating temprature at the reference condition, J is the current density of the DSSC.
#G=100mW/cm^2;
#yita_opt is the optical efficient of the glass cover;
#A_D is the effective area of the DSSC under illumination. G is the solar spectrum intensity;
#because of the effeciency of tempreture 

J_sc=q*fai*l*arfa/(1-l*l*arfa*arfa)*(-l*arfa+math.tanh(d/l)+l*arfa*math.exp(-d*arfa)/math.cosh(d/l))
V_func=lambda J : k_B*T_ref*M/q*math.log1p(l*(J_sc-J)/(q*D*n0*math.tanh(d/l)))
#find the max current density J_sc
N=100
J=[x*J_sc/N for x in range(N)]
V=list(map(V_func,J))
P_D_ref=list(map(lambda J,V : J*A_D*V,J,V))
yita_D_ref=list(map(lambda P_D_ref : P_D_ref/(G*A_D),P_D_ref))
T=340
lambda0=0.00506e-3
beta=0.0278
P_D=list(map(lambda P_D_ref : P_D_ref-lambda0*(T-T_ref),P_D_ref))
yita_D=list(map(lambda yita_D_ref : yita_D_ref*(1-beta*(T-T_ref)),yita_D_ref))
print('Do you want to print the curve?y/n')
select1=input()
if select1=='y' :
    plt.figure(figsize=(8,4))
#    plt.plot(J, P_D_ref, 'b*')#,label="$sin(x)$",color="red",linewidth=2)
    plt.plot(J, P_D_ref, 'r')
    plt.plot(J,P_D)
    plt.xlabel("current density (A/cm^2)")
    plt.ylabel("output power (W)")
    plt.ylim(0, max(P_D_ref)*1.1)
    plt.title('the output charactoristics of DSSC')
    plt.legend()
    plt.show()
select2=input()
print(yita_D)
if select2=='y' :
    plt.figure(figsize=(8,4))
#    plt.plot(J, yita_D_ref, 'b*')#,label="$sin(x)$",color="red",linewidth=2)
    plt.plot(J, yita_D_ref)
    plt.plot(J,yita_D)
    plt.xlabel("current density (A/cm^2)")
    plt.ylabel("output effeciency")
    plt.ylim(0, max(yita_D_ref)*1.1)
    plt.title('the output charactoristics of DSSC')
    plt.legend()
    plt.show()
#when lambda0=0.00506mWK^-1,it fits well experimental data
#beta=0.0278K^-1

