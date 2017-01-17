import simulation_DSSC as sD
import simulation_TTG2 as sT
import simulation_TTG1 as sTr
import math
h=0.28e-3 #W*cm^-2*K^-1
T_H=T=340
QL = h*(T-sD.T_ref)
Qfun = lambda J : sD.G * sD.miu - sD.OutputPower(J)[0] - QL 
N=100
J=[(x+1)/N*sD.J_sc for x in range(N)]
Q_heat_can_be_used=list(map(Qfun, J))
Q_min=min(Q_heat_can_be_used)
J_0=J[Q_heat_can_be_used.index(Q_min)]
#find the proper value of n
OutFunc=sT.Output(T_H)
n0=0
i_max=2*sT.Z*(T_H-330)
n1=Q_min/(sT.Kg*(T_H-330))
eps1=Q_min/1000
while(math.fabs(sT.Q1-Q_min) > eps1) :
    sT.n=1/2*(n0+n1)
    N=100
    i=[(x+1)/N*i_max for x in range(N)]
    P_TTG=[OutFunc(x)[0] for x in i]
    i_max=i[P_TTG.index(max(P_TTG))]
    OutFunc(i_max)
    if(sT.Q1 < Q_min) :
        n0=sT.n
    else :
        n1=sT.n
print(n0,n1)
sT.n=1/2*(n0+n1)       
#find the optimisive value of TTG's current 
i_max=2*sT.Z*(T_H-330)
OutFunc=sT.Output(T_H)
N=100
i=[(x+1)/N*i_max for x in range(N)]
P_TTG=[OutFunc(x)[0] for x in i]
i_max=i[P_TTG.index(max(P_TTG))]
######
sT.Output(T_H)(1e-10)
#Q_max=sT.Q1
n=sT.n
x_opt=sT.x_opt
print(Q_min==sT.Qa_min)
if (n==0) :
    print(n)
    exit()
#####
b=1/x_opt
Q1Func=lambda i : x_opt*n*(i*sT.Kg*T_H-sT.Kg/2*i**2/sT.Z+sT.Kg*(T_H-(i**2/(2*sT.Z)*(1+b)+T_H+b*sT.T_L)/(b*i-i+b+1)))
Q_max=Q1Func(i_max)
def GetTTGCurrent(Q_absorb) :
    if (Q_absorb >= Q_min and Q_absorb <= Q_max) :
        i0=0
        i1=i_max
        eps=i_max/1000
        while((i1-i0) > eps) : 
            if (Q1Func((i0+i1)/2) < Q_absorb) :
                i0=(i0+i1)/2
            else :
                i1=(i0+i1)/2
        return (i1+i0)/2
    else :
        return 0

P_DSSC_s=[sD.OutputPower(x)[0] for x in J]
P_TTG_s=[sTr.Output(n,x_opt,GetTTGCurrent(x))[0] for x in Q_heat_can_be_used]

sT.DrawGraphics(J,P_TTG_s)
