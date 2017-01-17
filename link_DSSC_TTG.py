import simulation_DSSC as sD
import simulation_TTG as sT
import simulation_TTG1 as sTr
h=0.28e-3 #W*cm^-2*K^-1
T_H=T=340
QL = h*(T-sD.T_ref)
Qfun = lambda J : sD.G * sD.miu - sD.OutputPower(J)[0] - QL 
N=100
J=[(x+1)/N*sD.J_sc for x in range(N)]
Q_heat_can_be_used=list(map(Qfun, J))
Q_min=min(Q_heat_can_be_used)
J_0=J[Q_heat_can_be_used.index(Q_min)]
#find the optimisive value of TTG's current 
i_max=2*sT.Z*(T_H-330)
OutFunc=sT.Output(Q_min,T_H)
N=100
i=[(x+1)/N*i_max for x in range(N)]
P_DSSC=[OutFunc(x)[0] for x in i]
i_max=i[P_DSSC.index(max(P_DSSC))]
######
sT.Output(Q_min,T_H)(i_max)
Q_max=sT.Q1
n=sT.n
x_opt=sT.x_opt
if (n==1) :
    print(n)
    exit()
#####
b=1/x_opt
Q1Func=lambda i : x_opt*n*(i*sT.Kg*T_H-sT.Kg/2*i**2/sT.Z+sT.Kg*(T_H-(i**2/(2*sT.Z)*(1+b)+T_H+b*sT.T_L)/(b*i-i+b+1)))
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


P_TTG=[sTr.Output(n,x_opt,GetTTGCurrent(x))[0] for x in Q_heat_can_be_used]
#sT.DrawGraphics(J,P_TTG)
