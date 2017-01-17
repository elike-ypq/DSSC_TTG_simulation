import simulation_TTG as sT
import link_DSSC_TTG as link
T_H=340
i_max=2*sT.Z*(T_H-330)
OutFunc=sT.Output(link.Q_min,T_H)
N=100
i=[(x+1)/N*i_max for x in range(N)]
P=[OutFunc(x)[0] for x in i]
print(max(P))
print(i_max,i[P.index(max(P))])
sT.DrawGraphics(i,P)