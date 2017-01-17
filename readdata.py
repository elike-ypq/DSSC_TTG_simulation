#can read data from file
result=[]
i=1
f=open('OutputPower.txt','r')
for line in f:
    if i == 1 :
	    i += 1
    else :
        result.append(list(map(float,line.split('\t'))))
#print(result)
J=[x[0] for x in result]
P_DSSC=[x[2] for x in result]
#print(J,P_DSSC)
f.close()
