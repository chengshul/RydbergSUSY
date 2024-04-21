import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse

r1_l = [0.85,0.86,0.87,0.88,
        0.881,0.882,0.883,0.884,0.885,
        0.886,0.887,0.888,0.889,0.89,
        0.892,0.894,0.896,0.898,0.9,
        0.91,0.92,0.93,0.94,0.95]
D_l = -np.array([2.42,2.24,2.07,1.9,
       1.887,1.872,1.856,1.841,1.826,
       1.811,1.796,1.781,1.766,1.751,
       1.721,1.692,1.662,1.633,1.605,
       1.46,1.32,1.18,1.05,0.91])

fig, ax = plt.subplots()

x = -np.array([1.3,1.4,1.5,1.6])
y = [0.86, 0.857, 0.87, 0.867]

for i in range(4):
    for j in range(4):
        plt.plot(x[i],y[j],'o',color=f'C{i%2}')
        if j in [1,2]:
            plt.plot(x[i],y[j],'o',color=f'C{i%2}',mfc='None',ms=12)

x = -np.array([2.6,2.7,2.8,2.9])
y = [0.86, 0.857]

for i in range(4):
    for j in range(2):
        plt.plot(x[i],y[j],'o',color=f'C{i%2}')

x = -np.array([1.8,1.9,2.0,2.1])
y = [0.904, 0.901]
theta = [np.linspace(-0.7,3.8,201),np.linspace(0.7,-3.8,201)]
rx = 0.025
ry = 0.00125

for i in range(4):
    for j in range(2):
        plt.plot(x[i],y[j],'o',color=f'C{i%2}')
        plt.plot(x[i]+rx*np.cos(theta[j]),y[j]+ry*np.sin(theta[j]),
                 linestyle=(0,(1,1)),color=f'C{i%2}')
    plt.plot(x[i]+rx*np.cos([4,0.65]),[y[0]+ry*np.sin(4),y[1]+ry*np.sin(0.65)],
             linestyle=(0,(1,1)),color=f'C{i%2}')
    plt.plot(x[i]+rx*np.cos([1,4]),[y[0]+ry*np.sin(-1),y[1]+ry*np.sin(-4)],
             linestyle=(0,(1,1)),color=f'C{i%2}')

plt.text(-1.9,0.862,"Ferro",color='k',size=20)
plt.text(-2.9,0.885,"Disordered",color='k',size=20)
plt.text(-2.05,0.882,"TCI",color='k',size=20)

i_c = 7
print(r1_l[i_c])

plt.plot(D_l[:i_c+1],r1_l[:i_c+1],'--',lw=3,color='C0')
plt.plot(D_l[i_c:],r1_l[i_c:],'-',lw=3,color='C0')
plt.plot(D_l[i_c],r1_l[i_c],'o',markersize=10,color='C0')
plt.xlim([-3,-1.2])
plt.ylim([0.85,0.91])
plt.xlabel("$\Delta$",fontsize=20)
plt.ylabel("$r_1$",fontsize=20)
plt.xticks([-1.3,-1.7,-2.1,-2.5,-2.9],fontsize=18)
plt.yticks([0.86,0.88,0.9],fontsize=18)
plt.tight_layout()
plt.savefig("plot/pd.pdf",format='pdf')

plt.show()
