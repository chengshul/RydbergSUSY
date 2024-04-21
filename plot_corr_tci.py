import numpy as np
import matplotlib.pyplot as plt
import sys

L = 40 #int(sys.argv[1])
r1 = 0.8838 #float(sys.argv[2])
r2 = 0.5
D = 1.8442 #float(sys.argv[3])
chi = 400 #int(sys.argv[4])

corr1 = np.fromfile(f"data0/corr1_2a_L_{L}_r1_{r1}_r2_{r2}_D_{D}_chi_{chi}.dat").reshape((2,L),order='F')
corr2 = np.fromfile(f"data0/corr2_2a_L_{L}_r1_{r1}_r2_{r2}_D_{D}_chi_{chi}.dat").reshape((2,L,2,L),order='F')
corrf = np.fromfile(f"data0/corrf_2a_L_{L}_r1_{r1}_r2_{r2}_D_{D}_chi_{chi}.dat").reshape((L-1,4),order='F')

c_z2z2 = corr2[0,0,0,:]+corr2[1,0,0,:]+corr2[0,0,1,:]+corr2[1,0,1,:]\
        -(corr1[0,0]+corr1[1,0])*(corr1[0,:]+corr1[1,:])
c_f = corrf[:,0]-corrf[:,1]-corrf[:,2]+corrf[:,3]

c_z2z2 = c_z2z2[1:]
c_f = c_f+c_f[-1::-1]

x = np.sin(np.arange(1,L)*np.pi/L)
p2 = np.poly1d(np.polyfit(np.log(x[2:-2]),np.log(c_z2z2[2:-2]),1))
p3 = np.poly1d(np.polyfit(np.log(x[2:-2]),np.log(c_f[2:-2]),1))

x1 = np.sin(np.arange(1,L)*np.pi/L)[:L//2]

plt.plot(x,c_z2z2,'o',label="$\epsilon$")
plt.plot(x1[2:],np.exp(p2(np.log(x1)))[2:],
         '--',color='grey',label=f"$\eta={-np.round(p2[1],3)}$")

plt.loglog(x,c_f,'o',label="$\psi$")
plt.plot(x1[2:],np.exp(p3(np.log(x1)))[2:],
         '--',color='grey',label=f"$\eta={-np.round(p3[1],3)}$")

plt.xlabel("$\sin(\pi l/L)$",fontsize=30)
plt.ylabel(r"$\langle O_0O_l\rangle_c$",fontsize=30)
plt.legend(fontsize=20)
plt.xticks([0.2,0.3,0.6,1],[0.2,0.3,0.6,1],fontsize=26)
plt.yticks([1e-2,0.1,1],fontsize=26)
ax = plt.gca()
ax.minorticks_off()
plt.tight_layout()
plt.savefig(f"plot/corr_L_{L}_r1_{r1}_D_{D}.pdf",format='pdf')

plt.show()
