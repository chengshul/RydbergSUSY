import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import sys

L = 40 #int(sys.argv[1])
r1 = 0.884 #float(sys.argv[2])
r2 = 0.5
D1 = 1.838 #float(sys.argv[3])
D2 = 1.843 #float(sys.argv[4])
chi = 200 #int(sys.argv[5])
D_l = np.linspace(D1,D2,21)
iD = 13 #int(sys.argv[6])

S = np.fromfile(f"data0/S_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.dat").reshape((len(D_l),2*L-1),order='F')
corr1 = np.fromfile(f"data0/corr1_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.dat").reshape((len(D_l),2,L),order='F')
corr2 = np.fromfile(f"data0/corr2_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.dat").reshape((len(D_l),2,L,2,L),order='F')
corrf = np.fromfile(f"data0/corrf_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.dat").reshape((len(D_l),L-1,4),order='F')

c_z2z2 = corr2[:,0,0,0,:]+corr2[:,1,0,0,:]+corr2[:,0,0,1,:]+corr2[:,1,0,1,:]\
        -(corr1[:,0,0]+corr1[:,1,0])[:,None]*(corr1[:,0,:]+corr1[:,1,:])
c_z2z2 = c_z2z2[:,1:]

c_f = corrf[:,:,0]-corrf[:,:,1]-corrf[:,:,2]+corrf[:,:,3]
c_f = c_f+c_f[:,-1::-1]

x2 = np.sin(np.arange(1,L)*np.pi/L)
eta1 = [-np.polyfit(np.log(x2[2:-2]),np.log(y[2:-2]),1)[0] for y in c_z2z2]
eta2 = [-np.polyfit(np.log(x2[2:-2]),np.log(y[2:-2]),1)[0] for y in c_f]

def func(x,a,b):
    return a*x/3.+b

cc = np.zeros((len(D_l)))

x = np.log(np.sin(np.linspace(np.pi/2/L,np.pi*(2*L-1)/2/L,2*L-1)))
fit0 = 5
fitP = 4
xi = x[fit0:2*L-1-fit0:fitP]

for i_D in range(len(D_l)):
    yi = S[i_D,fit0:2*L-1-fit0:fitP]
    popt,pcov = fit(func,xi,yi)
    cc[i_D] = popt[0]
    print(popt[0])
    yiloc = np.arange(1,x.size+1)[fit0:2*L-1-fit0:fitP]
    #plt.plot(S,'.-')
    #plt.figure(1)
    x1 = np.linspace(0,2*L,100)[1:-1]
    #if popt[0] > 0.1:
    if i_D == iD:
        plt.plot(yiloc/2,yi,'o')#, label=f'$L={L}$')
        plt.plot(x1/2,popt[0]*np.log(np.sin(x1*np.pi/2/L))/3+popt[1],'--',color='grey',label=f"$c={np.round(popt[0],3)}$")
        

np.save(f"data1/cc_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}",cc)
np.save(f"data1/eta_eps_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}",eta1)
np.save(f"data1/eta_psi_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}",eta2)

plt.xlabel("$j$",fontsize=30)
plt.ylabel("$S_j$",fontsize=30)
plt.legend(fontsize=30)
plt.xticks([0,20,40],fontsize=30)
plt.yticks([0,0.5,1],fontsize=30)
plt.tight_layout()
plt.savefig(f"plot/cc_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.pdf",format='pdf')

plt.figure()
plt.plot(D_l,cc,'.-',label="c")
plt.plot(D_l,np.abs(np.sum(corr1[:,0,:]-corr1[:,1,:],axis=1)/2/L),'.-',label="M")
plt.legend()
plt.show()

