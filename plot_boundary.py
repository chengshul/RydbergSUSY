import numpy as np
import matplotlib.pyplot as plt
import sys

L = int(sys.argv[1])
r2 = 0.5
chi = 200

para_l = np.array([[0.85,2.3,2.5],[0.86,2.1,2.3],[0.87,2.0,2.2],
                   [0.88,1.898,1.903],[0.881,1.883,1.888],[0.882,1.868,1.873],
                   [0.883,1.853,1.858],[0.884,1.838,1.843],[0.885,1.823,1.828],
                   [0.886,1.808,1.813],[0.887,1.793,1.798],[0.888,1.778,1.783],
                   [0.889,1.755,1.775],[0.89,1.74,1.76],[0.892,1.71,1.73],
                   [0.894,1.68,1.7],[0.896,1.65,1.67],[0.898,1.62,1.64],
                   [0.9,1.55,1.65],[0.91,1.4,1.6],[0.92,1.25,1.45],
                   [0.93,1.1,1.3],[0.94,1.0,1.2],[0.95,0.8,1.0]])
r1_l = para_l[:,0]
cc_l = []
cc_fit = []
D_fit = []
eta1_l = []
eta2_l = []

for para in para_l:
    r1, D1, D2 = para
    D = np.linspace(D1, D2, 21)
    cc = np.load(f"data1/cc_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.npy")
    cc_l.append(cc)
    #plt.plot(D, cc, '.-',label=f"$r_1={r1}$")
    
    if r1 > 0.8835 and L == 40:
        eta1_l.append(np.load(f"data1/eta_eps_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.npy"))
        eta2_l.append(np.load(f"data1/eta_psi_2a_L_{L}_r1_{r1}_r2_{r2}_D1_{D1}_D2_{D2}_chi_{chi}.npy"))

cc_l = np.array(cc_l)
eta1_l = np.array(eta1_l)
eta2_l = np.array(eta2_l)

np.save(f"data1/cc_all_L_{L}",[r1_l,np.max(cc_l,axis=1)])

plt.figure()
plt.plot(r1_l,np.max(cc_l,axis=1),'.-')
plt.plot(r1_l[[0,-1]],[0.5,0.5],'--',color='grey')
plt.plot(r1_l[[0,-1]],[0.7,0.7],'--',color='grey')
plt.xlabel("$r_1$")
plt.ylabel("$c$")

if L == 40:
    plt.figure()
    i_c = np.argmax(cc_l[7:],axis=1)

    plt.plot(r1_l[[7,-1]],[0.4,0.4],'--',color='grey')
    plt.plot(r1_l[[7,-1]],[1.4,1.4],'--',color='grey')

    plt.plot(r1_l[7:],[eta1_l[i,i_c[i]] for i in range(17)],
         color='C0',label="$\epsilon$")
    plt.plot(r1_l[7:],[eta2_l[i,i_c[i]] for i in range(17)],
         color='C1',label="$\psi$")
    plt.plot(0.8838,0.416,'o',color='C0')
    plt.plot(0.8838,1.424,'o',color='C1')
    plt.legend(fontsize=20)
    plt.xticks([0.89,0.91,0.93,0.95],fontsize=26)
    plt.yticks([0.5,1,1.5,2],fontsize=26)
    plt.xlabel("$r_1$",fontsize=30)
    plt.ylabel("$\eta$",fontsize=30)
    plt.tight_layout()
    plt.savefig("plot/eta1.pdf",format='pdf')
plt.show()
