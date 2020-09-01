import matplotlib.pyplot as plt
import numpy as np

bms = np.array([5, 4.5, 4, 3.75, 3.82])
bmc = 3.68423

hs = np.sqrt(bms/bmc) - 1

taus = np.array([18.5, 29.85, 54, 82, 66])

logh = np.log(hs)
logt = np.log(taus)




fit = np.polyfit(logh, logt, deg=1)
p = np.poly1d(fit)


plt.title(f"k = {fit[0]:.3}", fontsize=16)
plt.grid()
hplot = np.array([-4.7, -1.7])
plt.plot(hplot, p(hplot))
plt.scatter(logh, logt)

#plt.yscale("log")
#plt.xscale("log")

plt.ylabel(r"$\log{\ \tau}$", fontsize=16)
plt.xlabel(r"$\log{\ h},\  h = \sqrt{\frac{Bm}{Bm_c}} - 1$", fontsize=16)