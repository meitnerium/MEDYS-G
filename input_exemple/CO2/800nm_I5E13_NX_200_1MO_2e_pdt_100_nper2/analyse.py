import numpy as np
import matplotlib.pyplot as plt
champ = np.genfromtxt("champ.dat")
plt.plot(champ[:,0],champ[:,1])
plt.save("champ.dat")
