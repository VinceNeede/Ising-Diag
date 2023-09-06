import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


data=pd.read_csv('critic g/g_max.dat')
plt.plot(data['L'], data['g_max'])
plt.show()