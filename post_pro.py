import pandas as pd
import matplotlib.pyplot as plt
data=pd.read_csv("quenching.dat")
elles=[3,7,8,10,12]
for ell in elles:
    filtered_data=data[data['L']==ell]
    plt.plot(filtered_data['theta'], filtered_data['mL'])
plt.show()