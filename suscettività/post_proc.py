import pandas as pd
import matplotlib.pyplot as plt

data=pd.read_csv('suscettività/susc.dat')

for ispin in range(3,20):
    filtered_data = data[data['L'] == ispin]
    plt.plot(filtered_data['g'], filtered_data['chi'])

plt.show()

filtered_data = data[data['g'] == 1.0]
plt.plot(filtered_data['L'], filtered_data['chi'])
plt.show()