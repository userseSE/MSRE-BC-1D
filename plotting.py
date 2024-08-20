import matplotlib.pyplot as plt
import numpy as np
import os

# Change directory to "dataset"
os.chdir('dataset')

def plot_csv(filename, title, xlabel='', ylabel=''):
    data = np.loadtxt(filename, delimiter=',')
    plt.figure(figsize=(14, 6))
    if len(data.shape) == 1:
        plt.plot(data)
    else:
        print(title)
        print(data.shape)
        for i in range(data.shape[1]):
            plt.plot(data[:,i], label=f'Column {i+1}')
        plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

plot_csv('rho_matrix.csv', 'Reactivity')
plot_csv('phi_middle_matrix.csv', 'Neutron Flux in the Middle')
plot_csv('ci_middle_matrix.csv', 'Delayed Neutron Precursors in the Middle')
plot_csv('temperature_fuel_middle_matrix.csv', 'Temperature in the Core Middle')

plot_csv('phi.csv', 'Neutron Flux')
plot_csv('ci.csv', 'Delayed Neutron Precursors')
plot_csv('temperature_core.csv', 'Temperature in the core')
# plot_csv('temperature_graphite.csv', 'Temperature in the Graphite')
plot_csv('temperature_HX1.csv', 'Temperature in Heat Exchanger 1')
# plot_csv('Tss_HX1.csv', 'Temperature in Heat Exchanger 1 (Tss_HX1)')
plot_csv('temperature_HX2.csv', 'Temperature in Heat Exchanger 2')
# plot_csv('Tsss_HX2.csv', 'Temperature in Heat Exchanger 2 (Tsss_HX2)')
