import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def plot_csv(filename, title, xlabel='', ylabel=''):
    data = np.loadtxt(filename, delimiter=',')
    plt.figure(figsize=(14, 6))
    if len(data.shape) == 1:
        plt.plot(data)
    else:
        for i in range(data.shape[1]):
            plt.plot(data[:, i], label=f'Column {i+1}')
        plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def main():
    # Get the folder path from the command line argument
    if len(sys.argv) < 2:
        print("Usage: python plotting.py <folder_path>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    if not os.path.exists(folder_path):
        print(f"Error: The folder '{folder_path}' does not exist.")
        sys.exit(1)

    # Change directory to the specified folder
    os.chdir(folder_path)

    # Plot data from CSV files
    plot_csv('rho_matrix.csv', 'Reactivity (pcm, avg)')
    plot_csv('phi_middle_matrix.csv', 'Neutron Flux in the Middle (at L/2)')
    plot_csv('ci_middle_matrix.csv', 'Delayed Neutron Precursors in the Middle (at L/2)')
    plot_csv('temperature_fuel_middle_matrix.csv', 'Temperature in the Core Middle')

    plot_csv('phi.csv', 'Neutron Flux (n/(cm^2/s))')
    plot_csv('ci.csv', 'Delayed Neutron Precursors')
    plot_csv('rho.csv', 'Reactivity (pcm, avg)')
    plot_csv('temperature_core.csv', 'Temperature in the Core (Kelvin, col1-fuel, col2-graphite)')
    plot_csv('temperature_HX1.csv', 'Temperature in Heat Exchanger 1 (Kelvin, col1-fuel, col2-coolant1)')
    plot_csv('temperature_HX2.csv', 'Temperature in Heat Exchanger 2 (Kelvin)')

if __name__ == "__main__":
    main()