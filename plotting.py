import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def plot_multiple_sections(data, sections, title, xlabel='', ylabel='', labels=None):
    plt.figure(figsize=(14, 6))
    for i, (start, end) in enumerate(sections):
        label = labels[i] if labels and i < len(labels) else f'Section {i+1}'
        plt.plot(data[start:end], label=label)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()

def main():
    time_span = 2500
    # Get the folder path from the command line argument
    if len(sys.argv) < 2:
        print("Usage: python plotting.py <folder_path>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    if not os.path.exists(folder_path):
        print(f"Error: The folder '{folder_path}' does not exist.")
        sys.exit(1)

    # Load data from CSV files
    rho_matrix = np.loadtxt(os.path.join(folder_path, 'rho_matrix.csv'), delimiter=',')
    rho= np.loadtxt(os.path.join(folder_path, 'rho.csv'), delimiter=',')
    y_n = np.loadtxt(os.path.join(folder_path, 'y_n.csv'), delimiter=',')
    y_th = np.loadtxt(os.path.join(folder_path, 'y_th.csv'), delimiter=',')
    y_hx1 = np.loadtxt(os.path.join(folder_path, 'y_hx1.csv'), delimiter=',')
    y_hx2 = np.loadtxt(os.path.join(folder_path, 'y_hx2.csv'), delimiter=',')
    
    ci_middle_matrix = np.loadtxt(os.path.join(folder_path, 'ci_middle_matrix.csv'), delimiter=',')
    phi_middle_matrix = np.loadtxt(os.path.join(folder_path, 'phi_middle_matrix.csv'), delimiter=',')

    # Define the sections and labels for each plot
    plot_multiple_sections(rho_matrix, [(0,time_span)], 'avg Reactivity (pcm)')
    plot_multiple_sections(rho, [(0, 200)], 'Reactivity (pcm)')
    # Plot phi1 and phi2 together
    plot_multiple_sections(y_n, [(0, 200), (200, 400)], 'Neutron Flux φ1 and φ2 (n/(cm^2/s))',
                           ylabel='Flux', labels=['φ1', 'φ2'])

    # Plot all ci (c1 to c6) together
    ci_sections = [(400 + 200 * i, 600 + 200 * i) for i in range(6)]
    ci_labels = [f'c{i+1}' for i in range(6)]
    plot_multiple_sections(y_n, ci_sections, 'Delayed Neutron Precursors c1 to c6',
                           ylabel='Concentration', labels=ci_labels)

    # Plot fuel and graphite temperatures in core together
    plot_multiple_sections(y_th, [(0, 200), (200, 400)], 'Core Temperatures (Kelvin)',
                           ylabel='Temperature (K)', labels=['Fuel', 'Graphite'])

    # Plot fuel and coolant temperatures in Heat Exchanger 1 together
    plot_multiple_sections(y_hx1, [(0, 200), (200, 400)], 'Temperatures in Heat Exchanger 1 (Kelvin)',
                           ylabel='Temperature (K)', labels=['Fuel', 'Coolant'])

    # Plot fuel and coolant temperatures in Heat Exchanger 2 together
    plot_multiple_sections(y_hx2, [(0, 200), (200, 400)], 'Temperatures in Heat Exchanger 2 (Kelvin)',
                           ylabel='Temperature (K)', labels=['Fuel', 'Coolant'])
    
    ci_middle_sections = [(time_span * i, time_span * (i+1)) for i in range(6)]
    ci_middle_labels = [f'c{i+1}' for i in range(6)]
    plot_multiple_sections(ci_middle_matrix, ci_middle_sections, 'Delayed Neutron Precursors in the Middle (at L/2)',
                           ylabel='Concentration', labels=ci_middle_labels)
    
    phi_middle_sections = [(time_span * i, time_span * (i+1)) for i in range(2)]
    phi_middle_labels = ['φ1', 'φ2']
    plot_multiple_sections(phi_middle_matrix, phi_middle_sections, 'Neutron Flux in the Middle (at L/2)',
                           ylabel='Flux', labels=phi_middle_labels)

if __name__ == "__main__":
    main()
