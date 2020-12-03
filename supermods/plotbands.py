#!/usr/bin/env python

import numpy as np
import os
from os.path import expanduser 
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gs
import supermods.col_key as col
from pathlib import Path

path=str(Path(__file__).parent.absolute())+'/../saved_plots/'

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 14,
        }

def clear_plot(figure):
    ax = figure.add_subplot(111)
    ax.set_ylabel("Energy (eV)",fontsize=14)
    ax.set_yticks([])
    ax.set_xticks([])

def plotbands(bands_dat, bands_out, bands_in, bg, fermi, fermi_check, fermi_shift,
        ev_a, ev_b, figure, save_key, label, band_lw, only_bands, gs1):
    def find_path(fstring, ax): 
        f = open(fstring,'r')
        x = np.zeros(0)
        for i in f:
            if "high-symmetry" in i:
                x = np.append(x,float(i.split()[-1]))
        ticks=np.unique(x, axis=0)
        f.close()
        ax.set(xticks=ticks)
        return x
    # This function takes in bands.dat.gnu, the fermi level, bands.out, a subplot, and the label
    # It then extracts the band data, plots the bands, the fermi level, and q-points
    def bandplot(bands_dat_gnu, bands_out, bands_in, bg, fermi, fermi_check, fermi_shift, figure, save_key, label, band_lw, only_bands, gs1):
        if only_bands:
            ax = figure.add_subplot(111)
            ax.set_title(label, fontname='serif', fontsize=14)
        else:
            ax = figure.add_subplot(gs1[0, 0])
        ax.set_xlabel("K-Path", fontname='serif', fontsize=14, labelpad=10)
        ax.set_ylabel("Energy (eV)", fontname='serif', fontsize=14)
        if fermi_shift:
            fermi_shift=fermi
        z = np.loadtxt(bands_dat_gnu)
        # Evaluate the number of bands associated with each k-point, used for managing duplicate q-points
        kpt, counts = np.unique(z[:,0], return_counts=True)
        # List of all x coordinates for the structure's high-sym points, includes duplicates 
        q_points=find_path(bands_out, ax)
        k_minus_q=np.array(list(set(kpt).symmetric_difference(q_points)))
        kpts=np.sort(np.concatenate((k_minus_q, q_points)))
        nbnd = len(z[z[:,0]==kpt[1]])
        k_bnds=dict(zip(kpt, counts))
        axis = [min(kpts), max(kpts), ev_a, ev_b]
        bands = []
        for i in range(0, nbnd):
            bands.append(np.zeros([len(kpts), 2])) # Initialize empty array to store band data
        q_count=q_points.tolist()
        for i in range(0, len(kpts)):
            if int(k_bnds[kpts[i]]/nbnd)>1:
                sel = z[z[:,0] == kpts[i]][int(k_bnds[kpts[i]]/nbnd)-q_count.count(kpts[i])::int(k_bnds[kpts[i]]/nbnd)]
                q_count.remove(kpts[i])
            else:
                sel = z[z[:,0] == kpts[i]]
            for j in range(0, nbnd):# Separation of energy values into single bands
                bands[j][i][0] = kpts[i]
                bands[j][i][1] = sel[j][1]-fermi_shift
        # Band gap calculation, first we find the min/max energy values at each kpt
        if fermi_shift:
            bg_shift=0
        else:
            bg_shift=fermi
        lumo_eigs=np.array([(np.minimum.reduce(array, axis=0))[1]-bg_shift for array in bands])
        homo_eigs=np.array([(np.maximum.reduce(array, axis=0))[1]-bg_shift for array in bands])
        # Determination of absolute min/max energy with respect to the fermi level
        min_lumo_eig=lumo_eigs[np.where(lumo_eigs > 0, lumo_eigs, np.inf).argmin()]
        max_homo_eig=homo_eigs[np.where(lumo_eigs < 0, lumo_eigs, -np.inf).argmax()]
        bandgap=round(min_lumo_eig + abs(max_homo_eig), 4)
        # Retrieve all kpt values with the max/min energies and create an array of coordinate sets to label the plot
        gap_ticks=np.vstack([eig for array in bands for eig in array if eig[1]-bg_shift in [max_homo_eig, min_lumo_eig]])
        for i in bands: # Plot bands
            ax.plot(i[:,0], i[:,1], color="black")
        # Plot the band gap markers
        if bg:
            for eig in gap_ticks:
                print(eig, fermi)
                ax.plot([eig[0]], [eig[1]], marker='x',
                    markersize=7, color=['red' if eig[1]-fermi>0 else 'blue'][0])
            ax.text(0, ev_a - abs(abs(ev_b)-ev_a)/10, "Band gap: {} eV".format(bandgap), fontsize=10, color='red')
        ax.set_ylim([axis[2], axis[3]])
        ax.set_xlim([axis[0], axis[1]])
        for line in ax.get_lines():
            line.set_linewidth(band_lw)
        if fermi_check:
            ax.text(0, ev_b + abs(abs(ev_b) - ev_a)/40, "Fermi level: {} eV".format(fermi), fontsize=10, color='red')
            if not fermi_shift:
                ax.plot([min(kpt),max(kpt)],[fermi, fermi],color='red', linewidth=1.0)
        for j in np.unique(q_points, axis=0):
            x1 = [j, j]
            x2 = [axis[2], axis[3]]
            ax.plot(x1, x2,'--', linewidth=1.0, color='black', alpha=0.75)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
        ax.set_xticklabels(['X', '\u0393', 'Y', 'L', '\u0393', 'Z', 'N', '\u0393', 'M', 'R', '\u0393'], fontdict=font)
        plt.yticks(fontsize=12)

        if save_key:
            figure.savefig(path+label.replace(" ",".")+'.pdf', format='pdf', dpi=1200)
        print("Total bands: {}\nBangap: {} eV\nFermi level: {} eV".format(nbnd, bandgap, fermi))
    
    return bandplot(bands_dat, bands_out, bands_in, bg, fermi, fermi_check, fermi_shift, figure, save_key, label, band_lw, only_bands, gs1)

def dosplot(dos_files, figure, save_key, label, fermi, y_range, x_range, col_keys, fermi_check, fermi_shift, only_pdos, gs1, fill, dos_lw):
    if only_pdos:
        ax = figure.add_subplot(111)
        ax.set_title(label, fontname='serif', fontsize=14)
        ax.set_ylabel("Energy (eV)", fontname='serif', fontsize=14)
    else:
        ax = figure.add_subplot(gs1[0, 1])
        ax.set_yticks([])
    ax.set_xticks([])   
    keys = {}
    leg_labels={}
    for i in dos_files:
        atom=i[0]
        leg_labels.update({atom: [i[3], i[1]]})
        if atom not in keys:
            t = np.loadtxt(i[2])
            keys[atom] = t[:,[0,1]] # Only the energy and ldos are extracted
        else:
            t = np.loadtxt(i[2])
            t = t[:,[0,1]]
            keys[atom][:,1] = np.maximum(keys[atom][:,1],t[:,1]) # Since we already had a numpy array here, we take the maximum of each point
    if fermi_check and not fermi_shift:
        ax.plot(x_range,[fermi, fermi], color='red', linewidth=1.0)
    if fermi_shift:
        fermi_shift=fermi
    for i in keys:
        keys[i] = np.fliplr(keys[i]) # We flip the array for vertical plotting of the pdos
    for i in keys:
        ax.plot(keys[i][:,0], keys[i][:,1]-fermi_shift, leg_labels[i][1], label=i)
        if fill:
            ax.fill_betweenx(keys[i][:,1]-fermi_shift,keys[i][:,0], color=leg_labels[i][1], alpha=0.25) # Optional filling on PDOS
    ax.set_ylim(y_range)
    ax.set_xlim(x_range)
    for line in ax.get_lines():
        line.set_linewidth(dos_lw)

    ax.set_xlabel("PDOS (arb)", fontname='serif', fontsize=14, labelpad=10)
    ax.legend([v[0] for v in leg_labels.values()], bbox_to_anchor=(1, 1), loc=1, borderaxespad=0, fontsize=10)
    if save_key:
        figure.savefig(path+label.replace(" ",".")+'.pdf', format='pdf', dpi=1200)   

def bands_and_pdos(bands_dat, bands_out, bands_in, fermi, ev_a, ev_b, 
    figure, save_key, label, band_lw, dos_lw, bg, fermi_check, fermi_shift, dos_files, x_range, col_keys, fill):

    gs1 = gs.GridSpec(1,2, width_ratios=[2,1], figure=figure)
    plotbands(bands_dat, bands_out, bands_in, bg, fermi, fermi_check, fermi_shift, 
        ev_a, ev_b, figure, False, label, band_lw, False, gs1)
    dosplot(dos_files, figure, False, label, fermi, [ev_a, ev_b], 
        x_range, col_keys, fermi_check, fermi_shift, False, gs1, fill, dos_lw)
    gs1.update(wspace=0.0, hspace=0.0)
    figure.suptitle(label, fontname='serif', fontsize=14, y=0.92)
    if save_key:
        figure.savefig(path+label.replace(" ",".")+'.pdf', format='pdf', dpi=1200)
