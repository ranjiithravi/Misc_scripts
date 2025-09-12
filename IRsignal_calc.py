"""
A script to calculate integrated Planck radiation for a given wavelength band

Author: Ranjith Ravichandran
"""
import pandas as pd
import matplotlib.pyplot as mplt
import numpy as np
import scipy as sc

fs_labels = 10
fs_ticks = 10
mplt.rcParams["font.family"] = "Times New Roman"
mplt.rcParams["mathtext.fontset"] = "custom"
mplt.rcParams["mathtext.it"] = "STIXGeneral:italic"
mplt.rcParams["axes.linewidth"] = 0.75

# Inputs
# T_range = np.arange(1000, 3000, 50)
T_range = np.arange(500, 3500, 50)
# longPass_wl = 850
camera_sen_range = [900, 1700]
emissivity = 1

def Planck(wl, T, factor):
    """ Return the Planck function """
    c = sc.constants.c
    h = sc.constants.h
    k = sc.constants.k
    B = 2 * h * c ** 2 / (wl * 1e-9) ** 5 / (np.exp(h * c / (wl * 1e-9) / k / T) - 1)
    B = B * 1e-9
    return factor * B

def integrated_radiance(spRad, wl):
    """ Return integrated radiance in chosen wavelength window"""
    int_Planck = sc.integrate.trapezoid(spRad, wl)

    return int_Planck

def Planck_distribution():
    T_range = np.arange(500, 3500, 300)
    wl_range = np.arange(300, 5000, 10)
    # Sun_Temp = 5762 # Kelvin

    for T in T_range:
        Planck_rad_wide = Planck(wl_range, T, emissivity)
        mplt.plot(wl_range, Planck_rad_wide, label=str(T) + ' K', linewidth=0.5)
        wl_max = 2.898*1e6/T
        # mplt.vlines(x=wl_max, ymin=0, ymax=np.max(Planck_rad_wide), color='k', linewidth=0.5)
        mplt.plot(wl_max, np.max(Planck_rad_wide), marker='o', color='k', markersize=3)
        #print(Planck_rad_wide[20])
        #print(np.max(Planck_rad_wide))
        #print('% rad: ', Planck_rad_wide[20]/np.max(Planck_rad_wide))

    # blue_sky_rad = Planck(wl_range, Sun_Temp, 1)/(wl_range*1e-9)**4
    # mplt.plot(wl_range, blue_sky_rad/np.max(blue_sky_rad), label='Rayleigh sky')
    mplt.legend(fontsize=fs_ticks, loc='upper right', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Wavelength (nm)', fontsize=fs_labels, fontweight='bold')
    mplt.ylabel(r'Spectral radiance $\mathregular{(Wm^{-2}sr^{-1}nm^{-1})}$', fontsize=fs_labels, fontweight='bold')
    mplt.xlim(250, 5050)
    # mplt.ylim(0, 0.045)
    mplt.axvspan(1000, 5000, alpha=0.15, color='red')
    mplt.axvspan(3000, 5000, alpha=0.15, color='black')

    # mplt.axvspan(longPass_wl, 1700, alpha=0.15, color='blue')
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('Planck_curve_dist.pdf', format='pdf', bbox_inches='tight')
    mplt.show()


def intPlanck():
    T_range = np.arange(500, 3100, 50)

    IR_1to5_range = np.arange(1000, 5000, 1)
    IR_1p5to5p4_range = np.arange(1500, 5400, 1)
    IR_1p5to5p5_range = np.arange(1500, 5500, 1)
    IR_3to5_range = np.arange(3000, 5000, 1)
    IR_3to5p5_range = np.arange(3000, 5500, 1)
    IR_3p7to4p15_range = np.arange(3700, 4150, 1)

    IR_1to5_int_data = []
    IR_1p5to5p4_int_data = []
    IR_1p5to5p5_int_data = []
    IR_3to5_int_data = []
    IR_3to5p5_int_data = []
    IR_3p7to4p15_int_data =[]

    for T in T_range:
        Planck_rad_IR_1to5 = Planck(IR_1to5_range, T, emissivity)
        Planck_rad_IR_1p5to5p4 = Planck(IR_1p5to5p4_range, T, emissivity)
        Planck_rad_IR_1p5to5p5 = Planck(IR_1p5to5p5_range, T, emissivity)
        Planck_rad_IR_3to5 = Planck(IR_3to5_range, T, emissivity)
        Planck_rad_IR_3to5p5 = Planck(IR_3to5p5_range, T, emissivity)
        Planck_rad_IR_3p7to4p15 = Planck(IR_3p7to4p15_range, T, emissivity)

        IR_1to5_int = integrated_radiance(Planck_rad_IR_1to5, IR_1to5_range)
        IR_1p5to5p4_int = integrated_radiance(Planck_rad_IR_1p5to5p4, IR_1p5to5p4_range)
        IR_1p5to5p5_int = integrated_radiance(Planck_rad_IR_1p5to5p5, IR_1p5to5p5_range)
        IR_3to5_int = integrated_radiance(Planck_rad_IR_3to5, IR_3to5_range)
        IR_3to5p5_int = integrated_radiance(Planck_rad_IR_3to5p5, IR_3to5p5_range)
        IR_3p7to4p15_int = integrated_radiance(Planck_rad_IR_3p7to4p15, IR_3p7to4p15_range)

        IR_1to5_int_data.append(IR_1to5_int)
        IR_1p5to5p4_int_data.append(IR_1p5to5p4_int)
        IR_1p5to5p5_int_data.append(IR_1p5to5p5_int)
        IR_3to5_int_data.append(IR_3to5_int)
        IR_3to5p5_int_data.append(IR_3to5p5_int)
        IR_3p7to4p15_int_data.append(IR_3p7to4p15_int)

    mplt.figure(1)
    mplt.plot(T_range, IR_1to5_int_data, 'o-', label=r'FLIR A6780, 1.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, IR_1p5to5p4_int_data, 'o-', label=r'FAST M350, 1.5 $-$ 5.4 $\mu$m', markersize=3)
    mplt.plot(T_range, IR_1p5to5p5_int_data, 'o-', label=r'VELOX 327k SM, 1.5 $-$ 5.5 $\mu$m', markersize=3)
    mplt.plot(T_range, IR_3to5_int_data, 'o-', label=r'FLIR, FAST Lens, 3.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, IR_3to5p5_int_data, 'o-', label=r'VELOX Lens, 3.0 $-$ 5.5 $\mu$m', markersize=3)
    mplt.plot(T_range, IR_3p7to4p15_int_data, 'o-', label=r'IR6300, 3.7 $-$ 4.15 $\mu$m', markersize=3)
    mplt.legend(fontsize=fs_ticks, loc='lower right', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (K)', fontsize=fs_labels, fontweight='bold')
    mplt.ylabel(r'Radiance $\mathregular{(Wm^{-2}sr^{-1})}$', fontsize=fs_labels, fontweight='bold')
    mplt.yscale('log')
    # mplt.xlim(350, 1750)
    # mplt.ylim(0.3, 1)
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('int_signal.pdf', format='pdf', bbox_inches='tight')
    '''
    mplt.figure(2)
    mplt.plot(T_range, np.array(IR_1to5_int_data) / np.array(IR_3to5_int_data),
              'o-', label=r'FLIR A6780, 1.0 $-$ 5.0 $\mu$m / 3.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_1p5to5p5_int_data) / np.array(IR_3to5p5_int_data),
              'o-', label=r'VELOX 327k SM, 1.5 $-$ 5.5 $\mu$m / 3.0 $-$ 5.5 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_1p5to5p4_int_data) / np.array(IR_3to5_int_data),
              'o-', label=r'FAST M350, 1.5 $-$ 5.4 $\mu$m / 3.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.array(IR_3p7to4p15_int_data),
              'o-', label=r'IR6300, 3.0 $-$ 5.0 $\mu$m / 3.7 $-$ 4.15 $\mu$m', markersize=3)
    mplt.ylabel('Signal ratio (-)', fontsize=fs_labels, fontweight='bold')
    #mplt.yscale('log')
    mplt.legend(fontsize=fs_ticks, loc='upper left', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (K)', fontsize=fs_labels, fontweight='bold')
    # mplt.xlim(350, 1750)
    # mplt.ylim(0.3, 1)
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('fraction_int_signal.pdf', format='pdf', bbox_inches='tight')
    '''
    '''
    mplt.figure(3)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.array(IR_1to5_int_data),
              'o-', label=r'FLIR A6780, 3.0 $-$ 5.0 $\mu$m / 1.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5p5_int_data) / np.array(IR_1p5to5p5_int_data),
              'o-', label=r'VELOX 327k SM, 3.0 $-$ 5.5 $\mu$m / 1.5 $-$ 5.5 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.array(IR_1p5to5p4_int_data),
              'o-', label=r'FAST M350, 3.0 $-$ 5.0 $\mu$m / 1.5 $-$ 5.4 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3p7to4p15_int_data) / np.array(IR_3to5_int_data),
              'o-', label=r'IR6300, 3.7 $-$ 4.15 $\mu$m / 3.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.ylabel('Signal ratio (-)', fontsize=fs_labels, fontweight='bold')
    # mplt.yscale('log')
    mplt.legend(fontsize=fs_ticks, loc='upper left', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (K)', fontsize=fs_labels, fontweight='bold')
    mplt.xlim(500, 2000)
    mplt.ylim(1, 5)
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('fraction_int_signal_zoomed.pdf', format='pdf', bbox_inches='tight')
    '''
    mplt.figure(4)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.array(IR_1to5_int_data),
              'o-', label=r'FLIR A6780, 3.0 $-$ 5.0 $\mu$m / 1.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5p5_int_data) / np.array(IR_1p5to5p5_int_data),
              'o-', label=r'VELOX 327k SM, 3.0 $-$ 5.5 $\mu$m / 1.5 $-$ 5.5 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.array(IR_1p5to5p4_int_data),
              'o-', label=r'FAST M350, 3.0 $-$ 5.0 $\mu$m / 1.5 $-$ 5.4 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3p7to4p15_int_data) / np.array(IR_3to5_int_data),
              'o-', label=r'IR6300, 3.7 $-$ 4.15 $\mu$m / 3.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.ylabel('Signal ratio (-)', fontsize=fs_labels, fontweight='bold')
    # mplt.yscale('log')
    mplt.legend(fontsize=fs_ticks, loc='upper right', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (K)', fontsize=fs_labels, fontweight='bold')
    #mplt.xlim(500, 2000)
    mplt.ylim(0, 1)
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('fraction_int_signal.pdf', format='pdf', bbox_inches='tight')

    mplt.figure(5)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.array(IR_1to5_int_data),
              'o-', label=r'FLIR A6780, 3.0 $-$ 5.0 $\mu$m / 1.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5p5_int_data) / np.array(IR_1p5to5p5_int_data),
              'o-', label=r'VELOX 327k SM, 3.0 $-$ 5.5 $\mu$m / 1.5 $-$ 5.5 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.array(IR_1p5to5p4_int_data),
              'o-', label=r'FAST M350, 3.0 $-$ 5.0 $\mu$m / 1.5 $-$ 5.4 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3p7to4p15_int_data) / np.array(IR_3to5_int_data),
              'o-', label=r'IR6300, 3.7 $-$ 4.15 $\mu$m / 3.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.ylabel('Signal ratio (-)', fontsize=fs_labels, fontweight='bold')
    # mplt.yscale('log')
    mplt.legend(fontsize=fs_ticks, loc='upper right', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (K)', fontsize=fs_labels, fontweight='bold')
    mplt.xlim(500, 1800)
    mplt.ylim(0.2, 1)
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('fraction_int_signal_zoomed.pdf', format='pdf', bbox_inches='tight')
    '''
    mplt.figure(4)
    mplt.plot(T_range, np.array(IR_1to5_int_data) / np.min(np.array(IR_1to5_int_data)),
              'o-', label=r'normalised, 1.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_1p5to5p5_int_data) / np.min(np.array(IR_1p5to5p5_int_data)),
              'o-', label=r'normalised, 1.5 $-$ 5.5 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3to5_int_data) / np.min(np.array(IR_3to5_int_data)),
              'o-', label=r'normalised, 3.0 $-$ 5.0 $\mu$m', markersize=3)
    mplt.plot(T_range, np.array(IR_3p7to4p15_int_data) / np.min(np.array(IR_3p7to4p15_int_data)),
              'o-', label=r'normalised, 3.7 $-$ 4.15 $\mu$m', markersize=3)
    mplt.ylabel('Normalised signal (-)', fontsize=fs_labels, fontweight='bold')
    # mplt.yscale('log')
    mplt.legend(fontsize=fs_ticks, loc='upper left', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (K)', fontsize=fs_labels, fontweight='bold')
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('normalised_int_signal.pdf', format='pdf', bbox_inches='tight')
    '''
    mplt.show()

if __name__ == '__main__':
    Planck_distribution()
    intPlanck()
