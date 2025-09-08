"""
A script to calculate integrated Planck radiation for BlastOne project

Author: Ranjith Ravichandran
"""
import pandas as pd
import matplotlib.pyplot as mplt
import numpy as np
import scipy as sc
# import scipy.optimize as sciopt
# import datetime as dt
# import peakutils
# from specutils.spectra import Spectrum1D, SpectralRegion
# from specutils.fitting import fit_generic_continuum
# from specutils.fitting.continuum import fit_continuum
# import astropy.units as u
# import warnings

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
    wl_range = np.arange(300, 3000, 10)
    # Sun_Temp = 5762 # Kelvin

    for T in T_range:
        Planck_rad_wide = Planck(wl_range, T, emissivity)
        mplt.plot(wl_range, Planck_rad_wide, label=str(T) + ' K', linewidth=0.5)
        wl_max = 2.898*1e6/T
        # mplt.vlines(x=wl_max, ymin=0, ymax=np.max(Planck_rad_wide), color='k', linewidth=0.5)
        mplt.plot(wl_max, np.max(Planck_rad_wide), marker='o', color='k', markersize=3)
        print(Planck_rad_wide[20])
        print(np.max(Planck_rad_wide))
        print('% rad: ', Planck_rad_wide[20]/np.max(Planck_rad_wide))

    # blue_sky_rad = Planck(wl_range, Sun_Temp, 1)/(wl_range*1e-9)**4
    # mplt.plot(wl_range, blue_sky_rad/np.max(blue_sky_rad), label='Rayleigh sky')
    mplt.legend(fontsize=fs_ticks, loc='upper right', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Wavelength (nm)', fontsize=fs_labels, fontweight='bold')
    mplt.ylabel('Spectral radiance $\mathregular{(Wm^{-2}sr^{-1}nm^{-1})}$', fontsize=fs_labels, fontweight='bold')
    mplt.xlim(250, 3050)
    # mplt.ylim(0, 0.045)
    mplt.axvspan(1000, 2200, alpha=0.15, color='red')

    # mplt.axvspan(longPass_wl, 1700, alpha=0.15, color='blue')
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('Planck_curve_dist.pdf', format='pdf', bbox_inches='tight')
    mplt.show()


def intPlanck():
    T_range = np.arange(500, 3500, 100)
    cred1p9_range = np.arange(1000, 1900, 1)
    cred2p2_range = np.arange(1200, 2200, 1)

    cred1p9_int_data = []
    cred2p2_int_data = []

    for T in T_range:
        Planck_rad_cred1p9 = Planck(cred1p9_range, T, emissivity)
        Planck_rad_cred2p2 = Planck(cred2p2_range, T, emissivity)

        cred1p9_int = integrated_radiance(Planck_rad_cred1p9, cred1p9_range)
        cred2p2_int = integrated_radiance(Planck_rad_cred2p2, cred2p2_range)
        cred1p9_int_data.append(cred1p9_int)
        cred2p2_int_data.append(cred2p2_int)
        # mplt.plot(T, xenics_int, 'ro', label='xenics')
        # mplt.plot(T, peak_int, 'ko', label='peak')
        # mplt.plot(T-273, xenics_int/peak_int, 'bo')

    mplt.figure(1)
    mplt.plot(T_range, cred1p9_int_data, 'ro-', label='CRED 1.0 to 1.9 um', markersize=3)
    mplt.plot(T_range, cred2p2_int_data, 'ko-', label='CRED 1.2 to 2.2 um', markersize=3)
    mplt.legend(fontsize=fs_ticks, loc='upper left', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (K)', fontsize=fs_labels, fontweight='bold')
    mplt.ylabel('Radiance $\mathregular{(Wm^{-2}sr^{-1})}$', fontsize=fs_labels, fontweight='bold')
    mplt.yscale('log')
    # mplt.xlim(350, 1750)
    # mplt.ylim(0.3, 1)
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('int_signal.pdf', format='pdf', bbox_inches='tight')

    mplt.figure(2)
    mplt.plot(T_range-273, (np.array(cred2p2_int_data) - np.array(cred1p9_int_data)) / np.array(cred2p2_int_data) * 100, 'bo-', label='% of radiance', markersize=3)
    mplt.legend(fontsize=fs_ticks, loc='upper right', fancybox=False).get_frame().set_linewidth(0.25)
    mplt.xlabel('Temperature (deg C)', fontsize=fs_labels, fontweight='bold')
    mplt.ylabel('(2.2 um - 1.9 um)/2.2 um', fontsize=fs_labels, fontweight='bold')
    # mplt.xlim(350, 1750)
    # mplt.ylim(0.3, 1)
    mplt.xticks(fontsize=fs_ticks), mplt.yticks(fontsize=fs_ticks)
    mplt.grid(linestyle=':', linewidth=0.05)
    mplt.savefig('fraction_int_signal.pdf', format='pdf', bbox_inches='tight')
    mplt.show()


if __name__ == '__main__':
    Planck_distribution()
    intPlanck()
