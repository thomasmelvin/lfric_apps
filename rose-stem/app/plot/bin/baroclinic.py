#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Need to set a non-interactive backend for suites
'''
Plots horizontal slices and vertical cross sections specified at particular
longitude and latitude values.
'''
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
# Need to set a non-interactive backend for suites
matplotlib.use('Agg')  # noqa: E402
from six.moves import range
import sys
import iris

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from read_data import read_ugrid_data
# Use magma colormap
try:
    magma = plt.get_cmap("magma")
except ValueError:
    from matplotlib.colors import ListedColormap
    from magma import magma_data
    magma = ListedColormap(magma_data, name='magma')
    plt.register_cmap(name='magma', cmap=cmap)

# Size of regular grid
ny, nx = 360, 720

plot_lon = 360
plot_lat = 280


def make_figures(filein, plotpath, fields, vertical_spacing):

    # Compute the vertical grid
    lid = 30.
    n_full = 31

    if vertical_spacing == 'um':
        # UM L38 set
        lid = 40.
        n_full = 39
        zi_f = np.array([
            .0, .0005095,  .0020380,  .0045854,  .0081519,  .0127373,
            .0183417,  .0249651,  .0326074,  .0412688,  .0509491,
            .0616485,  .0733668,  .0861040,  .0998603,  .1146356,
            .1304298,  .1472430,  .1650752,  .1839264,  .2037966,
            .2246857,  .2465938,  .2695209,  .2934670,  .3184321,
            .3444162,  .3714396,  .3998142,  .4298913,  .4620737,
            .4968308,  .5347160,  .5763897,  .6230643,  .6772068,
            .7443435,  .8383348, 1.000000])*lid

    elif vertical_spacing == 'dcmip':
        # DCMIP Stretched grid
        mu = 15.
        zi_f = np.zeros([n_full])
        for k in range(n_full):
            zi_f[k] = lid * \
                      (np.sqrt(mu * (float(k)/float(n_full))**2 + 1.) - 1.) \
                      / (np.sqrt(mu+1.) - 1)
    else:
        # Assume uniform grid
        zmin = 0.0
        zmax = lid
        zi_f = np.linspace(zmin, zmax, n_full)

    zi_h = 0.5*(zi_f[1:] + zi_f[0:n_full-1])

    direction = 'xy'  # , 'yz', 'xz'

    Wcube = read_ugrid_data(filein, 'w_in_wth')
    w_time = Wcube.coord('time')
    if w_time.units.is_time_reference():
        w_start = w_time.units.num2date(0)
        w_seconds = np.array([
            (date - w_start).total_seconds()
            for date in w_time.units.num2date(w_time.points)
        ])
    else:
        w_seconds = w_time.units.convert(w_time.points, 'seconds')
    day_one = np.flatnonzero(np.isclose(w_seconds, 24*3600, rtol=0, atol=0.1))
    if day_one.size:
        wind_level = 1
        xw = np.around(Wcube.coord('longitude').points, decimals=5)
        yw = np.around(Wcube.coord('latitude').points, decimals=5)
        xiw = np.linspace(np.amin(xw), np.amax(xw), nx)
        yiw = np.linspace(np.amin(yw), np.amax(yw), ny)
        xfw, yfw = np.meshgrid(xiw, yiw)
        data = Wcube[int(day_one[0]), wind_level].data
        w_day_one = griddata((xw, yw), data, (xfw, yfw), method='linear')
        w_nearest = griddata((xw, yw), data, (xfw, yfw), method='nearest')
        missing = np.isnan(w_day_one)
        w_day_one[missing] = w_nearest[missing]
        w_day_one -= np.mean(w_day_one, axis=1, keepdims=True)

        #w_limit = max(float(np.max(np.abs(w_day_one))), np.finfo(float).eps)
        w_limit = 5.5e-6
        w_levels = np.linspace(-w_limit, w_limit, 21)
        fig, ax = plt.subplots(figsize=(6, 3.4))
        cf_w = ax.contourf(
            xfw, yfw, w_day_one, levels=w_levels, cmap="RdBu_r",
            extend="both",
        )
        ax.set_xlabel("Longitude [degrees E]")
        ax.set_ylabel("Latitude [degrees N]")
        ax.set_title("Vertical velocity anomaly, level 1, day 1")
        fig.colorbar(cf_w, ax=ax, pad=0.02, label="w - zonal mean(w) [m/s]")
        fig.tight_layout()
        out = '%s/baroclinic_plot-w-anomaly-level1-day1.0-%s.png' % \
            (plotpath, direction)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out}")
    else:
        print("  No vertical velocity output at day 1; skipping anomaly plot")

    for t in [-2, -1]:

        fig, axes = plt.subplots(
            1, 2, figsize=(10, 3.4), sharex=True, sharey=True)
        axes = axes[None, :]
        Tcube = read_ugrid_data(filein, 'theta')
        Pcube = read_ugrid_data(filein, 'exner')
        Ucube = read_ugrid_data(filein, 'u_in_w2h')
        Wcube = read_ugrid_data(filein, 'w_in_wth')

        T_levels = np.arange(220, 310 + 10, 10)
        ps_levels = np.arange(930, 1030 + 10, 10)
        #ps_levels = np.arange(916, 1012 + 8, 8) # Match Gungho Spherical paper

        Tn_levs = len(Tcube.coord(Tcube.dim_coords[-1].name()).points)
        Pn_levs = len(Pcube.coord(Pcube.dim_coords[-1].name()).points)
        Un_levs = len(Ucube.coord(Ucube.dim_coords[-1].name()).points)
        Wn_levs = len(Wcube.coord(Wcube.dim_coords[-1].name()).points)

        T_data = np.zeros((ny, nx, Tn_levs))
        P_data = np.zeros((ny, nx, Pn_levs))
        U_data = np.zeros((ny, nx))
        W_data = np.zeros((ny, nx))

        time = np.around(Tcube.coord('time').points, decimals=1)
        day = time[t]/(24*3600)

        # Compute the horizontal grid
        x = np.around(Tcube.coord('longitude').points, decimals=5)
        y = np.around(Tcube.coord('latitude').points, decimals=5)

        xmin = np.amin(x)
        xmax = np.amax(x)
        ymin = np.amin(y)
        ymax = np.amax(y)

        # Generate a regular grid to interpolate the data.
        xi = np.linspace(xmin, xmax, nx)
        yi = np.linspace(ymin, ymax, ny)

        xf, yf = np.meshgrid(xi, yi)

        # Choose the correct vertical level set
        ziT = zi_f
        ziP = zi_h

        # Interpolate using delaunay triangularization
        for p, level in enumerate(range(Tn_levs)):
            data = Tcube[t, level].data
            fi = griddata((x, y), data, (xf, yf), method='linear')
            fi_n = griddata((x, y), data, (xf, yf), method='nearest')
            fi[np.isnan(fi)] = fi_n[np.isnan(fi)]
            T_data[:, :, level] = fi

        for p, level in enumerate(range(Pn_levs)):
            data = Pcube[t, level].data
            fi = griddata((x, y), data, (xf, yf), method='linear')
            fi_n = griddata((x, y), data, (xf, yf), method='nearest')
            fi[np.isnan(fi)] = fi_n[np.isnan(fi)]
            # Convert to hPa
            rd = 287.05
            p0 = 100000.0
            kappa = rd/1005.0
            P_data[:, :, level] = 0.01*fi**(1.0/kappa) * p0

        lat, lon = np.meshgrid(yi, xi)
        # Extrapolate pressure to the surface
        P0 = P_data[:, :, 0] + (zi_f[0] - zi_h[0]) * \
            (P_data[:, :, 0] - P_data[:, :, 1]) \
            / (zi_h[0] - zi_h[1])


        # Compute T850
        exner850 = (850.0*100.0/p0)**(kappa)
        T850 = np.zeros([ny,nx])       
        for i in range(nx):
            for j in range(ny):
              f = interp1d(P_data[j,i,:], zi_h)
              z850 =  f(850.0)
              f = interp1d(zi_f, T_data[j,i,:])
              T850[j,i] = f(z850)*exner850

        # Plot Surface Pressure
        ax_ps = axes[0, 0]
        cf_ps = ax_ps.contourf(
            lon, lat, P0.T, levels=ps_levels, cmap="RdBu_r", extend="both",
        )
        cs_ps = ax_ps.contour(
            lon, lat, P0.T, levels=ps_levels, colors="k", linewidths=0.5,
        )
        ax_ps.clabel(cs_ps, fmt="%d", fontsize=6)
        ax_ps.set_ylabel("Lat [°N]")
        ax_ps.set_title(f"ps [hPa]")
        plt.colorbar(cf_ps, ax=ax_ps, pad=0.02, label="ps [hPa]")

        # Plot T850
        ax_T = axes[0, 1]
        cf_T = ax_T.contourf(
            lon, lat, T850.T, levels=T_levels, cmap="RdYlBu_r", extend="both",
        )
        cs_T = ax_T.contour(
            lon, lat, T850.T, levels=T_levels, colors="k", linewidths=0.5,
        )
        ax_T.clabel(cs_T, fmt="%d", fontsize=6)
        ax_T.set_title(f"T @ 850 hPa [K]")
        plt.colorbar(cf_T, ax=ax_T, pad=0.02, label="T [K]")

        for ax in axes[-1, :]:
            ax.set_xlabel("Lon [°E]")

        fig.suptitle(
            "Surface pressure and 850 hPa temperature (JW2006 comparison)",
            fontsize=12, y=1.005,
        )
        fig.tight_layout()
        out = '%s/baroclinic_plot-%s-day%s-%s.png' % \
            (plotpath, "p0-T850", day, direction)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out}")

        # Wind plot
        wind_level = 10

        u_levels = np.arange(15, 35+2, 2)
        w_levels = np.arange(-0.1, 0.1+0.01, 0.01)

        xh = np.around(Ucube.coord('longitude').points, decimals=5)
        yh = np.around(Ucube.coord('latitude').points, decimals=5)

        data = Ucube[t, wind_level].data
        fi = griddata((xh, yh), data, (xf, yf), method='linear')
        fi_n = griddata((xh, yh), data, (xf, yf), method='nearest')
        fi[np.isnan(fi)] = fi_n[np.isnan(fi)]
        U_data[:, :] = fi

        data = Wcube[t, wind_level].data
        fi = griddata((x, y), data, (xf, yf), method='linear')
        fi_n = griddata((x, y), data, (xf, yf), method='nearest')
        fi[np.isnan(fi)] = fi_n[np.isnan(fi)]
        W_data[:, :] = fi

        fig, axes = plt.subplots(
            1, 2, figsize=(10, 3.4), sharex=True, sharey=True)
        axes = axes[None, :]

        # Zonal Wind
        ax_u = axes[0, 0]
        cf_u = ax_u.contourf(
            lon, lat, U_data.T, levels=u_levels, cmap="RdBu_r", extend="both",
        )
        cs_u = ax_u.contour(
            lon, lat, U_data.T, levels=u_levels, colors="k", linewidths=0.5,
        )
        ax_u.clabel(cs_u, fmt="%d", fontsize=6)
        ax_u.set_ylabel("Lat [°N]")
        ax_u.set_title(f"u [m/s]")
        plt.colorbar(cf_u, ax=ax_u, pad=0.02, label="zonal wind [m/s]")
  
        # Vertical Wind
        ax_w = axes[0, 1]
        cf_w = ax_w.contourf(
            lon, lat, W_data.T, levels=w_levels, cmap="RdBu_r", extend="both",
        )
        cs_w = ax_w.contour(
            lon, lat, W_data.T, levels=w_levels, colors="k", linewidths=0.5,
        )
        ax_w.clabel(cs_w, fmt="%d", fontsize=6)
        ax_w.set_ylabel("Lat [°N]")
        ax_w.set_title(f"w [m/s]")
        plt.colorbar(cf_w, ax=ax_w, pad=0.02, label="zonal wind [m/s]")
      

        for ax in axes[-1, :]:
            ax.set_xlabel("Lon [°E]")

        fig.suptitle(
            "Zonal and vertical wind",
            fontsize=12, y=1.005,
        )
        fig.tight_layout()
        out = '%s/baroclinic_plot-%s-day%s-%s.png' % \
            (plotpath, "wind", day, direction)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out}")



if __name__ == "__main__":

    try:
        args = sys.argv[:]
        filein, plotpath, vertical_grid = args[1:4]
        field_list = None
        if len(args[:]) > 4:
            field_list = args[4].split(':')
    except ValueError:
        print("Usage: {0} <filein> <plotpath> <vertical_grid> [<fields_list>] "
              .format(sys.argv[0]))
        exit(1)

    make_figures(filein, plotpath, field_list, vertical_grid)
