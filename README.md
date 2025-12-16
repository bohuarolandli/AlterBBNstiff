# AlterBBNstiff

The original AlterBBN 2.2 code is modified to accommodate new ingredients/implementations of Beyond-Standad-Model (BSM) physics, which include:
1. A dark radiation component parametrized by its energy density ratio to photons at the initial temperature of BBN calculations ($T=27\times10^9\\mathrm{K}$), **$\kappa_{\mathrm{rad},i}$**.
This parameter is proportional to the effective number of relativistic species, **$\Delta N_\mathrm{eff}$**.
2. A stiff component (with equation of state $w=1$) parametrized by its energy density ratio to photons at the initial temperature, **$\kappa_{\mathrm{s},i}$**.

Some of the key nuclear reaction rates are also updated, which affect the primordial helium and deuterium abundances ($Y_\mathrm{P}$ and D/H).

We include a python interface, "alter_stiff.py", for applications of parameter inferences (e.g., to be incorporated into a Cobaya pipeline).

---------------

See the README_AlterBBN file for documatation of the original AlterBBN 2.2 code.
