"""datanal: A Python package for data analysis and visualization.

datanal is a Python package for data analysis and visualization. It is designed to be used in conjunction with matplotlib and seaborn to create publication-quality figures.

The package is developed and maintained by Zhaoyang Li. If you have any questions or suggestions, please feel free to contact me at zhaoyangli@stanford.edu.

Modules
-------
datanal.utils: A module for general utilities
datanal.enzyme_kinetics: A module for enzyme kinetics data analysis
datanal.mst: A module for microscale thermophoresis (MST) data analysis
datanal.spr: A module for surface plasmon resonance (SPR) data analysis
datanal.itc: A module for isothermal titration calorimetry (ITC) data analysis

Classes
-------
EnzymeKinetics: A class for enzyme kinetics data analysis
MST: A class for microscale thermophoresis (MST) data analysis
SPR: A class for surface plasmon resonance (SPR) data analysis
ITC: A class for isothermal titration calorimetry (ITC) data analysis
"""

__version__ = "1.2.0"
__author__ = "Zhaoyang Li"
__email__ = "zhaoyangli@stanford.edu"

from warnings import warn

from matplotlib import ticker
from matplotlib.axes import Axes


def auto_ticks(ax: Axes) -> None:

    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6, min_n_ticks=3))
    x_tick_num = len(ax.get_xticks())
    if x_tick_num >= 6:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    elif x_tick_num >= 5:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    else:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6, min_n_ticks=3))
    y_tick_num = len(ax.get_yticks())
    if y_tick_num >= 6:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    elif y_tick_num >= 5:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    else:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))


def auto_units(quant: float, old_unit: str = "M") -> tuple[str, float]:

    units = {
        "M": (1, 2e-1),
        "mM": (1e-3, 2e-4),
        "μM": (1e-6, 2e-7),
        "nM": (1e-9, 2e-10),
        "pM": (1e-12, 2e-13),
        "fM": (1e-15, 2e-16),
    }

    if old_unit not in units:
        warn(f"Unknown unit {old_unit} was detected.")
        return old_unit, 1

    old_scale = units[old_unit][0]
    quant *= old_scale

    for new_unit, (new_scale, bound) in units.items():
        if quant >= bound:
            factor = old_scale / new_scale
            return new_unit, factor

    return old_unit, 1
