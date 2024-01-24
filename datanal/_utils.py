import os
import json
from typing import Any
from warnings import warn

import numpy as np
import matplotlib.pyplot as plt  # For variable ax in scale_tick
from matplotlib import rcParams, ticker
import seaborn as sns  # For color palette in config file
from cycler import cycler  # For color cycle in config file


def set_fig(config_file: str) -> None:
    """Set figure style.

    Parameters
    ----------
    config_file : str, optional
        Path to the config file.

    Returns
    -------
    None.
    """

    if os.path.exists(config_file):
        with open(config_file) as cf:
            config = json.load(cf)
        new_params = {key: eval(value) if type(
            value) == str else value for key, value in config.items()}
        rcParams.update(new_params)
    else:
        warn(
            f'{config_file} does not exist. Default settings are used.')


def scale_tick(ax: Any, data: np.ndarray, axis: str, bound: float) -> None:
    """Automatically scale ticks in GraphPad Prism style.

    Parameters  
    ----------
    ax : Any
        Axes object.
    data : numpy.ndarray
        x or y data.
    axis : str
        'x' for x axis and 'y' for y axis.
    bound : float
        Minimum bound between the maximum (minimum) data and the maximum (minimum) tick.

    Returns
    -------
    None.
    """

    assert axis in ['x', 'y'], "Argument axis should be 'x' or 'y'."

    # Generate ticks
    D_min = np.min(data)
    D_max = np.max(data)
    d_locator = ticker.MaxNLocator(4)
    d_ticks = d_locator.tick_values(D_min, D_max)

    # Expand ticks to include the minimum and maximum values
    d_gap = d_ticks[1] - d_ticks[0]
    if np.min(d_ticks) + d_gap*bound >= D_min:
        d_ticks = np.insert(d_ticks, 0, np.min(d_ticks) - d_gap)
    if np.max(d_ticks) - d_gap*bound <= D_max:
        d_ticks = np.append(d_ticks, np.max(d_ticks) + d_gap)

    # Set ticks to ensure that the limits of axes match the major ticks
    if axis == 'x':
        ax.set_xticks(d_ticks)
        if len(d_ticks) >= 4:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        else:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.set_xlim(np.min(d_ticks), np.max(d_ticks))
    if axis == 'y':
        ax.set_yticks(d_ticks)
        if len(d_ticks) >= 4:
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        else:
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.set_ylim(np.min(d_ticks), np.max(d_ticks))


def determine_unit(const: float) -> tuple[str, float]:
    """Determine the unit by a constant.

    Parameters
    ----------
    const : float
        A constant.

    Returns
    -------
    unit : str
        Unit name.
    factor : float
        Factor to convert the constant to the unit.
    """

    units = [
        ('M', 1, 2e-1),
        ('mM', 1e3, 2e-4),
        ('μM', 1e6, 2e-7),
        ('nM', 1e9, 2e-10),
        ('pM', 1e12, 2e-13)
    ]

    for suffix, multiplier, threshold in units:
        if const > threshold:
            unit, factor = suffix, multiplier
            break
    else:
        unit, factor = 'fM', 1e15

    return unit, factor
