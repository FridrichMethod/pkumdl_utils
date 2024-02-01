"""Utility functions for the package.

Functions
---------
set_fig(config_file: str) -> None
    Set figure style.
scale_tick(ax: Any, data: np.ndarray, axis: str, bound: float) -> None
    Automatically scale ticks in GraphPad Prism style.
determine_unit(const: float) -> tuple[str, float]
    Determine the unit by a constant.
"""

import os
import json
from io import StringIO
from typing import Any
from warnings import warn

import numpy as np
import pandas as pd
from matplotlib import rcParams, ticker

units = {
    "M": (1, 2e-1),
    "mM": (1e-3, 2e-4),
    "μM": (1e-6, 2e-7),
    "nM": (1e-9, 2e-10),
    "pM": (1e-12, 2e-13),
}


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
        with open(config_file, encoding="utf-8") as cf:
            config = json.load(cf)
        rcParams.update(config)
    else:
        warn(f"{config_file} does not exist. Default settings are used.")


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

    if axis not in ("x", "y"):
        raise ValueError('axis must be "x" or "y".')

    # Generate ticks
    D_min = np.min(data)
    D_max = np.max(data)
    d_locator = ticker.MaxNLocator(4)
    d_ticks = d_locator.tick_values(D_min, D_max)

    # Expand ticks to include the minimum and maximum values
    d_gap = d_ticks[1] - d_ticks[0]
    if np.min(d_ticks) + d_gap * bound >= D_min:
        d_ticks = np.insert(d_ticks, 0, np.min(d_ticks) - d_gap)
    if np.max(d_ticks) - d_gap * bound <= D_max:
        d_ticks = np.append(d_ticks, np.max(d_ticks) + d_gap)

    # Set ticks to ensure that the limits of axes match the major ticks
    if axis == "x":
        ax.set_xticks(d_ticks)
        if len(d_ticks) >= 4:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        else:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.set_xlim(np.min(d_ticks), np.max(d_ticks))
    if axis == "y":
        ax.set_yticks(d_ticks)
        if len(d_ticks) >= 4:
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        else:
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.set_ylim(np.min(d_ticks), np.max(d_ticks))


def read_curves(input_file: str) -> pd.DataFrame:
    """Read curves from a file encoded in UTF-8 or ISO-8859-1.

    Parameters
    ----------
    input_file : str
        Path to the file.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the curves.
    """

    try:
        df = pd.read_csv(input_file, sep="\t", encoding="utf-8").dropna()
    except UnicodeDecodeError:
        # Some files are not encoded in UTF-8 but ISO-8859-1
        warn(f"{input_file} is not encoded in UTF-8.")
        # Replace the unknown characters to the replacement character
        with open(input_file, "r", encoding="iso-8859-1", errors="replace") as f:
            curves = f.read()
        curves_utf8 = (
            curves.encode("iso-8859-1")
            .decode("utf-8", errors="replace")
            .replace("\ufffd", "μ")
        )
        df = pd.read_csv(StringIO(curves_utf8), sep="\t").dropna()

    return df


def determine_unit(quant: float) -> tuple[str, float]:
    """Determine the unit by a constant.

    Parameters
    ----------
    quant : float
        A constant.

    Returns
    -------
    unit : str
        Unit name.
    factor : float
        Factor to convert the constant to the unit.
    """

    for unit, (scale, bound) in units.items():
        if quant >= bound:
            break
    else:
        unit, scale = "fM", 1e-15

    factor = 1 / scale

    return unit, factor


def transform_unit(quant_unit: str) -> str:
    """Transform the unit of a constant.

    Parameters
    ----------
    quant_unit : str
        A constant with old unit.

    Returns
    -------
    str
        A constant with new unit.
    """

    quant, old_unit = quant_unit.split()
    quant = float(quant)
    try:
        quant *= units[old_unit][0]
    except KeyError:
        print(f"Unit {old_unit} is not supported.")
    new_unit, factor = determine_unit(quant)
    quant *= factor

    return f"{quant:.2f} {new_unit}" if quant <= 2 else f"{quant:.1f} {new_unit}"
