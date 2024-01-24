import os
import sys
from warnings import warn

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.signal import medfilt

from ._utils import set_fig, scale_tick, determine_unit

global MODULE_PATH
global MODULE_DIR
global CONFIG

MODULE_PATH = os.path.abspath(__file__)
MODULE_DIR = os.path.dirname(MODULE_PATH)
CONFIG = os.path.join(MODULE_DIR, 'config/config.json')


def _langmuir(x: np.ndarray, R_m: float, K_d: float, offset: float) -> np.ndarray:
    """Define equation to fit affinity.

    Here we use the simplest form of the Langmuir 1:1 binding model.

    Parameters
    ----------
    x : numpy.ndarray
        x data.
    R_m : float
        Maximum response.
    K_d : float
        Dissociation constant.
    offset : float
        Offset.

    Returns 
    -------
    y : numpy.ndarray
        y data.
    """

    return R_m * x / (K_d + x) + offset


def _langmuir_fit(xdata: np.ndarray, ydata: np.ndarray) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """Fit data using the Langmuir 1:1 binding model 

    Fit with the initial guess and bounds.

    Parameters
    ----------
    xdata : numpy.ndarray
        x data.
    ydata : numpy.ndarray
        y data.

    Returns
    -------
    popt : numpy.ndarray
        Fitted parameters.
    pcov : numpy.ndarray
        Covariance matrix of the fitted parameters.
    """

    Y_max = np.max(ydata)
    Y_mid = np.median(ydata)

    # Initial guess
    R_m_initial = 0.75 * Y_max
    K_d_initial = xdata[np.argmin(np.abs(ydata - Y_mid))]
    offset_initial = 0
    initial_guess = [R_m_initial, K_d_initial, offset_initial]

    # Bounds
    lower_bounds = [0, 0, -np.inf]
    upper_bounds = [np.inf, np.inf, np.inf]

    # Fit affinity
    try:
        # Fit affinity
        popt, pcov = curve_fit(
            _langmuir, xdata, ydata, p0=initial_guess, bounds=(lower_bounds, upper_bounds))

        # Check the standard error of the fitted parameters
        perr = np.sqrt(np.diag(pcov))
        if perr[1] > popt[1] or np.sum(perr > popt) > 1:
            warn('Too large SEM!')
            return None, None
    except RuntimeError:
        # Some affinity data are not fitted successfully
        return None, None
    return popt, pcov


def _affinity_fit(affinity_file: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, str, float] | tuple[np.ndarray, np.ndarray, None, None, None, None, str, float]:
    """Fit affinity and generate fitted curve.

    Parameters
    ----------
    affinity_data : str
        Path to the affinity data.

    Returns
    -------
    xdata : numpy.ndarray
        Original x data of the points.
    ydata : numpy.ndarray
        Original y data of the points.
    xfit : numpy.ndarray
        x data of the fitted curve.
    yfit : numpy.ndarray
        y data of the fitted curve.
    K_d_mean : float
        Mean of the dissociation constant.
    K_d_std : float
        Standard deviation of the dissociation constant.
    unit : str
        Unit of the dissociation constant.
    factor : float
        Factor to convert the dissociation constant to the unit.
    """

    # Read the first two data frames (concentration and response) of the affinity data
    df = pd.read_csv(affinity_file, sep='\t').iloc[:, :2].dropna()
    xdata: np.ndarray = np.array(df.iloc[:, 0].values)
    ydata: np.ndarray = np.array(df.iloc[:, 1].values)
    unit, factor = determine_unit(xdata.mean())

    # Fit affinity
    popt, pcov = _langmuir_fit(xdata, ydata)
    if popt is None or pcov is None:
        warn(
            f'{os.path.basename(affinity_file)} is not fitted successfully.')
        return xdata*factor, ydata, None, None, None, None, unit, factor

    # Generate fitted curve
    xfit = np.linspace(np.min(xdata), np.max(xdata), 1000)
    yfit = _langmuir(xfit, *popt)
    K_d_mean = popt[1]
    K_d_std = np.sqrt(pcov[1, 1])

    # Check the dissociation constant
    if K_d_mean > 1e-4:
        warn(
            f'{os.path.basename(affinity_file)} has too large dissociation constant.')
    if K_d_mean < np.min(xdata) or K_d_mean > (np.max(xdata) + np.min(xdata)) / 2:
        warn(
            f'{os.path.basename(affinity_file)} has dissociation constant out of confident interval.')

    return xdata*factor, ydata, xfit*factor, yfit, K_d_mean*factor, K_d_std*factor, unit, factor


def affinity_plot(
    input_file: str, output_file: str, config_file: str = CONFIG,
    fig_size: tuple = (4.5, 3), dpi: float = 300, marker_color: str = '#000000', marker_size: float = 15, line_color: str = '#FF0000', line_width: float = 2, bound: float = 0.5
) -> None:
    """Plot affinity data.

    Parameters
    ----------
    input_file : str
        Path to the affinity data.
    output_file : str
        Path to the output figure.
    config_file : str, optional
        Path to the config file. The default is './config/config.json'.
    fig_size : tuple, optional
        Figure size. The default is (4.5, 3).
    dpi : float, optional
        Figure dpi. The default is 300.
    marker_color : tuple, optional
        Color of the data points. The default is (0, 0, 0).
    marker_size : float, optional
        Size of the data points. The default is 5.
    line_color : tuple, optional
        Color of the fitted curve. The default is (255, 0, 0).
    line_width : float, optional
        Line width of the fitted curve. The default is 2.
    bound : float, optional
        Minimum bound between the maximum (minimum) data and the maximum (minimum) tick. The default is 0.5.
    """

    # Figure Settings
    set_fig(config_file)

    # Initializing
    fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)
    xdata, ydata, xfit, yfit, K_d_mean, K_d_std, unit, factor = _affinity_fit(
        input_file)

    # Plot data points
    ax.scatter(xdata, ydata, color=marker_color, s=marker_size)

    # Auto scale ticks of x and y axes
    scale_tick(ax, xdata, 'x', bound)
    scale_tick(ax, ydata, 'y', bound)

    # Plot fitted curve
    if xfit is None or yfit is None or K_d_mean is None or K_d_std is None:  # If fitted unsuccessfully
        ax.text(0.95, 0.05, r'Fitted unsuccessfully',
                transform=ax.transAxes, va='bottom', ha='right')
    else:
        ax.plot(xfit, yfit, color=line_color, linewidth=line_width)
        ax.text(0.95, 0.05, r'$K_\mathrm{{d}} = {:.1f} \pm {:.1f}\ \mathrm{{{}}}$'.format(
            K_d_mean, K_d_std, unit), transform=ax.transAxes, va='bottom', ha='right')

    # Set labels
    ax.set_xlabel(f'Concentration/{unit}')
    ax.set_ylabel('Response/RU')

    # Save figure
    plt.savefig(output_file)
    print(f'Successfully saved {os.path.basename(output_file)}')
    plt.close(fig)


def curves_plot(
    input_file: str, output_file: str, config_file: str = CONFIG,
    affinity_file: None | str = None,
    min_filter_order: int = 7, max_filter_order: int = 31, peak_cutoff: float = 0.2,
    x_lower_bound: float = -50, x_upper_bound: float = 200, tick_num: int = 5, auto_scale_x: bool = False,
    fig_size: tuple = (4.5, 3), dpi: float = 300, custom_palette: str | list[str] = 'bright', line_width: float = 2, legend_dist: float = 0.2, bound: float = 0.5
) -> None:
    """Plot SPR curves.

    Parameters
    ----------
    input_file : str
        Path to the SPR curves data.
    output_file : str
        Path to the output figure.
    config_file : str, optional
        Path to the config file. The default is './config/config.json'.
    affinity_file : str, optional
        Path to the affinity data.
    min_filter_order : int, optional
        Minimum filter order. Should be an odd number. The default is 7.
    max_filter_order : int, optional
        Maximum filter order. Should be an odd number. The default is 31.
    peak_cutoff : float, optional
        The cutoff of the relative response of the sharp peaks to the maximum response. The default is 0.2.
    x_lower_bound : float, optional
        Lower bound of the x axis. The default is -50.
    x_upper_bound : float, optional
        Upper bound of the x axis. The default is 200.
    tick_num : int, optional
        Number of ticks. The default is 5.
    auto_scale_x : bool, optional
        Whether to auto scale the x axis. The default is False.
    fig_size : tuple, optional
        Figure size. The default is (4.5, 3).
    dpi : float, optional
        Figure dpi. The default is 300.
    custom_palette : str or list, optional
        Custom palette. The default is 'bright'.
    line_width : float, optional
        Line width of the SPR curves. The default is 2.
    legend_dist : float, optional
        Distance between the legend and plot. The default is 0.2.
    bound : float, optional
        Minimum bound between the maximum (minimum) data and the maximum (minimum) tick. The default is 0.5.

    Returns
    -------
    None.
    """

    # Figure Settings
    set_fig(config_file)

    # Initializing
    fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)
    try:
        df = pd.read_csv(input_file, sep='\t', encoding='utf-8').dropna()
    except UnicodeDecodeError:
        # Some files are not encoded in UTF-8
        df = pd.read_csv(input_file, sep='\t', encoding='ISO-8859-1').dropna()
    df_x = df.iloc[:, 0::2]
    df_y = df.iloc[:, 1::2]

    # Plot SPR curves
    assert max_filter_order % 2 == 1 and max_filter_order > min_filter_order, 'Argument max_filter_order should be an odd number and greater than min_filter_order.'
    conc: list[str] = [' '.join(i.split()[-2:])[:-2] for i in df_x.columns]
    try:
        # Filter the sharp peaks using median filter based on the moving median
        # Also try scipy.signal.butter and scipy.signal.filtfilt
        ydata_max = np.max(medfilt(df_y.values.T.flatten(), max_filter_order))
        # Used to auto scale ticks
        xdata_flat = np.array([])
        ydata_flat = np.array([])
        for i, c in enumerate(conc):
            xfiltered: np.ndarray = np.array(df_x.iloc[:, i].values)
            ytofilter: np.ndarray = np.array(df_y.iloc[:, i].values)
            yfiltered = medfilt(ytofilter, min_filter_order)

            # Determine the filter order
            for filter_order in range(min_filter_order+2, max_filter_order+2, 2):
                if np.max(np.abs(yfiltered[np.abs(xfiltered) < 1])) < peak_cutoff*ydata_max:
                    break
                yfiltered = medfilt(ytofilter, filter_order)
            else:
                warn(
                    f'Concentration {c} in {os.path.basename(output_file)} has reached the max filter order.')

            # Used to auto scale ticks
            xdata_flat = np.append(xdata_flat, xfiltered)
            ydata_flat = np.append(ydata_flat, yfiltered)

            # Plot one SPR curve each time
            ax.plot(xfiltered, yfiltered, label=c, linewidth=line_width)
    except ValueError:
        print(
            f'ValueError: {os.path.basename(output_file)} contains invalid concentration values. SPR curves data should be processed again.')
        return

    # Set lines color
    sns.set_palette(custom_palette)

    # Set legend and labels
    ax.legend(bbox_to_anchor=(1+legend_dist, 1))
    ax.set_xlabel('Time/s')
    ax.set_ylabel('Response/RU')

    # Scale ticks
    if auto_scale_x:
        scale_tick(ax, xdata_flat, 'x', bound)
        scale_tick(ax, ydata_flat, 'y', bound)
    else:
        # Scale y ticks
        scale_tick(ax, ydata_flat[(xdata_flat > x_lower_bound) & (
            xdata_flat < x_upper_bound)], 'y', bound=bound)

        # Scale x ticks manually
        ax.set_xlim(x_lower_bound, x_upper_bound)
        ax.set_xticks(np.linspace(x_lower_bound, x_upper_bound, tick_num + 1))
        if tick_num >= 4:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        else:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))

    # Set title
    if affinity_file is None:
        ax.set_title(f'Compound {os.path.basename(output_file)[:-5]}', y=1.05)
    else:
        K_d_mean, K_d_std, unit, factor = _affinity_fit(affinity_file)[4:]
        if K_d_mean is None or K_d_std is None:
            ax.set_title(
                f'Compound {os.path.basename(output_file)[:-5]}: Fitted unsuccessfully', y=1.05)
        else:
            ax.set_title(f'Compound {os.path.basename(output_file)[:-5]}: ' +
                         r'$K_\mathrm{{d}} = {:.1f} \pm {:.1f}\ \mathrm{{{}}}$'.format(K_d_mean, K_d_std, unit), y=1.05)

    # Save figure
    plt.savefig(output_file)
    print(f'Successfully saved {os.path.basename(output_file)}')
    plt.close(fig)


def affinity(input_dir: str, config_file: str = CONFIG, **kwargs) -> None:
    """Batch process affinity data.

    Parameters
    ----------
    input_dir : str
        Path to the folder containing the affinity data.
    config_file : str, optional
        Path to the config file. The default is './config/config.json'.
    **kwargs : dict
        Keyword arguments for affinity_plot.

    Returns
    -------
    None.
    """

    for file in os.listdir(input_dir):
        if file.endswith('a.txt'):
            affinity_plot(os.path.join(input_dir, file), os.path.join(
                input_dir, file[:-4] + '.png'), config_file, **kwargs)


def curves(input_dir: str, config_file: str = CONFIG, detect_affinity_file: bool = False, **kwargs) -> None:
    """Batch process SPR curves.

    Parameters
    ----------
    input_dir : str
        Path to the folder containing the SPR curves data.
    config_file : str, optional
        Path to the config file. The default is './config/config.json'.
    affinity_mode : bool, optional
        Whether to plot affinity. The default is False.
    """

    if detect_affinity_file == True:
        for file in os.listdir(input_dir):
            if file.endswith('c.txt'):
                if os.path.exists(os.path.join(input_dir, file[:-5] + 'a.txt')):
                    kwargs['affinity_file'] = os.path.join(
                        input_dir, file[:-5] + 'a.txt')  # Use affinity file with the same name
                else:
                    warn(
                        f'{os.path.basename(file)} does not have affinity file. Parsed affinity file will be used if exists.')
                try:
                    curves_plot(os.path.join(input_dir, file), os.path.join(
                        input_dir, file[:-4] + '.png'), config_file, **kwargs)
                except:
                    warn(
                        f'Unknown errors occur. {os.path.basename(file)} is not plotted successfully and skipped.')
    if detect_affinity_file == False:
        for file in os.listdir(input_dir):
            if file.endswith('c.txt'):
                try:
                    curves_plot(os.path.join(input_dir, file), os.path.join(
                        input_dir, file[:-4] + '.png'), config_file, **kwargs)
                except:
                    warn(
                        f'Unknown errors occur. {os.path.basename(file)} is not plotted successfully and skipped.')


if __name__ == '__main__':
    # Enter the path to the folder containing the affinity data and SPR curves.
    input_dir = sys.argv[1]
    affinity(input_dir)
    curves(input_dir)
