"""SPR data analysis and plotting.

Example
-------
>>> from datanal import spr
>>> spr.affinity(r"path/to/affinity/data")
>>> spr.curves(r"path/to/curves/data")

Notes
-----
The affinity data should be in the following format (directly export from Biacore T200):

    title (conc.)     title (resp.)
    Concentration1    Response1
    Concentration2    Response2
    Concentration3    Response3
    ...

The SPR curves data should be in the following format (directly export from Biacore T200):

    title1   title1       title2   title2       title3   title3       ...
    Time1    Response1    Time2    Response2    Time3    Response3    ...
    Time1    Response1    Time2    Response2    Time3    Response3    ...
    Time1    Response1    Time2    Response2    Time3    Response3    ...
    ...

The affinity file should end with "a.txt"; the curves file should end with "c.txt".
"""

import os
import argparse
from io import StringIO
from typing import Any
from warnings import warn

import numpy as np
import pandas as pd

# import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.signal import medfilt

# from scipy.signal import savgol_filter


from ._utils import auto_style, auto_ticks, auto_units


class SPR:
    """Surface plasmon resonance (SPR) data analysis.

    Methods
    -------
    langmuir
        The Langmuir isotherm equation.
    langmuir_fit
        Fit the Langmuir isotherm equation to the affinity data.
    """

    @staticmethod
    def langmuir(
        x: np.ndarray,
        R_m: float,
        K_d: float,
        offset: float,
    ) -> np.ndarray:
        """The Langmuir isotherm equation.

        Parameters
        ----------
        x : np.ndarray
            The concentration of the ligand.
        R_m : float
            The maximal response.
        K_d : float
            The dissociation constant.
        offset : float
            The offset.

        Returns
        -------
        np.ndarray
            The response.
        """

        return R_m * x / (K_d + x) + offset

    @staticmethod
    def langmuir_fit(
        x_data: np.ndarray, y_data: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
        """Fit the Langmuir isotherm equation to the affinity data.

        Parameters
        ----------
        x_data : np.ndarray
            The concentration of the ligand.
        y_data : np.ndarray
            The response.

        Returns
        -------
        tuple[np.ndarray, np.ndarray] | tuple[None, None]
            The fitted parameters and the covariance matrix.
        """

        y_max = np.max(y_data)
        y_mid = np.median(y_data)

        # Initial guess
        R_m_initial = 0.75 * y_max
        K_d_initial = x_data[np.argmin(np.abs(y_data - y_mid))]
        offset_initial = 0
        initial_guess = [R_m_initial, K_d_initial, offset_initial]

        # Bounds
        lower_bounds = [0, 0, -np.inf]
        upper_bounds = [np.inf, np.inf, np.inf]

        # Fit affinity
        try:
            popt, pcov, *_ = curve_fit(
                SPR.langmuir,
                x_data,
                y_data,
                p0=initial_guess,
                bounds=(lower_bounds, upper_bounds),
            )
        except RuntimeError:
            # Some affinity data are not fitted successfully
            return None, None

        # Check the standard error of the fitted parameters
        perr = np.sqrt(np.diag(pcov))
        if perr[1] > popt[1] or np.sum(perr > popt) > 1:
            warn("Too large SEM for affinity fitted!")
            return None, None

        return popt, pcov


class SPRAffinity:
    """Surface plasmon resonance (SPR) affinity data analysis.

    Attributes
    ----------
    affinity_file : str
        The affinity data file.

    Methods
    -------
    __repr__
        Return the name of the affinity file.
    __len__
        Return the number of data points in the affinity file.
    affinity_data
        Read the first two data frames (concentration and response) of the affinity file.
    affinity_fit
        Fit the affinity data to the Langmuir isotherm equation.
    affinity_plot
        Plot the affinity data and the fitted curve.
    """

    def __init__(
        self,
        affinity_file: str,
        *,
        rc_mplstyle: dict[str, Any] | None = None,
        fname_mplstyle: str | None = None,
        palette_snsstyle: str | list[str] | None = None,
    ) -> None:
        """Initialize the matplotlib and seaborn styles and the affinity file.

        Parameters
        ----------
        affinity_file : str
            The affinity data file.
        rc_mplstyle : dict[str, Any] | None
            The matplotlib style for rcParams.
        fname_mplstyle : str | None
            The matplotlib style for figure.
        palette_snsstyle : str | list[str] | None
            The seaborn style for palette.

        Returns
        -------
        None
        """

        auto_style(rc_mplstyle, fname_mplstyle, palette_snsstyle)

        if os.path.exists(affinity_file):
            self._affinity_file = affinity_file
        else:
            raise FileNotFoundError(f"{affinity_file} does not exist.")

    def __repr__(self) -> str:
        """Return the name of the affinity file."""

        return os.path.basename(self._affinity_file)

    def __len__(self) -> int:
        """Return the number of data points in the affinity file."""

        return len(self.affinity_data)

    @property
    def affinity_data(self) -> pd.DataFrame:
        """Read the first two data frames (concentration and response) of the affinity file."""

        return (
            pd.read_csv(self._affinity_file, sep="\t", encoding="utf-8")
            .iloc[:, :2]
            .dropna()
        )

    def affinity_fit(
        self,
    ) -> (
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, str, float]
        | tuple[np.ndarray, np.ndarray, None, None, None, None, str, float]
    ):
        """Fit the affinity data to the Langmuir isotherm equation.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, str, float] | tuple[np.ndarray, np.ndarray, None, None, None, None, str, float]
            The concentration, response, fitted concentration, fitted response, dissociation constant, standard error of the dissociation constant, unit, and factor.
        """

        # Read affinity data
        df = self.affinity_data
        x_data: np.ndarray = np.array(df.iloc[:, 0].values)
        y_data: np.ndarray = np.array(df.iloc[:, 1].values)
        unit, factor = auto_units(x_data.mean())

        # Fit affinity
        fit = SPR.langmuir_fit(x_data, y_data)
        if any(p is None for p in fit):
            warn(f"{os.path.basename(self._affinity_file)} is not fitted successfully.")
            return x_data * factor, y_data, None, None, None, None, unit, factor

        popt, pcov = fit

        # Generate fitted curve
        assert popt is not None and pcov is not None
        x_fit = np.linspace(np.min(x_data), np.max(x_data), 1000)
        y_fit = SPR.langmuir(x_fit, *popt)
        K_d_mean = popt[1]
        K_d_std = np.sqrt(pcov[1, 1])

        # Check the dissociation constant
        if K_d_mean > 1e-4:
            warn(
                f"{os.path.basename(self._affinity_file)} has too large dissociation constant."
            )
        if (
            K_d_mean < np.min(x_data)
            or K_d_mean > (np.max(x_data) + np.min(x_data)) / 2
        ):
            warn(
                f"{os.path.basename(self._affinity_file)} has dissociation constant out of confident interval."
            )

        return (
            x_data * factor,
            y_data,
            x_fit * factor,
            y_fit,
            K_d_mean * factor,
            K_d_std * factor,
            unit,
            factor,
        )

    def affinity_plot(
        self,
        *,
        fig_size: tuple[float, float] = (4.5, 3),
        dpi: float = 300,
        marker_color: str = "#000000",
        marker_size: float = 16,
        line_color: str = "#FF0000",
        line_width: float = 2,
    ) -> None:
        """Plot the affinity data and the fitted curve.

        Parameters
        ----------
        fig_size : tuple[float, float], optional
            The size of the figure.
        dpi : float, optional
            The resolution of the figure.
        marker_color : str, optional
            The color of the data points.
        marker_size : float, optional
            The size of the data points.
        line_color : str, optional
            The color of the fitted curve.
        line_width : float, optional
            The width of the fitted curve.

        Returns
        -------
        None
        """

        # Initializing
        fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)
        fit = self.affinity_fit()
        x_data, y_data, x_fit, y_fit, K_d_mean, K_d_std, unit, _ = fit

        # Auto scale ticks
        auto_ticks(ax)

        # Plot data points
        ax.scatter(x_data, y_data, c=marker_color, s=marker_size)

        # Plot fitted curve
        if any(p is None for p in fit):  # If fitted unsuccessfully
            ax.text(
                0.95,
                0.05,
                "Fitted unsuccessfully",
                transform=ax.transAxes,
                va="bottom",
                ha="right",
            )
        else:
            ax.plot(x_fit, y_fit, color=line_color, linewidth=line_width)
            ax.text(
                0.95,
                0.05,
                rf"$K_\mathrm{{d}} = {K_d_mean:.1f} \pm {K_d_std:.1f}\ \mathrm{{{unit}}}$",
                transform=ax.transAxes,
                va="bottom",
                ha="right",
            )

        # Set labels
        ax.set_xlabel(f"Concentration/{unit}")
        ax.set_ylabel("Response/RU")

        # Save figure
        output_file = os.path.splitext(self._affinity_file)[0] + ".png"
        plt.savefig(output_file)
        print(f"Successfully saved {os.path.basename(output_file)}")

        plt.close(fig)


class SPRCurves:
    """Surface plasmon resonance (SPR) curves data analysis.

    Attributes
    ----------
    curves_file : str
        The curves data file.

    Methods
    -------
    __repr__
        Return the name of the curves file.
    __len__
        Return the number of curves in the curves file.
    curves_data
        Read curves from the curves file and change the encoding to UTF-8.
    concentrations
        Extract the concentrations from the curves file.
    curves_plot
        Plot the SPR curves.
    """

    def __init__(
        self,
        curves_file: str,
        affinity_file: str | None = None,
        *,
        rc_mplstyle: dict[str, Any] | None = None,
        fname_mplstyle: str | None = None,
        palette_snsstyle: str | list[str] | None = None,
    ) -> None:
        """Initialize the matplotlib and seaborn styles and the curves file.

        Parameters
        ----------
        curves_file : str
            The curves data file.
        affinity_file : str | None
            The affinity data file.
        rc_mplstyle : dict[str, Any] | None
            The matplotlib style for rcParams.
        fname_mplstyle : str | None
            The matplotlib style for figure.
        palette_snsstyle : str | list[str] | None
            The seaborn style for palette.

        Returns
        -------
        None
        """

        if affinity_file is None:
            self.affinity_fit_ = None
        else:
            self.affinity_fit_ = SPRAffinity(affinity_file).affinity_fit()

        auto_style(rc_mplstyle, fname_mplstyle, palette_snsstyle)
        if os.path.exists(curves_file):
            self._curves_file = curves_file
        else:
            raise FileNotFoundError(f"{curves_file} does not exist.")

    def __repr__(self) -> str:
        """Return the name of the curves file."""

        return os.path.basename(self._curves_file)

    def __len__(self) -> int:
        """Return the number of curves in the curves file."""

        return len(self.curves_data)

    @property
    def curves_data(self) -> pd.DataFrame:
        """Read curves from the curves file and change the encoding to UTF-8."""

        try:
            df = pd.read_csv(self._curves_file, sep="\t", encoding="utf-8").dropna()
        except UnicodeDecodeError:
            # Some files are not encoded in UTF-8 but ISO-8859-1
            warn(f"{self._curves_file} is not encoded in UTF-8. Try ISO-8859-1...")
            # Replace the unknown characters to the replacement character
            with open(self._curves_file, "r", encoding="iso-8859-1") as f:
                curves_iso = f.read()
            curves_utf8 = (
                curves_iso.encode("iso-8859-1")
                .decode("utf-8", errors="replace")
                .replace("\ufffd", "μ")
            )
            df = pd.read_csv(StringIO(curves_utf8), sep="\t").dropna()

        return df

    @property
    def concentrations(self) -> list[tuple[float, str]]:
        """Extract the concentrations from the curves file."""

        def _extract_conc(title):
            quant, old_unit = title.strip()[:-2].split()[-2:]
            quant = float(quant)
            new_unit, factor = auto_units(quant, old_unit)
            quant *= factor
            return quant, new_unit

        df = self.curves_data
        conc = list(map(_extract_conc, df.iloc[:, 0::2].columns))

        return conc

    def curves_plot(
        self,
        *,
        x_lower_bound: float = -50,
        x_upper_bound: float = 200,
        min_filter_order: int = 7,
        max_filter_order: int = 31,
        peak_cutoff: float = 0.2,
        fig_size: tuple[float, float] = (4.5, 3),
        dpi: float = 300,
        line_width: float = 2,
        legend_distance: float = 0.33,
    ) -> None:
        """Plot the SPR curves.

        Parameters
        ----------
        x_lower_bound : float, optional
            The lower bound of the x-axis.
        x_upper_bound : float, optional
            The upper bound of the x-axis.
        min_filter_order : int, optional
            The minimum filter order for the signal filter.
        max_filter_order : int, optional
            The maximum filter order for the signal filter.
        peak_cutoff : float, optional
            The cutoff for the sharp peaks.
        fig_size : tuple[float, float], optional
            The size of the figure.
        dpi : float, optional
            The resolution of the figure.
        line_width : float, optional
            The width of the SPR curves.
        legend_distance : float, optional
            The distance of the legend from the upper right corner.

        Returns
        -------
        None
        """

        # Initializing
        fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)
        df = self.curves_data
        df_x = df.iloc[:, 0::2]
        df_y = df.iloc[:, 1::2]

        # Auto scale ticks
        auto_ticks(ax, x_lim=(x_lower_bound, x_upper_bound))

        # Plot SPR curves
        if max_filter_order % 2 == 0 or max_filter_order <= min_filter_order:
            raise ValueError(
                "max_filter_order should be an odd number and greater than or equal to min_filter_order."
            )
        conc = self.concentrations
        # Filter the sharp peaks using signal filter
        y_data_max = medfilt(df_y.values.T.flatten(), max_filter_order).max()

        for i, (quant, old_unit) in enumerate(conc):
            new_unit, factor = auto_units(quant, old_unit)
            quant *= factor
            x_data: np.ndarray = np.array(df_x.iloc[:, i].values)
            y_data: np.ndarray = np.array(df_y.iloc[:, i].values)
            mask = (x_data >= x_lower_bound) & (x_data <= x_upper_bound)
            x_between: np.ndarray = x_data[mask]
            y_between: np.ndarray = y_data[mask]
            y_filtered = medfilt(y_between, min_filter_order)
            # Determine the filter order
            for filter_order in range(min_filter_order + 2, max_filter_order + 2, 2):
                if (
                    np.max(np.abs(y_filtered[np.abs(x_between) < 1]))
                    < peak_cutoff * y_data_max
                ):
                    break
                y_filtered = medfilt(y_between, filter_order)
            else:
                warn(
                    f"Concentration {quant:.2f} {new_unit} in {os.path.basename(self._curves_file)} has reached the max filter order."
                )

            # Plot one SPR curve each time
            ax.plot(
                x_between,
                y_filtered,
                label=f"{quant:.2f} {new_unit}",
                linewidth=line_width,
            )

        # Set legend
        ax.legend(bbox_to_anchor=(1 + legend_distance, 1))

        # Set labels
        ax.set_xlabel("Time/s")
        ax.set_ylabel("Response/RU")

        # Set title
        fit = self.affinity_fit_
        if fit is not None:
            K_d_mean, K_d_std, unit, _ = fit[4:]

            compound_name = (
                os.path.splitext(os.path.basename(self._curves_file))[0][:-1]
                .replace("_", " ")
                .strip()
            )
            if any(p is None for p in fit):  # If fitted unsuccessfully
                ax.set_title(f"{compound_name}: Fitted unsuccessfully")
            else:
                ax.set_title(
                    rf"{compound_name}: $K_\mathrm{{d}} = {K_d_mean:.1f} \pm {K_d_std:.1f}\ \mathrm{{{unit}}}$"
                )

        # Save figure
        output_file = os.path.splitext(self._curves_file)[0] + ".png"
        plt.savefig(output_file)
        print(f"Successfully saved {os.path.basename(output_file)}")

        plt.close(fig)


def affinity(
    affinity_dir: str,
    *,
    rc_mplstyle: dict[str, Any] | None = None,
    fname_mplstyle: str | None = None,
    palette_snsstyle: str | list[str] | None = None,
    **kwargs: Any,
) -> None:
    """Plot the affinity data.

    Parameters
    ----------
    affinity_dir : str
        The directory containing the affinity data files.
    rc_mplstyle : dict[str, Any] | None
        The matplotlib style for rcParams.
    fname_mplstyle : str | None
        The matplotlib style for figure.
    palette_snsstyle : str | list[str] | None
        The seaborn style for palette.
    **kwargs : Any
        The keyword arguments for the affinity plot, see below.

    Keyword Arguments
    -----------------
    fig_size : tuple[float, float], optional
        The size of the figure.
    dpi : float, optional
        The resolution of the figure.
    marker_color : str, optional
        The color of the data points.
    marker_size : float, optional
        The size of the data points.
    line_color : str, optional
        The color of the fitted curve.
    line_width : float, optional
        The width of the fitted curve.

    Returns
    -------
    None
    """

    for file in os.listdir(affinity_dir):
        if file.endswith("a.txt"):
            spr_affinity = SPRAffinity(
                os.path.join(affinity_dir, file),
                rc_mplstyle=rc_mplstyle,
                fname_mplstyle=fname_mplstyle,
                palette_snsstyle=palette_snsstyle,
            )
            spr_affinity.affinity_plot(**kwargs)


def curves(
    curves_dir: str,
    *,
    rc_mplstyle: dict[str, Any] | None = None,
    fname_mplstyle: str | None = None,
    palette_snsstyle: str | list[str] | None = None,
    **kwargs: Any,
) -> None:
    """Plot the SPR curves.

    Parameters
    ----------
    curves_dir : str
        The directory containing the curves data files.
    rc_mplstyle : dict[str, Any] | None
        The matplotlib style for rcParams.
    fname_mplstyle : str | None
        The matplotlib style for figure.
    palette_snsstyle : str | list[str] | None
        The seaborn style for palette.
    **kwargs : Any
        The keyword arguments for the curves plot, see below.

    Keyword Arguments
    -----------------
    x_lower_bound : float, optional
        The lower bound of the x-axis.
    x_upper_bound : float, optional
        The upper bound of the x-axis.
    min_filter_order : int, optional
        The minimum filter order for the signal filter.
    max_filter_order : int, optional
        The maximum filter order for the signal filter.
    peak_cutoff : float, optional
        The cutoff for the sharp peaks.
    fig_size : tuple[float, float], optional
        The size of the figure.
    dpi : float, optional
        The resolution of the figure.
    line_width : float, optional
        The width of the SPR curves.
    legend_distance : float, optional
        The distance of the legend from the upper right corner.

    Returns
    -------
    None
    """

    for file in os.listdir(curves_dir):
        if file.endswith("c.txt"):
            affinity_file = os.path.join(curves_dir, file[:-5] + "a.txt")
            if os.path.exists(affinity_file):
                spr_curves = SPRCurves(
                    os.path.join(curves_dir, file),
                    affinity_file,
                    rc_mplstyle=rc_mplstyle,
                    fname_mplstyle=fname_mplstyle,
                    palette_snsstyle=palette_snsstyle,
                )
            else:
                spr_curves = SPRCurves(
                    os.path.join(curves_dir, file),
                    rc_mplstyle=rc_mplstyle,
                    fname_mplstyle=fname_mplstyle,
                    palette_snsstyle=palette_snsstyle,
                )
            spr_curves.curves_plot(**kwargs)


def main() -> None:
    """Parse the command line arguments and execute the corresponding function."""

    parser = argparse.ArgumentParser(
        description="Surface plasmon resonance (SPR) data analysis and plotting."
    )
    parser.add_argument(
        "-a",
        "--affinity",
        type=str,
        help="The directory containing the affinity data files.",
    )
    parser.add_argument(
        "-c",
        "--curves",
        type=str,
        help="The directory containing the curves data files.",
    )
    parser.add_argument(
        "-l",
        "--lower",
        type=float,
        help="The lower bound of the x-axis in curves plot.",
        default=-50,
    )
    parser.add_argument(
        "-u",
        "--upper",
        type=float,
        help="The upper bound of the x-axis in curves plot.",
        default=200,
    )
    parser.add_argument(
        "-r",
        "--rc",
        type=str,
        help="The matplotlib style file for rcParams.",
        default=None,
    )
    parser.add_argument(
        "-f",
        "--fname",
        type=str,
        help="The matplotlib style file for figure.",
        default=None,
    )
    parser.add_argument(
        "-p",
        "--palette",
        type=str,
        help="The seaborn style for palette.",
        default=None,
    )
    args = parser.parse_args()

    if args.affinity is not None:
        affinity(
            args.affinity,
            rc_mplstyle=args.rc,
            fname_mplstyle=args.fname,
            palette_snsstyle=args.palette,
        )
    if args.curves is not None:
        curves(
            args.curves,
            rc_mplstyle=args.rc,
            fname_mplstyle=args.fname,
            palette_snsstyle=args.palette,
            x_lower_bound=args.lower,
            x_upper_bound=args.upper,
        )


if __name__ == "__main__":
    main()
