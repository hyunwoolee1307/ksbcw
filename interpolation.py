import numpy as np
import scipy as sp
import pandas as pd

def depth_interp_spline(df: pd.DataFrame, variable: str = 'temperature') -> pd.Series:
    """
    Interpolate a variable to standard depth levels using a cubic spline.

    Extrapolates by holding the first and last observed values constant.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'depth' and the specified variable columns for a single profile.
    variable : str, optional
        The name of the column to be interpolated, by default 'temperature'.

    Returns
    -------
    pd.Series
        A pandas Series with the interpolated variable, indexed by standard depth levels.
    """
    # Ensure the dataframe is sorted by depth and has unique depth values
    df = df.sort_values(by='depth').drop_duplicates(subset='depth').reset_index(drop=True)

    if df.shape[0] < 2:
        raise ValueError("Input DataFrame must have at least 2 unique depth points for interpolation.")

    standard_levels = np.array([
        0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500
    ])

    # Initialize the final Series with the correct index and name
    result_series = pd.Series(index=standard_levels, dtype=float, name=variable)

    min_depth, max_depth = df['depth'].min(), df['depth'].max()
    first_val, last_val = df[variable].iloc[0], 0

    # --- Define masks for each region ---
    # 1. Standard levels requiring interpolation (within the data range)
    interp_mask = (standard_levels >= min_depth) & (standard_levels <= max_depth)
    interp_targets = standard_levels[interp_mask]

    # 2. Standard levels requiring extrapolation (shallower than data)
    extrap_shallow_mask = standard_levels < min_depth

    # 3. Standard levels requiring extrapolation (deeper than data)
    extrap_deep_mask = standard_levels > max_depth

    # --- Perform interpolation and extrapolation ---
    # Interpolate for levels within the observed depth range
    if interp_targets.size > 0:
        # CubicSpline requires at least 2 points
        cs = sp.interpolate.CubicSpline(df['depth'], df[variable])
        result_series.loc[interp_mask] = cs(interp_targets)

    # Extrapolate for shallower levels by holding the first value constant
    result_series.loc[extrap_shallow_mask] = first_val

    # Extrapolate for deeper levels by holding the last value constant
    result_series.loc[extrap_deep_mask] = last_val

    return result_series
