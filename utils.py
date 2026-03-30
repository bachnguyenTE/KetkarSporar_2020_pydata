"""
Utility functions for loading and processing calcium imaging data
from Ketkar & Sporar 2020 (Current Biology).

Ports the key MATLAB functions to Python:
- .mat file loading and struct unwrapping
- Summary.xlsx metadata parsing
- 10Hz interpolation (matching MATLAB's interp1)
- Aggregation functions for each figure type
- Per-fly averaging (mean_cat_full)
- ROI classification by correlation
- HDF5 storage
"""

import os
import re
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import mode as scipy_mode


# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent
DATA_ROOT = REPO_ROOT / "data" / "in vivo two-photon calcium imaging"
SUMMARY_PATH = (
    REPO_ROOT
    / "code"
    / "in vivo two photon calcium imaging"
    / "Common functions"
    / "Summary.xlsx.xlsx"
)

# Column name mapping: xlsx column → internal key
_COL_MAP = {
    "fname": "fname",
    " stimbouts": "stimbouts",
    "locnum": "locnum",
    "inverted cell flag": "inverted",
    "stimulus code": "stimcode",
    "qualitative goodness (0-10)": "quality",
    "Driver": "driver",
    "Moving": "moving",
    "Layer(2/4)": "layer",
    "Stimulus Wavelength": "wavelength",
    "flyID": "flyID",
    "responsivefly": "responsivefly",
}

IRATE = 10  # interpolation rate in Hz


# ---------------------------------------------------------------------------
# .mat file loading
# ---------------------------------------------------------------------------


def load_mat_struct(filepath):
    """Load a pData .mat file and return a flat dict of fields.

    Handles the nested struct unwrapping that scipy.io.loadmat requires.

    Returns
    -------
    dict with keys:
        dRatio   : (nROIs, nFrames) float64
        fstimval : (nFrames,) uint8
        avrstimval : (nFrames,) float64
        ch3      : (nFrames,) uint8
        framerate : float
        frame_nums : (N,) array
        fstimpos1 : (nFrames,) or None
        fstimpos2 : (nFrames,) or None
    """
    from scipy.io import loadmat

    mat = loadmat(filepath, squeeze_me=False)
    strct = mat["strct"]

    result = {}
    result["dRatio"] = strct["dRatio"][0, 0].astype(np.float64)
    result["fstimval"] = strct["fstimval"][0, 0].squeeze().astype(np.float64)
    result["avrstimval"] = strct["avrstimval"][0, 0].squeeze().astype(np.float64)
    result["ch3"] = strct["ch3"][0, 0].squeeze()
    result["frame_nums"] = strct["frame_nums"][0, 0].squeeze()

    # dataID identifies the biological fly (e.g. "160727ks_fly2")
    result["dataID"] = str(strct["dataID"][0, 0].flat[0])

    # Nested xml struct
    xml = strct["xml"][0, 0]
    result["framerate"] = float(xml["framerate"][0, 0].flat[0])

    # Optional fields
    if "fstimpos1" in strct.dtype.names:
        result["fstimpos1"] = strct["fstimpos1"][0, 0].squeeze().astype(np.float64)
    else:
        result["fstimpos1"] = None
    if "fstimpos2" in strct.dtype.names:
        result["fstimpos2"] = strct["fstimpos2"][0, 0].squeeze().astype(np.float64)
    else:
        result["fstimpos2"] = None

    return result


# ---------------------------------------------------------------------------
# Summary.xlsx parsing
# ---------------------------------------------------------------------------


def load_summary(path=None):
    """Load Summary.xlsx and return a DataFrame with standardized column names.

    Parameters
    ----------
    path : str or Path, optional
        Path to the xlsx file.  Defaults to the repository's copy.

    Returns
    -------
    pd.DataFrame with columns renamed via _COL_MAP.
    """
    if path is None:
        path = SUMMARY_PATH
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        df = pd.read_excel(path)
    df = df.rename(columns=_COL_MAP)
    # Strip whitespace from fname
    df["fname"] = df["fname"].str.strip()
    return df


def get_fly_id(filename, summary_df):
    """Look up the numeric flyID for a .mat filename from the summary table."""
    basename = os.path.basename(filename).strip()
    match = summary_df.loc[summary_df["fname"] == basename, "flyID"]
    if len(match) == 0:
        return None
    return int(match.iloc[0])


def parse_layer_region(layer_str, region="AT"):
    """Parse a Layer string like 'AT1:9DE1:6' to extract ROI indices for a region.

    MATLAB's str2num('1:9') evaluates to [1,2,...,9].

    Returns
    -------
    list of int (1-based ROI indices), or empty list if region not found.
    """
    if not isinstance(layer_str, str):
        return []
    idx = layer_str.find(region)
    if idx < 0:
        return []
    # Skip past the region tag
    start = idx + len(region)
    # Read until next letter or end
    end = start
    while end < len(layer_str) and not layer_str[end].isalpha():
        end += 1
    range_str = layer_str[start:end].strip()
    if not range_str:
        return []
    # Parse MATLAB-style range "1:9" or space-separated "1 3 5"
    if ":" in range_str:
        parts = range_str.split(":")
        return list(range(int(parts[0]), int(parts[1]) + 1))
    else:
        return [int(x) for x in range_str.split()]


# ---------------------------------------------------------------------------
# 10Hz interpolation
# ---------------------------------------------------------------------------


def interpolate_to_10hz(ratio, stim, fps, avrstimval=None, fstimpos1=None, fstimpos2=None):
    """Interpolate calcium and stimulus traces to 10 Hz.

    Matches MATLAB code exactly:
        t  = [1:N] / fps           (1-based)
        it = [0.5/fps : 0.1 : (N+0.5)/fps]
        istim  = interp1(t, stim, it, 'nearest', 'extrap')
        iratio = interp1(t, ratio, it, 'linear', 'extrap')

    Parameters
    ----------
    ratio : (nFrames,) array - calcium fluorescence ratio for one ROI
    stim  : (nFrames,) array - stimulus trace
    fps   : float - original frame rate

    Returns
    -------
    dict with keys: it, iratio, istim, iavrstim, ifstimpos1, ifstimpos2
    """
    n_frames = len(ratio)
    t = np.arange(1, n_frames + 1) / fps
    it = np.arange(0.5 / fps, (n_frames + 0.5) / fps + 0.1, 0.1)
    # Clip it to avoid floating-point overshoot
    it = it[it <= (n_frames + 0.5) / fps + 1e-10]

    # Nearest-neighbor for stimulus (preserves discrete values)
    f_stim = interp1d(t, stim, kind="nearest", fill_value="extrapolate")
    istim = f_stim(it)

    # Linear for calcium (smooth)
    f_ratio = interp1d(t, ratio, kind="linear", fill_value="extrapolate")
    iratio = f_ratio(it)

    result = {"it": it, "iratio": iratio, "istim": istim}

    if avrstimval is not None:
        f_avr = interp1d(t, avrstimval, kind="nearest", fill_value="extrapolate")
        result["iavrstim"] = f_avr(it)

    if fstimpos1 is not None:
        f_p1 = interp1d(t, fstimpos1, kind="nearest", fill_value="extrapolate")
        result["ifstimpos1"] = f_p1(it)

    if fstimpos2 is not None:
        f_p2 = interp1d(t, fstimpos2, kind="nearest", fill_value="extrapolate")
        result["ifstimpos2"] = f_p2(it)

    return result


# ---------------------------------------------------------------------------
# Data loading pipeline
# ---------------------------------------------------------------------------


def load_figure_data(data_dir, summary_df=None, region=None):
    """Load all .mat files from a directory and return per-ROI data at 10 Hz.

    Parameters
    ----------
    data_dir : str or Path
        Directory containing *_pData.mat files.
    summary_df : pd.DataFrame, optional
        Summary table for flyID lookup. Loaded automatically if None.
    region : str, optional
        If given (e.g. 'AT'), only include ROIs in that region.
        Requires summary_df with 'layer' column.

    Returns
    -------
    list of dicts, one per ROI, with keys:
        iratio, istim, flyID, name, cell, ifstimpos1, ifstimpos2
    """
    if summary_df is None:
        summary_df = load_summary()

    data_dir = Path(data_dir)
    mat_files = sorted(data_dir.glob("*_pData.mat"))

    all_rois = []

    for mat_path in mat_files:
        fname = mat_path.name
        mat_data = load_mat_struct(str(mat_path))

        # Use dataID from .mat file as fly grouping key
        fly_id = mat_data["dataID"]
        n_rois = mat_data["dRatio"].shape[0]
        fps = mat_data["framerate"]

        # Determine which ROIs to include
        if region is not None:
            layer_row = summary_df.loc[summary_df["fname"] == fname, "layer"]
            if len(layer_row) == 0:
                # File not in summary — include all ROIs
                roi_indices = list(range(1, n_rois + 1))
            else:
                layer_str = str(layer_row.iloc[0])
                roi_indices = parse_layer_region(layer_str, region)
        else:
            roi_indices = list(range(1, n_rois + 1))

        for cell in roi_indices:
            if cell < 1 or cell > n_rois:
                continue

            ratio = mat_data["dRatio"][cell - 1, :]  # 1-based → 0-based
            stim = mat_data["fstimval"]

            interp_data = interpolate_to_10hz(
                ratio,
                stim,
                fps,
                avrstimval=mat_data["avrstimval"],
                fstimpos1=mat_data["fstimpos1"],
                fstimpos2=mat_data["fstimpos2"],
            )

            roi_dict = {
                "iratio": interp_data["iratio"],
                "istim": interp_data["istim"],
                "it": interp_data["it"],
                "flyID": fly_id,
                "name": fname,
                "cell": cell,
            }
            if "ifstimpos1" in interp_data:
                roi_dict["ifstimpos1"] = interp_data["ifstimpos1"]
            if "ifstimpos2" in interp_data:
                roi_dict["ifstimpos2"] = interp_data["ifstimpos2"]

            all_rois.append(roi_dict)

    return all_rois


# ---------------------------------------------------------------------------
# Cantor pairing
# ---------------------------------------------------------------------------


def cantor(a, b):
    """Cantor pairing function: maps two non-negative integers to a unique integer."""
    a, b = int(a), int(b)
    return (a + b) * (a + b + 1) // 2 + b


def cantor_x(z):
    """Inverse cantor pairing: recover first element."""
    z = int(z)
    w = int(np.floor((np.sqrt(8 * z + 1) - 1) / 2))
    t = (w * w + w) // 2
    y = z - t
    return w - y


def cantor_y(z):
    """Inverse cantor pairing: recover second element."""
    z = int(z)
    w = int(np.floor((np.sqrt(8 * z + 1) - 1) / 2))
    t = (w * w + w) // 2
    return z - t


# ---------------------------------------------------------------------------
# Per-fly averaging
# ---------------------------------------------------------------------------


def mean_cat_full(matrix, fly_ids):
    """Group ROIs by fly, compute per-fly means, then grand mean and SEM.

    Parameters
    ----------
    matrix : (nROIs, ...) array
        Data matrix. First axis is ROIs.
    fly_ids : (nROIs,) array
        Fly identifier for each ROI.

    Returns
    -------
    fly_means : (nFlies, ...) array - mean trace per fly
    grand_mean : (...) array - mean across flies
    sem : (...) array - SEM across flies
    """
    unique_flies = np.unique(fly_ids)
    fly_means = []
    for fid in unique_flies:
        mask = np.array(fly_ids) == fid
        fly_mean = np.nanmean(matrix[mask], axis=0)
        fly_means.append(fly_mean)
    fly_means = np.array(fly_means)
    n_flies = len(unique_flies)
    grand_mean = np.nanmean(fly_means, axis=0)
    sem = np.nanstd(fly_means, axis=0, ddof=0) / np.sqrt(n_flies)
    return fly_means, grand_mean, sem


# ---------------------------------------------------------------------------
# Aggregation functions
# ---------------------------------------------------------------------------


def aggregate_fff(roi_data_list):
    """Aggregate full-field flash data (Figure 2B).

    Triggers on falling edges (dark onset), computes dF/F = iratio/mean(iratio) - 1,
    averages across stimulus repetitions per ROI.

    Parameters
    ----------
    roi_data_list : list of dicts from load_figure_data

    Returns
    -------
    dict with keys: rats (nROIs, dur), stims (nROIs, dur), flyID, name
    """
    # Determine epoch length from first ROI
    dxi_all = []
    for roi in roi_data_list:
        s = roi["istim"]
        xi = np.where(np.diff(s) < 0)[0] + 1  # trigger on going to dark
        if len(xi) > 1:
            dxi_all.append(np.mean(np.diff(xi)))
    epochlength = int(round(np.mean(dxi_all)))
    dur = int(epochlength * 1.3)

    n_rois = len(roi_data_list)
    rats = np.zeros((n_rois, dur))
    stims = np.zeros((n_rois, dur))
    fly_ids = []
    names = []

    for ii, roi in enumerate(roi_data_list):
        s = roi["istim"]
        iratio = roi["iratio"]

        # dF/F with whole-trace mean baseline
        iratio_df = iratio / np.mean(iratio) - 1

        # Trigger on falling edges (exclude those too close to end)
        xi = np.where(np.diff(s[: len(s) - dur]) < 0)[0] + 1

        if len(xi) == 0:
            continue

        mr = np.zeros((len(xi), dur))
        ms = np.zeros((len(xi), dur))

        for jj, trigger in enumerate(xi):
            mr[jj, :] = iratio_df[trigger : trigger + dur]
            ms[jj, :] = s[trigger : trigger + dur]

        rats[ii, :] = np.mean(mr[: len(xi)], axis=0)
        stims[ii, :] = np.mean(ms[: len(xi)], axis=0)
        fly_ids.append(roi["flyID"])
        names.append(roi["name"])

    fly_ids = np.array(fly_ids)
    return {"rats": rats, "stims": stims, "flyID": fly_ids, "name": names}


def aggregate_fff60s(roi_data_list):
    """Aggregate 60s full-field flash data (Figure 2C).

    Triggers on light-ON transitions, computes dF/F = iratio/mean(iratio) - 1.

    Parameters
    ----------
    roi_data_list : list of dicts from load_figure_data

    Returns
    -------
    dict with keys: rats (nROIs, dur), stims (nROIs, dur), flyID, name
    """
    # Determine epoch length
    dxi_all = []
    for roi in roi_data_list:
        s = roi["istim"]
        # Include position 0 as a trigger, then find falling edges
        xi_fall = np.where(np.diff(s) < 0)[0] + 1
        xi = np.concatenate([[0], xi_fall]) if len(xi_fall) > 0 else np.array([0])
        if len(xi) > 1:
            dxi_all.append(np.mean(np.diff(xi)))

    if len(dxi_all) == 0:
        epochlength = 600  # fallback: 60s at 10Hz
    else:
        epochlength = int(round(np.mean(dxi_all)))
    dur = epochlength

    n_rois = len(roi_data_list)
    rats = np.zeros((n_rois, dur))
    stims = np.zeros((n_rois, dur))
    fly_ids = []
    names = []

    for ii, roi in enumerate(roi_data_list):
        s = roi["istim"]
        iratio = roi["iratio"]

        # dF/F with whole-trace mean baseline
        iratio_df = iratio / np.mean(iratio) - 1

        # Trigger on falling edges
        xi_fall = np.where(np.diff(s[: len(s) - dur]) < 0)[0] + 1
        xi = np.concatenate([[0], xi_fall]) if len(xi_fall) > 0 else np.array([0])

        if len(xi) == 0:
            continue

        mr = np.zeros((len(xi), dur))
        ms = np.zeros((len(xi), dur))

        for jj, trigger in enumerate(xi):
            end_idx = trigger + dur
            if end_idx > len(iratio_df):
                continue
            mr[jj, :] = iratio_df[trigger:end_idx]
            ms[jj, :] = s[trigger:end_idx]

        rats[ii, :] = np.mean(mr[: len(xi)], axis=0)
        stims[ii, :] = np.mean(ms[: len(xi)], axis=0)
        fly_ids.append(roi["flyID"])
        names.append(roi["name"])

    fly_ids = np.array(fly_ids)
    return {"rats": rats, "stims": stims, "flyID": fly_ids, "name": names}


def aggregate_standing_stripe(roi_data_list):
    """Aggregate contrast step data with grey baseline (Figure 2D-G).

    10 contrast epochs, dF/F normalized to grey baseline.

    Parameters
    ----------
    roi_data_list : list of dicts from load_figure_data

    Returns
    -------
    dict with keys:
        rats (nROIs, n_epoch, full_dur), stimstruct (full_dur,),
        flyID, name
    """
    # Determine number of epochs and durations from first ROI
    s0 = roi_data_list[0]["istim"]
    unique_vals = np.unique(s0)
    n_epoch = int(np.sum(unique_vals > 0))

    # Compute epoch durations across all ROIs
    epoch_durs = []
    zero_durs = []
    for roi in roi_data_list:
        s = roi["istim"]
        ds = np.diff(s)
        # Falling edges: epoch→grey
        binds = np.where(ds < 0)[0] + 1
        # Rising edges: grey→epoch
        einds = np.where(ds > 0)[0] + 1
        length = min(len(binds), len(einds))
        if length < 2:
            continue
        # Epoch duration: from rising to falling
        dinds = binds[:length] - einds[:length]
        epoch_durs.append(np.mean(dinds))
        # Zero interleave: from falling to next rising
        zinds = einds[1:length] - binds[: length - 1]
        zero_durs.append(np.mean(zinds))

    epoch_dur = int(round(np.mean(epoch_durs)))
    zero_dur = int(round(np.mean(zero_durs)))
    full_dur = zero_dur + epoch_dur + zero_dur

    # Stimulus pattern template
    stimstruct = np.zeros(full_dur)
    stimstruct[zero_dur : zero_dur + epoch_dur] = 1.0

    n_rois = len(roi_data_list)
    rats = np.zeros((n_rois, n_epoch, full_dur))
    fly_ids = []
    names = []

    for ii, roi in enumerate(roi_data_list):
        iratio = roi["iratio"]
        s = roi["istim"]

        # Compute grey baseline: average of grey (stimulus==0) periods, last half
        is_grey = s == 0
        begin = np.diff(is_grey.astype(float))
        b_inds_grey = np.where(begin == 1)[0] + 1
        e_inds_grey = b_inds_grey + zero_dur - 1

        grey_vals = []
        for kk in range(len(b_inds_grey)):
            if e_inds_grey[kk] >= len(iratio):
                continue
            grey_vals.append(iratio[b_inds_grey[kk] : e_inds_grey[kk] + 1])
        if len(grey_vals) > 0:
            grey_all = np.stack(grey_vals)
            grey_all[grey_all == 0] = np.nan
            grey_mean = np.nanmean(grey_all, axis=0)
            baseline = np.nanmean(grey_mean[len(grey_mean) // 2 :])
        else:
            baseline = np.nanmean(iratio)

        if baseline == 0:
            baseline = 1.0  # avoid division by zero
        d_iratio = iratio / baseline - 1

        # Find epoch starts (rising edges), skip first and last
        ds = np.diff(s)
        inds = np.where(ds > 0)[0] + 1
        if len(inds) > 2:
            inds = inds[1:-1]  # clip first and last
        epochs = s[inds].astype(int)

        for jj in range(1, n_epoch + 1):
            I = np.where(epochs == jj)[0]
            if len(I) == 0:
                continue
            bind = inds[I]

            all_epochs = np.zeros((len(bind), full_dur))
            valid = 0
            for mm in range(len(bind)):
                i_begin = bind[mm] - zero_dur
                i_end = bind[mm] + epoch_dur + zero_dur
                if i_begin < 0 or i_end > len(d_iratio):
                    continue
                all_epochs[valid, :] = d_iratio[i_begin:i_end]
                valid += 1

            if valid > 0:
                temp = np.mean(all_epochs[:valid], axis=0)
                # Subtract pre-epoch baseline
                temp = temp - np.mean(temp[:zero_dur])
                rats[ii, jj - 1, :] = temp

        fly_ids.append(roi["flyID"])
        names.append(roi["name"])

    fly_ids = np.array(fly_ids)
    return {
        "rats": rats,
        "stimstruct": stimstruct,
        "flyID": fly_ids,
        "name": names,
        "n_epoch": n_epoch,
        "epoch_dur": epoch_dur,
        "zero_dur": zero_dur,
    }


def aggregate_adapting_contrast(roi_data_list, off_stimulus=True):
    """Aggregate adapting contrast step data (Figure 3A-E).

    Groups OFF-OFF epoch combinations using cantor pairing.
    Each epoch pair = 360 timepoints (300 + 60 at 10Hz = 36s).

    Parameters
    ----------
    roi_data_list : list of dicts from load_figure_data
    off_stimulus : bool
        True for OFF-OFF stimulus (default), False for ON-ON.

    Returns
    -------
    dict with keys:
        rats (nROIs, 7, 360), flyID, name, stims
    """
    N_EPOCHS = 7
    realepochlength1 = 300
    realepochlength2 = 60
    total_len = realepochlength1 + realepochlength2  # 360

    n_rois = len(roi_data_list)
    rats = np.full((n_rois, N_EPOCHS, total_len), np.nan)
    fly_ids = []
    names = []

    # Build the cantor→epoch mapping for OFF stimulus
    if off_stimulus:
        epoch_map = {
            cantor(1, 0): 0,  # 100% OFF → 50% OFF
            cantor(2, 0): 1,  # 100% OFF → GREY
            cantor(3, 0): 2,  # 100% OFF → 50% ON
            cantor(4, 0): 3,  # 100% OFF → 100% ON
            cantor(5, 0): 4,  # 50% OFF → 100% OFF
            cantor(6, 0): 5,  # 50% OFF → GREY
            cantor(7, 0): 6,  # 50% OFF → 50% ON
        }
    else:
        # ON-ON stimulus mapping (not currently used for figures of interest)
        epoch_map = {}

    for ii, roi in enumerate(roi_data_list):
        s = roi["istim"]
        iratio = roi["iratio"]

        # dF/F: normalized to IQR of lowest-luminance values
        is_dark = s == 0
        lowest_vals = iratio[is_dark]
        if len(lowest_vals) > 0:
            lowest_sorted = np.sort(lowest_vals)
            q25 = int(round(0.25 * len(lowest_sorted)))
            q75 = int(round(0.75 * len(lowest_sorted)))
            if q25 < q75:
                baseline_vals = lowest_sorted[q25:q75]
            else:
                baseline_vals = lowest_sorted
            baseline = np.mean(baseline_vals)
        else:
            baseline = np.mean(iratio)

        if baseline == 0:
            baseline = 1.0
        d_iratio = iratio / baseline - 1

        # Find epoch transitions
        ds = np.diff(s)
        b_inds_raw = np.where(ds != 0)[0] + 1

        if len(b_inds_raw) < 2:
            fly_ids.append(roi["flyID"])
            names.append(roi["name"])
            continue

        # Compute epoch lengths for this ROI
        epochlengths1 = []
        epochlengths2 = []
        for kk in range(0, len(b_inds_raw) - 1, 2):
            if kk + 1 < len(b_inds_raw):
                epochlengths2.append(b_inds_raw[kk + 1] - b_inds_raw[kk])
            if kk + 2 < len(b_inds_raw):
                epochlengths1.append(b_inds_raw[kk + 2] - b_inds_raw[kk + 1])

        # Shift b_inds back by 149 to start before the transition (matching MATLAB)
        b_inds = np.where(ds != 0)[0] - 149

        # Compute end indices
        e_inds = np.zeros(len(b_inds), dtype=int)
        for kk in range(len(b_inds)):
            half_kk = int(round(0.5 * (kk + 1)))  # MATLAB uses round(0.5*kk) with 1-based
            if half_kk > 0 and half_kk <= len(epochlengths1):
                el1 = epochlengths1[min(half_kk - 1, len(epochlengths1) - 1)]
            else:
                el1 = realepochlength1
            if half_kk > 0 and half_kk <= len(epochlengths2):
                el2 = epochlengths2[min(half_kk - 1, len(epochlengths2) - 1)]
            else:
                el2 = realepochlength2
            e_inds[kk] = b_inds[kk] + el1 + el2 - 1

        # Identify pairs using cantor
        pairs = np.full(len(b_inds), -1, dtype=int)
        for kk in range(len(b_inds)):
            if e_inds[kk] >= len(s) or b_inds[kk] + 150 >= len(s):
                continue
            stim_at_begin = int(s[b_inds[kk] + 150])
            stim_at_end = int(s[e_inds[kk]])
            pairs[kk] = cantor(stim_at_begin, stim_at_end)

        # Accumulate traces by epoch type
        num_pairs = max(1, int(round(len(s) / total_len)))
        temp = np.zeros((num_pairs, N_EPOCHS, total_len))

        for jj in range(len(pairs)):
            if pairs[jj] < 0:
                continue
            if jj >= len(e_inds) or e_inds[jj] >= len(d_iratio):
                continue
            if pairs[jj] not in epoch_map:
                continue

            epoch_idx = epoch_map[pairs[jj]]
            b_ind = b_inds[jj]
            e_ind = e_inds[jj]

            # Adjust to exactly total_len
            trace_len = e_ind - b_ind + 1
            if trace_len > total_len:
                e_ind -= trace_len - total_len
            elif trace_len < total_len:
                e_ind += total_len - trace_len

            if b_ind < 0 or e_ind >= len(d_iratio) or e_ind - b_ind + 1 != total_len:
                continue

            if jj < num_pairs:
                temp[jj, epoch_idx, :] = d_iratio[b_ind : b_ind + total_len]

        # Replace zeros with NaN and average
        temp[temp == 0] = np.nan
        rats[ii, :, :] = np.nanmean(temp, axis=0)

        fly_ids.append(roi["flyID"])
        names.append(roi["name"])

    fly_ids = np.array(fly_ids)
    return {"rats": rats, "flyID": fly_ids, "name": names}


def aggregate_05steps(roi_data_list):
    """Aggregate random 5-level luminance step data (Figure 3F-I).

    20 possible transitions between 5 luminance levels, each held ~10s.
    dF/F normalized to 100%ON baseline.

    Parameters
    ----------
    roi_data_list : list of dicts from load_figure_data

    Returns
    -------
    dict with keys:
        rats (nROIs, 20, 2*epochlength), flyID, name
    """
    N_EPOCHS = 20

    # Determine epoch length using mode of inter-transition durations
    epochlengths = []
    for roi in roi_data_list:
        s = roi["istim"]
        ds = np.diff(s)
        b_inds = np.where(ds != 0)[0] + 1
        if len(b_inds) < 2:
            continue
        durations = np.diff(b_inds)
        epochlengths.append(scipy_mode(durations, keepdims=False).mode)
    epochlength = int(round(np.mean(epochlengths)))

    n_rois = len(roi_data_list)
    rats = np.full((n_rois, N_EPOCHS, 2 * epochlength), np.nan)
    fly_ids = []
    names = []

    # Build cantor→epoch mapping for all 20 non-self transitions
    # Levels: 0=100%OFF, 1=50%OFF, 2=GREY, 3=50%ON, 4=100%ON
    transitions = [
        (0, 1), (0, 2), (0, 3), (0, 4),  # from 100%OFF
        (1, 0), (1, 2), (1, 3), (1, 4),  # from 50%OFF
        (2, 0), (2, 1), (2, 3), (2, 4),  # from GREY
        (3, 0), (3, 1), (3, 2), (3, 4),  # from 50%ON
        (4, 0), (4, 1), (4, 2), (4, 3),  # from 100%ON
    ]
    epoch_map = {cantor(a, b): idx for idx, (a, b) in enumerate(transitions)}

    for ii, roi in enumerate(roi_data_list):
        s = roi["istim"]
        iratio = roi["iratio"]

        # dF/F baseline: 100%ON epochs (stim == 4), last 20%
        is_100on = s == 4
        begin_100 = np.diff(is_100on.astype(float))
        b_inds_100 = np.where(begin_100 == 1)[0] + 1
        e_inds_100 = b_inds_100 + epochlength - 1

        on_100_vals = []
        for kk in range(len(b_inds_100)):
            if e_inds_100[kk] >= len(iratio):
                continue
            on_100_vals.append(iratio[b_inds_100[kk] : e_inds_100[kk] + 1])

        if len(on_100_vals) > 0:
            on_100_all = np.stack(on_100_vals)
            on_100_all[on_100_all == 0] = np.nan
            on_100_mean = np.nanmean(on_100_all, axis=0)
            baseline = np.nanmean(on_100_mean[int(round(0.8 * len(on_100_mean))) :])
        else:
            baseline = np.nanmean(iratio)

        if baseline == 0:
            baseline = 1.0
        d_iratio = iratio / baseline - 1

        # Find transitions
        ds = np.diff(s)
        b_inds = np.where(ds != 0)[0] + 1
        e_inds = b_inds + 2 * epochlength - 1

        # Identify pairs
        pairs = np.full(len(b_inds), -1, dtype=int)
        for kk in range(len(b_inds)):
            if e_inds[kk] >= len(s):
                continue
            pairs[kk] = cantor(int(s[b_inds[kk]]), int(s[e_inds[kk]]))

        # Accumulate traces
        num_pairs = max(1, int(round(len(s) / epochlength * 2)))
        temp = np.zeros((num_pairs, N_EPOCHS, 2 * epochlength))

        for jj in range(len(pairs)):
            if pairs[jj] < 0 or pairs[jj] not in epoch_map:
                continue
            if jj >= len(e_inds) or e_inds[jj] >= len(d_iratio):
                continue

            epoch_idx = epoch_map[pairs[jj]]
            b_ind = b_inds[jj]
            e_ind = e_inds[jj]

            if e_ind >= len(d_iratio):
                continue

            if jj < num_pairs:
                temp[jj, epoch_idx, :] = d_iratio[b_ind : e_ind + 1]

        temp[temp == 0] = np.nan
        rats[ii, :, :] = np.nanmean(temp, axis=0)

        fly_ids.append(roi["flyID"])
        names.append(roi["name"])

    fly_ids = np.array(fly_ids)

    # Build transition labels
    level_names = ["100%OFF", "50%OFF", "GREY", "50%ON", "100%ON"]
    transition_labels = [f"{level_names[a]}->{level_names[b]}" for a, b in transitions]

    return {
        "rats": rats,
        "flyID": fly_ids,
        "name": names,
        "epochlength": epochlength,
        "transition_labels": transition_labels,
    }


# ---------------------------------------------------------------------------
# ROI classification
# ---------------------------------------------------------------------------


def classify_rois_by_correlation(rats, stims, threshold=0.5):
    """Classify ROIs as positively or negatively correlated with stimulus.

    Used for Figure 2B, 2C.

    Parameters
    ----------
    rats : (nROIs, nTime) array
    stims : (nROIs, nTime) array
    threshold : float

    Returns
    -------
    neg_mask : (nROIs,) bool - True for negatively correlated ROIs
    pos_mask : (nROIs,) bool - True for positively correlated ROIs
    """
    mean_stim = np.nanmean(stims, axis=0)

    # Remove zero-only and NaN rows
    valid = (np.sum(rats, axis=1) != 0) & ~np.any(np.isnan(rats), axis=1)
    Q = np.array([np.corrcoef(mean_stim, rats[i])[0, 1] if valid[i] else 0
                  for i in range(len(rats))])

    neg_mask = Q < -threshold
    pos_mask = Q > threshold
    return neg_mask, pos_mask


def classify_rois_05steps(rats, n_epochs=20):
    """Classify ROIs for the 5-step ON-OFF stimulus (Figure 3F).

    Uses the mean ON-step vs OFF-step response to determine correlation.

    Parameters
    ----------
    rats : (nROIs, 20, nTime) array

    Returns
    -------
    neg_mask : (nROIs,) bool
    pos_mask : (nROIs,) bool
    """
    # Transition directions: from→to
    from_levels = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
    to_levels = [1, 2, 3, 4, 0, 2, 3, 4, 0, 1, 3, 4, 0, 1, 2, 4, 0, 1, 2, 3]
    step = np.array(to_levels) - np.array(from_levels)
    on_step = step > 0

    n_rois = rats.shape[0]
    neg_mask = np.zeros(n_rois, dtype=bool)
    pos_mask = np.zeros(n_rois, dtype=bool)

    for ii in range(n_rois):
        roi_rats = rats[ii]
        mean_on = np.nanmean(roi_rats[on_step], axis=0)
        mean_off = np.nanmean(roi_rats[~on_step], axis=0)

        n_time = roi_rats.shape[1]
        # Compare response window (after transition) to pre-transition baseline
        t_start = int(round(0.5 * n_time))  # halfway = transition point
        reaction_on = np.nanmean(mean_on[t_start : t_start + 20]) - np.nanmean(
            mean_on[t_start - 20 : t_start]
        )
        reaction_off = np.nanmean(mean_off[t_start : t_start + 20]) - np.nanmean(
            mean_off[t_start - 20 : t_start]
        )

        if reaction_off > 0 and reaction_on < 0:
            neg_mask[ii] = True
        elif reaction_off < 0 and reaction_on > 0:
            pos_mask[ii] = True

    return neg_mask, pos_mask


def classify_rois_adapting(rats, fly_ids):
    """Classify ROIs for adapting contrast (Figure 3A).

    Checks if response in epoch window is higher than baseline,
    matching crosscorrelation_negative.m logic.

    Parameters
    ----------
    rats : (nROIs, 7, 360) array
    fly_ids : (nROIs,) array

    Returns
    -------
    neg_mask : (nROIs,) bool
    """
    n_rois = rats.shape[0]
    neg_mask = np.zeros(n_rois, dtype=bool)

    for ii in range(n_rois):
        # Check response in the B-step window (timepoints 160-210 = 16-21s)
        # vs baseline (timepoints 1-140 = 0-14s)
        trace = rats[ii, -1, :]  # last epoch combination
        if np.all(np.isnan(trace)):
            continue

        baseline = np.nanmean(trace[:140])
        response = np.nanmax(trace[160:210])

        if response > baseline:
            neg_mask[ii] = True

    return neg_mask


# ---------------------------------------------------------------------------
# HDF5 storage
# ---------------------------------------------------------------------------


def save_figure_to_hdf5(
    h5file,
    group_name,
    rats,
    fly_ids,
    time_axis,
    stimulus,
    cell_type,
    stimulus_type,
    fly_means=None,
    grand_mean=None,
    sem=None,
    fly_ids_unique=None,
    roi_names=None,
    epoch_labels=None,
    extra_attrs=None,
):
    """Save processed figure data to an HDF5 group.

    Parameters
    ----------
    h5file : h5py.File (open for writing)
    group_name : str - e.g. 'fig2b_L3'
    rats : ndarray - per-ROI traces
    fly_ids : (nROIs,) array
    time_axis : (nTime,) array in seconds
    stimulus : stimulus trace or pattern
    cell_type : str - 'L2' or 'L3'
    stimulus_type : str - description
    fly_means, grand_mean, sem : optional precomputed per-fly stats
    fly_ids_unique : (nFlies,) array
    roi_names : list of str
    epoch_labels : list of str for multi-epoch data
    extra_attrs : dict of additional attributes
    """
    grp = h5file.require_group(group_name)

    # Attributes
    grp.attrs["cell_type"] = cell_type
    grp.attrs["stimulus_type"] = stimulus_type
    grp.attrs["sampling_rate"] = IRATE
    grp.attrs["n_rois"] = rats.shape[0]
    if fly_ids_unique is not None:
        grp.attrs["n_flies"] = len(fly_ids_unique)
    if extra_attrs:
        for k, v in extra_attrs.items():
            grp.attrs[k] = v

    # Per-ROI data
    roi_grp = grp.require_group("per_roi")
    _write_dataset(roi_grp, "traces", rats)
    fly_id_arr = np.array(fly_ids)
    if fly_id_arr.dtype.kind in ("U", "O"):  # string fly IDs (dataID)
        _write_dataset(roi_grp, "fly_id", fly_id_arr.astype(object), dtype=h5py.string_dtype())
    else:
        _write_dataset(roi_grp, "fly_id", fly_id_arr)
    if roi_names is not None:
        dt = h5py.string_dtype()
        _write_dataset(roi_grp, "roi_name", np.array(roi_names, dtype=object), dtype=dt)

    # Per-fly data
    if fly_means is not None:
        fly_grp = grp.require_group("per_fly")
        _write_dataset(fly_grp, "fly_means", fly_means)
        _write_dataset(fly_grp, "grand_mean", grand_mean)
        _write_dataset(fly_grp, "sem", sem)
        if fly_ids_unique is not None:
            if fly_ids_unique.dtype.kind in ("U", "O"):
                _write_dataset(fly_grp, "fly_ids", fly_ids_unique.astype(object),
                               dtype=h5py.string_dtype())
            else:
                _write_dataset(fly_grp, "fly_ids", fly_ids_unique)

    # Time axis
    _write_dataset(grp, "time", time_axis)

    # Stimulus
    _write_dataset(grp, "stimulus_trace", np.array(stimulus))

    # Epoch labels
    if epoch_labels is not None:
        dt = h5py.string_dtype()
        _write_dataset(grp, "epoch_labels", np.array(epoch_labels, dtype=object), dtype=dt)


def _write_dataset(group, name, data, dtype=None):
    """Write a dataset, replacing if it already exists."""
    if name in group:
        del group[name]
    if dtype is not None:
        group.create_dataset(name, data=data, dtype=dtype)
    else:
        group.create_dataset(name, data=data)


def load_figure_from_hdf5(h5file, group_name):
    """Load processed figure data from an HDF5 group.

    Returns
    -------
    dict with keys matching the HDF5 structure.
    """
    grp = h5file[group_name]
    result = {
        "attrs": dict(grp.attrs),
        "time": grp["time"][:],
        "stimulus_trace": grp["stimulus_trace"][:],
    }

    if "per_roi" in grp:
        roi = grp["per_roi"]
        result["traces"] = roi["traces"][:]
        result["fly_id"] = roi["fly_id"][:]
        if "roi_name" in roi:
            result["roi_name"] = [s.decode() if isinstance(s, bytes) else s
                                  for s in roi["roi_name"][:]]

    if "per_fly" in grp:
        fly = grp["per_fly"]
        result["fly_means"] = fly["fly_means"][:]
        result["grand_mean"] = fly["grand_mean"][:]
        result["sem"] = fly["sem"][:]
        if "fly_ids" in fly:
            result["fly_ids"] = fly["fly_ids"][:]

    if "epoch_labels" in grp:
        result["epoch_labels"] = [s.decode() if isinstance(s, bytes) else s
                                  for s in grp["epoch_labels"][:]]

    return result
