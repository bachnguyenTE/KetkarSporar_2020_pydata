#!/usr/bin/env python3
"""
Convert Ketkar & Sporar 2020 calcium imaging data from .mat to HDF5.

Processes Figures 2B, 2C, 2D-E, 3A, 3F for L2 and L3 neurons.
Stores per-ROI traces, per-fly means, grand mean + SEM, and stimulus,
all aligned at 10 Hz.

Usage:
    python convert_to_hdf5.py [output_path]

Default output: processed_data.h5
"""

import sys
from pathlib import Path

import h5py
import numpy as np

import utils


DATA_ROOT = Path("data/in vivo two-photon calcium imaging")


def process_and_save(h5f, group_name, rats, fly_ids, time_axis, stimulus,
                     cell_type, stimulus_type, roi_names=None,
                     epoch_labels=None, extra_attrs=None):
    """Clean, compute per-fly stats, save to HDF5."""
    if rats.ndim == 2:
        valid = (np.nansum(np.abs(rats), axis=1) > 0) & ~np.all(np.isnan(rats), axis=1)
    else:
        valid = np.nansum(np.abs(rats), axis=tuple(range(1, rats.ndim))) > 0

    rats_clean = rats[valid]
    fids_clean = fly_ids[valid]
    names_clean = ([roi_names[i] for i in range(len(roi_names)) if valid[i]]
                   if roi_names else None)

    if rats_clean.ndim == 2:
        fly_means, grand_mean, sem = utils.mean_cat_full(rats_clean, fids_clean)
    else:
        orig_shape = rats_clean.shape
        flat = rats_clean.reshape(orig_shape[0], -1)
        fm_flat, gm_flat, sem_flat = utils.mean_cat_full(flat, fids_clean)
        fly_means = fm_flat.reshape(fm_flat.shape[0], *orig_shape[1:])
        grand_mean = gm_flat.reshape(orig_shape[1:])
        sem = sem_flat.reshape(orig_shape[1:])

    fly_ids_unique = np.unique(fids_clean)

    utils.save_figure_to_hdf5(
        h5f, group_name, rats_clean, fids_clean, time_axis, stimulus,
        cell_type, stimulus_type,
        fly_means=fly_means, grand_mean=grand_mean, sem=sem,
        fly_ids_unique=fly_ids_unique, roi_names=names_clean,
        epoch_labels=epoch_labels, extra_attrs=extra_attrs,
    )
    print(f"  {group_name}: {rats_clean.shape[0]} ROIs, "
          f"{len(fly_ids_unique)} flies, shape={rats_clean.shape}")


def main():
    output_path = sys.argv[1] if len(sys.argv) > 1 else "processed_data.h5"
    summary = utils.load_summary()

    with h5py.File(output_path, "w") as h5f:
        # ---- Figure 2B: 5s FFF ----
        print("Figure 2B: 5s full-field flash")
        for cell_type in ["L2", "L3"]:
            data_dir = DATA_ROOT / "Figure2" / "Figure2B" / f"{cell_type}_pData"
            roi_data = utils.load_figure_data(data_dir, summary, region="AT")
            agg = utils.aggregate_fff(roi_data)
            neg, _ = utils.classify_rois_by_correlation(agg["rats"], agg["stims"])
            dur = agg["rats"].shape[1]
            t = np.arange(dur) / utils.IRATE
            stim_trace = np.nanmean(agg["stims"][neg], axis=0)
            process_and_save(
                h5f, f"fig2b_{cell_type}", agg["rats"][neg], agg["flyID"][neg],
                t, stim_trace, cell_type, "5s_full_field_flash",
                roi_names=[agg["name"][i] for i in range(len(agg["name"])) if neg[i]])

        # ---- Figure 2C: 60s FFF ----
        print("Figure 2C: 60s full-field flash")
        data_dir = DATA_ROOT / "Figure2" / "Figure2C" / "L3_pData"
        roi_data = utils.load_figure_data(data_dir, summary, region="AT")
        agg = utils.aggregate_fff60s(roi_data)
        neg, _ = utils.classify_rois_by_correlation(agg["rats"], agg["stims"])
        dur = agg["rats"].shape[1]
        t = np.arange(dur) / utils.IRATE
        stim_trace = np.nanmean(agg["stims"][neg], axis=0)
        process_and_save(
            h5f, "fig2c_L3", agg["rats"][neg], agg["flyID"][neg],
            t, stim_trace, "L3", "60s_full_field_flash",
            roi_names=[agg["name"][i] for i in range(len(agg["name"])) if neg[i]])

        # ---- Figure 2D-E: Contrast steps ----
        print("Figure 2D-E: Contrast steps from grey")
        contrast_labels = [
            "-100% OFF", "-80% OFF", "-60% OFF", "-40% OFF", "-20% OFF",
            "20% ON", "40% ON", "60% ON", "80% ON", "100% ON",
        ]
        for cell_type in ["L2", "L3"]:
            data_dir = DATA_ROOT / "Figure2" / "Figure2D_G" / f"{cell_type}_pData"
            roi_data = utils.load_figure_data(data_dir, summary, region="AT")
            agg = utils.aggregate_standing_stripe(roi_data)
            full_dur = agg["rats"].shape[2]
            t = np.arange(full_dur) / utils.IRATE
            process_and_save(
                h5f, f"fig2dg_{cell_type}", agg["rats"], agg["flyID"],
                t, agg["stimstruct"], cell_type, "contrast_steps_grey_baseline",
                roi_names=agg["name"], epoch_labels=contrast_labels,
                extra_attrs={"contrast_values": [-100, -80, -60, -40, -20, 20, 40, 60, 80, 100]})

        # ---- Figure 3A: Adapting contrast ----
        print("Figure 3A: Adapting contrast steps")
        adapt_labels = [
            "100%OFF->50%OFF", "100%OFF->GREY", "100%OFF->50%ON",
            "100%OFF->100%ON", "50%OFF->100%OFF", "50%OFF->GREY",
            "50%OFF->50%ON",
        ]
        for cell_type in ["L2", "L3"]:
            data_dir = DATA_ROOT / "Figure3" / "Figure3A_E" / f"{cell_type}_pData"
            roi_data = utils.load_figure_data(data_dir, summary)
            agg = utils.aggregate_adapting_contrast(roi_data, off_stimulus=True)
            neg_mask = utils.classify_rois_adapting(agg["rats"], agg["flyID"])
            total_len = agg["rats"].shape[2]
            t = np.arange(total_len) / utils.IRATE
            stim_pattern = np.zeros(total_len)
            stim_pattern[300:] = 1.0
            process_and_save(
                h5f, f"fig3a_{cell_type}", agg["rats"][neg_mask], agg["flyID"][neg_mask],
                t, stim_pattern, cell_type, "adapting_contrast_OFF_OFF",
                roi_names=[agg["name"][i] for i in range(len(agg["name"])) if neg_mask[i]],
                epoch_labels=adapt_labels,
                extra_attrs={"luminance_values_cdm2": [75, 64, 54, 43, 32, 21]})

        # ---- Figure 3F: Random 5-level steps ----
        print("Figure 3F: Random luminance steps")
        for cell_type in ["L2", "L3"]:
            data_dir = DATA_ROOT / "Figure3" / "Figure3F_H_I" / f"{cell_type}_pData"
            roi_data = utils.load_figure_data(data_dir, summary)
            agg = utils.aggregate_05steps(roi_data)
            neg_mask, _ = utils.classify_rois_05steps(agg["rats"])
            n_time = agg["rats"].shape[2]
            t = np.arange(n_time) / utils.IRATE
            stim_pattern = np.zeros(n_time)
            stim_pattern[n_time // 2:] = 1.0
            process_and_save(
                h5f, f"fig3f_{cell_type}", agg["rats"][neg_mask], agg["flyID"][neg_mask],
                t, stim_pattern, cell_type, "random_5level_luminance_steps",
                roi_names=[agg["name"][i] for i in range(len(agg["name"])) if neg_mask[i]],
                epoch_labels=agg["transition_labels"],
                extra_attrs={"epochlength_samples": agg["epochlength"]})

    print(f"\nDone. Output: {output_path} "
          f"({Path(output_path).stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
