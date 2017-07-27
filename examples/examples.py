#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from dimspy.workflow import process_scans
from dimspy.workflow import replicate_filter
from dimspy.workflow import align_samples
from dimspy.workflow import blank_filter
from dimspy.workflow import sample_filter
from dimspy.portals.hdf5_portal import save_peaklists_as_hdf5, save_peak_matrix_as_hdf5
from dimspy.portals.hdf5_portal import load_peaklists_from_hdf5, load_peak_matrix_from_hdf5


def main():

    # Example 1 - mzML files (zip file)
    #source = os.path.join("..", "tests", "data", "MTBLS79_subset", "MTBLS79_mzml_triplicates.zip")
    #fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_mzml_triplicates.txt")
    output = os.path.join("..", "tests", "test_results")

    # source = os.path.join("..", "tests", "data", "FracResus_DIMS")
    source = os.path.join("..", "tests", "data", "A08_daph_PlateC.raw")

    fse = {"include": [[50.0000, 100.0000, "sim"],
                       [90.0000, 140.0000, "sim"],
                       [130.0000, 180.0000, "sim"],
                       [170.0000, 220.0000, "sim"],
                       [210.0000, 260.0000, "sim"],
                       [250.0000, 300.0000, "sim"],
                       [290.0000, 340.0000, "sim"],
                       [330.0000, 380.0000, "sim"],
                       [370.0000, 420.0000, "sim"],
                       [410.0000, 460.0000, "sim"],
                       [450.0000, 500.0000, "sim"],
                       [490.0000, 540.0000, "sim"],
                       [530.0000, 580.0000, "sim"],
                       [570.0000, 620.0000, "sim"],
                       [610.0000, 660.0000, "sim"],
                       [650.0000, 700.0000, "sim"],
                       [690.0000, 740.0000, "sim"],
                       [730.0000, 780.0000, "sim"],
                       [770.0000, 820.0000, "sim"],
                       [810.0000, 860.0000, "sim"],
                       [850.0000, 900.0000, "sim"],
                       [890.0000, 940.0000, "sim"],
                       [930.0000, 980.0000, "sim"]]}

    print "Process Scans....."
    pls = process_scans(source, min_scans=1, function_noise="median",
                        snr_thres=1.0, ppm=2.0, min_fraction=None, rsd_thres=None,
                        filelist=None, remove_mz_range=[], filter_scan_events=fse, block_size=2000, ncpus=None)
    print "Finished"
    print

    for pl in pls:
        print pl.ID, pl.shape
        with open(os.path.join(output, pl.ID + ".txt"), "w") as out: out.write(pl.to_str("\t"))
    """
    save_peaklists_as_hdf5(pls, os.path.join(output, "MTBLS79_mzml_triplicates_multiple_scans.hdf5"))

    pls = load_peaklists_from_hdf5(os.path.join(output, "MTBLS79_mzml_triplicates_multiple_scans.hdf5"))

    print
    print "Replicate Filter....."
    pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=20.0)
    print "Finished"
    print

    print "Align Samples...."
    pm = align_samples(pls_rf, ppm=2.0)
    print "Finished", pm.shape
    print
    print "Blank Filter"
    pm_bf = blank_filter(pm, "blank", min_fraction=1.0, min_fold_change=10.0, function="mean", rm_samples=False)
    print "Finished", pm_bf.shape
    print

    # Example 2 - RAW files (Directory)
    source = os.path.join("..", "tests", "data", "MTBLS79_subset", "raw")
    fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_raw_triplicates.txt")
    print
    print "Process Scans....."
    pls = process_scans(source, min_scans=1, function_noise="noise_packets",
                        snr_thres=3.0, ppm=2.0, min_fraction=None, rsd_thres=None,
                        filelist=fn_filelist, remove_mz_range=[], filter_scan_events={}, block_size=2000, ncpus=None)
    print
    print "Finished"
    save_peaklists_as_hdf5(pls, os.path.join(output, "MTBLS79_raw_triplicates_single_scan.hdf5"))
    pls = load_peaklists_from_hdf5(os.path.join(output, "MTBLS79_raw_triplicates_single_scan.hdf5"))

    print "Align Samples...."
    pm = align_samples(pls, ppm=2.0)
    print "Finished", pm.shape

    save_peak_matrix_as_hdf5(pm, os.path.join(output, "MTBLS79_raw_peak_matrix_as.hdf5"))
    pm = load_peak_matrix_from_hdf5(os.path.join(output, "MTBLS79_raw_peak_matrix_as.hdf5"))

    print "Sample Filter"
    pm_sf = sample_filter(pm, 0.5, within=False)
    print "Finished", pm_sf.shape
    print "Sample Filter"
    pm_sf = sample_filter(pm, 0.8, within=False)
    print "Finished", pm_sf.shape

    # Example 3 - Subset m/z ranges
    source = os.path.join("..", "tests", "data", "MTBLS79_subset", "raw")
    fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_raw_triplicates.txt")

    print
    print "Process Scans....."
    pls = process_scans(source, min_scans=1, function_noise="noise_packets",
                        snr_thres=3.0, ppm=2.0, filelist=fn_filelist, filter_scan_events={"include": [[70.0, 170.0, "sim"]]})
    print
    print
    """
if __name__ == '__main__':
    main()
