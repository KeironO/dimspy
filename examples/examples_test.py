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
    source = os.path.join("E:\\GitHub\\dimspy\\tests\\test_results")
    fn_filelist = os.path.join("E:\\GitHub\\dimspy\\tests\\test_results\\filelist_test.txt")
    output = os.path.join("E:\\GitHub\\dimspy\\tests\\test_results")

    print "Process Scans....."
    pls = process_scans(source, min_scans=1, function_noise="median",
                        snr_thres=1.0, ppm=2.0, min_fraction=None, rsd_thres=None,
                        filelist=fn_filelist, remove_mz_range=[], filter_scan_events={"exclude":[["50", "620", "full"]]}, block_size=2000, ncpus=None)
    print "Finished"
    print

    for pl in pls:
        print pl.ID, pl.shape
        with open(os.path.join(output, pl.ID + ".txt"), "w") as out: out.write(pl.to_str("\t"))

    save_peaklists_as_hdf5(pls, os.path.join(output, "MTBLS79_mzml_triplicates_multiple_scans.hdf5"))
    pls = load_peaklists_from_hdf5(os.path.join(output, "MTBLS79_mzml_triplicates_multiple_scans.hdf5"))

    print
    print "Replicate Filter....."
    pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None)

    for pl in pls_rf:
        print pl.ID, pl.shape
        with open(os.path.join(output, pl.ID + ".txt"), "w") as out: out.write(pl.to_str("\t"))

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

if __name__ == '__main__':
    main()