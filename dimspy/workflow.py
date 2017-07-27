#!/usr/bin/python
# -*- coding: utf-8 -*-


import os
import logging
import operator
import collections
import h5py
import numpy as np
from itertools import combinations
import zipfile
from models.peaklist import PeakList
from models.peak_matrix import PeakMatrix
from portals import hdf5_portal
from portals import txt_portal
from portals.paths import check_paths
from experiment import check_metadata
from experiment import update_class_labels
from experiment import update_metadata
from experiment import idxs_reps_from_filelist
from experiment import mz_range_from_header
from experiment import interpret_experiment
from process.peak_alignment import align_peaks
from process.peak_filters import filter_fraction
from process.peak_filters import filter_blank_peaks
from process.peak_filters import filter_rsd
from process.peak_filters import filter_mz_ranges
from process.peak_filters import filter_attr
from process.peak_filters import filter_ringing
from process.scan_processing import average_replicate_scans
from process.scan_processing import join_peaklists
from process.scan_processing import read_scans
from process.scan_processing import remove_edges


def process_scans(source, function_noise, snr_thres, ppm, min_fraction=None, rsd_thres=None, min_scans=1, filelist=None,
                  skip_stitching=False, remove_mz_range=None, ringing_thres=None, filter_scan_events=None, block_size=2000, ncpus=None):

    if filter_scan_events is None:
        filter_scan_events = {}
    if remove_mz_range is None:
        remove_mz_range = []
    filenames = check_paths(filelist, source)
    if len([fn for fn in filenames if not fn.lower().endswith(".mzml") or not fn.lower().endswith(".raw")]) == 0:
        raise IOError("Incorrect file format. Provide .mzML and .raw files")

    if filelist is not None:
        fl = check_metadata(filelist)
    else:
        fl = collections.OrderedDict()

    pls = []
    for i in range(len(filenames)):

        print
        print os.path.basename(filenames[i])

        if type(source) is not str:
            source = ""

        print "Reading scans...."
        pls_scans = read_scans(filenames[i], source, function_noise, min_scans, filter_scan_events)

        if type(remove_mz_range) == list and len(remove_mz_range) > 0:
            print "Removing m/z ranges....."
            for h in pls_scans:
                pls_scans[h] = [filter_mz_ranges(pl, remove_mz_range) for pl in pls_scans[h] if len(pl.mz) > 0]

        if not skip_stitching:
            # Temp fix lock mass (1)
            mz_ranges = []
            for h in pls_scans:
                mzr = mz_range_from_header(h)
                if mzr not in mz_ranges:
                    mz_ranges.append(mzr)
                else:
                    del pls_scans[h]
                    print mzr, h, "exclude from SIM-Stitching"
            # End temp fix  lock mass(1)

            exp = interpret_experiment(mz_ranges)
            if exp == "overlapping":
                print "Removing 'edges' from SIM windows....."
                pls_scans = remove_edges(pls_scans)

        if ringing_thres is not None and float(ringing_thres) > 0.0:
            print "Removing ringing artifacts....."
            for h in pls_scans:
                pls_scans[h] = [filter_ringing(pl, threshold=ringing_thres, bin_size=1.0)
                                for pl in pls_scans[h] if len(pl.mz) > 0]

        print "Removing noise....."
        for h in pls_scans:

            # Temp fix  lock mass (3)
            pls_scans[h] = [filter_attr(pl, "snr", min_threshold=snr_thres) for pl in pls_scans[h][0:1] if len(pl.mz) > 0]
            # End temp fix  lock mass (3)

        print "Aligning, averaging and filtering peaks....."
        pls_avg = []
        print "--------------------------------------------"
        print "event\tpeaks\tmedian_rsd"
        for h in pls_scans:
            if len(pls_scans[h]) >= 1:
                if sum(pl.shape[0] for pl in pls_scans[h]) == 0:
                    logging.warning("No scan data available for {}".format(h))
                else:
                    pl_avg = average_replicate_scans(h, pls_scans[h], ppm, min_fraction, rsd_thres, block_size, ncpus)
                    pls_avg.append(pl_avg)
                    print "{}\t{}\t{}".format(h, pl_avg.shape[0], np.nanmedian(pl_avg.rsd))
            else:
                logging.warning("No scan data available for {}".format(h))
                print "{}\t{}\t{}".format(h, 0, "na")
        
        if len(pls_avg) == 0:
            raise IOError("No peaks remaining after filtering. Remove MS data file / sample.")

        if not skip_stitching:
            pl = join_peaklists(os.path.basename(filenames[i]), pls_avg)
            if "class" in fl:
                pl.tags.add_tags(class_label=fl["class"][i])
            for k in fl.keys():
                pl.metadata[k] = fl[k][i]
            pls.append(pl)
        else:
            for pl in pls_avg:
                pl = join_peaklists("{}#{}".format(os.path.basename(filenames[i]), pl.metadata["header"][0]), [pl])
                for k in fl.keys():
                    pl.metadata[k] = fl[k][i]
                pls.append(pl)
    return pls


# placeholder (synonym)
def stitch(source, function_noise, snr_thres, ppm, min_fraction=None, rsd_thres=None, min_scans=1, filelist=None,
           skip_stitching=False, remove_mz_range=None, ringing_thres=None, filter_scan_events=None, block_size=2000, ncpus=None):

    if filter_scan_events is None:
        filter_scan_events = {}
    if remove_mz_range is None:
        remove_mz_range = []

    return process_scans(source, function_noise, snr_thres, min_scans, ppm, min_fraction, rsd_thres, filelist,
                         skip_stitching, remove_mz_range, ringing_thres, filter_scan_events, block_size, ncpus)


def replicate_filter(source, ppm, replicates, min_peaks, rsd_thres=None, filelist=None, block_size=2000, ncpus=None):

    if replicates < min_peaks:
        raise IOError("Provide realistic values for the number of replicates and minimum number of peaks present (min_peaks)")

    filenames = check_paths(filelist, source)
    if len(filenames) == 0:
        raise IOError("Provide a filelist that list all the text files (columnname:filename) and assign replicate numbers to each filename/sample (columnname:replicate)")
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata(peaklists, fl)

    if not hasattr(peaklists[0].metadata, "replicate"):
        raise IOError("Provide a filelist and assign replicate numbers (columnname:replicate) to each filename/sample")

    idxs_peaklists = idxs_reps_from_filelist([pl.metadata.replicate for pl in peaklists])
    unique, counts = np.unique([pl.metadata.replicate for pl in peaklists], return_counts=True)

    if len(counts) <= 1:
        raise ValueError("No technical replicates available (single) - Skip 'replicate filter'")
    if max(unique) < replicates:
        raise ValueError("Replicates incorrectly labeled")
    if sum(counts) != len(peaklists):
        raise ValueError("Replicates incorrectly labeled")

    if len([True for idxs_pls in idxs_peaklists if len(idxs_pls) >= replicates]) == len(idxs_peaklists):
        print
        print "All combinations (n={}) for each each set of replicates will be " \
              "processed to calculate the most reproducible set".format(replicates)
        print
        print "rank\tID\tpeaks\tmedian_rsd({}/{})".format(replicates, replicates)
    else:
        raise ValueError("Not enough (technical) replicates available for each sample {}.")

    pls_rep_filt = []

    for idxs_pls in idxs_peaklists:

        temp = []

        for j, pls_comb in enumerate(combinations(peaklists[idxs_pls[0]:idxs_pls[-1] + 1], replicates)):

            pm = align_peaks(pls_comb, ppm, block_size, ncpus=ncpus)

            prefix = os.path.commonprefix([p.ID for p in pls_comb])
            merged_id = "{}{}".format(prefix, "_".join(map(str, [p.ID.replace(prefix, "").split(".")[0] for p in pls_comb])))

            pl = pm.to_peaklist(ID=merged_id)
            if "snr" in pm.attributes:
                pl.add_attribute("snr", pm.attr_mean_vector("snr"), on_index=2)

            pl.tags.add_tags(*pls_comb[0].tags.tag_of(None), **{t: pls_comb[0].tags.tag_of(t) for t in pls_comb[0].tags.tag_types})
            pl.add_attribute("present_flag", pm.present >= min_peaks, is_flag=True)

            if rsd_thres is not None:
                rsd_flag = map(lambda x: not np.isnan(x) and x < rsd_thres, pm.rsd)
                pl.add_attribute("rsd_flag", rsd_flag, flagged_only=False, is_flag=True)

            pl_filt = filter_attr(pl.copy(), attr_name="present", min_threshold=replicates, flag_name="pres_rsd")
            temp.append([pl, pl_filt.shape[0], np.median(pl_filt.rsd)])

        temp.sort(key=operator.itemgetter(2, 1))
        pls_rep_filt.append(temp[0][0])  # Most reproducible set of replicates

        for p in range(0, len(temp)):
            print "{}\t{}\t{}\t{}\t".format(p+1, temp[p][0].ID, temp[p][1], temp[p][2])
        print

    return pls_rep_filt


def align_samples(source, ppm, filelist=None, block_size=2000, ncpus=None):

    filenames = check_paths(filelist, source)
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata(peaklists, fl)

    return align_peaks(peaklists, ppm=ppm, block_size=block_size, ncpus=ncpus)


def blank_filter(peak_matrix, blank_label, min_fraction=1.0, min_fold_change=1.0, function="mean", rm_samples=True, class_labels=None):

    if min_fraction < 0.0 or min_fraction > 1.0:
        raise ValueError("Provide a value between 0. and 1.")
    if min_fold_change < 0:
        raise ValueError("Provide a value larger than zero.")
    if function not in ("mean", "median", "max"):
        raise ValueError("Mean, median or max intensity")

    if not isinstance(peak_matrix, PeakMatrix):
        if h5py.is_hdf5(peak_matrix):
            peak_matrix = hdf5_portal.load_peak_matrix_from_hdf5(peak_matrix)
        else:
            peak_matrix = txt_portal.load_peak_matrix_from_txt(peak_matrix)

    if class_labels is not None:
        peak_matrix = update_class_labels(peak_matrix, class_labels)

    if blank_label not in peak_matrix.peaklist_tag_values:
        raise IOError("Blank label ({}) does not exist".format(blank_label))

    return filter_blank_peaks(peak_matrix, blank_label, min_fraction, min_fold_change, function, rm_samples)


def sample_filter(peak_matrix, min_fraction, within=False, rsd=None, qc_label=None, class_labels=None):

    if not isinstance(peak_matrix, PeakMatrix):
        if h5py.is_hdf5(peak_matrix):
            peak_matrix = hdf5_portal.load_peak_matrix_from_hdf5(peak_matrix)
        else:
            peak_matrix = txt_portal.load_peak_matrix_from_txt(peak_matrix)

    if class_labels is not None:
        if not os.path.isfile(class_labels):
            raise IOError("{} does not exist".format(class_labels))
        peak_matrix = update_class_labels(peak_matrix, class_labels)

    if qc_label is not None:
        if qc_label not in peak_matrix.peaklist_tag_values:
            raise IOError("QC label ({}) does not exist".format(qc_label))

    peak_matrix = filter_fraction(peak_matrix, min_fraction, within_classes=within, class_tag_type=None)

    if rsd is not None:
        peak_matrix = filter_rsd(peak_matrix, rsd, qc_label)
    return peak_matrix


def hdf5_to_txt(fname, path_out, attr_name="intensity", separator="\t", transpose=False, comprehensive=False):
    assert os.path.isfile(fname), 'HDF5 database [%s] not exists' % fname
    assert h5py.is_hdf5(fname), 'input file [%s] is not a valid HDF5 database' % fname
    seps = {"comma": ",", "tab": "\t"}
    if separator in seps:
        separator = seps[separator]
    assert separator in [",", "\t"], "Incorrect separator ('tab', 'comma', ',', '\t')"
    f = h5py.File(fname, 'r')

    if "mz" in f:
        obj = hdf5_portal.load_peak_matrix_from_hdf5(fname)
        assert isinstance(obj, PeakMatrix)
        with open(os.path.join(path_out), "w") as pk_out:
            pk_out.write(obj.to_str(attr_name=attr_name, delimiter=separator, transpose=transpose, comprehensive=comprehensive))
    else:
        assert os.path.isdir(path_out), "File or Directory does not exist:".format(path_out)
        obj = hdf5_portal.load_peaklists_from_hdf5(fname)
        assert isinstance(obj[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
        if "#" in obj[0].ID:
            fns = set([pl.ID.split("#")[0] for pl in obj])
            sub_ids = [pl.ID.split("#")[1] for pl in obj]
            for fn in fns:
                with open(os.path.join(path_out, os.path.splitext(fn)[0] + ".txt"), "w") as pk_out:
                    for i, pl in enumerate(obj):
                        if fn in pl.ID:
                            pl.add_attribute("window", pl.full_shape[0] * [sub_ids[i]], flagged_only=False, on_index=3)
                            pk_out.write(pl.to_str(delimiter=separator))
                            pl.drop_attribute("window")
        else:
            for pl in obj:
                with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                    pk_out.write(pl.to_str(delimiter=separator))
    return


def merge_peaklists(source, filelist=None):

    if not isinstance(source, list):
        raise IOError("Incorrect input: list of lists of peaklists, list of peak matrix objects or list of HDF5 files expected.")

    pls_merged = []
    for s in source:
        if isinstance(s, list) or isinstance(s, tuple):
            if isinstance(s[0], PeakList):
                pls_merged.extend(s)
            else:
                raise IOError("Incorrect Object in list. Peaklist Object expected.")
        elif isinstance(s, PeakMatrix):
            pls = s.get_peaklists()
            pls_merged.extend(pls)
        elif h5py.is_hdf5(s):
            f = h5py.File(s, 'r')
            if "mz" in f:
                pm = txt_portal.load_peak_matrix_from_txt(s)
                pls = pm.get_peaklists()
            else:
                pls = hdf5_portal.load_peaklists_from_hdf5(s)
            f.close()
            pls_merged.extend(pls)
        else:
            raise IOError("Incorrect input: list of lists of peaklists, list of peak matrix objects or list of HDF5 files expected.")

    if filelist is not None:
        fl = check_metadata(filelist)
        pls_merged = update_metadata(pls_merged, fl)

    return pls_merged


def load_peaklists(source):

    if type(source) == str:
        source = source.encode('string-escape')
        if h5py.is_hdf5(source):
            peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
        elif zipfile.is_zipfile(source):
            zf = zipfile.ZipFile(source)
            filenames = zf.namelist()
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [txt_portal.load_peaklist_from_txt(zf.open(fn), ID=os.path.basename(fn), has_flag_col=True) for fn in filenames]
        elif os.path.isdir(source):
            filenames = os.listdir(source)
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [txt_portal.load_peaklist_from_txt(os.path.join(source, fn), ID=os.path.basename(fn), delimiter="\t", has_flag_col=False) for fn in filenames]
        else:
            raise TypeError("Incorrect format. Process .mzML and .raw files first using the 'process scans' function")
    elif type(source) == list or type(source) == tuple:
        if not isinstance(source[0], PeakList):
            raise TypeError("List has incorrect format. PeakList objects required.")
        else:
            peaklists = source
    else:
        raise IOError("Inccorrect input: list with peaklist objects or path")

    return peaklists
