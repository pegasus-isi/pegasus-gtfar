#!/usr/bin/env python

import os
import sys

import argparse

import pylab as pl

import matplotlib.pyplot as plt

from collections import defaultdict as dd


def load_bar_data(infofile):
    INIT = True
    k = 0
    vals = []
    FEATURE_DICT = dd(int)
    FAMILY_DICT = dd(int)
    for line in infofile:
        lp = line.split()
        if lp[0] == 'FeatureTypes':
            for i in range(1, len(lp), 5):
                #print FEATURE_DICT[lp[i]]

                if lp[i] == "IJXN" or lp[i] == "ITRN":
                    FEATURE_DICT["INTRONIC"] += int(float(lp[i + 1]))

                elif lp[i] == "EXON":
                    FEATURE_DICT["EXONIC"] += int(float(lp[i + 1]))

                elif lp[i] == "KJXN":
                    FEATURE_DICT["KNOWN-JXN"] += int(float(lp[i + 1]))
                elif lp[i] == "KJXN" or lp[i] == "NJXN":
                    FEATURE_DICT["NOVEL-JXN"] += int(float(lp[i + 1]))
                else:
                    FEATURE_DICT[lp[i].upper()] = int(float(lp[i + 1]))
        elif lp[0] == 'FamilyTypes':
            for i in range(6, len(lp), 5):
                if lp[i] == "ribosome":
                    FEATURE_DICT["RIBOSOMAL"] = int(float(lp[i + 1]))
                elif lp[i] == "mitochondria":
                    FEATURE_DICT["MITOCHONDRIAL"] = int(float(lp[i + 1]))
                elif lp[i].upper() == "PROCESSED_TRANSCRIPT" or lp[i].upper() == "MISC":
                    FAMILY_DICT["NONCODING"] = int(float(lp[i + 1]))
                else:
                    FAMILY_DICT[lp[i].upper()] = int(float(lp[i + 1]))
        else:
            continue

    return FEATURE_DICT, FAMILY_DICT


def make_adjacent_bars(feats, fams, output_file):
    default_colors = ['blue', 'Crimson', 'black', 'cyan', 'magenta', 'green']
    more_colors = ['Pink', 'Bisque', 'Brown', 'CadetBlue', 'Chartreuse', 'Chocolate', 'DarkGray', 'DarkKhaki']
    even_more_colors = ['DarkOrange', 'DarkSalmon', 'DarkSlateGray', 'DeepPink', 'DarkTurquoise', 'ForestGreen',
                        'GoldenRod', 'Olive', 'Orange', 'Aquamarine', 'MediumSpringGreen']
    final_colors = ['PapayaWhip', 'Tomato', 'Thistle', 'Salmon', 'Plum', 'Lime', 'Fuchsia', 'FireBrick', 'black',
                    'white']
    colors = default_colors + more_colors + even_more_colors + final_colors
    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8, 4))
    fig.subplots_adjust(bottom=0.15)
    feat_tuples = []
    fam_tuples = []
    for f in feats:
        feat_tuples.append([feats[f], f])
    for f in fams:
        fam_tuples.append([fams[f], f])

    feat_tuples.sort(reverse=True)
    fam_tuples.sort(reverse=True)
    feature_vals, feature_labels, fam_vals, fam_labels, bar_labels = [], [], [], [], []
    WIDTH = 0.5
    k = 0
    for f in feat_tuples:
        if f[1] == "FILTER": continue

        feature_vals.append(f[0])
        feature_labels.append(f[1])
        bar_labels.append(ax0.bar(k * WIDTH, f[0], WIDTH, bottom=0, color=colors[k]))
        k += 1
    ax0.set_xticks([k * WIDTH + (WIDTH / 2.5) for k in range(len(bar_labels))])
    ax0.set_xticklabels(feature_labels, rotation=70, size=12)
    ax0.set_title("Alignment Feature Summary")
    k = 0
    WIDTH = 0.5
    for f in fam_tuples:
        print f[0]
        print f[1]
        fam_labels.append(f[1])
        bar_labels.append(ax1.bar(k * WIDTH, f[0], WIDTH, bottom=0, color=even_more_colors[k]))
        k += 1
    ax1.set_xticks([k * WIDTH + (WIDTH / 2.5) for k in range(len(bar_labels))])
    ax1.set_xticklabels(fam_labels, rotation=70, size=12)
    ax1.set_xlim([0, len(fam_labels) / 2.0])
    ax1.set_title("Alignment Function Summary")

    #pl.show()
    pl.savefig(output_file)


def main():
    parser = argparse.ArgumentParser(description='GTFAR status update utility')

    parser.add_argument('-o', '--output-file', required=True, help='Output file name')
    parser.add_argument('input_file', help='Input file name')

    args = parser.parse_args(sys.argv[1:])

    output_file = args.output_file

    FILE = open(args.input_file)
    feature_dict, family_dict = load_bar_data(FILE)
    make_adjacent_bars(feature_dict, family_dict, output_file)


if __name__ == '__main__':
    main()
