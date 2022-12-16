#!/usr/bin/env python

# DESCRIPTION: Extracts DSP counts data for specified TMA and Ab and plots Levey-Jennings chart. It will also execute
# Westgard multirule QC and reject samples that break the ruleset
# USAGE: python westgard.py <filepath> <TMA-name> <Ab-name>

# By Benson Chong

import os
import numpy as np
import pandas as pd
import sys
import re
import matplotlib.pyplot as plt
import openpyxl
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from datetime import date

VERSION = '0.3.0'

#report_path = f'{sys.argv[2]}_{sys.argv[3]}_{date.today().strftime("%Y%m%d")}'
#report_path = re.sub(r'[^\w]', '', report_path).replace(" ", "") + ".pdf"
report_path = sys.argv[4]
# fix -- IF path is not properly formatted *.pdf, then stop script

c = canvas.Canvas(report_path, pagesize=letter)

# RegEx keyword search for TMA cell line name
def cellsearch(tma_input: str):
    cell_line_names = ['Tonsil', 'OvCar8+PARPi', 'BT474', 'MCF7', '468', 'Jurkat', 'Raji', 'THP1', 'PBMC+PHA',
                       '231', '453', 'Spleen', 'Liver', 'T47D']

    tma_oi = None

    for name in cell_line_names:
        if re.search(re.sub(r'[^\w]', '', tma_input), re.sub(r'[^\w]', '', name), re.IGNORECASE) or \
                re.search(re.sub(r'[^\w]', '', name), re.sub(r'[^\w]', '', tma_input), re.IGNORECASE):
            tma_oi = name
            break

    if tma_oi is None:
        raise NameError('Please enter valid TMA cell line name.')

    return tma_oi


# Extract count data from each initial batch dataset
def parse_batches(abcount_path, tma_oi, ab_oi):
    df = pd.read_csv(abcount_path, header=0)
    df = df[(df['name'] == tma_oi) & (df['ProbeName'] == ab_oi)].reset_index()

    # Take mean and sd from table rather than calculate
    #mean, sigma = df['Mean'][0], df['StDev'][0]

    df = df[['name_ProbeName', 'abundance', 'year', 'monthday', 'batch']]
    df = df.sort_values(['year', 'monthday']).reset_index()

    # Following format of previous version
    counts = {}
    ids ={}

    for i, row in df.iterrows():
        id_name = f'Batch_{i+1}_'+str(row['batch'])
        counts[row['batch']] = row['abundance']
        ids[id_name] = row['abundance']

    return [counts, ids]


# For getting annotations, e.g., 'Reference'
def get_key(my_dict, val):
    for k, v in my_dict.items():

        if val == v:
            return k
    return "Key doesn't exist"


# Obtain count values from dictionary
def get_counts(counts):
    x = []
    if not isinstance(counts, dict):
        raise TypeError('Please provide a dict')

    for k, v in counts.items():
        x.append(v)

    x = pd.Series(x)

    return x


# Creates LJ chart
def chart(y, i, segment, plotname, counts):
    x = np.arange(len(y))

    mu, sigma = np.mean(y), np.std(y)
    pos_1s = mu + sigma
    pos_2s = mu + (2 * sigma)
    pos_3s = mu + (3 * sigma)
    neg_1s = mu + (-1 * sigma)
    neg_2s = mu + (-2 * sigma)
    neg_3s = mu + (-3 * sigma)

    fig, ax = plt.subplots()
    plt.yticks([neg_3s, neg_2s, neg_1s, mu, pos_1s, pos_2s, pos_3s])
    ax.set_xticks(range(len(x)))
    ax.set_yticklabels(['-3s', '-2s', '-1s', 'Mean', '1s', '2s', '3s'])
    ax.set_xticklabels(counts.keys())
    plt.xticks(fontsize=5, rotation=45, ha='right')

    plt.axhline(y=mu, color='blue')
    plt.axhline(y=pos_1s, color='green', linestyle=':')
    plt.axhline(y=pos_2s, color='gold', linestyle='--')
    plt.axhline(y=pos_3s, color='red', linestyle='-')
    plt.axhline(y=neg_1s, color='green', linestyle=':')
    plt.axhline(y=neg_2s, color='gold', linestyle='--')
    plt.axhline(y=neg_3s, color='red', linestyle='-')

    plt.figure(i)
    figure = plt.gcf()
    figure.set_size_inches(8,8)
    plt.title("Rule: " + plotname)
    plt.plot(x, y, color='blue', marker='o', markersize=3, linestyle='--')
    plt.plot(x[segment], y[segment], color='red', linewidth=3, marker='o', markersize=6)

    image = re.sub(r'[^\w]', '', (str(i) + "_Rule:" + plotname)).replace(" ", "")
    image += ".png"

    plt.savefig(image, dpi=500)

    c.drawImage(image, 30, 250, width=500, height=500)
    c.showPage()


# multirule QC
def westgard_qc(counts: dict, ids: dict, tma_oi: str, ab_oi: str):
    x = get_counts(counts)

    # refine target name
    ab_oi = ab_oi.replace(r"/", "_")
    ab_oi = ab_oi.replace(" ", "")

    # fix -- redundant section
    mu, sigma = np.mean(x), np.std(x)
    pos_1s = mu + sigma
    pos_2s = mu + (2 * sigma)
    pos_3s = mu + (3 * sigma)
    neg_1s = mu + (-1 * sigma)
    neg_2s = mu + (-2 * sigma)
    neg_3s = mu + (-3 * sigma)

    # Initialize ruleset dict with 0 values
    ruleset = {
        '1_3s': [],
        '2_2s': [],
        'R_4s': [],
        '2of3_2s': [],
        '3_1s': [],
        '6x': [],
        '7T': []
    }

    j = 1
    plotname = '_'+ tma_oi + '_' + ab_oi

    for i in range(len(x)):

        # 1_3s rule
        if x[i] >= pos_3s or x[i] <= neg_3s:
            ruleset['1_3s'].append(list(ids.keys())[list(ids.values()).index(x[i])])
            chart(x, j, [i], '1_3s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
            j += 1

        # Rules spanning multiple data points
        if i > 0:

            # 2_2s rule
            if x[i - 1] > pos_2s and x[i] > pos_2s:
                ruleset['2_2s'].append(list(ids.keys())[list(ids.values()).index(x[i - 1])])
                chart(x, j, range((i - 1), (i + 1)), '2_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            elif x[i - 1] < neg_2s and x[i] < neg_2s:
                ruleset['2_2s'].append(list(ids.keys())[list(ids.values()).index(x[i - 1])])
                chart(x, j, range((i - 1), (i + 1)), '2_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            # R_4s rule
            if x[i - 1] >= pos_2s and x[i] <= neg_2s:
                ruleset['R_4s'].append(list(ids.keys())[list(ids.values()).index(x[i - 1])])
                chart(x, j, range((i - 1), (i + 1)), 'R_4s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            elif x[i] >= pos_2s and x[i - 1] <= neg_2s:
                ruleset['R_4s'].append(list(ids.keys())[list(ids.values()).index(x[i - 1])])
                chart(x, j, range((i - 1), (i + 1)), 'R_4s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

        # If rule spans 3+ data points
        if i >= 3:

            # 3_1s rule
            if x[i - 2] >= pos_1s and x[i - 1] >= pos_1s and x[i] >= pos_1s:
                ruleset['3_1s'].append(list(ids.keys())[list(ids.values()).index(x[i-2])])
                chart(x, j, range((i - 2), (i + 1)), '3_1s | '+str(list(ids.keys())[list(ids.values()).index(x[i])])+plotname, ids)
                j += 1

            elif x[i - 2] <= neg_1s and x[i - 1] <= neg_1s and x[i] <= neg_1s:
                ruleset['3_1s'].append(list(ids.keys())[list(ids.values()).index(x[i-2])])
                chart(x, j, range((i - 2), (i + 1)), '3_1s | '+str(list(ids.keys())[list(ids.values()).index(x[i])])+plotname, ids)
                j += 1

            # 2of3_2s rule
            if (x[i - 2] > mu and x[i - 1] > mu and x[i] > mu) and x[i - 2] > pos_2s and x[i - 1] > pos_2s:
                ruleset['2of3_2s'].append(list(ids.keys())[list(ids.values()).index(x[i-2])])
                chart(x, j, range((i - 2), (i + 1)), '2of3_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            elif (x[i - 2] > mu and x[i - 1] and x[i] > mu > mu) and x[i - 2] > pos_2s and x[i] > pos_2s:
                ruleset['2of3_2s'].append(list(ids.keys())[list(ids.values()).index(x[i-2])])
                chart(x, j, range((i - 2), (i + 1)), '2of3_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            elif (x[i - 2] > mu and x[i - 1] > mu and x[i] > mu) and x[i - 1] > pos_2s and x[i] > pos_2s:
                ruleset['2of3_2s'].append(list(ids.keys())[list(ids.values()).index(x[i-2])])
                chart(x, j, range((i - 2), (i + 1)), '2of3_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            elif (x[i - 2] > mu and x[i - 1] > mu and x[i] > mu) and x[i - 2] > pos_2s and x[i - 1] < neg_2s:
                ruleset['2of3_2s'].append(list(ids.keys())[list(ids.values()).index(x[i-2])])
                chart(x, j, range((i - 2), (i + 1)), '2of3_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            elif (x[i - 2] > mu and x[i - 1] > mu and x[i] > mu) and x[i - 2] > pos_2s and x[i] < neg_2s:
                ruleset['2of3_2s'].append(list(ids.keys())[list(ids.values()).index(x[i - 2])])
                chart(x, j, range((i - 2), (i + 1)), '2of3_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

            elif (x[i - 2] > mu and x[i - 1] > mu and x[i] > mu) and x[i - 1] > pos_2s and x[i] < neg_2s:
                ruleset['2of3_2s'].append(list(ids.keys())[list(ids.values()).index(x[i-2])])
                chart(x, j, range((i - 2), (i + 1)), '2of3_2s | '+ list(ids.keys())[list(ids.values()).index(x[i])]+plotname, ids)
                j += 1

        # If rule spans 6+ data points
        if i >= 6:

            # 6x rule
            if x[i - 5] > mu and x[i - 4] > mu and x[i - 3] > mu and x[i - 2] > mu and x[i - 1] > mu and x[i] > mu:
                ruleset['6x'].append(list(ids.keys())[list(ids.values()).index(x[i - 5])])
                chart(x, j, range((i - 5), (i + 1)), '6x | '+ list(ids.keys())[list(ids.values()).index(x[i - 5])]+plotname, ids)
                j += 1

            elif x[i - 5] < mu and x[i - 4] < mu and x[i - 3] < mu and x[i - 2] < mu and x[i - 1] < mu and x[i] < mu:
                ruleset['6x'].append(list(ids.keys())[list(ids.values()).index(x[i - 5])])
                chart(x, j, range((i - 5), (i + 1)), '6x | '+ list(ids.keys())[list(ids.values()).index(x[i - 5])]+plotname, ids)
                j += 1

        if i >= 7:

            # 7T rule -- alternative: if equal to or greater/less than
            if x[i - 6] < x[i - 5] < x[i - 4] < x[i - 3] < x[i - 2] < x[i - 1] < x[i]:
                ruleset['7T'].append(list(ids.keys())[list(ids.values()).index(x[i -6])])
                chart(x, j, range((i - 6), (i + 1)), '7T | '+ list(ids.keys())[list(ids.values()).index(x[i - 5])]+plotname, ids)
                j += 1

            elif x[i - 6] > x[i - 5] > x[i - 4] > x[i - 3] > x[i - 2] > x[i - 1] > x[i]:
                ruleset['7T'].append(list(ids.keys())[list(ids.values()).index(x[i-6])])
                chart(x, j, range((i - 6), (i + 1)), '7T | '+ list(ids.keys())[list(ids.values()).index(x[i - 5])]+plotname, ids)
                j += 1

    y = 750 # starting y-position for PDF writing
    rule_broke = False
    for k, v in ruleset.items():
        if v:
            rule_broke = True
            print('Rule', k, ' is broken starting at batch(es): ', v, 'for ', tma_oi, '|', ab_oi)
            y -= 15

    return rule_broke


def main():
    tma_oi = cellsearch(sys.argv[2])
    batches = parse_batches(sys.argv[1], tma_oi, str(sys.argv[3]))
    counts, ids = batches[0], batches[1]
    print('Running Westgard QC...')
    rule_broke = westgard_qc(counts, ids, tma_oi, str(sys.argv[3]))
    if rule_broke:
        c.save()
    else:
        print('No rules broken!')

if __name__ == '__main__':
    main()