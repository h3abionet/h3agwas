#!/usr/bin/env python2.7

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inSNP_acc", help="")
parser.add_argument("--report_acc", help="")
parser.add_argument("--group", help="")
args = parser.parse_args()


def acc_by_maf(inSNP_acc, outSNP_acc, group=''):
    """
    :return:
    """
    datas = {}

    outSNP_acc_out = open(outSNP_acc, 'w')
    outSNP_acc_out.writelines('\\t'.join(
            [group, '(0,0.001]', '(0.001,0.01]', '(0.01,0.02]', '(0.02,0.05]', '(0.05,0.2]', '(0.2,0.5]']) + '\\n')
    outWell_imputed_out_1 = open(outSNP_acc + "_summary.tsv", 'w')
    outWell_imputed_out_1.writelines('\\t'.join(
            [group, '(0,0.001]', '(0.001,0.01]', '(0.01,0.02]', '(0.02,0.05]', '(0.05,0.2]', '(0.2,0.5]',
             'TOTAL']) + '\\n')
    info_datas = open(inSNP_acc).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if "SNP" not in line and dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['extreme_rare'] = []
            datas[dataset]['moderate_rare'] = []
            datas[dataset]['rare'] = []
            datas[dataset]['moderate'] = []
            datas[dataset]['common'] = []
            datas[dataset]['extreme_common'] = []
            datas[dataset]['total'] = 0
        if "SNP" in line and "Rsq" in line:
            idx_exp_freq_a1 = data.index('MAF')
            # idx_conc = data.index("concord_type0")
            idx_conc = data.index("LooRsq")
        else:
            maf = float(data[idx_exp_freq_a1])
            acc = float(data[idx_conc])
            datas[dataset]['total'] += 1
            if maf >= 0.5:
                maf = 1 - maf
            if acc > 0 and acc <= 0.001:
                datas[dataset]['extreme_rare'].append(acc)
                datas[dataset]['total'] += 1
            elif acc > 0.001 and acc <= 0.01:
                datas[dataset]['moderate_rare'].append(acc)
                datas[dataset]['total'] += 1
            elif acc > 0.01 and acc <= 0.02:
                datas[dataset]['rare'].append(acc)
                datas[dataset]['total'] += 1
            elif acc > 0.02 and acc <= 0.05:
                datas[dataset]['moderate'].append(acc)
                datas[dataset]['total'] += 1
            elif acc > 0.05 and acc <= 0.2:
                datas[dataset]['common'].append(acc)
                datas[dataset]['total'] += 1
            elif acc > 0.2 and acc <= 0.5:
                datas[dataset]['extreme_common'].append(acc)
                datas[dataset]['total'] += 1
    try:
        datasets = [str(it) for it in sorted([int(it) for it in datas])]
    except:
        datasets = sorted(datas)
    for dataset in datasets:
        tot = datas[dataset]['total']
        if len(datas[dataset]['extreme_rare']) == 0:
            extreme_rare = 0.0
        else:
            extreme_rare = sum(datas[dataset]['extreme_rare']) / float(len(datas[dataset]['extreme_rare']))
        if len(datas[dataset]['moderate_rare']) == 0:
            moderate_rare = 0.0
        else:
            moderate_rare = sum(datas[dataset]['moderate_rare']) / float(len(datas[dataset]['moderate_rare']))
        if len(datas[dataset]['rare']) == 0:
            rare = 0.0
        else:
            rare = sum(datas[dataset]['rare']) / float(len(datas[dataset]['rare']))
        if len(datas[dataset]['moderate']) == 0:
            moderate = 0.0
        else:
            moderate = sum(datas[dataset]['moderate']) / float(len(datas[dataset]['moderate']))
        if len(datas[dataset]['common']) == 0:
            common = 0.0
        else:
            common = sum(datas[dataset]['common']) / float(len(datas[dataset]['common']))
        if len(datas[dataset]['extreme_common']) == 0:
            extreme_common = 0.0
        else:
            extreme_common = sum(datas[dataset]['extreme_common']) / float(len(datas[dataset]['extreme_common']))

        outSNP_acc_out.write(
            "{}\\t{:3.3f}\\t{:3.3f}\\t{:3.3f}\\t{:3.3f}\\t{:3.3f}\\t{:3.3f}\\n".format(dataset, extreme_rare,
                                                                                       moderate_rare, rare, moderate,
                                                                                       common, extreme_common))

    outSNP_acc_out.close()


args.inSNP_acc = "${inSNP_acc}"
args.outSNP_acc = "${outSNP_acc}"
args.group = "${group}"
if args.inSNP_acc and args.outSNP_acc:
    acc_by_maf(args.inSNP_acc, args.outSNP_acc, args.group)
