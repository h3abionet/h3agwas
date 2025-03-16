#!/usr/bin/env python2.7

import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--infoFiles", help="")
parser.add_argument("--outWell_imputed", help="")
parser.add_argument("--outSNP_acc", help="")
parser.add_argument("--infoCutoff", help="")
parser.add_argument("--inWell_imputed", help="")
parser.add_argument("--inSNP_acc", help="")
parser.add_argument("--report_acc", help="")
parser.add_argument("--ldFiles", help="")
parser.add_argument("--report_ld", help="")

args = parser.parse_args()

def filter_info(infoFiles, infoCutoff, outWell_imputed, outSNP_acc):
    """
    Return:
        well_imputed: certainy >= 1
        SNP_concordance: concord_type0 != -1
    """
    well_imputed = {}
    SNP_concordance = {}
    count = 0
    infoFiles = infoFiles.split(',')
    datas = {}
    header = []
    outWell_imputed_out = open(outWell_imputed, 'w')
    outWell_imputed_snp_out = open(outWell_imputed+"_snp", 'w')
    outSNP_accuracy_out = open(outSNP_acc, 'w')
    for infoFile in infoFiles:
        infoFile = infoFile.strip().split('==')
        dataset = infoFile[0]
        info = infoFile[1]
        well_imputed[dataset] = []
        SNP_concordance[dataset] = []
        # outWell_imputed_out_dataset = open(dataset+'_'+outWell_imputed, 'w')
        for line in open(info):
            data = line.strip().split()
            if "snp_id" in line and "info" in line:
                if len(header) == 0:
                    header = data
                    info_idx = header.index("info")
                    conc_idx = header.index("concord_type0")
                    outWell_imputed_out.writelines(' '.join([dataset]+data)+'\n')
                    # outWell_imputed_out_dataset.writelines(' '.join([dataset]+data)+'\n')
                    outWell_imputed_snp_out.writelines(data[1]+'\n')
                    outSNP_accuracy_out.writelines(' '.join([dataset]+data)+'\n')
            else:
                if float(data[info_idx]) >= float(infoCutoff):
                    outWell_imputed_out.writelines(' '.join([dataset]+data)+'\n')
                    # outWell_imputed_out_dataset.writelines(' '.join([dataset]+data)+'\n')
                    outWell_imputed_snp_out.writelines(data[1]+'\n')
                if float(data[conc_idx]) != float(-1):
                    outSNP_accuracy_out.writelines(' '.join([dataset]+data)+'\n')
                count += 1
        # outWell_imputed_out_dataset.close()
    outWell_imputed_out.close()
    outWell_imputed_snp_out.close()
    outSNP_accuracy_out.close()

def well_imputed_by_maf(inWell_imputed, outWell_imputed):
    """

    :return:
    """
    datas = {}
    outWell_imputed_out = open(outWell_imputed, 'w')
    outWell_imputed_out.writelines('\t'.join(['CHRM', 'MAF>=5%', 'MAF>=1%', 'MAF 1-5%', 'TOTAL'])+'\n')
    info_datas = open(inWell_imputed).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['1_5'] = []
            datas[dataset]['1'] = []
            datas[dataset]['5'] = []
            datas[dataset]['total'] = 0
        if "snp_id" in line and "info" in line:
            idx_exp_freq_a1 = data.index('exp_freq_a1')
        else:
            maf = float(data[idx_exp_freq_a1])
            datas[dataset]['total'] += 1
            if maf >= 0.5:
                maf = 1-maf
            if maf >= 0.01:
                datas[dataset]['1'].append(maf)
                if maf >= 0.05:
                    datas[dataset]['5'].append(maf)
                if maf <= 0.05 and maf >= 0.01:
                    datas[dataset]['1_5'].append(maf)
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        if tot == 0:
            outWell_imputed_out.write("dataset {} is empty (tot=0)".format(dataset))
        else:
            outWell_imputed_out.writelines('\t'.join([dataset, str(format(len(datas[dataset]['5'])/1000000., '0,.1f'))+'M ('+str(len(datas[dataset]['5']) * 100/tot)+'%)', str(format(len(datas[dataset]['1'])/1000000., '0,.1f'))+'M ('+str(len(datas[dataset]['1']) * 100/tot)+'%)', str(format(len(datas[dataset]['1_5'])/1000000., '0,.1f'))+'M ('+str(len(datas[dataset]['1_5']) * 100/tot)+'%)', str(format(tot, '0,.0f'))])+'\n')
    outWell_imputed_out.close()

def acc_by_maf(inSNP_acc, outSNP_acc):
    """
    :return:
    """
    datas = {}

    outSNP_acc_out = open(outSNP_acc, 'w')
    outSNP_acc_out.writelines('\t'.join(['CHRM', 'MAF>=5%', 'MAF>=1%', 'MAF 1-5%', 'TOTAL'])+'\n')
    info_datas = open(inSNP_acc).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['1_5'] = []
            datas[dataset]['1'] = []
            datas[dataset]['5'] = []
            datas[dataset]['total'] = 0
        if "snp_id" in line and "info" in line:
            idx_exp_freq_a1 = data.index('exp_freq_a1')
            # idx_conc = data.index("concord_type0")
            idx_conc = data.index("r2_type0")
        else:
            maf = float(data[idx_exp_freq_a1])
            acc = float(data[idx_conc])
            datas[dataset]['total'] += 1
            if maf >= 0.5:
                maf = 1-maf
            if maf >= 0.01:
                datas[dataset]['1'].append(acc)
                if maf >= 0.05:
                    datas[dataset]['5'].append(acc)
                if maf <= 0.05 and maf >= 0.01:
                    datas[dataset]['1_5'].append(acc)
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        print tot

        if len(datas[dataset]['5'])==0:
            maf_5 = 0.0
        else:
            maf_5 = sum(datas[dataset]['5'])/float(len(datas[dataset]['5']))
        if len(datas[dataset]['1'])==0:
            maf_1 = 0.0
        else:
            maf_1 = sum(datas[dataset]['1'])/float(len(datas[dataset]['1']))
        if len(datas[dataset]['1_5'])==0:
            maf_1_5 = 0.0
        else:
            maf_1_5 = sum(datas[dataset]['1_5'])/float(len(datas[dataset]['1_5']))

        outSNP_acc_out.write("{0:.3f}\t{0:.3f}\t{0:.3f}\n".format(maf_5,maf_1, maf_1_5))
        #outSNP_acc_out.writelines('\t'.join([dataset, str(format(maf_5), '0,.3f')), str(format(maf_1), '0,.3f')), str(format(maf_1_5), '0,.3f')), str(format(tot, '0,.0f'))])+'\n')

    outSNP_acc_out.close()

def ld_by_maf(ldFiles, report_ld, inWell_imputed, infoCutoff):
    """

    :return:
    """
    maf_data = {}
    print 'Reading', inWell_imputed
    # inWell_imputed_data = open(inWell_imputed).readlines()
    for line in open(inWell_imputed):
        data = line.split(' ')
        dataset = data[0]
        if "snp_id" in line and "info" in line:
            idx_exp_freq_a1 = data.index('exp_freq_a1')
            rs_id = data.index('rs_id')
        if dataset not in maf_data:
            maf_data[dataset] = {}
        maf_data[dataset][data[rs_id]] = data[idx_exp_freq_a1]
    datas = {}
    ld_data = {}
    report_ld_out = open(report_ld, 'w')
    report_ld_out.writelines('\t'.join(['Population', 'MAF>=5%', 'MAF>=1%', 'MAF 1-5%', 'TOTAL'])+'\n')
    ldFiles = ldFiles.split(',')
    not_ = 0
    in_ = 0
    header = []
    infoCutoff = float(infoCutoff)
    for ldFile in ldFiles:
        ldFile = ldFile.strip().split('==')
        dataset =ldFile[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['1_5'] = 0
            datas[dataset]['1'] = 0
            datas[dataset]['5'] = 0
            datas[dataset]['ALL'] = set()
            datas[dataset]['total'] = len(maf_data[dataset])
        ld = ldFile[1]
        print 'Reading', ld
        ld_data[dataset] = []
        for line in open(ld):
            data = line.strip().split()
            if "SNP_A" in line and "SNP_B" in line:
                if len(header) == 0:
                    header = data
                    snpA_idx = header.index("SNP_A")
                    snpB_idx = header.index("SNP_B")
                    r2_idx = header.index("R2")
            else:
                if float(data[r2_idx]) >= infoCutoff:
                    snpA = data[snpA_idx]
                    snpB = data[snpB_idx]
                    for snp in [snpA, snpB]:
                        if snp not in datas[dataset]['ALL']:
                            try:
                                datas[dataset]['ALL'].add(snp)
                                maf = float(maf_data[dataset][snp])
                                if maf >= 0.5:
                                    maf = 1-maf
                                if maf >= 0.01:
                                    datas[dataset]['1'] += 1
                                    if maf >= 0.05:
                                        datas[dataset]['5'] += 1
                                    if maf <= 0.05 and maf >= 0.01:
                                        datas[dataset]['1_5'] += 1
                                # in_ += 1
                            except:
                                continue
                                # not_ += 1
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        # tot = len(datas[dataset]['ALL'])
        # print len(datas[dataset]['ALL']), tot
        if tot == 0:
            report_ld_out.write("dataset {} is empty (tot=0)".format(dataset))
        else:
            report_ld_out.writelines('\t'.join([dataset, str(format(datas[dataset]['5']/1000000., '0,.1f'))+'M ('+str(datas[dataset]['5'] * 100/tot)+'%)', str(format(datas[dataset]['1']/1000000., '0,.1f'))+'M ('+str(datas[dataset]['1'] * 100/tot)+'%)', str(format(datas[dataset]['1_5']/1000000., '0,.1f'))+'M ('+str(datas[dataset]['1_5'] * 100/tot)+'%)', str(format(tot, '0,.0f'))])+'\n')
    report_ld_out.close()

if args.infoFiles and args.infoCutoff:
    filter_info(args.infoFiles, args.infoCutoff, args.outWell_imputed, args.outSNP_acc)
if args.inWell_imputed and args.outWell_imputed:
    well_imputed_by_maf(args.inWell_imputed, args.outWell_imputed)
if args.inSNP_acc and args.outSNP_acc:
    acc_by_maf(args.inSNP_acc, args.outSNP_acc)
if args.report_ld and args.ldFiles:
    ld_by_maf(args.ldFiles, args.report_ld, args.inWell_imputed, args.infoCutoff)
