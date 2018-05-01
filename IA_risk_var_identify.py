'''
Damaging variants filtering pipeline used in IA study.
The control used in this analysis were come from
the public databases.
The method we used for normolized the allele count in to the 
family unit to control the influence of relateness to the case- 
control analysis was also included in the code; 
The input file is the variants vcf file with annotation information
'''

import sys
import re
import tempfile
import configparser 
import subprocess

sample_config = sys.argv[1]

config = configparser.RawConfigParser()
config.read(sample_config)

input = config.get('MAF1000G', 'input')
pop_1000G = config.get('MAF1000G', 'pop')
cutoff_1000G = config.getfloat('MAF1000G', 'cutoff')
pop_exac = config.get('exacAF', 'pop')
cutoff_exac = config.getfloat('exacAF', 'cutoff')
ped = config.get('pValue', 'ped_file')
control = config.get('pValue', 'control_file')
cutoff_cadd = config.get('caddPrediction', 'cutoff')

def get_1000G_AF(pop, info_str):
    maf = re.findall(';T{}=(\d+\.\d+);'.format(pop), info_str)
    if not maf:
        maf = re.findall(';T{}=(\d+\.\d+),'.format(pop), info_str)
        if not maf:
            return 0
    else:
        return float(maf[0])

def get_exac_af(pop, info_str):
    if pop == "MAF":
        ac = re.findall(';EX_AC=(\d+);',info_str)
        an = re.findall(';EX_AN=(\d+);',info_str)
        #print(ac, an)
        if not an:
            return 0, 0
        elif not ac:
            ac = re.findall(';EX_AC=(\d+),',info_str)
            if not ac:
                return 0, int(an[0])
            else:
                return int(ac[0]), int(an[0])
        else:
            return int(ac[0]), int(an[0])
    else:
        ac = re.findall(';EX_{}=(\d+);'.format(pop), info_str)
        an = re.findall(';EX_{}C=(\d+);'.format(pop), info_str)
        #print(ac, an)
        if not an:
            return 0, 0
        elif not ac:
            ac = re.findall(';EX_{}=(\d+),'.format(pop), info_str)
            #print(ac, an)
            if not ac:
                return 0, int(an[0])
            else:
                return int(ac[0]), int(an[0])
        else:
            return int(ac[0]), int(an[0])

def read_pedigree(ped_file):
    ped_hash = {}
    uniq_families = []
    for line in open(ped_file, "r"):
        data = line.rstrip().split('\t')
        sample = data[0]
        family = int(data[1])
        ped_hash[sample] = family
        if family not in uniq_families:
            uniq_families.append(family)
    num_families = int(len(uniq_families))
    return ped_hash, num_families


if pop_1000G != "MAF" and pop_1000G != "EUR" and \
    pop_1000G != "ASN" and pop_1000G != "AFR" and pop_1000G != "AMR":
    print("Wong 1000G population")
    exit(0)

if pop_exac != "MAF" and pop_exac !="AMR" and pop_exac != "AFR" and pop_exac != "EAS" \
    and pop_exac != "FIN" and pop_exac != "NFE" and pop_exac != "SAS":
    print("Wong exac population")
    exit(0)

if float(cutoff_1000G) > 1.0:
    print("Wrong 1000G MAF cutoff")
    exit(0)
if float(cutoff_exac) > 1.0:
    print("Wrong exac population")
    exit(0)

with open('lessCommonIn1000G.vcf', 'wt') as maf_filter:
    for line in open(input, 'r'):
        if line[0] ==  "#":
            #print(line.strip())
            maf_filter.write(line)
        else:
            data = line.strip().split('\t')
            maf_1000G = get_1000G_AF(pop_1000G, data[7])
            exac_ac, exac_an = get_exac_af(pop_exac, data[7])
            if maf_1000G < float(cutoff_1000G) and exac_an == 0:
                maf_filter.write(line)
            elif maf_1000G < float(cutoff_1000G) and float(exac_ac/exac_an) < float(cutoff_exac):
                #print(maf_all)
                #print(line.strip())
                maf_filter.write(line)

maf_filter.close()

def get_1000G_AN(control_file):
    tabix_command1 = 'tabix -H {} | grep -v ^##'.format(control_file)
    proc = subprocess.Popen(tabix_command1, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    proc.wait()
    header = proc.stdout.read()
    line = header.strip().split("\t")
    #print(line)
    an_control = int(len(line) - 9)
    return an_control*2


def get_1000G_AC(region1, region2, region3, control_file):
    region = region1 + ":" + region2 + "-" + region3
    #print(region)
    tabix_command = 'tabix {} {}'.format(control_file, region)
    proc = subprocess.Popen(tabix_command, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    proc.wait()
    out = proc.stdout.read()
    if out != '':
        #print(out)
        info = out.strip().split('\t')
    #    print(info[7])
        ac_control = re.findall(';AC=(\d+);', info[7])
        if not ac_control:
            ac_control = re.findall(';AC=(\d+),', info[7])
            if not ac_control:
                return 0
            else:
                return int(ac_control[0])
        else:
            return int(ac_control[0])
    else:
        return 0

def count_samples(vcf_file, ped_hash, num_families):
    num_sample = len(ped_hash)
    col_ped_hash = {}
    for line in open(vcf_file, "r"):
        if line[:2] == '##':
            continue

        data = line.rstrip().split('\t')
        if line[0] == '#':
            # read sample name 
            for i in range(9, len(data)):
                col_ped_hash[i] = ped_hash[data[i]]
            print(line.strip(), "family_AC", "family_AN", "AC", "AN", "1000G_AC", "1000G_AN", sep = "\t") 
            continue
        if len(data) != 9 + num_sample:
            print("*ERROR*: missing column in line:\n{}".format(line.rstrip()))
            exit(0)
        
        #get 1000G control AC
        no_control_ac = get_1000G_AC(data[0], data[1], data[1], control_file)
        #print(no_control_ac)

        ac_list = []
        sample_num_list = [[] for i in range(num_families)]
        for i in range(9, len(data)):
            gt = data[i].strip().split(':')[0]
            #print(gt)
            if gt[0] == '0' and gt[-1] == '0':
                count = 0
            elif gt[0] == '.' and gt[-1] == '.':
                count = 0
            elif gt[0] =='1' and gt[-1] == '1':
                count = 2
            else:
                count = 1
            #print(count)
            sample_num_list[col_ped_hash[i]].append(count)
            ac_list.append(count)
            #if count == 2:
                #ac_list.append(1)
            #else:
                #ac_list.append(count)
        total_ac = sum(ac_list)
        an = 2*int(len(data)-9) - total_ac
        total_num_case = 0
        rewight_per = 0
        rewight_all = 0
        n1 = 0
        #total_families = 0
        reweight_result = []
        for i in range(len(sample_num_list)):
            num_family_carry = sum(sample_num_list[i])
            #print(str(num_family_carry))
            total_num_case = len(sample_num_list[i])
            #print(str(total_num_case))
            #total_families += int(len(sample_num_list[i])
            rewight_per = num_family_carry/(total_num_case*2.0)
            reweight_result.append(rewight_per)

        reweight_all = sum(reweight_result)*2
        n1 = num_families*2 - reweight_all
        an_control_left = no_an_control - no_control_ac 
        print(line.strip(), '\t{0:4.2f}\t{1:4.2f}\t{2:4d}\t{3:4d}\t{4:4d}\t{5:4d}'.format(
            reweight_all, n1,  total_ac, an, no_control_ac, an_control_left))

ped_file = ped
vcf_file = "lessCommonIn1000G.vcf"
control_file = control

no_an_control = get_1000G_AN(control_file)
ped_hash, num_families = read_pedigree(ped_file)
count_samples(vcf_file, ped_hash, num_families)
 
