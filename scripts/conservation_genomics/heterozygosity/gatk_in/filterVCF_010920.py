'''
Custom filtering for VCF files

NOTE: This script was developed for personal use and is ever-evolving. It is not
guaranteed to work on any and all possible VCF files. Use with discretion.

Input: compressed raw VCF
Output: filtered VCF (prints to screen)
- Sites failing filters are marked as FAIL_? in the 7th column
- Sites that pass site-level filtering go on to genotype filtering
- Filtered genotypes are changed to './.', all others reported as is

Possible usage:
SCRIPT=filterVCF.py
python ${SCRIPT} myfile.vcf.gz | bgzip > myfile_filtered.vcf.gz
tabix -p vcf myfile_filtered.vcf.gz

'''

import sys
import gzip
import re

vcf_file = sys.argv[1]
VCF = gzip.open(vcf_file, 'rt')

# Minimum genotype quality score
minGQ=20.0

# Minimum individual depth (1/3x mean)
minD = int(sys.argv[2])

# Maximum individual depth (2x mean)
maxD = int(sys.argv[3])

# Individual genotype filtering function
#  - Genotypes failing filters are set to missing (./.)
#  - Applies individual min and max depth filters (thresholds set above)
#  - Applies minimum GQ filter (minGQ threshold set above)
#  - Filters heterozygotes if the allele balance (REF/DP) is <20% or >80%
#  - 'sample' is the sample name
#  - 'GT_entry' is the entire genotype entry for that individual (typically GT:AD:DP:GQ)
#  - 'ADpos' is the position of the AD field in FORMAT (as determined below)
#  - 'DPpos' is the position of the DP field in FORMAT (as determined below)
#  - 'GQpos' is the position of the GQ field in FORMAT (as determined below)

def GTfilter(sample, GT_entry, ADpos, DPpos, GQpos):
    if GT_entry[:1]=='.' : return GT_entry
    else:
        gt=GT_entry.split(':')
        if gt[0] in ('0/0','0/1','1/1') and gt[DPpos]!='.' and gt[GQpos]!='.':
            DP=int(gt[DPpos])
            GQ=float(gt[GQpos])
            if GQ>=minGQ and minD<=DP<=maxD:
                if gt[0]=='0/1':
                    REF=float(gt[ADpos].split(',')[0])
                    AB=float(REF/DP)
                    if 0.2<=AB<=0.8: return GT_entry
                    else: return './.:' + ':'.join(gt[1:])
                else: return GT_entry
            else: return './.:' + ':'.join(gt[1:])
        else: return './.:' + ':'.join(gt[1:])


# Get list of samples in VCF file
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break

# Return to the beginning of file
VCF.seek(0)

# Write pre-existing header lines & add new lines describing filters being applied to FILTER column
for line0 in VCF:
    if line0.startswith('#'):
        if line0.startswith('##FORMAT'):
            sys.stdout.write('##FILTER=<ID=FAIL_badRef,Description="Reference allele not one of [A,C,G,T].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_badAlt,Description="Alternate allele not one of [A,C,G,T,.].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noINFO,Description="No INFO present.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noADi,Description="AD not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noDPi,Description="DP not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noGQi,Description="GQ not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noGT,Description="No called genotypes remain after filtering.">\n')
            sys.stdout.write(line0)
            break
        else: sys.stdout.write(line0)


# Go through VCF file line by line to apply filters
for line0 in VCF:
    if line0.startswith('#'):
        sys.stdout.write(line0); continue
    line=line0.strip().split('\t')

### Site filtering:
### Keep any filters that have already been applied
    filter=[]
    if line[6] not in ('.', 'PASS'):
        filter.append(line[6])

### Check REF allele
    if line[3] not in ('A','C','G','T'):       
        filter.append('FAIL_badRef') 

### Check ALT allele
    if line[4] not in ('A','C','G','T','.'):
        filter.append('FAIL_badAlt') 

### Access INFO field annotations
    if ';' in line[7]:
        INFO=line[7].split(';')
        d=dict(x.split('=') for x in INFO)
    else:
        INFO=line[7]
        if '=' in INFO:
            d={INFO.split('=')[0]:INFO.split('=')[1]}
        else: filter.append('FAIL_noINFO')

### Get the position of AD, DP, GQ in genotype fields
    if 'AD' in line[8]:
        ADpos=line[8].split(':').index('AD')
    else: filter.append('FAIL_noADi')

    if 'DP' in line[8]:
        DPpos=line[8].split(':').index('DP')
    else: filter.append('FAIL_noDPi')

    if 'GQ' in line[8]:
        ff=line[8].split(':')
        GQpos=[ff.index(x) for x in ff if 'GQ' in x][0]
    else: filter.append('FAIL_noGQi')

### If any filters failed, write out line and 'continue' (i.e. skip genotype filtering)
    if filter!=[]:
        sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:])) ) ; continue

### Genotype filtering:
    GT_list=[]
    for i in range(0,len(samples)):
        GT=GTfilter(samples[i],line[i+9],ADpos,DPpos,GQpos)
        GT_list.append(GT)

### Recalculate AC, AN, AF for INFO
    REF=2*[x[:3] for x in GT_list].count('0/0') + [x[:3] for x in GT_list].count('0/1')
    ALT=2*[x[:3] for x in GT_list].count('1/1') + [x[:3] for x in GT_list].count('0/1')
    if REF+ALT==0:
        filter.append('FAIL_noGT')
        sys.stdout.write('%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:9]), '\t'.join(GT_list)) ); continue    
    d['AC']=ALT
    d['AN']=REF+ALT
    d['AF']=round(float(ALT)/(float(REF)+float(ALT)), 4)

### Write out new line
    if filter==[]:
        filter.append('PASS')
    sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )


# Close files and exit
VCF.close()
exit()

