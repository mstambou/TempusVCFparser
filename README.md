# tempusVCFparser
Github repository for tempus coding challenge for the Bioinformatics Scientist position.
There is one main python3 script that takes a VCF file as input pased by the -i option as well as an output directory where it will store the parsed output in table format, specified by the -o option. Example to run this script:

```python3
python3 parseVCF.py -i Challenge_data_\(1\).vcf -o .
```

The python packages imported by this script are listed here:

```
    import sys, getopt
    import random
    import pandas as pd
    import requests
    import json
    import os
    from collections import Counter
```

    

```python

def get_info(info):
    """
    simple function that will parse the info filed into dictionary maps
    """
    return {item.split('=')[0]:item.split('=')[1] for item in info.split(';')}

def get_mostDeleterious(varType, AC, AF):
    """
    function that will return the most deleterious variant effect if more than one is
    found. In case there are both insertions and deletions it will return the one 
    that is most frequently found, in case both insertions and deletions are mentioned
    as frequent it will return the one that has more allele count, and in case they 
    have the same allele count then it will randomly choose between an insertion
    or a deletion with equal probabilities.
    """        
    ranks = ['del', 'ins', 'complex', 'mnp', 'snp']
    variants = varType.split(',')
    allele_counts = [int(item) for item in AC.split(',')]
    allele_frequencies = [float(item) for item in AF.split(',')]

    varIdx_dic = {}
    for i, var in enumerate(variants):
        if var not in varIdx_dic:
            varIdx_dic[var] = [i]
        else:
            varIdx_dic[var].append(i)

    if 'del' in variants and 'ins' in variants:
        if variants.count('del') > variants.count('ins'):
            mostDelVarType = 'del'
        elif variants.count('del') < variants.count('ins'):
            mostDelVarType = 'ins'
        elif variants.count('del') == variants.count('ins'):
            max_del_AC =  max([item for item in [allele_counts[i] for i in varIdx_dic['del']] ])
            max_ins_AC =  max([item for item in [allele_counts[i] for i in varIdx_dic['del']] ])
            if max_del_AC > max_ins_AC:
                mostDelVarType = 'del'
            if max_ins_AC > max_del_AC:
                mostDelVarType = 'ins'
            if max_ins_AC == max_del_AC:
                mostDelVarType = random.choice(['ins', 'del'])
    else:
        for var in ranks:
            if var in variants:
                mostDelVarType = var
                break

    mostDelVarType_idx = varIdx_dic[mostDelVarType]
    if len(mostDelVarType_idx) == 1:
        idx = mostDelVarType_idx[0]
        allele_count, allele_frequency = allele_counts[idx], allele_frequencies[idx]
    else:
        max_count = -1
        idx = -1
        for i in mostDelVarType_idx:
            if allele_counts[i] > max_count:
                max_count = allele_counts[i]
                idx = i
        allele_count, allele_frequency = allele_counts[idx], allele_frequencies[idx]

    return mostDelVarType, allele_count, allele_frequency, idx


def get_ExACDB_info(exac_col):
    """
    function that will connect to the ExAC database via their REST API,
    and will request metadata for the variants that are identified in the 
    given VCF, for more details concerning the ExAC REST API please refer
    to: http://exac.hms.harvard.edu/#rest-bulk-variant
    """
    data = {item:item for item in exac_col}
    url = 'http://exac.hms.harvard.edu/rest/bulk/variant'

    print('retrieving query requests from ExAC db REST API, pleaese be patient ...')
    query = requests.post(url, json = data)
    exac_dic = json.loads(query.text)

    consequence, alleleFreq, ENSGs, ENSTs = [], [], [], []
    for k in exac_dic:
        if exac_dic[k]['consequence'] != None:
            consequence.append(','.join(list(exac_dic[k]['consequence'].keys())))        
        else:
            consequence.append('NA')        
        if 'allele_freq' in exac_dic[k]['variant']:
            alleleFreq.append(exac_dic[k]['variant']['allele_freq'])        
        else:
            alleleFreq.append('NA')        
        if 'genes' in exac_dic[k]['variant']:
            ENSGs.append(','.join(exac_dic[k]['variant']['genes']))        
        else:
            ENSGs.append('NA')    
        if 'transcripts' in exac_dic[k]['variant']:
            ENSTs.append(','.join(exac_dic[k]['variant']['transcripts']))        
        else:
            ENSTs.append('NA')
            
    return (consequence, alleleFreq, ENSGs, ENSTs)
    

def parseVCF(vcf_f):
    """
    main function for the VCF parser that will read the VCF input file line by line
    and will parse the necesarry fields, will call the necesarry functions predefined
    at the begining of this script and eventually will write the parsed VCF file 
    into a tab separated table file at the specified output directory.
    """
    
    with open(vcf_f, 'r') as in_f:
        out_df = pd.DataFrame(columns = ['CHROM', 'POS', 'refAllele', 'alternativeAllele', 'variantType', 'depthOfCoverage', 'nReadsVariant', 'variantReads_%', 'ExAC_col'])
        r_count = 0
        print('reading the vcf file ...')
        for line in in_f:
            if line.startswith('#') == False:
                l = line.strip().split('\t')
                
                r_count += 1
                chrom, pos, ref_allele, alt_allele = l[0], l[1], l[3], l[4]
                info = l[7]
                info_dic = get_info(info)
                varType = info_dic['TYPE']
                depthOfCov = int(info_dic['DP'])
                AC = info_dic['AC']
                AF = info_dic['AF']

                FORMAT, ref, alt = l[8], l[9], l[10]

                ref_format_dic = {k:v for k,v in zip(FORMAT.split(':'), ref.split(':'))}
                alt_format_dic = {k:v for k,v in zip(FORMAT.split(':'), alt.split(':'))}

                if len(varType.split(',')) == 1:
                    varType = varType.split(',')
                    mostDelVarType, allele_count, allele_frequency = varType[0], int(AC), float(AF)*100
                    nAlleleReads = int(alt_format_dic['AO'])
                    #nAlleleReadPercent = (int(alt_format_dic['AO'])/float(alt_format_dic['DP']))*100
                    nAlleleReadPercent = round((nAlleleReads/float(depthOfCov))*100, 3)
                else:
                    mostDelVarType, allele_count, allele_frequency, mostDelVarType_idx = get_mostDeleterious(varType, AC, AF)
                    nAlleleReads = int(alt_format_dic['AO'].split(',')[mostDelVarType_idx] )
                    #nAlleleReadPercent = (nAlleleReads/float(alt_format_dic['DP']))*100
                    nAlleleReadPercent = round((nAlleleReads/float(depthOfCov))*100, 3)
                allele_frequency = allele_frequency*100
                exac_col = f'{chrom}-{pos}-{ref_allele}-{alt_allele}'
                out_df.loc[r_count] = [chrom, pos, ref_allele, alt_allele, mostDelVarType, depthOfCov, nAlleleReads, nAlleleReadPercent, exac_col]
    exac_col = list(out_df['ExAC_col'])
    
    extra_cols = get_ExACDB_info(exac_col)
    out_df['consequence'] = extra_cols[0]
    out_df['alleleFreq'] = extra_cols[1]
    out_df['ensembl_gid'] = extra_cols[2]
    out_df['ensembl_tid'] = extra_cols[3]
    print('writing parsed VCF into table file ...')
    
    return out_df
    
    
if __name__ == "__main__":
    
    import sys, getopt
    import random
    import pandas as pd
    import requests
    import json
    import os
    from collections import Counter


    vcf_f = '/data/mstambou/tempus/Challenge_data_(1).vcf'
    out_dir = '/data/mstambou/tempus/'
    fname = vcf_f.rsplit('/', 1)[-1].rsplit('.', 1)[0]
    
    if not out_dir.endswith('/'):
        out_dir = out_dir + '/'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
        
    out_df = parseVCF(vcf_f)
    varCounts = Counter(out_df['variantType'])
    total = sum(varCounts.values())
    print(f'There were a total of {total}, variants identified, composed of {len(varCounts)}, types')
    print( '\n'.join(['\t:\t'.join([i[0],str(i[1])]) for i in list(varCounts.items()) ]) )
    
    
    import plotly.express as px
    variantsSummary_df = pd.DataFrame(columns = ['variant', 'counts'])
    varCountTuples = list(varCounts.items())
    variantsSummary_df['variant'] = [i[0] for i in varCountTuples]
    variantsSummary_df['counts'] = [i[1] for i in varCountTuples]
    fig = px.pie(variantsSummary_df, values='counts', names='variant', title=f'Summary of Variant compositions in the sample file {out_fname}')
    fig.show()
```

    reading the vcf file ...
    retrieving query requests from ExAC db REST API, pleaese be patient ...
    writing parsed VCF into table file ...
    There were a total of 6977, variants identified, composed of 5, types
    snp	:	5377
    complex	:	112
    del	:	919
    ins	:	531
    mnp	:	38




