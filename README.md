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



<div>
        
        
            <div id="407fb868-8c7a-4e1d-a35e-fcb5172ced9a" class="plotly-graph-div" style="height:525px; width:100%;"></div>
            <script type="text/javascript">
                require(["plotly"], function(Plotly) {
                    window.PLOTLYENV=window.PLOTLYENV || {};
                    
                if (document.getElementById("407fb868-8c7a-4e1d-a35e-fcb5172ced9a")) {
                    Plotly.newPlot(
                        '407fb868-8c7a-4e1d-a35e-fcb5172ced9a',
                        [{"domain": {"x": [0.0, 1.0], "y": [0.0, 1.0]}, "hoverlabel": {"namelength": 0}, "hovertemplate": "variant=%{label}<br>counts=%{value}", "labels": ["snp", "complex", "del", "ins", "mnp"], "legendgroup": "", "name": "", "showlegend": true, "type": "pie", "values": [5377, 112, 919, 531, 38]}],
                        {"legend": {"tracegroupgap": 0}, "template": {"data": {"bar": [{"error_x": {"color": "#2a3f5f"}, "error_y": {"color": "#2a3f5f"}, "marker": {"line": {"color": "#E5ECF6", "width": 0.5}}, "type": "bar"}], "barpolar": [{"marker": {"line": {"color": "#E5ECF6", "width": 0.5}}, "type": "barpolar"}], "carpet": [{"aaxis": {"endlinecolor": "#2a3f5f", "gridcolor": "white", "linecolor": "white", "minorgridcolor": "white", "startlinecolor": "#2a3f5f"}, "baxis": {"endlinecolor": "#2a3f5f", "gridcolor": "white", "linecolor": "white", "minorgridcolor": "white", "startlinecolor": "#2a3f5f"}, "type": "carpet"}], "choropleth": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "type": "choropleth"}], "contour": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "contour"}], "contourcarpet": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "type": "contourcarpet"}], "heatmap": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "heatmap"}], "heatmapgl": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "heatmapgl"}], "histogram": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "histogram"}], "histogram2d": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "histogram2d"}], "histogram2dcontour": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "histogram2dcontour"}], "mesh3d": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "type": "mesh3d"}], "parcoords": [{"line": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "parcoords"}], "pie": [{"automargin": true, "type": "pie"}], "scatter": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatter"}], "scatter3d": [{"line": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatter3d"}], "scattercarpet": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattercarpet"}], "scattergeo": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattergeo"}], "scattergl": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattergl"}], "scattermapbox": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattermapbox"}], "scatterpolar": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatterpolar"}], "scatterpolargl": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatterpolargl"}], "scatterternary": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatterternary"}], "surface": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "surface"}], "table": [{"cells": {"fill": {"color": "#EBF0F8"}, "line": {"color": "white"}}, "header": {"fill": {"color": "#C8D4E3"}, "line": {"color": "white"}}, "type": "table"}]}, "layout": {"annotationdefaults": {"arrowcolor": "#2a3f5f", "arrowhead": 0, "arrowwidth": 1}, "coloraxis": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "colorscale": {"diverging": [[0, "#8e0152"], [0.1, "#c51b7d"], [0.2, "#de77ae"], [0.3, "#f1b6da"], [0.4, "#fde0ef"], [0.5, "#f7f7f7"], [0.6, "#e6f5d0"], [0.7, "#b8e186"], [0.8, "#7fbc41"], [0.9, "#4d9221"], [1, "#276419"]], "sequential": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "sequentialminus": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]]}, "colorway": ["#636efa", "#EF553B", "#00cc96", "#ab63fa", "#FFA15A", "#19d3f3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"], "font": {"color": "#2a3f5f"}, "geo": {"bgcolor": "white", "lakecolor": "white", "landcolor": "#E5ECF6", "showlakes": true, "showland": true, "subunitcolor": "white"}, "hoverlabel": {"align": "left"}, "hovermode": "closest", "mapbox": {"style": "light"}, "paper_bgcolor": "white", "plot_bgcolor": "#E5ECF6", "polar": {"angularaxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}, "bgcolor": "#E5ECF6", "radialaxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}}, "scene": {"xaxis": {"backgroundcolor": "#E5ECF6", "gridcolor": "white", "gridwidth": 2, "linecolor": "white", "showbackground": true, "ticks": "", "zerolinecolor": "white"}, "yaxis": {"backgroundcolor": "#E5ECF6", "gridcolor": "white", "gridwidth": 2, "linecolor": "white", "showbackground": true, "ticks": "", "zerolinecolor": "white"}, "zaxis": {"backgroundcolor": "#E5ECF6", "gridcolor": "white", "gridwidth": 2, "linecolor": "white", "showbackground": true, "ticks": "", "zerolinecolor": "white"}}, "shapedefaults": {"line": {"color": "#2a3f5f"}}, "ternary": {"aaxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}, "baxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}, "bgcolor": "#E5ECF6", "caxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}}, "title": {"x": 0.05}, "xaxis": {"automargin": true, "gridcolor": "white", "linecolor": "white", "ticks": "", "title": {"standoff": 15}, "zerolinecolor": "white", "zerolinewidth": 2}, "yaxis": {"automargin": true, "gridcolor": "white", "linecolor": "white", "ticks": "", "title": {"standoff": 15}, "zerolinecolor": "white", "zerolinewidth": 2}}}, "title": {"text": "Summary of Variant compositions in the sample file /data/mstambou/tempus/Challenge_data_(1)parsedVCF.tsv"}},
                        {"responsive": true}
                    ).then(function(){
                            
var gd = document.getElementById('407fb868-8c7a-4e1d-a35e-fcb5172ced9a');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })
                };
                });
            </script>
        </div>

