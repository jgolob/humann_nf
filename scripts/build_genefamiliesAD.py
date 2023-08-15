#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import anndata as ad
import scipy

BASEDIR = '.'
raw_data = {
    'genefamilies.tsv.gz': [],
}

for f_i, fn in enumerate(sorted(os.listdir(
        os.path.join(
            BASEDIR,
            'per_specimen'
        )
    ))):
    specimen = "_".join(fn.split("_")[:-1])
    filetype = fn.split('_')[-1]
    if filetype != 'genefamilies.tsv.gz':
        continue

    print(
        f_i+1,
        fn,
        specimen,
        filetype
    )
    raw_data[filetype] += [
        {
            'specimen': specimen,
            'feature': r.feature,
            'abundance': float(r.abundance)
        } for(i, r) in pd.read_csv(
            os.path.join(BASEDIR, 'per_specimen', fn),
            comment='#',
            names=['feature', 'abundance'],
            sep='\t'
        ).iterrows()
        if float(r.abundance) > 0
    ]

# Make sets of specimens and features
gf_sp_set = {
    gf['specimen'] for gf in raw_data['genefamilies.tsv.gz']
}
gf_f_set = {
    gf['feature'] for gf in raw_data['genefamilies.tsv.gz']
}
# Create index 
gf_obs_names = sorted(gf_sp_set)
gf_var_names = sorted(gf_f_set)
# Create index lookup
gf_obs_ilu = {
    o: i
    for (i, o) in enumerate(
        gf_obs_names
    )
}
gf_var_ilu = {
    o: i
    for (i, o) in enumerate(
        gf_var_names
    )
}
# Make sparse matrix:
# Anndata
gf_ad = ad.AnnData(
    scipy.sparse.coo_matrix(
        (
            [v['abundance'] for v in raw_data['genefamilies.tsv.gz']], # Values
            (
                [gf_obs_ilu.get(v['specimen']) for v in raw_data['genefamilies.tsv.gz']],
                [gf_var_ilu.get(v['feature']) for v in raw_data['genefamilies.tsv.gz']],
            )

        ),
        shape=(
            len(gf_obs_names), len(gf_var_names),
        ),
        dtype=np.float32
    ).tocsr(),
    obs=pd.DataFrame(
        index=gf_obs_names
    ),
    var=pd.DataFrame(
        index=gf_var_names
    ),
    filename=os.path.join(
            BASEDIR,
            'genefamilies.h5ad'
        ),
    filemode='w',
)

# var
gf_ad.var['Uniref'] = gf_ad.var.index.map(lambda f: f.split('|')[0])
gf_ad.var['Taxon'] = gf_ad.var.index.map(lambda f: "".join(f.split('|')[1:]))
# Output

gf_ad.write()