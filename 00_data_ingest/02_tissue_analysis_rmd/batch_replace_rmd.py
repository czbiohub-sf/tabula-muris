# coding: utf-8

# How to use:
# type 'ipython' at the terminal to get to the Interactive Python (IPython) environments
# Type: 
# %run batch_replace_rmd.py

import re

methods = 'facs', 'droplet'

for method in methods:

    rmds = get_ipython().getoutput(f'ls *{method}.Rmd')

    for rmd in rmds:
        backup = rmd + '.backup'
        get_ipython().system(' cp $rmd $backup')
        replaced = rmd + '.replaced'
        print(replaced)
        batch_name = 'channel' if method == 'droplet' else 'plate.barcode'
        with open(rmd) as f:
            content = f.read()
            content += \
f'''
Write the cell ontology and free annotations to CSV.

```{{r}}
filename = here('00_data_ingest', '03_tissue_annotation_csv', 
                    paste0(tissue_of_interest, "_{method}_annotation.csv"))
write.csv(tiss@meta.data[,c('{batch_name},'cell_ontology_class','cell_ontology_id', 'free.annotation')], file=filename)
```
'''
            # print(content)
    #         content = content.replace("""tissue_of_interest = "Bladder"

    # tiss <- ScaleData(object = tiss, vars.to.regress = c("nUMI", "percent.ribo","Rn45s"))""", """tiss <- NormalizeData(object = tiss)
    # tiss <- ScaleData(object = tiss)""")
            with open(replaced, 'w') as g:
                g.write(content)
           
    # get_ipython().system(f'diff Bladder_{method}.Rmd Bladder_{method}.Rmd.replaced')
    get_ipython().system(f'diff Lung_{method}.Rmd Lung_{method}.Rmd.replaced')

            
# replaceds = get_ipython().getoutput('ls *.Rmd.replaced')
# for replaced in replaceds:
#     rmd = replaced.split('.replaced')[0]
#     print(rmd)
#     get_ipython().system(' mv $replaced $rmd')
    
