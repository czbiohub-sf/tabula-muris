# coding: utf-8
import re

methods = 'facs', 'droplet'

for method in methods:

    rmds = get_ipython().getoutput(f'ls *{method}.Rmd')

    for rmd in rmds:
        backup = rmd + '.backup'
        get_ipython().system(' cp $rmd $backup')
        replaced = rmd + '.replaced'
        print(replaced)
        with open(rmd) as f:
            content = f.read()
            content = re.sub('''Enter the directory of the maca folder on your drive and the name of the tissue you want to analyze.

(```{r}
tissue_of_interest = \"[\w_-]+\")(.*)(PCHeatmap\(object = tiss)''', 
r'''Specify the tissue of interest, run the boilerplate code which sets up the functions and environment, load the tissue object.

\1

library(here)
source(here("00_data_ingest", "02_tissue_analysis_rmd", "boilerplate.R"))
''' 
+ f'load_tissue_{method}(tissue_of_interest)' 
+ r'''
```

Visualize top genes in principal components

```{r, echo=FALSE, fig.height=4, fig.width=8}
\3''', content, flags=re.DOTALL)
            # print(content)
    #         content = content.replace("""tissue_of_interest = "Bladder"

    # tiss <- ScaleData(object = tiss, vars.to.regress = c("nUMI", "percent.ribo","Rn45s"))""", """tiss <- NormalizeData(object = tiss)
    # tiss <- ScaleData(object = tiss)""")
            with open(replaced, 'w') as g:
                g.write(content)
           
    get_ipython().system(f'diff Bladder_{method}.Rmd Bladder_{method}.Rmd.replaced')

            
replaceds = get_ipython().getoutput('ls *.Rmd.replaced')
for replaced in replaceds:
    rmd = replaced.split('.replaced')[0]
    print(rmd)
    get_ipython().system(' mv $replaced $rmd')
    
