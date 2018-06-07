

Format:

- [ ] free annotation `relationship` cell ontology class

Relationships could be either:

- `is_a` - means that the free annotation label is a new "cell type" to add to the cell ontology
- `synonym` - means that the free annotation label is not a new cell type, but another term for the same thing that already exists in the cell ontology

I wasn't always sure which relationship to use so please comment below and I will correct it.

Examples:

- [ ] Inhibitory neuron `is_a` [neuron](http://purl.obolibrary.org/obo/CL_0000540)
- [ ] alveolar epithelial type 2 cells `synonym` [type II pneumocyte](http://purl.obolibrary.org/obo/CL_0002063)


FACS and droplet free annotation tables in the "Details" below.



<details>

## FACS free annotation table

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>tissue</th>
      <th>cell_ontology_class</th>
      <th>free_annotation</th>
      <th>n_cells</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Aorta</td>
      <td>endothelial cell</td>
      <td>endothelial cells_1</td>
      <td>74</td>
    </tr>
    <tr>
      <td>Aorta</td>
      <td>endothelial cell</td>
      <td>endothelial cells_2</td>
      <td>70</td>
    </tr>
    <tr>
      <td>Aorta</td>
      <td>endothelial cell</td>
      <td>endothelial cells_3</td>
      <td>44</td>
    </tr>
    <tr>
      <td>Bladder</td>
      <td>bladder cell</td>
      <td>Bladder mesenchymal cell</td>
      <td>695</td>
    </tr>
    <tr>
      <td>Bladder</td>
      <td>bladder urothelial cell</td>
      <td>Basal bladder epithelial cell</td>
      <td>101</td>
    </tr>
    <tr>
      <td>Bladder</td>
      <td>bladder urothelial cell</td>
      <td>Luminal bladder epithelial cell</td>
      <td>582</td>
    </tr>
    <tr>
      <td>Brain_Non-Myeloid</td>
      <td>neuron</td>
      <td>excitatory neurons and some neuronal stem cells</td>
      <td>194</td>
    </tr>
    <tr>
      <td>Brain_Non-Myeloid</td>
      <td>neuron</td>
      <td>inhibitory neurons</td>
      <td>87</td>
    </tr>
    <tr>
      <td>Fat</td>
      <td>mesenchymal stem cell of adipose</td>
      <td>mesenchymal progenitor</td>
      <td>2107</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>Brush cell of epithelium proper of large intes...</td>
      <td>Tuft cell</td>
      <td>63</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>enterocyte of epithelium of large intestine</td>
      <td>Enterocyte (Distal)</td>
      <td>191</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>enterocyte of epithelium of large intestine</td>
      <td>Enterocyte (Proximal)</td>
      <td>773</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>enteroendocrine cell</td>
      <td>Chromaffin Cell</td>
      <td>59</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>epithelial cell of large intestine</td>
      <td>Lgr5+ amplifying undifferentiated cell (Distal)</td>
      <td>106</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>epithelial cell of large intestine</td>
      <td>Lgr5+ amplifying undifferentiated cell (Proximal)</td>
      <td>172</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>epithelial cell of large intestine</td>
      <td>Lgr5+ undifferentiated cell (Distal)</td>
      <td>343</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>epithelial cell of large intestine</td>
      <td>Lgr5+ undifferentiated cell (Proximal)</td>
      <td>528</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>epithelial cell of large intestine</td>
      <td>Lgr5- amplifying undifferentiated cell</td>
      <td>576</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>epithelial cell of large intestine</td>
      <td>Lgr5- undifferentiated cell</td>
      <td>294</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>large intestine goblet cell</td>
      <td>Goblet cell (Distal)</td>
      <td>471</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>large intestine goblet cell</td>
      <td>Goblet cell (Proximal)</td>
      <td>264</td>
    </tr>
    <tr>
      <td>Large_Intestine</td>
      <td>large intestine goblet cell</td>
      <td>Goblet cell, top of crypt (Distal)</td>
      <td>98</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>endothelial cell of hepatic sinusoid</td>
      <td>endothelial cell</td>
      <td>182</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>natural killer cell</td>
      <td>NK/NKT cells</td>
      <td>39</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>ciliated columnar cell of tracheobronchial tree</td>
      <td>multiciliated cells</td>
      <td>25</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>classical monocyte</td>
      <td>invading monocytes</td>
      <td>90</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>epithelial cell of lung</td>
      <td>alveolar epithelial type 1 cells, alveolar epi...</td>
      <td>113</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>leukocyte</td>
      <td>mast cells and unknown immune cells</td>
      <td>35</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>monocyte</td>
      <td>circulating monocytes</td>
      <td>65</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>myeloid cell</td>
      <td>dendritic cells, alveolar macrophages, and int...</td>
      <td>85</td>
    </tr>
    <tr>
      <td>Mammary_Gland</td>
      <td>luminal epithelial cell of mammary gland</td>
      <td>luminal progenitor</td>
      <td>411</td>
    </tr>
    <tr>
      <td>Mammary_Gland</td>
      <td>luminal epithelial cell of mammary gland</td>
      <td>mature luminal cell</td>
      <td>167</td>
    </tr>
    <tr>
      <td>Marrow</td>
      <td>B cell</td>
      <td>Cd3e+ Klrb1+ B cell</td>
      <td>44</td>
    </tr>
    <tr>
      <td>Marrow</td>
      <td>late pro-B cell</td>
      <td>Dntt+ late pro-B cell</td>
      <td>117</td>
    </tr>
    <tr>
      <td>Marrow</td>
      <td>late pro-B cell</td>
      <td>Dntt- late pro-B cell</td>
      <td>189</td>
    </tr>
    <tr>
      <td>Marrow</td>
      <td>precursor B cell</td>
      <td>pre-B cell (Philadelphia nomenclature)</td>
      <td>517</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>endothelial cell</td>
      <td>endothelial cell</td>
      <td>66</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>leukocyte</td>
      <td>immune cell</td>
      <td>54</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>pancreatic A cell</td>
      <td>pancreatic A cell</td>
      <td>390</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>pancreatic D cell</td>
      <td>pancreatic D cell</td>
      <td>140</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>pancreatic PP cell</td>
      <td>pancreatic PP cell</td>
      <td>73</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>pancreatic acinar cell</td>
      <td>acinar cell</td>
      <td>182</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>pancreatic ductal cell</td>
      <td>ductal cell</td>
      <td>161</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>pancreatic stellate cell</td>
      <td>stellate cell</td>
      <td>49</td>
    </tr>
    <tr>
      <td>Pancreas</td>
      <td>type B pancreatic cell</td>
      <td>beta cell</td>
      <td>449</td>
    </tr>
    <tr>
      <td>Skin</td>
      <td>basal cell of epidermis</td>
      <td>Basal IFE</td>
      <td>562</td>
    </tr>
    <tr>
      <td>Skin</td>
      <td>epidermal cell</td>
      <td>Intermediate IFE</td>
      <td>276</td>
    </tr>
    <tr>
      <td>Skin</td>
      <td>keratinocyte stem cell</td>
      <td>Inner Bulge</td>
      <td>573</td>
    </tr>
    <tr>
      <td>Skin</td>
      <td>keratinocyte stem cell</td>
      <td>Outer Bulge</td>
      <td>831</td>
    </tr>
    <tr>
      <td>Skin</td>
      <td>leukocyte</td>
      <td>Leukocyte</td>
      <td>15</td>
    </tr>
    <tr>
      <td>Skin</td>
      <td>stem cell of epidermis</td>
      <td>Replicating Basal IFE</td>
      <td>53</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>DN1 thymic pro-T cell</td>
      <td>DN1 thymocytes</td>
      <td>32</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69 negative rapidly div...</td>
      <td>135</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69 negative thymocytes</td>
      <td>540</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69 negative thyomcytes</td>
      <td>39</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69 positive thymocytes</td>
      <td>563</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>professional antigen presenting cell</td>
      <td>antigen presenting cell</td>
      <td>40</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>basal cell of epidermis</td>
      <td>basal cells</td>
      <td>729</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>basal cell of epidermis</td>
      <td>differentiating basal cells</td>
      <td>161</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>basal cell of epidermis</td>
      <td>proliferating</td>
      <td>196</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>keratinocyte</td>
      <td>differentiated keratinocyte</td>
      <td>182</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>keratinocyte</td>
      <td>filiform differentiated keratinocytes</td>
      <td>27</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>keratinocyte</td>
      <td>suprabasal differentiating keratinocytes</td>
      <td>121</td>
    </tr>
  </tbody>
</table>


## Droplet free annotations

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>tissue</th>
      <th>cell_ontology_class</th>
      <th>free_annotation</th>
      <th>n_cells</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Bladder</td>
      <td>bladder cell</td>
      <td>Bladder mesenchymal cell</td>
      <td>1203</td>
    </tr>
    <tr>
      <td>Bladder</td>
      <td>bladder urothelial cell</td>
      <td>Basal bladder epithelial cell</td>
      <td>498</td>
    </tr>
    <tr>
      <td>Bladder</td>
      <td>bladder urothelial cell</td>
      <td>Luminal bladder epithelial cell</td>
      <td>669</td>
    </tr>
    <tr>
      <td>Bladder</td>
      <td>leukocyte</td>
      <td>Monocyte/Macrophage</td>
      <td>73</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>duct epithelial cell</td>
      <td>bile duct epithelial cells</td>
      <td>27</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>endothelial cell of hepatic sinusoid</td>
      <td>endothelial cells</td>
      <td>28</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>hepatocyte</td>
      <td>midlobular female</td>
      <td>234</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>hepatocyte</td>
      <td>midlobular male</td>
      <td>526</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>hepatocyte</td>
      <td>pericentral female</td>
      <td>274</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>hepatocyte</td>
      <td>pericentral male</td>
      <td>159</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>hepatocyte</td>
      <td>periportal female</td>
      <td>308</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>hepatocyte</td>
      <td>periportal male</td>
      <td>263</td>
    </tr>
    <tr>
      <td>Liver</td>
      <td>leukocyte</td>
      <td>immune cells</td>
      <td>26</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>ciliated columnar cell of tracheobronchial tree</td>
      <td>multiciliated cells</td>
      <td>49</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>classical monocyte</td>
      <td>invading monocytes</td>
      <td>161</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>myeloid cell</td>
      <td>dendritic cells and interstital macrophages</td>
      <td>87</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>non-classical monocyte</td>
      <td>circulating monocytes</td>
      <td>220</td>
    </tr>
    <tr>
      <td>Lung</td>
      <td>type II pneumocyte</td>
      <td>alveolar epithelial type 2 cells</td>
      <td>89</td>
    </tr>
    <tr>
      <td>Mammary_Gland</td>
      <td>luminal epithelial cell of mammary gland</td>
      <td>luminal progenitor cell</td>
      <td>243</td>
    </tr>
    <tr>
      <td>Mammary_Gland</td>
      <td>luminal epithelial cell of mammary gland</td>
      <td>mature luminal cell</td>
      <td>216</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>DN1 thymic pro-T cell</td>
      <td>DN1 thymocytes</td>
      <td>44</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN3-DN4 thymocytes</td>
      <td>33</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69 negative rapidly div...</td>
      <td>219</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69 negative thymocytes</td>
      <td>71</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69 positive thymocytes</td>
      <td>187</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP in transition Cd69_negative thymocytes</td>
      <td>334</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP transition Cd69 low thymocytes</td>
      <td>227</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>immature T cell</td>
      <td>DN4-DP transition Cd69 negative rapidly dividi...</td>
      <td>283</td>
    </tr>
    <tr>
      <td>Thymus</td>
      <td>professional antigen presenting cell</td>
      <td>antigen presenting cells</td>
      <td>31</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>basal cell of epidermis</td>
      <td>basal cells</td>
      <td>3345</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>basal cell of epidermis</td>
      <td>proliferating</td>
      <td>1079</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>keratinocyte</td>
      <td>differentiated</td>
      <td>1095</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>keratinocyte</td>
      <td>filiform differentiated</td>
      <td>199</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>keratinocyte</td>
      <td>suprabasal differentiated</td>
      <td>967</td>
    </tr>
    <tr>
      <td>Tongue</td>
      <td>keratinocyte</td>
      <td>suprabasal differentiating</td>
      <td>815</td>
    </tr>
  </tbody>
</table>
</details>

## Aorta

### FACS

(None)

## Bladder

@akershner @bkiss

### FACS *and* Droplet

- [ ] Bladder mesenchymal cell `is_a` [bladder cell](http://purl.obolibrary.org/obo/CL_1001319)
- [ ] Basal bladder epithelial cell `is_a` [bladder urothelial cell](http://purl.obolibrary.org/obo/CL_1001428)
- [ ] Luminal bladder epithelial cell`is_a` [bladder urothelial cell](http://purl.obolibrary.org/obo/CL_1001428)



## Brain Myeloid

### FACS

(None)

## Brain Non-Myeloid

@Taliram

### FACS

- [ ] Inhibitory neuron `is_a` [neuron](http://purl.obolibrary.org/obo/CL_0000540)
- [ ] Excitatory neuron `is_a` [neuron](http://purl.obolibrary.org/obo/CL_0000540)

## Fat

@nschaum

### FACS

- [ ] mesenchymal progenitor `is_a` [mesenchymal stem cell of adipose](http://purl.obolibrary.org/obo/CL_0002570)

## Heart

@guangli0817

### FACS

- [ ] Conduction cell `is_a` [heart cell](http://www.ontobee.org/ontology/CL?iri=http://purl.obolibrary.org/obo/CL_0002494)

### Droplet

(None)

## Kidney

### FACS

(None)

### Droplet

(None)

## Large Intestine

@tisobe1

### FACS

- [ ] tuft cell `is_a` [Brush cell of epithelium proper of large intestine](http://www.ontobee.org/ontology/CL?iri=http://purl.obolibrary.org/obo/CL_0002203)
- [ ] enterocyte (distal) `is_a` [enterocyte of epithelium of large intestine](http://purl.obolibrary.org/obo/CL_0002071)
- [ ] enterocyte (proximal) `is_a` [enterocyte of epithelium of large intestine](http://purl.obolibrary.org/obo/CL_0002071)
- [ ] chromaffin cell `is_a` [enteroendocrine cell](http://purl.obolibrary.org/obo/CL_0000164)
- [ ] Lgr5+ amplifying undifferentiated cell (Distal) `is_a` [epithelial cell of large intestine](http://purl.obolibrary.org/obo/CL_0002253)
- [ ] Lgr5+ amplifying undifferentiated cell (Proximal) `is_a` [epithelial cell of large intestine](http://purl.obolibrary.org/obo/CL_0002253)
- [ ] Lgr5+ undifferentiated cell (Distal) `is_a` [epithelial cell of large intestine](http://purl.obolibrary.org/obo/CL_0002253)
- [ ] Lgr5+ undifferentiated cell (Proximal) `is_a` [epithelial cell of large intestine](http://purl.obolibrary.org/obo/CL_0002253)
- [ ] Goblet cell (Distal) `is_a` [large intestine goblet cell](http://purl.obolibrary.org/obo/CL_1000320)
- [ ] Goblet cell (Proximal) `is_a` [large intestine goblet cell](http://purl.obolibrary.org/obo/CL_1000320)
- [ ] Goblet cell, top of crypt (Distal) `is_a`[large intestine goblet cell](http://purl.obolibrary.org/obo/CL_1000320)

## Limb Muscle

### FACS

(None)

### Droplet

(None)

## Liver

@patsika @batson

### FACS

(None)

### Droplet

- [ ] bile duct epithelial cells `is_a` [duct epithelial cell](http://purl.obolibrary.org/obo/CL_0000068)
- [ ] Periportal hepatocyte `is_a` [hepatocyte](http://purl.obolibrary.org/obo/CL_0000182)
- [ ] Midlobular hepatocyte `is_a` [hepatocyte](http://purl.obolibrary.org/obo/CL_0000182)
- [ ] Pericentral hepatocyte `is_a` [hepatocyte](http://purl.obolibrary.org/obo/CL_0000182)

## Lung

@ktravaglini

### FACS *and* Droplet

- [ ] multiciliated cells `is_a` [ciliated columnar cell of tracheobronchial tree](http://purl.obolibrary.org/obo/CL_0002145)
- [ ] invading monocytes `is_a` [classical monocyte](http://purl.obolibrary.org/obo/CL_0000860)
- [ ] circulating monocytes `is_a` [monocyte](http://purl.obolibrary.org/obo/CL_0000576)

### Droplet only

- [ ] circulating monocytes `is_a` [non-classical monocyte](http://purl.obolibrary.org/obo/CL_0000576)
- [ ] alveolar epithelial type 2 cells `synonym` [type II pneumocyte](http://purl.obolibrary.org/obo/CL_0002063)

## Mammary Gland

@sikandars

### FACS *and* Droplet

- [ ] Luminal progenitor `is_a` [luminal epithelial cell of mammary gland](http://purl.obolibrary.org/obo/CL_0002326)
- [ ] Mature luminal cell `is_a` [luminal epithelial cell of mammary gland](http://purl.obolibrary.org/obo/CL_0002326)

## Marrow

@transcriptomics

### FACS

- [ ] Cd3e+ Klrb1+ B cell `is_a` [B cell](http://purl.obolibrary.org/obo/CL_0000236)
- [ ] Dntt+ late pro-B cell `is_a` [late pro-B cell](http://purl.obolibrary.org/obo/CL_0002048)
- [ ] Dntt- late pro-B cell `is_a` [late pro-B cell](http://purl.obolibrary.org/obo/CL_0002048)

### Droplet

(None)


## Pancreas

@YtheRookie

### FACS

(None)

"alpha cell" and such as synonmyms for e.g. [pancreatic A cell](http://purl.obolibrary.org/obo/CL_0000171) but weren't in the list of terms that we used.

## Skin

@nschaum

### FACS

- [ ] basal skin interfollicular epidermis cell `is_a` [basal cell of epidermis](http://purl.obolibrary.org/obo/CL_0002187)
- [ ] intermediate skin interfollicular epidermis cell `is_a` [epidermal cell](http://purl.obolibrary.org/obo/CL_0000362)
- [ ] inner bulge cell `is_a` [keratinocyte stem cell](http://purl.obolibrary.org/obo/CL_0002337)
- [ ] replicating basal interfollicular epidermis cell `is_a` [stem cell of epidermis](http://purl.obolibrary.org/obo/CL_1000428)


## Spleen

### FACS

(None)

### Droplet

(None)

## Thymus

@pknguyen1

### FACS *and* Droplet

- [ ] double negative 1 (DN1) thymocyte `synonym` [DN1 thymic pro-T cell](http://purl.obolibrary.org/obo/CL_0000894)
- [ ] DN4-DP in transition Cd69 negative rapidly dividing thymocyte `is_a` [immature T cell](http://purl.obolibrary.org/obo/CL_0002420)
- [ ] DN4-DP in transition Cd69 positive rapidly dividing thymocyte `is_a` [immature T cell](http://purl.obolibrary.org/obo/CL_0002420)
- [ ] DN4-DP in transition Cd69 negative `is_a` [immature T cell](http://purl.obolibrary.org/obo/CL_0002420)
- [ ] DN4-DP in transition Cd69 positive `is_a` [immature T cell](http://purl.obolibrary.org/obo/CL_0002420)


### Droplet only

- [ ] DN3-DN4 thymocyte `is_a` [immature T cell](http://purl.obolibrary.org/obo/CL_0002420)
- [ ] DN4-DP transition Cd69 low thymocyte `is_a` [immature T cell](http://purl.obolibrary.org/obo/CL_0002420)
- [ ] DN4-DP transition Cd69 negative rapidly dividing thymocyte `is_a` [immature T cell](http://purl.obolibrary.org/obo/CL_0002420)

## Tongue

@melocactus

### FACS

- [ ] differentiating basal cells `is_a` [basal cell of epidermis](http://purl.obolibrary.org/obo/CL_0002187)
- [ ] differentiated keratinocyte `is_a` [keratinocyte](http://purl.obolibrary.org/obo/CL_0000312)
- [ ] filiform differentiated keratinocyte `is_a` [keratinocyte](http://purl.obolibrary.org/obo/CL_0000312)
- [ ] suprabasal differentiating keratinocyte `is_a` [keratinocyte](http://purl.obolibrary.org/obo/CL_0000312)

### Droplet

(None)

## Trachea

### FACS

(None)

### Droplet

(None)