Hello!
This discussion is to propose the `free_annotation` terms to the Cell Ontology so that future researchers may use your specific terms rather than the generic ones that you had to use. Below are  your `free_annotation` terms and what I (@olgabot) thought the relationship was to the existing cell ontology term, whether it was brand new term and cell type, in which case I used `is_a`, or it was merely a synonym for a term that already existed in the cell ontology, in which case I used `synonym`.

Please comment below with what is the correct relationship between terms, or correct the term itself and I'll deal with the submission to the Cell Ontology group.

Thank you for your help!
Warmest,
Olga



For the annotations below, here is the format:

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


### FACS free annotation table

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
  </tbody>
</table>


### Droplet free annotation table

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
  </tbody>
</table>



</details>

## Proposed annotations

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
CL?iri=http://purl.obolibrary.org/obo/CL_0002494)