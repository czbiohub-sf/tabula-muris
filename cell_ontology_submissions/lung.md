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
      <td>Lung</td>
      <td>NA</td>
      <td>lung neuroendocrine cells and unknown cells</td>
      <td>40</td>
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
  </tbody>
</table>


</details>

## Proposed annotations

@ktravaglini

### FACS *and* Droplet

- [ ] multiciliated cells `is_a` [ciliated columnar cell of tracheobronchial tree](http://purl.obolibrary.org/obo/CL_0002145)
- [ ] invading monocytes `is_a` [classical monocyte](http://purl.obolibrary.org/obo/CL_0000860)
- [ ] circulating monocytes `is_a` [monocyte](http://purl.obolibrary.org/obo/CL_0000576)

### Droplet only

- [ ] circulating monocytes `is_a` [non-classical monocyte](http://purl.obolibrary.org/obo/CL_0000576)
- [ ] alveolar epithelial type 2 cells `synonym` [type II pneumocyte](http://purl.obolibrary.org/obo/CL_0002063)
