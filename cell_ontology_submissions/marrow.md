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

@transcriptomics

### FACS

- [ ] Cd3e+ Klrb1+ B cell `is_a` [B cell](http://purl.obolibrary.org/obo/CL_0000236)
- [ ] Dntt+ late pro-B cell `is_a` [late pro-B cell](http://purl.obolibrary.org/obo/CL_0002048)
- [ ] Dntt- late pro-B cell `is_a` [late pro-B cell](http://purl.obolibrary.org/obo/CL_0002048)
