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

## Proposed annotations

@melocactus

### FACS

- [ ] differentiating basal cells `is_a` [basal cell of epidermis](http://purl.obolibrary.org/obo/CL_0002187)
- [ ] differentiated keratinocyte `is_a` [keratinocyte](http://purl.obolibrary.org/obo/CL_0000312)
- [ ] filiform differentiated keratinocyte `is_a` [keratinocyte](http://purl.obolibrary.org/obo/CL_0000312)
- [ ] suprabasal differentiating keratinocyte `is_a` [keratinocyte](http://purl.obolibrary.org/obo/CL_0000312)
