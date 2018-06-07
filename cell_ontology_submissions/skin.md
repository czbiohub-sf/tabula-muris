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

@nschaum

### FACS

- [ ] basal skin interfollicular epidermis cell `is_a` [basal cell of epidermis](http://purl.obolibrary.org/obo/CL_0002187)
- [ ] intermediate skin interfollicular epidermis cell `is_a` [epidermal cell](http://purl.obolibrary.org/obo/CL_0000362)
- [ ] inner bulge cell `is_a` [keratinocyte stem cell](http://purl.obolibrary.org/obo/CL_0002337)
- [ ] replicating basal interfollicular epidermis cell `is_a` [stem cell of epidermis](http://purl.obolibrary.org/obo/CL_1000428)
