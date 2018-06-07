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
  </tbody>
</table>

</details>

## Proposed annotations

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
