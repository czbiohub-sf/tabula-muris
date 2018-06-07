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
  </tbody>
</table>




</details>

## Proposed annotations

@patsika @batson

### FACS

(None)

### Droplet

- [ ] bile duct epithelial cells `is_a` [duct epithelial cell](http://purl.obolibrary.org/obo/CL_0000068)
- [ ] Periportal hepatocyte `is_a` [hepatocyte](http://purl.obolibrary.org/obo/CL_0000182)
- [ ] Midlobular hepatocyte `is_a` [hepatocyte](http://purl.obolibrary.org/obo/CL_0000182)
- [ ] Pericentral hepatocyte `is_a` [hepatocyte](http://purl.obolibrary.org/obo/CL_0000182)