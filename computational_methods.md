## Clustering

Standard procedures for cell selection, variable gene selection, dimensionality reduction, and clustering were performed using the Seurat package. For each tissue and each sequencing method (FACS and microfluidic droplet), we separately performed the following steps.

0. Exclude cells with fewer than 500 genes detected (a gene is detected if it has a count of at least 1). Exclude cells with fewer than 50,000 reads (FACS) or 1000 UMIs (microfluidic droplet).

1. Load a raw data matrix $A$ with rows genes, columns cells, and entries counts (reads or UMI).

2. Create a log normalized data matrix
$$N_{ij} = \log \left (1 + M \frac{A_{ij}}{\sum_{j^\prime} A_{i j^\prime}} \right ),$$
where we set $M = 10^6$ for FACS and $M = 10^4$ for droplets. The log is base $e$.

3. Create a scaled data matrix, rescaling each gene to mean zero and variance one.
$$X_{ij} = (N_{ij} - \mu_i)/\sigma_i,$$
where $\mu_i$ is the mean of $N_{ij}$ and $\sigma_i$ is the standard deviation of $N_{ij}$.

4. Variable genes were selected using a threshold for standardized log dispersion (FindVariableGenes in Seurat). More precisely, the log dispersion $d_i$ of a gene $i$ is
$$d_i = \log(v_i/m_i),$$
where $v_i$ is the variance of the raw counts $A_{ij}$ and $m_i$ is the mean of the raw counts $A_{ij}$. We bin the genes into 20 equal-spaced bins according to their log means $\log(m_i)$, then compute the mean and standard deviation of $d_i$ within each bin. The standardized log dispersion $\bar{d}_i$ is the dispersion $d_i$ shifted by the mean and rescaled by the standard deviation within its bin. We retain genes with $\bar{d}_i > 0.5$.

5.	The variable genes were projected onto a low-dimensional subspace using principal component analysis on the appropriate subset of $X_{ij}$. The number of principal components was selected based on inspection of the (elbow) plot of variance.

6.	A shared-nearest-neighbors graph was constructed based on the Euclidean distance in the low-dimensional subspace. Let $\mathcal{N}(j)$ denote the k-nearest neighborhood (k = 30) for a cell $j$. The shared-nearest-neighbor graph G has a vertex for each cell and an edge of weight
$$w_{jk} = \frac{
\vert \mathcal{N}(j) \cap \mathcal{N}(k) \vert}{\vert \mathcal{N}(j) \cup \mathcal{N}(k) \vert}$$
between cells $j$ and $k$.
Cells were clustered using a modified version of the Louvain method for modularity maximization. [cite https://journals.aps.org/pre/abstract/10.1103/PhysRevE.74.016110]. The modularity has a resolution parameter $\gamma$,

$$Q = \sum_{jk} \left (A_{jk} - \gamma \frac{k_i k_j}{2m} \right ) \delta(c_i, c_j),$$

where $A_{jk}$ is the weighted adjacency matrix, $k_i$ and $k_j$ are the weighted degrees, $m$ is the total weight of edges, $c_i$ denotes cluster membership, and $\delta(c_i,c_j)$ is $1$ if $j$ and $k$ are in the same cluster, and $0$ otherwise.
The resolution is a tuneable parameter in this analysis: larger $\gamma$ favors smaller clusters.

7.	Cells were visualized using a 2-dimensional t-distributed Stochastic Neighbor Embedding on the PCA embedding.

8.	Cell types were assigned to each cluster using the expression of known marker genes. Plots showing the expression of the markers for each tissue appear in the extended data.

9.	When clusters appeared to be mixtures of cell types, they were refined either by increasing the resolution parameter for Louvain clustering, or subsetting the data and rerunning steps 3-7.

A similar analysis was done globally for all FACS processed cells and for all microfluidic droplet processed cells (Supp. Fig. 4).

We elected not to use a regression step, as is commonly done against the number of counts per cell or the percent of mitochondrial genes. Because those metrics are often correlated with biology, including cell type, regression risks suppressing important signals (Sup Fig X).
