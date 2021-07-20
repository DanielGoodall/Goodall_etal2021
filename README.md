# Goodall_etal2021
- Scripts and data from Goodall et al 2021
- **A fork trap in the chromosomal termination area is highly conserved across all Escherichia coli phylogenetic groups**

---


# Install R
https://www.r-project.org/

# Install Rstudio
https://www.rstudio.com/

# Install Bioconductor
https://www.bioconductor.org/

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.13")
```

## Packages
- Biostrings (manipulation of genetic strings)
- DECIPHER (analysing and manipulating sequences)
- msa (alignment algorithm wrapper package)
- Rbowtie2
- Rsamtools

# GenoPlotR
https://genoplotr.r-forge.r-project.org/
genoPlotR is a R package to produce reproducible, publication-grade graphics of gene and genome maps. It allows the user to read from usual format such as protein table files and blast results, as well as home-made tabular files.

```{r}
install.packages("genoPlotR")
```

# BLAST
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

stand-alone set up for Windows PC
https://www.ncbi.nlm.nih.gov/books/NBK52637/







