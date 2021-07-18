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

```{r}
install.packages("genoPlotR")
```





