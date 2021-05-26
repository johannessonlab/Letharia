# Scripts and pipelines used in the *Letharia* study of the Johannesson's Lab

Here you'll find the scripts and pipelines associated to the study [Ament-Velásquez, Tuovinen, et al. (2020). The plot thickens: Haploid and triploid-like thalli, hybridization, and biased mating type ratios in Letharia. Front. Fungal Biol. 2:656386](https://www.frontiersin.org/articles/10.3389/ffunb.2021.656386/full).

There are four main folders:

- **MAT annotation:** gene annotation of the MAT locus of the genus *Letharia*, designed to run in UPPMAX (Uppsala Multidisciplinary Center for Advanced Computational Science).
- **Lichen Ploidy:** a Snakemake pipeline to infer the ploidy (*sensu lato*, since it would include heterokaryosis) from a (meta)genome, as well as additional analysis to show that the triploid-like samples look like hybrids.
- **Lichen Twisst:** a snakemake pipeline where I make little trees of SNP windows along the main contigs of the *L. lupina* pure culture. The format is inspired in Twisst ([Martin & Van Belleghem 2017 Genetics 206:429–438](https://www.genetics.org/content/206/1/429)) for plotting but it's not at all as complicated.
- **OtherFigures** script and data necessary for other figures, mostly related to the repeat content.