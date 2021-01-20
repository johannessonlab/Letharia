# Other figures

The scipts here are simple R code where I explore the data and I produce figures S5, S6 and S7. I also made other plots that didn't make it to the paper, specially for the McKenzie genome of *L. columbiana*.

I provide the necessary data in the folder `data`. I produced it with the output of RepeatMasker using repeat library from [McKenzie et al. 2020](https://www.sciencedirect.com/science/article/abs/pii/S0888754320304614). I processed the gff file produced by RepeatMasker and the genome assemblies with my script `totalcovergff.py` that can be ran like so:

    $ python totalcovergff.py RepeatMasker/Lecol_assembly_V1.0.fa.out.gff -f McKenzieData/Lecol_assembly_V1.0.fasta > RepeatContent_Lecol_assembly_V1.0.txt

- `RepeatsLelup.R` produces figures S5, S6 and S7
