# Letharia

All scripts related to genome annotation of the genus *Letharia*, designed to be run in [UPPMAX](https://www.uppmax.uu.se/) (Uppsala Multidisciplinary Center for Advanced Computational Science) as part of the lichen projects in [Johannesson Lab](http://www.iob.uu.se/research/systematic-biology/johannesson/) (Uppsala University).

In the directory `Training_Files` you can find the training files obtained for the *ab initio* gene prediction programs SNAP ([Korf, 2004](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-59)), Augustus ([Stanke and Waack, 2003](http://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/19/suppl_2/10.1093/bioinformatics/btg1080/2/btg1080.pdf?Expires=1486485018&Signature=cP369aPiM8h0yRrtGGMbqNUVCuffyrWzGjiS5CItfQGN27Gp5i1MLYP69u4tRrDgbRV-e13zj769V9uVO6vKaJX8946e1a6U0QhQ5RIK4rYqXRKTDnC92h3wubW2LgCHMY4xjw4oFvOfAhbsEiwyoYMhSFEpfuS5m7PHLW9sgnMIIB6JOELjZ6lSetEi8k8rQdIned~yI4Fb39LV5FQViT8uneLGL4aug3f3w6M9XvpkSFIveLc5keewO1iNAGcQSnrx1rjfE7Jtgpp178CP5jZh4DHxL5WHSn6IS~K4uVoRK5YMkDjGrg4bRUqX04nSsFiTd0w7yFbes0jrHyCRmg__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q)) and GeneMark-ES ([Lomsadze et al., 2005](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1298918/); [Ter-hovhannisyan et al., 2008](http://genome.cshlp.org/content/18/12/1979.long)). All training was performed on the SPAdes assembly of a pure culture of *Letharia lupina* (see Tuovinen et al. (in prep.); coming soon, I promise!).

The following notes described how I ran my UPPMAX scripts, but noticed I skipped the call to `sbatch` for all scripts, for clarity. So in effect the scripts might be called like such:

    $ sbatch ./myscript

See the Milou guide about how to run UPPMAX jobs [here](http://www.uppmax.uu.se/support-sv/user-guides/milou-user-guide/).

## Training of *ab initio* gene predictors

For all training I used an assembly multifasta file that excluded all contigs smaller than 50000bp (called `L.lupinapure.scaffolds_ed.fa`).

### Creating a SNAP HMM file

Following the instructions from the [MAKER tutorial](http://gmod.org/wiki/MAKER_Tutorial#MAKER.27s_Output), a very useful [blog](https://biowize.wordpress.com/2012/06/01/training-the-snap-ab-initio-gene-predictor/) written by Daniel Standage, and the [SNAP](http://korflab.ucdavis.edu/software.html) `00README` file, I prepared a script to create a SNAP HMM file for *Letharia* data. 

    $ ./Le_SNAPhmm.sh L.lupinapure.scaffolds_ed.fa lupina

### Creating an GeneMark HMM file

The script to run and self-train GeneMark (called `gmes_petap.pl` and provided with the program) has some perl dependencies that are now installed in Uppmax. The instructions looks straight forward, so I only checked [this](https://wiki.gacrc.uga.edu/wiki/GeneMark) website, and the help menu. I prepared a script for it and run it as such:

    $ ./GeneMarkhmm.sh L.lupinapure.scaffolds_ed.fa

I'm running it with options `--ES --fungus --max_intron 3000 --min_gene_prediction 120`. I'm not giving it any hints. It took less than half an hour.

I was very happy with GeneMark-ES results: they matched the Cufflinks transcripts very well.

### Creating an Augustus HMM file

The pipeline to create an HMM training file is based on the file `README.autoAug` from the Augustus distribution. Basically the script `autoAug.pl` does everything for you. However I ran into MANY problems with Uppmax due my lack of rights to write on Augustus configuration folder. I ended up making a local copy of Augustus directory with symlinks to everything in the original installation in Uppmax, and I tell Augustus where my new local configuration directory is. On top of all that, I had to copy locally and modify a couple of scripts to make it run.

To train Augustus I tried a number of strategies, including using the protein models from the *Xanthoria parietina* draft genome in JGI (not so much info on how that was produced that I could find...) and using Trinity transcripts form the RNAseq data of *L. lupina*. However the resulting models were awful. In the end, I created a bunch of protein models using SNAP with the training file from above and the SPAdes genome assembly. I selected 1000 proteins at random and used those to train Augustus like this:

    $ ./Le_Augustushmm.sh L.lupinapure.scaffolds_ed.fa lupina


### Producing *ab initio* gene models

I made a master script to run all three *ab initio* gene predictors:

    $ ./AbInitio.sh L.lupinapure.scaffolds_ed.fa lupina > AbInitio.log

### Creating transcript models with STAR and Cufflinks

My script calls Douglas Scofield's script `mergePileupColumns` that can be obtained from his GitHub:

    $ wget https://raw.githubusercontent.com/douglasgscofield/bioinfo/master/scripts/mergePileupColumns
    $ chmod a+x mergePileupColumns

Now produce the transcripts:

    $ ./mapRNA_STAR.sh -f L.lupinapure.scaffolds.fasta -a LethariaRNA_R1.postQC.1.fq.gz -b LethariaRNA_R2.postQC.2.fq.gz -s lupina

That is, the files `LethariaRNA_R*.postQC.*.fq.gz` are post-quality control RNAseq reads.

I had trouble running it for *L. lucida* but it got fixed when I added the flag `--genomeSAIndexNbases 3` during the indexing of STAR.

### Running TransDecoder on Cufflink transcript models

[TransDecoder](http://transdecoder.github.io/) is used to identify likely coding regions and to prepare a GFF3 file corresponding to the corresponding protein-coding gene structures used later by EVM. Like explained [here](http://transdecoder.github.io/): "TransDecoder identifies candidate coding regions within transcript sequences, such as those generated by *de novo* RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq alignments to the genome using Tophat and Cufflinks".

I installed it locally in Uppmax and prepared a script based on the instructions in the website. It is possible to run [homology searches](http://transdecoder.github.io/#incl_homology) to protein databases (like Swissprot or Uniref90), but I decided for now to run it without them.

Within the script I specify that we have stranded RNAseq data.

    $ ./runTransDecoder.sh $SOMEPATH/lupina/transcripts.gtf L.lupinapure.scaffolds.fasta lupina

Where `$SOMEPATH` is related to where did you run my script `./mapRNA_STAR.sh`.

### Aligning the MAT idiomorph proteins with Exonerate

I used Exonerate ([Slater & Birney, 2005](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-31)) to make splice aware aligners of orthologs of SLA2 and APN2 from *Aspergillus oryzae* (protein sequence) that I took from the OrthoMCL-DB. I also took the MAT1-1-1 (or alpha) and MAT1-2-1 (HMG) genes of *Xanthoria polycarpa* (Accession numbers AJ884598 and AJ884599). All these protein sequences were in the file `Ref_MATgenes_aa.fa` and I only aligned it to the relevant scaffolds for each species (in this case `L.lupinapure_Spades_MAT-scaffold.fas`).

    $ ./runExonerate.sh Ref_MATgenes_aa.fa L.lupinapure_Spades_MAT-scaffold.fas

### Running EVidenceModeler

I ran EVidenceModeler (EVM) v.1.1.1 following the instructions in the [documentation](https://evidencemodeler.github.io/), all in a single script.

I put all the gff3 files from the previous steps in the directory `InputEvidence`. The EVM weights file was also there. The weights file (prepared based on several trials) of *L. lupina*, *L. columbiana*, and *L. rugosa*:

    $ cat lupina_weights.txt
    ABINITIO_PREDICTION Augustus    1
    ABINITIO_PREDICTION SNAP    1
    ABINITIO_PREDICTION GeneMark.hmm3   1
    PROTEIN exonerate:protein2genome:local  5
    TRANSCRIPT  Cufflinks   8
    OTHER_PREDICTION    transdecoder    10

For *L. vulpina* a different set of weights was better to recover proper models:

    ABINITIO_PREDICTION Augustus    1
    ABINITIO_PREDICTION SNAP    1
    ABINITIO_PREDICTION GeneMark.hmm3   2
    PROTEIN exonerate:protein2genome:local  20
    TRANSCRIPT  Cufflinks   10
    OTHER_PREDICTION    transdecoder    1

Finally, run EVM:

    $ cd ..
    $ ./EVM.sh -g L.lupinapure_Spades_MAT-scaffold.fas -w InputEvidence/lupina_weights.txt -p InputEvidence/lupina_predicted.gff3 -a InputEvidence/lupina_exonerate.gff3 -r InputEvidence/lupina_transcripts.gff3 -s lupina

The output looks like this:

    $ ll
    total 160
    drwxrwsr-x 2 sandral 2048 Jan  2 08:00 InputEvidence
    -rw-rw-r-- 1 sandral 1566 Jan  2 08:00 lupina_commands.list
    -rw-rw-rw- 1 sandral  546 Jan  2 08:00 lupina_partitions_list.out
    -rw-rw-r-- 1 sandral 1576 Jan  2 08:01 lupina_run.log
    drwxrwsrwx 4 sandral 2048 Jan  2 08:01 NODE_87_length_133277_cov_84.8955

And

    $ ll NODE_87_length_133277_cov_84.8955/
    total 800
    -rw-rw-rw- 1 sandral 135534 Jan  2 08:00 L.lupinapure_Spades_MAT-scaffold.fas
    -rw-rw-r-- 1 sandral  66041 Jan  2 08:01 lupina_evm.out
    -rw-rw-r-- 1 sandral  81501 Jan  2 08:01 lupina_evm.out.gff3
    -rw-rw-rw- 1 sandral   2971 Jan  2 08:00 lupina_exonerate.gff3_EVM
    -rw-rw-rw- 1 sandral 299545 Jan  2 08:00 lupina_predicted.gff3
    -rw-rw-rw- 1 sandral  21847 Jan  2 08:00 lupina_transcripts.gff3_EVM
    drwxrwsrwx 2 sandral   2048 Jan  2 08:00 NODE_87_length_133277_cov_84.8955_1-100000
    drwxrwsrwx 2 sandral   2048 Jan  2 08:01 NODE_87_length_133277_cov_84.8955_90001-133277

Our precious almost final annotation file is `lupina_evm.out.gff3`. I then compared all the species evidence, fixed a few details of intron-exon boundaries in Artemis ([Rutherford et al. 2000](https://www.ncbi.nlm.nih.gov/pubmed/11120685)), and finally got the final annotations. 