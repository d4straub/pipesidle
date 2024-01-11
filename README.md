# pipesidle

Multipe region amplicon sequencing analysis with SMURF via Sidle in QIIME2.

The pipeline requires output from [nf-core/ampliseq](https://nf-co.re/ampliseq) and applies [Sidle (SMURF Implementation Done to acceLerate Efficiency)](https://github.com/jwdebelius/q2-sidle) within [QIIME2](https://qiime2.org/) to scaffolding multiple regions along a reference for improved resolution over a single region. 

For example, multiple variable regions of the 16S rRNA gene were sequenced with various primers and need to be unified. This leads to one unified abundance and taxonomy profile over all variable regions.

## Software requirement

- Java
- nextflow
- singularity / apptainer / docker (tested with singularity)

### Input data

Information to computed abundances, sequences, primers, and expected length, for each sequenced region. Pipesidle was made for output of [nf-core/ampliseq](https://nf-co.re/ampliseq), but should also accept any other pipeline output in correct format.
The information has to be provided in a tab-separated sample sheet via `--input`, one row per sequenced region.

| Parameter     | Description                      |
| ------------- | -------------------------------- |
| sampleID      | Unique ID*                       |
| asv_table     | Path to ASV count table (tsv)    |
| asv_seq       | Path to ASV sequences (fasta)    |
| primerfw      | Forward primer sequence (string) |
| primerrv      | Reverse primer sequence (string) |
| region_length | Minimal region length (integer)  |

> **Warning:**
> alpha-numerically sorted sampleID's must reflect sequence on the reference sequence, e.g. region with sampleID `1` must be located before region `2` in the reference database sequences.

Example sample sheet:

```csv
sampleID	asv_table	asv_seq	primerfw	primerrv	region_length
1	nfcoreampliseq/results_TGGCGAACGGGTGAGTAA_CCGTGTCTCAGTCCCARTG/dada2/ASV_table.tsv	nfcoreampliseq/results_TGGCGAACGGGTGAGTAA_CCGTGTCTCAGTCCCARTG/dada2/ASV_seqs.fasta	TGGCGAACGGGTGAGTAA	CCGTGTCTCAGTCCCARTG	145
2	nfcoreampliseq/results_ACTCCTACGGGAGGCAGC_GTATTACCGCGGCTGCTG/dada2/ASV_table.tsv	nfcoreampliseq/results_ACTCCTACGGGAGGCAGC_GTATTACCGCGGCTGCTG/dada2/ASV_seqs.fasta	ACTCCTACGGGAGGCAGC	GTATTACCGCGGCTGCTG	135
3	nfcoreampliseq/results_GTGTAGCGGTGRAATGCG_CCCGTCAATTCMTTTGAGTT/dada2/ASV_table.tsv	nfcoreampliseq/results_GTGTAGCGGTGRAATGCG_CCCGTCAATTCMTTTGAGTT/dada2/ASV_seqs.fasta	GTGTAGCGGTGRAATGCG	CCCGTCAATTCMTTTGAGTT	200
4	nfcoreampliseq/results_GGAGCATGTGGWTTAATTCGA_CGTTGCGGGACTTAACCC/dada2/ASV_table.tsv	nfcoreampliseq/results_GGAGCATGTGGWTTAATTCGA_CGTTGCGGGACTTAACCC/dada2/ASV_seqs.fasta	GGAGCATGTGGWTTAATTCGA	CGTTGCGGGACTTAACCC	115
5	nfcoreampliseq/results_GGAGGAAGGTGGGGATGAC_AAGGCCCGGGAACGTATT/dada2/ASV_table.tsv	nfcoreampliseq/results_GGAGGAAGGTGGGGATGAC_AAGGCCCGGGAACGTATT/dada2/ASV_seqs.fasta	GGAGGAAGGTGGGGATGAC	AAGGCCCGGGAACGTATT	150
```

### Taxonomic reference database

Three files are required:
- `fasta` file with reference sequences
- file with aligned reference sequences
- file with taxonomies of reference sequences

Here code for the use of greengenes 13_8:

```bash
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
tar -zxvf gg_13_8_otus.tar.gz gg_13_8_otus/rep_set/99_otus.fasta
mv gg_13_8_otus/rep_set/99_otus.fasta gg_13_8_otus_rep_set_99_otus.fasta
tar -zxvf gg_13_8_otus.tar.gz gg_13_8_otus/rep_set_aligned/99_otus.fasta
mv gg_13_8_otus/rep_set_aligned/99_otus.fasta gg_13_8_otus_rep_set_aligned_99_otus.fasta
tar -zxvf gg_13_8_otus.tar.gz gg_13_8_otus/taxonomy/99_otu_taxonomy.txt
mv gg_13_8_otus/taxonomy/99_otu_taxonomy.txt gg_13_8_otus_taxonomy_99_otu_taxonomy.txt
rm -r gg_13_8_otus
```

Alternatively, SILVA 128 has also all required files.

### Run the pipeline

As example with the above greengenes 13_8 reference database.

```bash
NXF_VER=23.10.0 nextflow run d4straub/pipesidle \
  -profile singularity
  --input 'samplesheet.tsv' \
  --db_sequences 'gg_13_8_otus_rep_set_99_otus.fasta' \
  --db_alignedsequences 'gg_13_8_otus_rep_set_aligned_99_otus.fasta' \
  --db_taxonomy 'gg_13_8_otus_taxonomy_99_otu_taxonomy.txt' \
  --db_tree 'https://data.qiime2.org/2021.4/common/sepp-refs-gg-13-8.qza' \
  --outdir result
```

The `--db_tree` parameter is optional but recommended.

Additional parameters are
- `--sidle_save_intermediate`: Save a bunch of intermediate files (default: false)
- `--max_memory`: maximum RAM allowed per process, to scale to computational resources at hand (default: `1992.GB`)
- `--max_cpus`: number of maximum cpus per process, to scale to computational resources at hand (default: `128`)
- `--max_time`: maximum time allowed per process, to scale to computational resources at hand  (default: `168.h`)

## Output

This section describes only the most important output files.

- `reconstructed/`
  - `reconstructed_feature-table.biom`: Unified abundance table in biom format
  - `reconstructed_feature-table.tsv`: Tab-separated unified abundance table
  - `reconstructed_taxonomy.tsv`: Tab-separated unified taxonomy table 
  - `reconstructed_merged.tsv`: Tab-separated unified table with merged abundance and taxonomy information
  - `reconstructed_tree.nwk`: Phylogenetic tree

- `barplot/`
  - `index.html`: Interactive barplot, open in your browser (i.e. double-click it)

## Credits

This pipeline was originally written by Daniel Straub ([@d4straub](https://github.com/d4straub)) for use at the [Quantitative Biology Center (QBiC)](http://www.qbic.life).

Code was inspired by [nf-core/ampliseq](https://nf-co.re/ampliseq) (doi: 10.5281/zenodo.1493841) ([Straub et al., 2020](https://doi.org/10.3389/fmicb.2020.550420)) of the [nf-core](https://nf-co.re) collection of workflows ([Ewels et al., 2020](https://dx.doi.org/10.1038/s41587-020-0439-x)).

## Citations

Please cite all tools that you used, i.e. all listed here, but only one containerization tool (here Singularity or Docker):

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [q2-sidle](https://doi.org/10.1101/2021.03.23.436606)

  > Debelius, J.W.; Robeson, M.; Lhugerth, L.W.; Boulund, F.; Ye, W.; Engstrand, L. "A comparison of approaches to scaffolding multiple regions along the 16S rRNA gene for improved resolution." Preprint in BioRxiv. doi: 10.1101/2021.03.23.436606

- [SMURF](https://doi.org/10.1186/s40168-017-0396-x)

  > Fuks, G.; Elgart, M.; Amir, A.; Zeisel, A.; Turnbaugh, P.J., Soen, Y.; and Shental, N. (2018). "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." Microbiome. 6: 17. doi: 10.1186/s40168-017-0396-x

- [RESCRIPt](https://doi.org/10.1371/journal.pcbi.1009581)

  > Robeson MS 2nd, O'Rourke DR, Kaehler BD, Ziemski M, Dillon MR, Foster JT, Bokulich NA. RESCRIPt: Reproducible sequence taxonomy reference database management. PLoS Comput Biol. 2021 Nov 8;17(11):e1009581. doi: 10.1371/journal.pcbi.1009581. PMID: 34748542; PMCID: PMC8601625.

- [SEPP](https://doi.org/10.1128/msystems.00021-18)

  > Janssen S, McDonald D, Gonzalez A, Navas-Molina JA, Jiang L, Xu ZZ, Winker K, Kado DM, Orwoll E, Manary M, Mirarab S, Knight R. Phylogenetic Placement of Exact Amplicon Sequences Improves Associations with Clinical Information. mSystems. 2018 Apr 17;3(3):e00021-18. doi: 10.1128/mSystems.00021-18. PMID: 29719869; PMCID: PMC5904434.

- [QIIME2](https://pubmed.ncbi.nlm.nih.gov/31341288/)

  > Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS 2nd, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, Caporaso JG. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol. 2019 Aug;37(8):852-857. doi: 10.1038/s41587-019-0209-9. Erratum in: Nat Biotechnol. 2019 Sep;37(9):1091. PMID: 31341288; PMCID: PMC7015180.
