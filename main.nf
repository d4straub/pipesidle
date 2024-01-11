//DEFAULT PARAMS

params.input = false // --input "samplesheet.tsv"
params.db_sequences = false // --db_sequences "gg_13_8_otus_rep_set_99_otus.fasta"
params.db_taxonomy = false // --db_taxonomy "gg_13_8_otus_taxonomy_99_otu_taxonomy.txt"
params.db_alignedsequences = false // --db_alignedsequences "gg_13_8_otus_rep_set_aligned_99_otus.fasta"
params.db_tree = false // --db_tree "https://data.qiime2.org/2021.4/common/sepp-refs-gg-13-8.qza"
if(!params.input){exit 1, log.info "Specify input samplesheet. Example: --input \"samplesheet.tsv\""}
if(!params.db_sequences){exit 1, log.info "Specfify database sequence file. Example: --db_sequences \"gg_13_8_otus_rep_set_99_otus.fasta\""}
if(!params.db_taxonomy){exit 1, log.info "Specfify database taxonomy file. Example: --db_taxonomy \"gg_13_8_otus_taxonomy_99_otu_taxonomy.txt\""}
if(!params.db_alignedsequences){exit 1, log.info "Specfify file with aligned sequences of database. Example: --db_alignedsequences \"gg_13_8_otus_rep_set_aligned_99_otus.fasta\""}

//CHANNELS

ch_db_sequences = file(params.db_sequences)
ch_db_taxonomy = file(params.db_taxonomy)
ch_db_alignedsequences = file(params.db_alignedsequences)
ch_db_tree = params.db_tree ? file(params.db_tree) : Channel.from([])

// PARSE SAMPLESHEET

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def parse_samplesheet(LinkedHashMap row) {
    //Check if samplesheet contains column sampleID  & asv_table & primerfw & primerrv
    if (row.sampleID == null || row.asv_table == null || row.asv_seq == null || row.primerfw == null || row.primerrv == null || row.region_length == null) {
        error("ERROR: Please check input samplesheet -> Column 'sampleID', 'asv_table', 'asv_seq', 'primerfw', 'primerrv', and 'region_length' are required but not detected.")
    }
    //read meta info
    def meta = [:]
    meta.id         = row.sampleID
    meta.primerfw   = row.primerfw
    meta.primerrv   = row.primerrv
    meta.region_length = row.region_length
    //read data info
    def array = []
    array = [ meta, file(row.asv_table), file(row.asv_seq) ]
    return array
}

// Sample sheet input
tsvFile = file(params.input).getName()
// extracts read files from TSV and distribute into channels
Channel
    .fromPath(params.input)
    .ifEmpty { error("Cannot find path file ${tsvFile}") }
    .splitCsv(header:true, sep:'\t')
    .map { parse_samplesheet(it) }
    .set { ch_samplesheet }

//PROCESSES

process QIIME2_IN {
    tag "$meta.id"
    label 'process_single'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(table), path(seq)

    output:
    tuple val(meta), path("*_table.qza"), path("*_rep-seqs.qza"), emit: table_seq
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # seq
    qiime tools import \\
        --input-path "$seq" \\
        --type 'FeatureData[Sequence]' \\
        --output-path ${prefix}_rep-seqs.qza

    # table
    biom convert -i "$table" -o table.biom --table-type="OTU table" --to-hdf5
    qiime tools import \\
        --input-path table.biom \\
        --type 'FeatureTable[Frequency]' \\
        --input-format BIOMV210Format \\
        --output-path ${prefix}_table.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

process QIIME2_INDB {
    label 'process_single'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(seq)
    path(tax)

    output:
    path("db_sequences.qza"), emit: seq
    path("db_taxonomy.qza") , emit: tax
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # db_seq
    qiime tools import \\
        --input-path $seq \\
        --output-path db_sequences.qza \\
        --type 'FeatureData[Sequence]'

    # db_tax
    qiime tools import \\
        --input-path $tax \\
        --output-path db_taxonomy.qza \\
        --type 'FeatureData[Taxonomy]' \\
        --input-format HeaderlessTSVTaxonomyFormat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

process QIIME2_INDBALIGNED {
    label 'process_single'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(seq)

    output:
    path("db_alignedsequences.qza"), emit: seq
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # db_seq
    qiime tools import \\
        --input-path $seq \\
        --output-path db_alignedsequences.qza \\
        --type 'FeatureData[AlignedSequence]'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

process DB_PREFILTERING {
    //tag "$meta.id"
    label 'process_low'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(seq)
    path(tax)

    output:
    path("db_filtered_sequences.qza")     , emit: seq
    path("db_filtered_sequences_tax.qza") , emit: tax
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    # https://q2-sidle.readthedocs.io/en/latest/database_preparation.html#filtering-the-database
    #pre-filtering should be very permissive!
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # authors of SMURF recommend "--p-num-degenerates 3" for greengenes 13_8 database at 99%
    # the RESCRIPt formatted Silva 128 database is filtered to exclude sequences with more than 5 degenerates [3], [4]
    qiime rescript cull-seqs \\
        --p-n-jobs $task.cpus \\
        --i-sequences $seq \\
        $args \\
        --o-clean-sequences db_filtered_sequences.qza
    
    #filtering a greengenes database for features missing a phylum (p__;) or kingdom(k__;) designation.
    #CPU=1
    qiime taxa filter-seqs \\
        --i-sequences db_filtered_sequences.qza \\
        --i-taxonomy $tax \\
        $args2 \\
        --o-filtered-sequences db_filtered_sequences_tax.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 rescript sidle\$( qiime rescript --version | sed 's/ (.*//' | sed 's/.*version //' )
    END_VERSIONS
    """
}

process DB_EXTRACTION {
    tag "$meta.id,$meta.region_length"
    label 'process_medium'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(table), path(seq), path(db_seq), path(db_tax)

    output:
    tuple val(meta), path("db_*_kmers.qza"), emit: kmers
    tuple val(meta), path("db_*_map.qza")  , emit: map
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def primerfw = "${meta.primerfw}"
    def primerrv = "${meta.primerrv}"
    def length = "${meta.region_length}"
    """
    # https://q2-sidle.readthedocs.io/en/latest/database_preparation.html#prepare-a-regional-database-for-each-primer-set
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"
    
    #extract sequences
    qiime feature-classifier extract-reads \\
        --p-n-jobs $task.cpus \\
        --i-sequences $db_seq \\
        $args \\
        --p-f-primer $primerfw \\
        --p-r-primer $primerrv \\
        --o-reads db_${prefix}.qza
    
    #prepare to be used in alignment
    qiime sidle prepare-extracted-region \\
        --p-n-workers $task.cpus \\
        --i-sequences db_${prefix}.qza \\
        --p-region "${prefix}" \\
        --p-fwd-primer $primerfw \\
        --p-rev-primer $primerrv \\
        --p-trim-length $length \\
        --o-collapsed-kmers db_${prefix}_${length}_kmers.qza \\
        --o-kmer-map db_${prefix}_${length}_map.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle\$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process ASV_TRIM {
    tag "$meta.id,$meta.region_length"
    label 'process_single'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(table), path(seq)

    output:
    tuple val(meta), path("*_table.qza")    , emit: table
    tuple val(meta), path("*_rep-seqs.qza") , emit: seq
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def primerfw = "${meta.primerfw}"
    def primerrv = "${meta.primerrv}"
    def length = "${meta.region_length}"
    """
    # https://q2-sidle.readthedocs.io/en/latest/read_preparation.html#dada2
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    #CPU=1
    qiime sidle trim-dada2-posthoc \
        --i-table ${table} \
        --i-representative-sequences ${seq} \
        --p-trim-length $length \
        --o-trimmed-table ${prefix}_${length}_table.qza \
        --o-trimmed-representative-sequences ${prefix}_${length}_rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle\$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process REGIONAL_ALIGNMENT {
    tag "$meta.id"
    label 'process_medium'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(kmers), path(seq)

    output:
    tuple val(meta), path("*rep-seqs_align-map.qza"), emit: aligned_map
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def primerfw = "${meta.primerfw}"
    def primerrv = "${meta.primerrv}"
    """
    # https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#regional-alignment
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"
    
    qiime sidle align-regional-kmers \\
        --p-n-workers $task.cpus \\
        --i-kmers ${kmers} \\
        --i-rep-seq ${seq} \\
        --p-region ${meta.id} \\
        $args \\
        --o-regional-alignment ${prefix}_rep-seqs_align-map.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle\$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process DB_RECONSTRUCTION {
    label 'process_medium'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    val(metaid)
    path(map)
    path(aligned_map)

    output:
    path("reconstruction_map.qza")    , emit: reconstruction_map
    path("reconstruction_summary.qza"), emit: reconstruction_summary
    path("reconstruction_summary/*")  , emit: visualisation
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_input = ""
    // sort the input so that the regions are sorted by sequence
    def df = [metaid, map, aligned_map].transpose().sort{ it[0] }
    for (i in df) {
        db_input += " --p-region "+i[0]+" --i-kmer-map "+i[1]+" --i-regional-alignment "+i[2]
    }
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#database-reconstruction
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime sidle reconstruct-database \\
        --p-n-workers $task.cpus \\
        $db_input \\
        $args \\
        --o-database-map reconstruction_map.qza \\
        --o-database-summary reconstruction_summary.qza

    #database summary can be used to evaluate the quality of the reconstruction; see Fuks, C; Elgart, M; Amir, A; et al (2018) “Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling.” Microbiome. 6:17. doi: 10.1186/s40168-017-0396-x
    qiime metadata tabulate \\
        --m-input-file reconstruction_summary.qza \\
        --o-visualization reconstruction_summary.qzv
    qiime tools export \\
        --input-path reconstruction_summary.qzv \\
        --output-path "reconstruction_summary"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle\$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process TABLE_RECONSTRUCTION {
    label 'process_medium'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    val(metaid)
    path(table)
    path(aligned_map)
    path(reconstruction_map)
    path(reconstruction_summary)

    output:
    path("reconstruction_table.qza")        , emit: qza
    path("reconstruction_table/*")          , emit: exported
    path("reconstructed_feature-table.biom"), emit: biom
    path("reconstructed_feature-table.tsv") , emit: tsv
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def region_input = ""
    // sort the input so that the regions are sorted by sequence
    def df = [metaid, aligned_map, table].transpose().sort{ it[0] }
    for (i in df) {
        region_input += " --p-region "+i[0]+" --i-regional-alignment "+i[1]+" --i-regional-table "+i[2]
    }
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#table-reconstruction
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime sidle reconstruct-counts \\
        --p-n-workers $task.cpus \\
        $region_input \\
        --i-database-map $reconstruction_map \\
        --i-database-summary $reconstruction_summary \\
        $args \\
        --o-reconstructed-table reconstruction_table.qza

    #export visualisation
    qiime feature-table summarize \\
        --i-table reconstruction_table.qza \\
        --o-visualization reconstruction_table.qzv
    qiime tools export \\
        --input-path reconstruction_table.qzv \\
        --output-path "reconstruction_table"

    #export feature table in biom and tsv format
    qiime tools export \\
        --input-path reconstruction_table.qza \\
        --output-path exported
    biom convert \\
        -i exported/feature-table.biom \\
        -o reconstructed_feature-table.tsv \\
        --to-tsv
    cp exported/feature-table.biom reconstructed_feature-table.biom

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle\$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process TAXONOMIC_RECONSTRUCTION {
    label 'process_single'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(reconstruction_map)
    path(tax)

    output:
    path("reconstruction_taxonomy.qza"), emit: qza
    path("reconstruction_taxonomy/*")  , emit: visualisation
    path("reconstruction_taxonomy.tsv"), emit: tsv
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#taxonomic-reconstruction
    #https://forum.qiime2.org/t/sidle-reconstruct-database/25439
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    #CPU=1
    qiime sidle reconstruct-taxonomy \\
        --i-reconstruction-map ${reconstruction_map} \\
        --i-taxonomy ${tax} \\
        $args \\
        --o-reconstructed-taxonomy reconstruction_taxonomy.qza

    #export visualisation
    qiime metadata tabulate \\
        --m-input-file reconstruction_taxonomy.qza \\
        --o-visualization reconstruction_taxonomy.qzv
    qiime tools export \\
        --input-path reconstruction_taxonomy.qzv \\
        --output-path "reconstruction_taxonomy"

    #export taxonomic tsv
    qiime tools export \\
        --input-path reconstruction_taxonomy.qza \\
        --output-path exported
    cp exported/taxonomy.tsv reconstruction_taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle\$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process FILTER_TAXONOMY {
    label 'process_single'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(table_tofilter)
    path(table_ref)

    output:
    path("reconstructed_taxonomy.tsv"), emit: filtered
    path("reconstructed_merged.tsv")  , emit: merged
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    df_tofilter <- read.table("$table_tofilter", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(df_tofilter)[1] <- "ID"

    df_ref <- read.table("$table_ref", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1, comment.char = "")
    colnames(df_ref)[1] <- "ID"

    df_merged <- merge(df_tofilter, df_ref, by="ID", all.x=FALSE, all.y=TRUE)
    write.table(df_merged, file = "reconstructed_merged.tsv", row.names=FALSE, sep="\t")

    df_filtered <- subset(df_tofilter, df_tofilter\$ID %in% df_ref\$ID)
    write.table(df_filtered, file = "reconstructed_taxonomy.tsv", row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")) ), "versions.yml")
    """
}

process FRAGMENT_RECONSTRUCTION {
    label 'process_single'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(reconstruction_map)
    path(reconstruction_summary)
    path(db_aligned_sequences)

    output:
    path("reconstruction_fragments.qza") , emit: qza
    path("reconstruction_fragments/*")   , emit: visualisation
    path("reconstructed_fragments.fasta"), emit: fasta
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#reconstructing-the-phylogenetic-tree
    #https://forum.qiime2.org/t/sidle-tutorial-missing-aligned-sequence-file/20604/4 for db_aligned_sequences
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    #CPU=1
    qiime sidle reconstruct-fragment-rep-seqs \\
        --i-reconstruction-map ${reconstruction_map} \\
        --i-reconstruction-summary ${reconstruction_summary} \\
        --i-aligned-sequences ${db_aligned_sequences} \\
        --o-representative-fragments reconstruction_fragments.qza

    #export visualisation
    qiime metadata tabulate \\
        --m-input-file reconstruction_fragments.qza \\
        --o-visualization reconstruction_fragments.qzv
    qiime tools export \\
        --input-path reconstruction_fragments.qzv \\
        --output-path "reconstruction_fragments"

    #export fasta file
    qiime tools export \\
        --input-path reconstruction_fragments.qza \\
        --output-path exported
    cp exported/dna-sequences.fasta reconstructed_fragments.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle\$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process TREE_RECONSTRUCTION {
    label 'process_medium'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(reconstruction_fragments)
    path(ref_db_tree)

    output:
    path("reconstructed_tree.qza")       , emit: qza
    path("reconstruction_placements.qza"), emit: qza_placements
    path("reconstructed_tree.nwk")       , emit: nwk
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#reconstructing-the-phylogenetic-tree
    # required: SEPP file https://forum.qiime2.org/t/sidle-tutorial-missing-aligned-sequence-file/20604/8
    # SEPP file only available for Greengenes 13_8 or SILVE 128 (not 138!): https://forum.qiime2.org/t/error-in-reconstructing-the-phylogenetic-tree/23757/8
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime fragment-insertion sepp \\
        --p-threads $task.cpus \\
        --i-representative-sequences $reconstruction_fragments \\
        --i-reference-database $ref_db_tree \\
        --o-tree reconstructed_tree.qza \\
        --o-placements reconstruction_placements.qza

    #export tree file
    qiime tools export \\
        --input-path reconstructed_tree.qza \\
        --output-path exported
    cp exported/tree.nwk reconstructed_tree.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        q2-fragment-insertion: \$( qiime fragment-insertion --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

// just to get a fast pretty picture
process QIIME2_BARPLOT {
    label 'process_low'

    container 'd4straub/pipesidle:0.1.0-beta'

    input:
    path(metadata)
    path(table)
    path(taxonomy)
    val(setting)

    output:
    path("barplot${suffix}/*"), emit: folder
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    suffix = setting ? "_${table.baseName}" : ""
    def metadata_cmd = metadata ? "--m-metadata-file ${metadata}": ""
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime taxa barplot  \\
        --i-table ${table}  \\
        --i-taxonomy ${taxonomy}  \\
        ${metadata_cmd}  \\
        --o-visualization taxa-bar-plots.qzv  \\
        --verbose
    qiime tools export \\
        --input-path taxa-bar-plots.qzv  \\
        --output-path barplot${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

//PIPELINE

workflow {
    // DB
    QIIME2_INDB ( ch_db_sequences, ch_db_taxonomy )
    QIIME2_INDBALIGNED ( ch_db_alignedsequences )
    DB_PREFILTERING ( QIIME2_INDB.out.seq, QIIME2_INDB.out.tax )

    // ASV
    QIIME2_IN ( ch_samplesheet )
    ASV_TRIM ( QIIME2_IN.out.table_seq )

    // combine & reconstruct
    DB_EXTRACTION (
        QIIME2_IN.out.table_seq
            .combine( DB_PREFILTERING.out.seq )
            .combine( DB_PREFILTERING.out.tax ) )
    REGIONAL_ALIGNMENT ( DB_EXTRACTION.out.kmers.join(ASV_TRIM.out.seq).dump(tag: 'into_REGIONAL_ALIGNMENT') )

    DB_EXTRACTION.out.map
        .join(REGIONAL_ALIGNMENT.out.aligned_map)
        .multiMap { meta, map, aligned_map ->
            sampleid: meta.id
            map: map
            aligned_map: aligned_map
        }
        .set { ch_db_reconstruction }

    DB_RECONSTRUCTION (
        ch_db_reconstruction.sampleid.collect(),
        ch_db_reconstruction.map.collect(),
        ch_db_reconstruction.aligned_map.collect() )

    ASV_TRIM.out.table
        .join(REGIONAL_ALIGNMENT.out.aligned_map)
        .multiMap { meta, table, aligned_map ->
            sampleid: meta.id
            table: table
            aligned_map: aligned_map
        }
        .set { ch_table_reconstruction }

    // Abundance table
    TABLE_RECONSTRUCTION (
        ch_table_reconstruction.sampleid.collect(),
        ch_table_reconstruction.table.collect(),
        ch_table_reconstruction.aligned_map.collect(),
        DB_RECONSTRUCTION.out.reconstruction_map,
        DB_RECONSTRUCTION.out.reconstruction_summary )

    // Taxonomic classification
    TAXONOMIC_RECONSTRUCTION (
        DB_RECONSTRUCTION.out.reconstruction_map,
        QIIME2_INDB.out.tax )
    FILTER_TAXONOMY ( TAXONOMIC_RECONSTRUCTION.out.tsv, TABLE_RECONSTRUCTION.out.tsv )

    // Reconstruct sequences/fragments
    // required: aligned sequences file: https://forum.qiime2.org/t/finding-alignment-files-for-sidle/23773/2
    FRAGMENT_RECONSTRUCTION (
        DB_RECONSTRUCTION.out.reconstruction_map,
        DB_RECONSTRUCTION.out.reconstruction_summary,
        QIIME2_INDBALIGNED.out.seq )
    // "The output of reconstruct-fragment-rep-seqs provides consensus sequences only if a reference sequence can't be resolved (ids that have a | symbol in them.) It's designed specifically to integrate with the fragment insertion and makes some downstream assumptions, including that you have the same database and insertion tree version.", see https://forum.qiime2.org/t/how-to-merge-q2-sidle-output-with-other-results/22823/2

    // Reconstruct phylogenetic tree
    TREE_RECONSTRUCTION (
        FRAGMENT_RECONSTRUCTION.out.qza,
        ch_db_tree )

    // just to get a fast pretty picture
    QIIME2_BARPLOT ( 
        [],
        TABLE_RECONSTRUCTION.out.qza,
        TAXONOMIC_RECONSTRUCTION.out.qza,
        []
    )
}
