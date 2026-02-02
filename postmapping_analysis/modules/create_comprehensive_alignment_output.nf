process CREATE_COMPREHENSIVE_ALIGNMENT_OUTPUT {
    tag "comprehensive_aln_${sample_name}"
    publishDir "${params.outdir}/comprehensive_alignments", mode: 'copy'

    input:
    path curated_dir                     // Directory with curated confident reads BED files
    path paf_b73                         // B73 PAF file with segments
    path paf_mo17                        // Mo17 PAF file with segments
    val sample_name                      // Sample name for output files

    output:
    path "${sample_name}_comprehensive_alignments.tsv", emit: comprehensive_tsv
    path "${sample_name}_confident_alignments_detailed.tsv", emit: detailed_tsv
    path "${sample_name}_confident_*_premapping_readcoords.bed", emit: premapping_bed_files
    path "${sample_name}_confident_*_postmapping_refcoords.bed", emit: postmapping_bed_files

    script:
    """
    Rscript ${projectDir}/scripts/19-create_comprehensive_alignment_output.R \\
        --curated_dir ${curated_dir} \\
        --paf_b73 ${paf_b73} \\
        --paf_mo17 ${paf_mo17} \\
        --sample_name ${sample_name} \\
        --output_dir .
    """
}
