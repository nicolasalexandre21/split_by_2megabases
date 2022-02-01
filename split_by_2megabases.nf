//source activate nextflow
//module load bedtools java
//nextflow run -params-file params.yml split_by_2megabases.nf -c nextflow.config -dump-channels

nextflow.enable.dsl = 2

params.jarfile = '/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar' 

workflow {

    reference_ch = file( params.ref_index, checkIfExists: true ) 
    MAKE_WINDOWS ( reference_ch )
    MAKE_WINDOWS.out.bed  // Replace channel name with process output stream
        .splitCsv( header:false, sep:'\t', strip:true )
        .map { it -> [ it[0], 1+(it[1] as Integer), it[2] ] } /* convert bed +0 to interval +1 */
        .set { interval_ch } // Intervals

    Channel.fromPath( params.vcfs, checkIfExists: true )
        // .into{ vcf_for_chromosomes; vcf_for_files }  // Channel forking is now implicit
        .set { vcf_ch }
    // Group the intervals with the matching VCF files, from the same chromosome
    vcf_ch
        .map { tuple(it.simpleName, it ) } // Convert to [ id, file(id_vcf) ] 
        // Get the chromosome from each per-chromosome VCF by reading the VCF file and taking the first column (the chromosome) in the first line
        // .cross( interval_ch ) // Wrong operation I think.
        .combine( interval_ch )            //  Builds tuple [ id, file(id_vcf), chr, start, stop ]
        .set { vcf_with_matched_intervals_ch }
    IMPUTE_VCF ( vcf_with_matched_intervals_ch )
    GATHER_VCF ( 
        IMPUTE_VCF.out.vcf.groupTuple(),  // Inputs are [ id, [vcf1, vcf2, vcf3 ] ]
        file( params.jarfile, checkIfExists: true )  // Implicitly converted to value type ch.  
    )
}

process MAKE_WINDOWS { // Name changed to uppercase only for readability purposes

    input:
    path ref_index // from reference_ch                

    output:
    path( "windows.bed" ), emit: bed //  into windows_bed        

    script:
    """
    bedtools makewindows -g $ref_index -w 2000000 > windows.bed    
    """
}

process IMPUTE_VCF {

    cpus 12

    input:
    tuple val(vcf_basename), path(vcf), val(chrom), val(start), val(end) // from vcf_with_matched_intervals_ch 

    output:
    tuple val(vcf_basename), path("${output}.vcf.gz"), emit: vcf // into vcf_out

    script:
    output = "${vcf_basename}__${chrom}__${start}-${end}"
    """
    ./loimpute -i ${vcf} -ne 80000 -range ${start} ${end} -h 
refpanel.vcf.gz -o ${output}
    """
}


process GATHER_VCF {

    cpus 12

    input:
    tuple val (vcf_basename), path (vcfs) // from vcf_out.groupTuple()
    path jarfile

    output:    
    path("${vcf_basename}.vcf") // into vcf_merged 

    script:
    """
    cat << EOF > tmp.list
    ${vcfs.join("\n")}
    EOF
    java -Xmx32G -jar $jarfile \\
        GatherVcfs I=tmp.list O=${vcf_basename}
    """
}

process GZIP_VCF { // Probably simpler to merge this command into the one above. 

    cpus 12

    input:
    tuple val(vcf_basename)

    output:
    path("${vcf_basename}.vcf.gz") // into vcf_zipped 

    script:
    output = "${vcf_basename}.vcf.gz"
    """
    gzip ${vcf_basename}.vcf
    """
}
