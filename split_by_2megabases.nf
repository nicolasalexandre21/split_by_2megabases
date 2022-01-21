//source activate nextflow
//module load bedtools java
//nextflow run -params-file params.yml split_by_2megabases.nf -c nextflow.config -dump-channels

nextflow.enable.dsl=2

reference_ch = file( params.ref_index, checkIfExists: true )  

Channel.fromPath( params.vcfs, checkIfExists: true )
    .dump(tag:"vcf_original")
    .into{ vcf_for_chromosomes;vcf_for_files }   

process makeWindows {
    input:
    path ref_index from reference_ch                

    output:
    path( "windows.bed" )  into windows_bed        
    
    script:
    """
    bedtools makewindows -g $ref_index -w 2000000 > windows.bed    
    """
}

windows_bed
    .splitCsv( header:false, sep:'\t', strip:true )
    .map { it -> [ it[0], 1+(it[1] as Integer), it[2] ] } /* convert bed 
+0 to interval +1 */
    .set { interval_ch } // Intervals

vcf_for_chromosomes.splitCsv(sep:"\t",limit:1)
    .dump(tag:"vcf_splitcsv")
    .map { it[0]}
    .set{ vcf_chromosomes_ch}

// Group the intervals with the matching VCF files, from the same chromosome
vcf_for_files
    // Get the chromosome from each per-chromosome VCF by reading the VCF file and taking the first column (the chromosome) in the first line
    .dump(tag:"vcf_before_map")    
    .map { tuple(it.baseName, it ) }
    .dump(tag:"vcf_tuple")    
    .cross( vcf_chromosomes_ch)
    .dump(tag:"vcf_tuple_chromosomes")
    .join( interval_ch, by:0, remainder:true )
    .dump(tag:"vcf_tuple_join")    
    .set { vcf_with_matched_intervals_ch }

process ImputeVcf {
    cpus 12

    input:
    tuple val(chrom), val(vcf_basename), val(vcf), val(start), val(end) 
from vcf_with_matched_intervals_ch 
    
    output:
    tuple val(vcf_basename), val("${chrom}:${start}-${end}"), 
path("${output}.txt")  into vcf_out

    script:
    output = "${vcf_basename}__${chrom}__${start}-${end}"
    """
    ./loimpute -i ${vcf} -ne 80000 -range ${start} ${end} -h 
refpanel.vcf.gz -o ${output}.txt
    """
}

process GatherVcf {
    cpus 12

    input:

    output:    

    script:
    output = "${vcf_basename}"
    """
    java -Xmx32G -jar 
/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar 
GatherVcfs I=?? O=${vcf_basename}.vcf
    """
}

process GzipVcf {
    cpus 12

    input:
    output:

    script:
    output = "${vcf_basename}.vcf"
    """
    gzip ${vcf_basename}.vcf
    """
}
