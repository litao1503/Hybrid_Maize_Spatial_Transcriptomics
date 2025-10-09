version 1.0

import "spaseASE_tasks.wdl" as tasks

struct SampleInfo_Visium {
    String sample_name
    String fastq1
    String fastq2
    String image
    String slide
    String area
    String slidefile
}

struct SampleInfo_RNA {
    String sample_name
    String fastq1
    String fastq2
}

struct SampleInfo_longReads {
    String sample_name
    String fastq1
}


workflow spaseASE {
    input {
        String sample_sheet
        String analysis_dir

        String nmasked_star_index = ""
        String snp_file = ""
        String STAR_path = ""
        String SNPsplit_path = ""

        # Options for 10x spaceranger
        String transcriptome = ""
        String spaceranger = ""

        # Option for Bulk RNAseq and longReads RNAseq
        String gtf = ""

        # Option for longReads RNAseq
        String hapA_only_kmers = ""
        String hapB_only_kmers = ""
        String hapA_prefix = ""
        String hapB_prefix = ""
        String reference_fa = ""
        String minimap2_path = ""

        Boolean run_SNPsplit_scrna = false
        Boolean run_SNPsplit_rna = false
        Boolean run_SNPsplit_longReads = false

        String pipefail_check = "true"
    }

    Array[String] hap_type_array = ["genome1", "genome2"]
    Array[String] hap_type_array_long = [hapA_prefix, hapB_prefix, "unclassified"]

    Array[Array[String]] sample_sheet_df = read_tsv(sample_sheet)
    scatter (row in sample_sheet_df) {

        if (run_SNPsplit_scrna) {
            # init SampleInfo
            SampleInfo_Visium sample = {
                "sample_name": row[0],
                "fastq1": row[1],
                "fastq2": row[2],
                "image": row[3],
                "slide": row[4],
                "area": row[5],
                "slidefile": row[6],
            }

            call tasks.STARmapping as STARmapping_scrna {
                input:
                    sample_name = sample.sample_name,
                    nmasked_star_index = nmasked_star_index,
                    fastq1 = sample.fastq2,
                    fastq2 = "",
                    analysis_dir = analysis_dir,
                    pipefail_check = pipefail_check,
                    threads = 60,
                    STAR_path = STAR_path,
            }

            call tasks.BamProcess as BamProcess_scrna {
                input:
                    sample_name = sample.sample_name,
                    mapping_res = STARmapping_scrna.star_bam,
                    filter_expression = "[NH] == 1",
                    result_dir = analysis_dir + "/" + sample.sample_name + "/STARmapping/",
                    pipefail_check = pipefail_check,
                    threads = 4,
            }

            call tasks.SNPsplit as SNPsplit_scrna {
                input:
                    sample_name = sample.sample_name,
                    bamfile = BamProcess_scrna.bam,
                    seqtype = "--single_end",
                    snp_file = snp_file,
                    analysis_dir = analysis_dir,
                    pipefail_check = pipefail_check,
                    SNPsplit_path = SNPsplit_path,
            }

            scatter (hap_type in hap_type_array) {
                call tasks.extractFastq as extractFastq_scrna {
                    input:
                        sample_name = sample.sample_name,
                        SNPbam = SNPsplit_scrna.bam_prefix + hap_type + ".bam",
                        hap_type = hap_type,
                        fastq1 = sample.fastq1,
                        fastq2 = sample.fastq2,
                        analysis_dir = analysis_dir,
                        pipefail_check = pipefail_check,
                }

                call tasks.spaceRnager_for_ASE as spaceRnager_for_ASE {
                    input:
                        sample_name = sample.sample_name,
                        hap_type = hap_type,
                        slide = sample.slide,
                        area = sample.area,
                        image = sample.image,
                        slidefile = sample.slidefile,
                        transcriptome = transcriptome,
                        fastq1 = extractFastq_scrna.SNPsplit_fq1,
                        fastq2 = extractFastq_scrna.SNPsplit_fq2,
                        analysis_dir = analysis_dir,
                        pipefail_check = pipefail_check,
                        threads = 60,
                        spaceranger = spaceranger,
                }
            }

            call tasks.spaceRnager_for_ASE as spaceRnager_for_full {
                input:
                    sample_name = sample.sample_name,
                    hap_type = "full",
                    slide = sample.slide,
                    area = sample.area,
                    image = sample.image,
                    slidefile = sample.slidefile,
                    transcriptome = transcriptome,
                    fastq1 = sample.fastq1,
                    fastq2 = sample.fastq2,
                    analysis_dir = analysis_dir,
                    pipefail_check = pipefail_check,
                    threads = 60,
                    spaceranger = spaceranger,
            }
        }

        if (run_SNPsplit_rna) {
            # init SampleInfo
            SampleInfo_RNA sample_rna = {
                "sample_name": row[0],
                "fastq1": row[1],
                "fastq2": row[2],
            }

            call tasks.fastp as fastp {
                input:
                    sample_name = sample_rna.sample_name,
                    fastq1 = sample_rna.fastq1,
                    fastq2 = sample_rna.fastq2,
                    analysis_dir = analysis_dir,
                    pipefail_check = pipefail_check,
            }

            call tasks.STARmapping as STARmapping_rna {
                input:
                    sample_name = sample_rna.sample_name,
                    nmasked_star_index = nmasked_star_index,
                    fastq1 = sample_rna.fastq1,
                    fastq2 = sample_rna.fastq2,
                    analysis_dir = analysis_dir,
                    pipefail_check = pipefail_check,
                    threads = 40,
                    STAR_path = STAR_path,
            }

            call tasks.BamProcess as BamProcess_rna {
                input:
                    sample_name = sample_rna.sample_name,
                    mapping_res = STARmapping_rna.star_bam,
                    filter_expression = "[NH] == 1 and proper_pair",
                    result_dir = analysis_dir + "/STARmapping/",
                    pipefail_check = pipefail_check,
                    threads = 4,
            }

            call tasks.SNPsplit as SNPsplit_rna {
                input:
                    sample_name = sample_rna.sample_name,
                    bamfile = BamProcess_rna.bam,
                    seqtype = "--paired",
                    snp_file = snp_file,
                    analysis_dir = analysis_dir,
                    pipefail_check = pipefail_check,
                    SNPsplit_path = SNPsplit_path,
            }

            scatter (hap_type in hap_type_array) {
                call tasks.featureCounts_for_ASE as featureCounts_for_ASE {
                    input:
                        sample_name = sample_rna.sample_name,
                        hap_type = hap_type,
                        cleanbam = SNPsplit_rna.bam_prefix + hap_type + ".bam",
                        gtf = gtf,
                        analysis_dir = analysis_dir,
                        pipefail_check = pipefail_check,
                        subdir = hap_type,
                        long_reads = false,
                }
            }
        }

        if (run_SNPsplit_longReads) {
            # init SampleInfo
            SampleInfo_longReads sample_longreads = {
                "sample_name": row[0],
                "fastq1": row[1],
            }

            call tasks.trio_binning_classify as trio_binning_classify {
                input:
                    long_reads_fq = sample_longreads.fastq1,
                    hapA_only_kmers = hapA_only_kmers,
                    hapB_only_kmers = hapB_only_kmers,
                    hapA_prefix = hapA_prefix,
                    hapB_prefix = hapB_prefix,
                    sample_name = sample_longreads.sample_name,
                    analysis_dir = analysis_dir,
                    pipefail_check = pipefail_check,
            }

            scatter (hap_type in hap_type_array_long) {
                call tasks.minimap2 as minimap2_pacbio {
                    input:
                        sample_name = sample_longreads.sample_name,
                        long_reads_fq = trio_binning_classify.prefix + hap_type + ".fastq.gz",
                        reference_fa = reference_fa,
                        minimap2_path = minimap2_path,
                        analysis_dir = analysis_dir,
                        pipefail_check = pipefail_check,
                        subdir = hap_type,
                }

                call tasks.sam2bam as sam2bam {
                    input:
                        sam = minimap2_pacbio.sam,
                        sample_name = sample_longreads.sample_name,
                        analysis_dir = analysis_dir,
                        pipefail_check = pipefail_check,
                        subdir = hap_type,
                }

                call tasks.featureCounts_for_ASE as featureCounts_for_long_reads {
                    input:
                        sample_name = sample_longreads.sample_name,
                        hap_type = hap_type,
                        cleanbam = sam2bam.bam_2308,
                        gtf = gtf,
                        subdir = hap_type,
                        long_reads = true,
                        analysis_dir = analysis_dir,
                        pipefail_check = pipefail_check,
                }
            }
        }
    }
}
