version 1.0

import "wgbs_tasks.wdl" as tasks

struct SampleInfo {
    String sample_name
    String fastq1
    String fastq2
}


workflow WGBS2SNPsplit {
    input {
        String sample_sheet
        String analysis_dir

        String nmasked_bismark_index
        String snp_file
        String SNPsplit_path

        String pipefail_check = "true"
    }

    Array[String] hap_type_array = ["genome1", "genome2"]

    Array[Array[String]] sample_sheet_df = read_tsv(sample_sheet)
    scatter (row in sample_sheet_df) {
        SampleInfo sample = {
            "sample_name": row[0],
            "fastq1": row[1],
            "fastq2": row[2],
        }

        call tasks.fastp as fastp {
            input:
                sample_name = sample.sample_name,
                fastq1 = sample.fastq1,
                fastq2 = sample.fastq2,
                analysis_dir = analysis_dir,
                pipefail_check = pipefail_check,
        }

        call tasks.bismark_alignment as bismark_alignment {
            input:
                clean_fq1 = fastp.fastp_fq1,
                clean_fq2 = fastp.fastp_fq2,
                ref_bismark_index = nmasked_bismark_index,
                non_directional = false,
                analysis_dir = analysis_dir,
                sample_name = sample.sample_name,
                pipefail_check = pipefail_check,
        }

        call tasks.deduplicated as bismark_deduplicated {
            input:
                bismark_bam = bismark_alignment.bismark_bam,
                analysis_dir = analysis_dir,
                sample_name = sample.sample_name,
                pipefail_check = pipefail_check,
        }

        call tasks.SNPsplit as SNPsplit_bisulfit {
            input:
                sample_name = sample.sample_name,
                bamfile = bismark_deduplicated.bismark_deduplicated_bam,
                seqtype = "--paired",
                snp_file = snp_file,
                analysis_dir = analysis_dir,
                SNPsplit_path = SNPsplit_path,
                pipefail_check = pipefail_check,
        }

        scatter (hap_type in hap_type_array) {
            call tasks.bismark_methylation_extractor as bismark_methylation_extractor {
                input:
                    bam = SNPsplit_bisulfit.bam_prefix + hap_type + ".bam",
                    ref_bismark_index = nmasked_bismark_index,
                    sub_dir = hap_type,
                    analysis_dir = analysis_dir,
                    sample_name = sample.sample_name,
                    pipefail_check = pipefail_check,
            }
        }

        call tasks.bismark_methylation_extractor as bismark_methylation_extractor2 {
            input:
                bam = bismark_deduplicated.bismark_deduplicated_bam,
                ref_bismark_index = nmasked_bismark_index,
                sub_dir = "all",
                analysis_dir = analysis_dir,
                sample_name = sample.sample_name,
                pipefail_check = pipefail_check,
        }
    }
}

