version 1.0

task fastp {
    input {
        String sample_name
        String fastq1
        String fastq2

        String analysis_dir
        String pipefail_check
    }

    command <<<
        TASK_NAME="fastp"

        if [[ "~{pipefail_check}" == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/${TASK_NAME}/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/${TASK_NAME}/
        fi

        cd ~{analysis_dir}/~{sample_name}/${TASK_NAME}/


        if [ ! -f ${TASK_NAME}.SUCCESS ]
        then
            # Parse a comma-separated fq path string
            # and merge it into a file to pass to fastp.
            # Will choose to call "zcat" or "cat" depending on the
            # file name suffix before the first comma.
            if [[ "~{fastq1}" == *","* && "~{fastq2}" == *","* ]]
            then
                IFS=',' read -r -a fastq1_array <<< "~{fastq1}"
                IFS=',' read -r -a fastq2_array <<< "~{fastq2}"

                > auto.concat.log

                # An ugly implementation to avoid WDL syntax errors
                len1=0
                for i in "${fastq1_array[@]}"; do
                len1=$((len1 + 1))
                done

                echo ${len1} >> auto.concat.log

                len2=0
                for i in "${fastq2_array[@]}"; do
                len2=$((len2 + 1))
                done

                echo ${len2} >> auto.concat.log

                if [ "$len1" -eq "$len2" ]
                then
                    if [[ "${fastq1_array[0]}" == *.gz ]]
                    then
                        > concat.R1.fq.gz
                        > concat.R2.fq.gz
                        input_fq1=concat.R1.fq.gz
                        input_fq2=concat.R2.fq.gz
                    else
                        > concat.R1.fq
                        > concat.R2.fq
                        input_fq1=concat.R1.fq
                        input_fq2=concat.R2.fq
                    fi

                    echo "Start concat fastq1 list" >> auto.concat.log
                    for fq1 in "${fastq1_array[@]}"
                    do
                        echo ${fq1} >> auto.concat.log
                        cat ${fq1} >> ${input_fq1}
                    done
                    echo "Start concat fastq2 list" >> auto.concat.log
                    for fq2 in "${fastq2_array[@]}"
                    do
                        echo ${fq2} >> auto.concat.log
                        cat ${fq2} >> ${input_fq2}
                    done
                else
                    echo "[ERROR] Fastq1 and Fastq2 are not equal in length" >> auto.concat.log
                    exit 22
                fi
            else
                input_fq1=~{fastq1}
                input_fq2=~{fastq2}
            fi

            fastp \
                -i ${input_fq1} \
                -I ${input_fq2} \
                -o ~{sample_name}.fastp.1.fq.gz \
                -O ~{sample_name}.fastp.2.fq.gz \
                -h ~{sample_name}.html \
                -j ~{sample_name}.json \
                --thread 10 \
                1>${TASK_NAME}.stdout 2>${TASK_NAME}.stderr \
                && touch ${TASK_NAME}.SUCCESS

            if [[ -f concat.R1.fq ]]; then
                echo "remove concat.R1.fq" >> auto.concat.log
                rm concat.R1.fq
            fi
            if [[ -f concat.R2.fq ]]; then
                echo "remove concat.R2.fq" >> auto.concat.log
                rm concat.R2.fq
            fi

            if [[ -f concat.R1.fq.gz ]]; then
                echo "remove concat.R1.fq.gz" >> auto.concat.log
                rm concat.R1.fq.gz
            fi
            if [[ -f concat.R2.fq.gz ]]; then
                echo "remove concat.R2.fq.gz" >> auto.concat.log
                rm concat.R2.fq.gz
            fi
        fi

        if [[ "~{pipefail_check}" == "false" ]]
        then
            exit 0
        fi

    >>>

    runtime {
        docker: "lt/fastp:0.23.4"
        cpu: 10
    }

    output {
        String return_code = "true"
        String fastp_fq1 = "~{analysis_dir}/~{sample_name}/fastp/~{sample_name}.fastp.1.fq.gz"
        String fastp_fq2 = "~{analysis_dir}/~{sample_name}/fastp/~{sample_name}.fastp.2.fq.gz"
    }
}


task bismark_alignment {
    input {
        String clean_fq1
        String clean_fq2
        String ref_bismark_index
        Boolean non_directional = false

        String analysis_dir
        String sample_name
        String pipefail_check
    }

    command <<<
        TASK_NAME="bismark_alignment"

        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/${TASK_NAME}/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/${TASK_NAME}/
        fi

        cd ~{analysis_dir}/~{sample_name}/${TASK_NAME}/


        if [ ! -f ${TASK_NAME}.SUCCESS ]
        then
            bismark \
                --genome ~{ref_bismark_index} \
                -1 ~{clean_fq1} \
                -2 ~{clean_fq2} \
                -p 8 \
                ~{true="--non_directional" false="" non_directional} \
                --output_dir ./ \
                1>${TASK_NAME}.stdout 2>${TASK_NAME}.stderr \
                && touch ${TASK_NAME}.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "josousa/bismark:0.24.2"
        cpu: 40
    }

    output {
        String return_code = "true"
        String bismark_bam = "~{analysis_dir}/~{sample_name}/bismark_alignment/~{sample_name}.fastp.1_bismark_bt2_pe.bam"
    }
}


task deduplicated {
    input {
        String bismark_bam

        String analysis_dir
        String sample_name
        String pipefail_check
    }

    command <<<
        TASK_NAME="deduplicated"

        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/bismark_alignment/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/bismark_alignment/
        fi

        cd ~{analysis_dir}/~{sample_name}/bismark_alignment/


        if [ ! -f ${TASK_NAME}.SUCCESS ]
        then
            deduplicate_bismark \
                --outfile ~{sample_name} \
                --output_dir ./ \
                ~{bismark_bam} \
                1>${TASK_NAME}.stdout 2>${TASK_NAME}.stderr \
                && touch ${TASK_NAME}.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "josousa/bismark:0.24.2"
        cpu: 2
    }

    output {
        String return_code = "true"
        String bismark_deduplicated_bam = "~{analysis_dir}/~{sample_name}/bismark_alignment/~{sample_name}.deduplicated.bam"
    }
}


task SNPsplit {
    input {
        String sample_name
        String bamfile
        # --single_end or --paired
        String seqtype
        String snp_file
        String analysis_dir
        String SNPsplit_path

        String pipefail_check = "true"
    }

    command <<<
        TASK_NAME="SNPsplit_bisulfit"

        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/${TASK_NAME}/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/${TASK_NAME}/
        fi

        cd ~{analysis_dir}/~{sample_name}/${TASK_NAME}/


        if [ ! -f ${TASK_NAME}.SUCCESS ]
        then
             ~{SNPsplit_path} \
                ~{seqtype} \
                --bisulfit \
                --conflicting \
                --snp_file ~{snp_file} \
                -o ./ \
                ~{bamfile} \
                1>${TASK_NAME}.stdout 2>${TASK_NAME}.stderr \
                && touch ${TASK_NAME}.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        cpu: 2
    }

    output {
        String return_code = "true"
        String bam_prefix = "~{analysis_dir}/~{sample_name}/SNPsplit_bisulfit/~{sample_name}.deduplicated."
        String paternal_bam = "~{analysis_dir}/~{sample_name}/SNPsplit_bisulfit/~{sample_name}.deduplicated.genome1.bam"
        String maternal_bam = "~{analysis_dir}/~{sample_name}/SNPsplit_bisulfit/~{sample_name}.deduplicated.genome2.bam"
        String conflicting_bam = "~{analysis_dir}/~{sample_name}/SNPsplit_bisulfit/~{sample_name}.deduplicated.conflicting.bam"
        String unassigned_bam = "~{analysis_dir}/~{sample_name}/SNPsplit_bisulfit/~{sample_name}.deduplicated.unassigned.bam"
    }
}


task bismark_methylation_extractor {
    input {
        String bam
        String ref_bismark_index
        String sub_dir = ""

        String analysis_dir
        String sample_name
        String pipefail_check
    }

    command <<<
        TASK_NAME="bismark_methylation_extractor"

        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/${TASK_NAME}/~{sub_dir} ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/${TASK_NAME}/~{sub_dir}
        fi

        cd ~{analysis_dir}/~{sample_name}/${TASK_NAME}/~{sub_dir}


        if [ ! -f ${TASK_NAME}.SUCCESS ]
        then
            bismark_methylation_extractor \
                --comprehensive \
                --output_dir ./ \
                --cytosine_report \
                --CX_context \
                --genome_folder ~{ref_bismark_index} \
                --split_by_chromosome \
                ~{bam} \
                1>${TASK_NAME}.stdout 2>${TASK_NAME}.stderr \
                && touch ${TASK_NAME}.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "josousa/bismark:0.24.2"
        cpu: 2
    }

    output {
        String return_code = "true"
    }
}