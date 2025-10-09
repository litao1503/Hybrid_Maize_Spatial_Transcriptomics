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


# Trio-binning
task find_unique_kmer {
    input {
        String paternal_fq
        String maternal_fq

        String sample_name
        String analysis_dir

        String kmer_size = "21"
        String pipefail_check = "true"
        Int threads = 60
    }

    command <<<

        TASK_NAME="find_unique_kmer"

        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/unique_kmer_~{kmer_size}/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/unique_kmer_~{kmer_size}/
        fi

        cd ~{analysis_dir}/~{sample_name}/unique_kmer_~{kmer_size}/


        if [ ! -f find_unique_kmer.SUCCESS ]
        then
            /tools/trio_binning/bin/find-unique-kmers \
                -k ~{kmer_size} \
                -p ~{threads} \
                -o ~{analysis_dir}/~{sample_name}/unique_kmer_~{kmer_size}/ \
                ~{paternal_fq} ~{maternal_fq} \
                1>find_unique_kmer.stdout 2>find_unique_kmer.stderr \
                && touch find_unique_kmer.SUCCESS

            echo "paternal_fq ~{paternal_fq}" > unique_kmer.log
            echo "maternal_fq ~{maternal_fq}" >> unique_kmer.log
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "lt/asecount:1.0.0"
        cpu: 60
    }

    output {
        String return_code = "true"
    }
}


task trio_binning_classify {
    input {
        String long_reads_fq
        String hapA_only_kmers
        String hapB_only_kmers
        String hapA_prefix
        String hapB_prefix

        String sample_name
        String analysis_dir

        String pipefail_check = "true"
    }

    command <<<

        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/trio_binning_classify/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/trio_binning_classify/
        fi

        cd ~{analysis_dir}/~{sample_name}/trio_binning_classify/


        if [ ! -f trio_binning_classify.SUCCESS ]
        then
            python /tools/trio_binning/src/trio_binning/classify_by_kmers.py \
                ~{long_reads_fq} \
                ~{hapA_only_kmers} \
                ~{hapB_only_kmers} \
                --haplotype-a-out-prefix  ~{analysis_dir}/~{sample_name}/trio_binning_classify/~{hapA_prefix} \
                --haplotype-b-out-prefix ~{analysis_dir}/~{sample_name}/trio_binning_classify/~{hapB_prefix} \
                --unclassified-out-prefix ./unclassified > ./classify-by-kmers.log \
                && touch trio_binning_classify.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "lt/asecount:1.0.0"
        cpu: 2
    }

    output {
        String return_code = "true"
        String prefix =  "~{analysis_dir}/~{sample_name}/trio_binning_classify/"
    }
}


task minimap2 {
    input {
        String sample_name
        String long_reads_fq
        String reference_fa
        String minimap2_path
        String subdir = ""

        String analysis_dir
        String pipefail_check
    }

    command <<<
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/minimap2/~{subdir} ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/minimap2/~{subdir}
        fi

        cd ~{analysis_dir}/~{sample_name}/minimap2/~{subdir}


        if [ ! -f minimap2.SUCCESS ]
        then
           ~{minimap2_path} \
                -t 20 \
                -ax splice:hq -uf --eqx \
                ~{reference_fa} \
                ~{long_reads_fq} \
                -o ~{sample_name}.sam \
                1>minimap2.stdout 2>minimap2.stderr \
                && touch minimap2.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "lt/asecount:1.0.0"
        cpu: 20
    }

    output {
        String return_code = "true"
        String sam = "~{analysis_dir}/~{sample_name}/minimap2/~{subdir}/~{sample_name}.sam"
    }
}


task sam2bam {
    input {
        String sam
        String sample_name
        String subdir = ""

        String analysis_dir
        String pipefail_check
    }

    command <<<
         if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/minimap2/~{subdir} ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/minimap2/~{subdir}
        fi

        cd ~{analysis_dir}/~{sample_name}/minimap2/~{subdir}

        #-F 2308. Discard unmapped/supplementary_alignment/not_primary_alignment reads.
        if [ ! -f samtools.view.SUCCESS ]
        then
            samtools view -S -b -F 2308 -o ~{sample_name}.2308.bam ~{sam} \
            1>samtools.view.stdout 2>samtools.view.stderr \
            && touch samtools.view.SUCCESS
        fi

        if [ ! -f samtools.sort.SUCCESS ]
        then
            samtools sort -o ~{sample_name}.2308.sorted.bam ~{sample_name}.2308.bam\
            1>samtools.sort.stdout 2>samtools.sort.stderr \
            && touch samtools.sort.SUCCESS
        fi

        if [ ! -f samtools.index.SUCCESS ]
        then
            samtools index ~{sample_name}.2308.sorted.bam \
            1>samtools.index.stdout 2>samtools.index.stderr \
            && touch samtools.index.SUCCESS
        fi

        if [ -f samtools.index.SUCCESS ] && [ -f ~{sample_name}.2308.bam ]
        then
            rm ~{sample_name}.2308.bam
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi

    >>>

    runtime {
        docker: "lt/asecount:1.0.0"
        cpu: 1
    }

    output {
        String return_code = "true"
        String bam_2308 = "~{analysis_dir}/~{sample_name}/minimap2/~{subdir}/~{sample_name}.2308.sorted.bam"
    }
}


# SNPsplit
task STARmapping {
    input {
        String sample_name
        String nmasked_star_index
        String fastq1
        String fastq2 = ""
        String analysis_dir
        String STAR_path

        String pipefail_check = "true"
        Int threads = 40
    }

    command <<<
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/STARmapping/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/STARmapping/
        fi

        cd ~{analysis_dir}/~{sample_name}/STARmapping/

        if [[ "~{fastq2}" == "" ]]
        then
            readFilesIn=~{fastq1}
        else
            readFilesIn="~{fastq1} ~{fastq2}"
        fi

        if [ ! -f STARmapping.SUCCESS ]
        then
            ~{STAR_path} \
                --genomeDir ~{nmasked_star_index} \
                --runThreadN ~{threads} \
                --readFilesIn ${readFilesIn} \
                --readFilesCommand zcat \
                --outFileNamePrefix ~{sample_name}. \
                --alignEndsType EndToEnd \
                --outSAMattributes NH HI NM MD \
                --outSAMtype BAM Unsorted \
                --outBAMsortingThreadN 10 \
                1>STARmapping.stdout 2>STARmapping.stderr \
                && touch STARmapping.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "lt/asecount:1.0.0"
        cpu: threads
    }

    output {
        String return_code = "true"
        String star_bam = "~{analysis_dir}/~{sample_name}/STARmapping/~{sample_name}.Aligned.out.bam"
    }
}


#TODO not implemented
task SNPprepare {
    input {
        String genome1
        String genome2

        String pipefail_check = "true"
    }

    command <<<
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        echo ~{genome1}
        echo ~{genome2}
    >>>

    runtime {
        docker: "lt/asecount:1.0.0"
        cpu: 1
    }

    output {
        String return_code = "true"
    }
}


task BamProcess {
    input {
        String sample_name
        String mapping_res
        String filter_expression = "[NH] == 1 and proper_pair"

        String result_dir
        String pipefail_check
        Int threads = 4

    }

    command <<<
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{result_dir}/ ]]
        then
            mkdir -p ~{result_dir}/
        fi

        cd ~{result_dir}/

        if [ ! -f BamProcess.SUCCESS ]
        then
            if [[ ~{mapping_res} == *.sam ]]
            then
                sambamba view -h -t ~{threads} -S -f bam \
                -F "~{filter_expression}" ~{mapping_res} \
                | sambamba sort -o ~{sample_name}.sorted.bam /dev/stdin \
                && touch BamProcess.SUCCESS
            else
                sambamba view -h -t ~{threads} -f bam \
                -F "~{filter_expression}" ~{mapping_res} \
                | sambamba sort -o ~{sample_name}.sorted.bam /dev/stdin \
                && touch BamProcess.SUCCESS
            fi
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "lt/rnaseq:1.0.0"
        cpu: 4
    }

    output {
        String return_code = "true"
        String bam = "~{result_dir}/~{sample_name}.sorted.bam"
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
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/SNPsplit/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/SNPsplit/
        fi

        cd ~{analysis_dir}/~{sample_name}/SNPsplit/

        if [ ! -f SNPsplit.SUCCESS ]
        then
            ~{SNPsplit_path} \
                ~{seqtype} \
                --conflicting \
                --snp_file ~{snp_file} \
                ~{bamfile} \
                1>SNPsplit.stdout 2>SNPsplit.stderr \
                && touch SNPsplit.SUCCESS
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
        String bam_prefix = "~{analysis_dir}/~{sample_name}/SNPsplit/~{sample_name}.sorted."
        String paternal_bam = "~{analysis_dir}/~{sample_name}/SNPsplit/~{sample_name}.sorted.genome1.bam"
        String maternal_bam = "~{analysis_dir}/~{sample_name}/SNPsplit/~{sample_name}.sorted.genome2.bam"
        String conflicting_bam = "~{analysis_dir}/~{sample_name}/SNPsplit/~{sample_name}.sorted.conflicting.bam"
        String unassigned_bam = "~{analysis_dir}/~{sample_name}/SNPsplit/~{sample_name}.sorted.unassigned.bam"
    }
}


task extractFastq {
    input {
        String sample_name
        String SNPbam
        # genome1 genome2 unassigned
        String hap_type
        String fastq1
        String fastq2 = "NA"
        String analysis_dir

        String pipefail_check = "true"
    }

    command <<<
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/extractFastq/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/extractFastq/
        fi

        cd ~{analysis_dir}/~{sample_name}/extractFastq/

        if [ ! -f extractFastq.~{hap_type}.SUCCESS ]
        then
            # extract reads name
            samtools view ~{SNPbam} | awk '{print $1}' > read_names_~{hap_type}.txt

            if [ ~{fastq2} == "NA" ]
            then
                seqtk subseq ~{fastq1} read_names_~{hap_type}.txt | gzip > ~{sample_name}.~{hap_type}.SNPsplit.R1.fq.gz
            else
                seqtk subseq ~{fastq1} read_names_~{hap_type}.txt | gzip > ~{sample_name}.~{hap_type}.SNPsplit.R1.fq.gz
                seqtk subseq ~{fastq2} read_names_~{hap_type}.txt | gzip > ~{sample_name}.~{hap_type}.SNPsplit.R2.fq.gz
            fi

            touch extractFastq.~{hap_type}.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "lt/asecount:1.0.0"
        cpu: 2
    }

    output {
        String return_code = "true"
        String SNPsplit_fq1 = "~{analysis_dir}/~{sample_name}/extractFastq/~{sample_name}.~{hap_type}.SNPsplit.R1.fq.gz"
        String SNPsplit_fq2 = "~{analysis_dir}/~{sample_name}/extractFastq/~{sample_name}.~{hap_type}.SNPsplit.R2.fq.gz"
    }
}


task featureCounts_for_ASE {
    input {
        String sample_name
        # genome1 genome2 unassigned
        String hap_type
        String cleanbam
        String gtf
        String subdir = ""
        Boolean long_reads

        String analysis_dir
        String pipefail_check = "true"
    }

    command <<<
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/featureCounts_for_ASE/~{subdir} ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/featureCounts_for_ASE/~{subdir}
        fi

        cd ~{analysis_dir}/~{sample_name}/featureCounts_for_ASE/~{subdir}

        if [ ! -f featureCounts_for_ASE.~{hap_type}.SUCCESS ]
        then
            featureCounts \
                -a ~{gtf} \
                -R CORE -t exon -g gene_id \
                ~{true="-L" false="-p" long_reads} \
                -o ~{sample_name}.~{hap_type}.featureCounts \
                ~{cleanbam} \
                1>featureCounts.~{hap_type}.stdout \
                2>featureCounts.~{hap_type}.stderr \
                && touch featureCounts_for_ASE.~{hap_type}.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        docker: "lt/rnaseq:1.0.0"
        cpu: 2
    }

    output {
        String return_code = "true"
    }
}


task spaceRnager_for_ASE {
    input {
        String sample_name
        # genome1 genome2 unassigned
        String hap_type
        String slide = "V13L17-403"
        String area = "C1"
        String image
        String transcriptome
        String fastq1
        String fastq2
        String spaceranger
        String slidefile

        String analysis_dir
        String pipefail_check
        Int threads = 20
    }

    command <<<
        if [[ ~{pipefail_check} == "true" ]]
        then
            set -euo pipefail
        fi

        if [[ ! -d ~{analysis_dir}/~{sample_name}/spaceRnager_~{hap_type}/ ]]
        then
            mkdir -p ~{analysis_dir}/~{sample_name}/spaceRnager_~{hap_type}/input
        fi

        cd ~{analysis_dir}/~{sample_name}/spaceRnager_~{hap_type}/

        if [ ! -f spaceRnager_for_ASE.~{hap_type}.SUCCESS ]
        then
            ln -s ~{fastq1} ./input/~{sample_name}_S1_L001_R1_001.fastq.gz
            ln -s ~{fastq2} ./input/~{sample_name}_S1_L001_R2_001.fastq.gz

            ~{spaceranger} count \
                --id ~{sample_name}_~{hap_type} \
                --description ~{sample_name}_~{hap_type}_10xVisium \
                --slide ~{slide} \
                --area ~{area} \
                --image ~{image} \
                --slidefile ~{slidefile} \
                --transcriptome ~{transcriptome} \
                --fastqs ~{analysis_dir}/~{sample_name}/spaceRnager_~{hap_type}/input \
                --localcores ~{threads} \
                1>spaceRnager_for_ASE.~{hap_type}.stdout \
                2>spaceRnager_for_ASE.~{hap_type}.stderr \
                && touch spaceRnager_for_ASE.~{hap_type}.SUCCESS
        fi

        if [[ ~{pipefail_check} == "false" ]]
        then
            exit 0
        fi
    >>>

    runtime {
        cpu: threads
    }

    output {
        String return_code = "true"
    }
}

