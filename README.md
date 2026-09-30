# Hybrid Maize Spatial Transcriptomics

Analysis code and workflow resources associated with the manuscript:

> **Spatial and single-nucleus transcriptomics reveal the complexity of genomic imprinting in maize**

Tao Li<sup>1</sup>, Yi Jiang<sup>1</sup>, Mingyue Zhang<sup>2</sup>, Kesen Zhu<sup>1</sup>,
Siqi Jiang<sup>1</sup>, Xuerong Yang<sup>2</sup>, Hongjun Liu<sup>2,3</sup>,
Jiechen Wang<sup>4</sup>, Xiaomei Dong<sup>5</sup>, and Junpeng Shi<sup>1,#</sup>

## Overview

This repository contains custom workflows used to investigate genomic imprinting in
hybrid maize using spatial, single-nucleus, bulk, and long-read transcriptomic data.
The principal workflow, **spaseASE**, resolves allele-specific expression (ASE) by
combining SNP-based read assignment with transcriptomic mapping and quantification.
The workflows are implemented in Workflow Description Language (WDL) 1.0.

<p align="center">
  <img src="./spaseASE.png" alt="Overview of the spaseASE allele-specific expression workflow" width="50%">
</p>


## Repository structure

| Path | Description |
| --- | --- |
| `spaseASE/spaseASE.wdl` | Main spaseASE workflow for spatial, bulk, and long-read transcriptomic data. |
| `spaseASE/spaseASE_tasks.wdl` | WDL task definitions used by the spaseASE workflow. |
| `scripts/wgbs2snpsplit.wdl` | Workflow for allele-specific processing of whole-genome bisulfite sequencing data. |
| `scripts/wgbs_tasks.wdl` | WDL task definitions used by the WGBS workflow. |
| `docker/Dockerfile.spaseASE` | Dockerfile describing the principal spaseASE software environment. |
| `data/` | Compact study-specific reference or annotation files, where provided. |

## Workflow components

The `spaseASE` workflow supports three analysis modes. Each mode is enabled by its
corresponding Boolean input.

| Analysis mode | Workflow input | Main processing steps |
| --- | --- | --- |
| Spatial or single-nucleus RNA sequencing | `run_SNPsplit_scrna` | STAR alignment, BAM filtering, SNPsplit read assignment, FASTQ extraction, and Space Ranger analysis. |
| Bulk RNA sequencing | `run_SNPsplit_rna` | Read preprocessing, STAR alignment, BAM filtering, SNPsplit read assignment, and featureCounts quantification. |
| Long-read RNA sequencing | `run_SNPsplit_longReads` | Trio-binning classification, minimap2 alignment, BAM processing, and featureCounts quantification. |

The WGBS workflow is provided separately in `scripts/wgbs2snpsplit.wdl`.

## Requirements

- A WDL 1.0-compatible execution engine, such as Cromwell or miniwdl.
- A container runtime supported by the selected WDL execution engine.
- Access to the reference genome, annotation files, indexes, and sequencing data
  required for the selected analysis mode.
- Access to the container images specified in the WDL task runtime blocks.

The workflows currently expect input paths and `analysis_dir` to be visible inside
their execution containers. Shared storage and container path mappings should
therefore be configured for the local computing environment.


## Input sample sheets

Sample sheets are tab-separated files without a header row. The required columns
depend on the selected analysis mode.

### Spatial or single-nucleus RNA sequencing

| Column | Description |
| --- | --- |
| 1 | Sample name. |
| 2 | Read 1 FASTQ path. |
| 3 | Read 2 FASTQ path. |
| 4 | Histology image path. |
| 5 | Visium slide identifier. |
| 6 | Visium capture area. |
| 7 | CytAssist slide file path. |

### Bulk RNA sequencing

| Column | Description |
| --- | --- |
| 1 | Sample name. |
| 2 | Read 1 FASTQ path. |
| 3 | Read 2 FASTQ path. |

### Long-read RNA sequencing

| Column | Description |
| --- | --- |
| 1 | Sample name. |
| 2 | Long-read FASTQ path. |

## Running the workflows

Workflow inputs should be supplied in a JSON file using namespaced WDL input names.
For example:

```json
{
  "spaseASE.sample_sheet": "/path/to/sample_sheet.tsv",
  "spaseASE.analysis_dir": "/path/to/analysis",
  "spaseASE.run_SNPsplit_rna": true,
  "spaseASE.nmasked_star_index": "/path/to/star_index",
  "spaseASE.snp_file": "/path/to/snps.txt.gz",
  "spaseASE.STAR_path": "/path/to/STAR",
  "spaseASE.SNPsplit_path": "/path/to/SNPsplit",
  "spaseASE.gtf": "/path/to/annotation.gtf"
}
```

Run the main workflow with an appropriately configured WDL engine. For example:


```bash
java -jar cromwell.jar run spaseASE/spaseASE.wdl --inputs inputs.json
```

Results are written beneath `analysis_dir`, with task-specific subdirectories for
each sample.


## Citation

If you use the workflows or other materials in this repository, please cite the
associated manuscript:

Li, T., Jiang, Y., Zhang, M. *et al.* Spatial and single-nucleus transcriptomics reveal the complexity of genomic imprinting in maize. *Genome Biol* (2026). [https://doi.org/10.1186/s13059-026-04291-9](https://doi.org/10.1186/s13059-026-04291-9)


## License

The code in this repository is distributed under the [MIT License](LICENSE).
Third-party software, reference data, and containerized dependencies remain subject
to their respective licenses and terms of use.

## Contact

For questions about the repository or workflows, contact Tao Li at
[litao1503@outlook.com](mailto:litao1503@outlook.com) or open a GitHub issue.
