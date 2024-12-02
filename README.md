<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-marsseq_logo_dark.png">
    <img alt="nf-core/marsseq" src="docs/images/nf-core-marsseq_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/marsseq/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/marsseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/marsseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/marsseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/marsseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.8063539-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.8063539)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/marsseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23marsseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/marsseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/marsseq** is a bioinformatics single-cell preprocessing pipeline for MARS-seq v2.0 experiments. MARS-seq is a plate-based technique that can be combined with FACS in order to study rare populations of cells. On top of the pre-existing pipeline, we have developed an RNA velocity workflow that can be used to study cell dynamics using `StarSolo`. We do so by converting the raw FASTQ reads into 10X v2 format.

![Workflow](docs/images/workflow.png)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

To run the pipeline you have create experiment metadata files:

- [amp_batches.xlsx](assets/amp_batches.xlsx)
- [wells_cells.xlsx](assets/wells_cells.xlsx)
- [seq_batches.xslx](assets/seq_batches.xlsx)

and samplesheet (`samplesheet.csv`). We provide test example [here](assets/samplesheet.csv).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/marsseq \
  -profile <docker/singularity/.../institute> \
  --fasta https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz \
  --gtf https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/marsseq/usage) and the [parameter documentation](https://nf-co.re/marsseq/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/marsseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/marsseq/output).

## Credits

nf-core/marsseq was originally written by [Martin Proks](https://github.com/matq007).

We thank the following people for their extensive assistance in the development of this pipeline:

- Jose Alejandro Romero Herrera ([@joseale2310](https://github.com/joseale2310))
- Maxime Garcia ([@maxulysse](https://github.com/maxulysse))

Keren-Shaul, H., Kenigsberg, E., Jaitin, D.A. et al. MARS-seq2.0: an experimental and analytical pipeline for indexed sorting combined with single-cell RNA sequencing. Nat Protoc 14, 1841–1862 (2019). https://doi.org/10.1038/s41596-019-0164-4

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#marsseq` channel](https://nfcore.slack.com/channels/marsseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/marsseq for your analysis, please cite it using the following doi: [10.5281/zenodo.8063539](https://doi.org/10.5281/zenodo.8063539)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
