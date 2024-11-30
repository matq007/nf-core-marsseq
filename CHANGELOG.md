# nf-core/marsseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### 🚀 Features

- Merge with template and full refactor
- Support for compressed genome references

### 🐛 Bug Fixes

- Linting issues
- Test config
- Check that provided xlsx files are not empty after subset
- Missing openpyxl in conda specification
- Missing openpyxl conda dependency
- Replace ubuntu image with conda packages
- Perl not anymore in bioconda channel
- Move uncompressed genome and annotation to reference folder
- Small and full AWS test

### 💼 Other

- Merge template 3.0.2 to dev

### 📚 Documentation

- Linting for CITATIONS

### ⚙️ Miscellaneous Tasks

- Merge pull request #22 from nf-core/nf-core-template-merge-3.0.2

## v1.0.3 - [11-07-2023]

### `Fixed`

- AWS testing

## v1.0.2 - [11-07-2023]

### `Dependencies`

- [#6](https://github.com/nf-core/marsseq/pull/6) sync with template 2.8 -> 2.9
- [#4](https://github.com/nf-core/marsseq/issues/4) - Fix AWS testing with s3 bucket

## v1.0.1 - [26-06-2023]

### `Fixed`

- [#4](https://github.com/nf-core/marsseq/issues/4) - Fix AWS testing with s3 bucket
- missing zenodo in `lib/`
- added credit and provided short MARS-seq description

## v1.0.0 - [21-06-2023]

### `Dependencies`

- Bump minimal Nextflow version to 23.04.0
- sync with template 2.8

## v1.0dev - [date]

Initial release of nf-core/marsseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
