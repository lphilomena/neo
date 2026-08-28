# ASCAT WGS hg38 reference files

Source: ASCAT official WGS reference files from Zenodo record 14008443, linked by `VanLoo-lab/ascat/ReferenceFiles/WGS`.

Downloaded archives:

- `downloads/G1000_loci_WGS_hg38.zip`
- `downloads/G1000_alleles_WGS_hg38.zip`

Archive SHA256:

```text
dc554d727d509862bbd1b5b9cc558fcf76fedc1e3e8fbd8404dcd16aa17b93e7  G1000_loci_WGS_hg38.zip
fe2244969b0bbefb4e8dcb33820a0dd99cb2068e932e10a2046dac40b1535bfd  G1000_alleles_WGS_hg38.zip
```

Official extracted prefixes:

```r
loci.prefix <- "data/ascat/reference/WGS_hg38/files/G1000_loci_hg38_chr"
alleles.prefix <- "data/ascat/reference/WGS_hg38/files/G1000_alleles_hg38_chr"
chrom_names <- c(1:22, "X")
genomeVersion <- "hg38"
```

For BAMs whose headers use `chr1`, `chr2`, ..., `chrX`, use the chr-prefixed loci copy below. The coordinates are from the official loci files; only the chromosome labels were changed from `1` to `chr1` etc. The alleles files remain the official extracted files.

```r
loci.prefix <- "data/ascat/reference/WGS_hg38/files_chr_bam/G1000_loci_hg38_"
alleles.prefix <- "data/ascat/reference/WGS_hg38/files/G1000_alleles_hg38_"
chrom_names <- paste0("chr", c(1:22, "X"))
genomeVersion <- "hg38"
```

These files are external reference data and should not be bundled into the lightweight online release.

For ASCAT v3 `prepareHTS` with `chr`-style BAMs, the most compatible prefixes are:

```r
loci.prefix <- "data/ascat/reference/WGS_hg38/files_chr_bam_nochr_name/G1000_loci_hg38_"
alleles.prefix <- "data/ascat/reference/WGS_hg38/files_nochr_name/G1000_alleles_hg38_"
chrom_names <- as.character(c(1:22, "X"))
genomeVersion <- "hg38"
```

Here loci file contents use `chr1`, `chr2`, ... for `alleleCounter`, while file names and ASCAT chromosome names use `1`, `2`, ... so ASCAT row names match after it strips `chr` from alleleCounter output.
