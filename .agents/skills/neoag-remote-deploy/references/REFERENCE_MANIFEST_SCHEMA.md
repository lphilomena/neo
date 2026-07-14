# Reference manifest guidance

Recommended local path: `configs/local/reference_manifest.yaml`.

Example:

```yaml
genome_build: GRCh38
reference_fasta: /data/refs/GRCh38/GRCh38.fa
gencode_gtf: /data/refs/gencode/gencode.v44.annotation.gtf
vep_cache: /data/refs/vep_cache/homo_sapiens/115_GRCh38
normal_proteome: /data/refs/neoag/normal_proteome.fa
normal_ligandome: /data/refs/neoag/normal_ligandome.tsv
normal_junctions: /data/refs/neoag/normal_junctions.tsv
facets_snp_vcf: /data/refs/facets/common_snp.hg38.vcf.gz
hla_reference: /data/refs/hla/
```

Rules:

- Use target-machine paths only.
- Do not invent fake paths to silence Doctor.
- Do not commit local manifests with private paths.
