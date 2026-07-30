# APPM Variant Evidence Chain

The APPM and immune-escape layers preserve a stable event-level chain from a
called presentation-pathway variant to its platform evidence and clonality.

## Outputs

`appm/appm_variant_evidence.tsv` contains one row per called APPM gene variant.
Rows retain the event ID, genomic coordinates, consequence, damaging status,
source VAF/depth/ALT count, WGS and WES VAF/depth/ALT count, normal support,
cross-platform status, and CCF fields joined by `event_id`.

`appm/appm_gene_status.tsv` retains `source_event_ids`, variant counts, maximum
source/WGS/WES VAF, maximum CCF, and an explicit variant evidence status. Genes
without a called variant are marked `NO_CALLED_VARIANT`; this does not claim
that every base was callable.

`immune_escape/immune_escape_events.tsv` carries the source event IDs and the
corresponding maximum source/WGS/WES VAF and CCF for each material escape
mechanism. The event IDs link back to `appm_variant_evidence.tsv` for exact
variant-level review.

## Completeness states

- `COMPLETE_CROSS_PLATFORM`: WGS and WES VAF/depth plus CCF are available.
- `COMPLETE_SOURCE_PLATFORM`: source VAF/depth plus CCF are available, but both
  WGS and WES measurements are not present.
- `VAF_DEPTH_NO_CCF`: source/platform frequency and depth exist without CCF.
- `PARTIAL`: the event is retained but one or more core evidence layers are
  unavailable.

CCF 2.1 now runs before APPM in the Python and Nextflow production paths, and
its output is passed to APPM explicitly. Existing output columns remain
available for backward compatibility.
