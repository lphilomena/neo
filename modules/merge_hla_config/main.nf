/*
 * MERGE_HLA_CONFIG — Inject computed HLA alleles into a TOML run config.
 *
 * Takes the original run TOML + an hla_alleles.txt file (one allele per line)
 * and produces a merged TOML with hla_alleles set in [inputs].  This mirrors
 * the behaviour of `neoag-v03 run-full` which reads hla_alleles from the TOML.
 */
process MERGE_HLA_CONFIG {
  tag "$sample_id"
  label 'small'
  publishDir "${params.outdir}/config", mode: 'copy'

  input:
    val sample_id
    path run_config
    path hla_alleles_txt

  output:
    path "run_merged.toml", emit: merged_config

  script:
  """
  python3 -c "
import sys
from pathlib import Path

# Read original TOML
toml = Path('${run_config}').read_text()

# Read HLA alleles from OptiType output
hla_raw = Path('${hla_alleles_txt}').read_text().strip()
hla_list = [a.strip() for a in hla_raw.split(chr(10)) if a.strip()]

# Format as TOML array
hla_toml = 'hla_alleles = [' + ', '.join(f'\"{a}\"' for a in hla_list) + ']'

# Replace or insert
import re
if re.search(r'hla_alleles\s*=\s*\[', toml):
    toml = re.sub(r'hla_alleles\s*=\s*\[[^\]]*\]', hla_toml, toml)
elif '[inputs]' in toml:
    toml = toml.replace('[inputs]', '[inputs]\n' + hla_toml)
else:
    toml += '\n[inputs]\n' + hla_toml + '\n'

Path('run_merged.toml').write_text(toml)
print(f'Merged HLA alleles into config: {hla_list}', flush=True)
  "
  """
}
