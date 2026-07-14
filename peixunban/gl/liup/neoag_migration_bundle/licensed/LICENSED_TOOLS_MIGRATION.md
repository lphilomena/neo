# Licensed or restricted tools migration notes

These tools may be license-sensitive or restricted. Do not redistribute them publicly. Copy or mount them only for authorized users.

| Tool | Current staged location | Environment variable / manifest path | Doctor behavior |
|---|---|---|---|
| NetMHCpan | /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/predictors/netMHCpan | NETMHCPAN_HOME, NEOAG_NETMHCPAN_BIN, test_env_tools/bin/netMHCpan | Path exists; license reminder is non-blocking |
| NetMHCstabpan | /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/predictors/netMHCstabpan | NETMHCSTABPAN_HOME, test_env_tools/bin/NetMHCstabpan | Path exists; license reminder is non-blocking |
| PRIME | /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/predictors/prime | PRIME_HOME, NEOAG_PRIME_BIN, test_env_tools/bin/PRIME wrapper | Wrapper provides stable help and delegates real runs |
| MixMHCpred | /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/predictors/mixMHCpred_install | MIXMHCPRED_HOME, MIXMHCPRED_BIN, test_env_tools/bin/MixMHCpred | Path exists |
| Novoalign license | /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/lohhla/novoalign.lic | NOVOALIGN_LICENSE_FILE | Required by LOHHLA/Polysolver paths when used |
| Polysolver / LOHHLA | /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/lohhla/polysolver | POLYSOLVER_HOME, lohhla_reference manifest | Reference/tool path is declared |

Activation helper:

```bash
source /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/activate_neoag_production_refs.sh
```

This helper sources neodata4git variables, then resets NEOAG_TOOLS_ROOT back to the fresh test_env_tools root so it does not clone or depend on old conda environments.
