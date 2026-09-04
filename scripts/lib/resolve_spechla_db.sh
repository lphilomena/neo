#!/usr/bin/env bash
# Portable SpecHLA database discovery.
# Do not hard-code host mounts. Accept the canonical upstream layout
# (data/hla/spechla/db) and the legacy alias (data/hla/spechla_db).

neoag_spechla_db_is_valid() {
  local dir="${1:-}"
  [[ -n "$dir" && -d "$dir" ]] || return 1
  [[ -s "$dir/ref/hla.ref.extend.fa" || -s "$dir/HLA/hla.ref.extend.fa" ]]
}

neoag_resolve_spechla_db() {
  local extra_root="${1:-}"
  local cand root rel
  local roots=()
  local rels=(
    "data/hla/spechla/db"
    "data/hla/spechla_db"
    "hla/spechla/db"
    "tools/SpecHLA/db"
  )

  if neoag_spechla_db_is_valid "${SPECHLA_DB:-}"; then
    (cd "$SPECHLA_DB" && pwd -P)
    return 0
  fi
  if neoag_spechla_db_is_valid "${SPECHLA_HOME:-}/db"; then
    (cd "$SPECHLA_HOME/db" && pwd -P)
    return 0
  fi

  for root in \
    "$extra_root" \
    "${NEOAG_REF_BUNDLE:-}" \
    "${NEOAG_REFERENCE_ROOT:-}" \
    "${OPEN_NEO_REFERENCE_ROOT:-}" \
    "${REFERENCE_ROOT:-}" \
    "${NEOAG_TOOLS_ROOT:-}" \
    "${OPEN_NEO_TOOLS_ROOT:-}" \
    "${NEOAG_PROJECT_ROOT:-}"; do
    [[ -n "$root" && -d "$root" ]] || continue
    roots+=("$(cd "$root" && pwd)")
  done

  for root in "${roots[@]}"; do
    for rel in "${rels[@]}"; do
      cand="$root/$rel"
      if neoag_spechla_db_is_valid "$cand"; then
        (cd "$cand" && pwd -P)
        return 0
      fi
    done
  done
  return 1
}
