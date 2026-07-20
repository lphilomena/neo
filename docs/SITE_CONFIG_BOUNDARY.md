# Site Config And Private Example Boundary

The lightweight release keeps portable examples in Git and keeps real deployment
state outside Git.

## Commit These

- `conf/site.config.example`
- `conf/run.private.example.toml`
- `conf/tools.env.local.example.sh`
- Small fixture configs that use relative paths under `data/fixtures*`

## Keep These Local

- `conf/site.config`
- `conf/private/*.toml`
- `conf/*.private.toml`
- `conf/*.local.toml`
- Real tool roots, cluster queues, reference bundle paths, patient/sample
  identifiers, and controlled-access data locations

## Why

Site files often contain absolute paths, licensed tool locations, sample
identifiers, and private storage conventions. Keeping them outside Git makes the
release portable and avoids leaking local deployment details.

## Recommended Pattern

1. Copy `conf/site.config.example` to `conf/site.config`.
2. Copy `conf/run.private.example.toml` to `conf/private/<project>.toml`.
3. Replace placeholder paths with local paths.
4. Run `neoag check-tools` before production runs.
