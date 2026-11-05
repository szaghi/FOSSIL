---
title: Contributing
---

# Contributing

This is a FOSS project — anyone interested in using, developing, or contributing
is welcome. The project follows a KISS (Keep It Simple and Stupid) philosophy.

## Reporting Issues

- Open a ticket on the repository **GitHub Issues** page
- Clearly describe the problem, including steps to reproduce for bugs
- Note the earliest version you know has the issue

## Pull Requests

1. Fork the repository on GitHub
2. Create a topic branch from `master`:
   ```bash
   git checkout -b fix/master/my_contribution master
   ```
3. Test your changes with `fobis build --mode tests-gnu && bash scripts/run_tests.sh` (the regression suite, see
   [Tests](#tests)); if they touch a heavy path (benches, large meshes, sweeps), run the validation suite too
4. Check for unnecessary whitespace: `git diff --check`
5. Submit a pull request with a clear commit message

## Tests

The tests are in two tiers (issue #45):

| Tier | Where | Built by | Run by | When |
|---|---|---|---|---|
| regression | `src/tests/*.f90` | `fobis build --mode tests-gnu` (CI: `tests-gnu-debug --coverage`) | `scripts/run_tests.sh` | every push |
| validation | `src/tests/validation/*.f90` | `fobis build --mode validation-gnu` (-O2, `exe-validation/`) | `scripts/run_validation.sh` | nightly (`.github/workflows/validation.yml`) and before a release |

A regression test must run in **2 s or less** in the debug coverage build (the CI build: `-O0 -fcheck=all`,
coverage instrumented; `fossil_test_boolean_coplanar` is the one approved exception, 3 s). Keep it fast by the choice
of input, never by a looser tolerance: `dragon.stl` (dirty, 0.14 s to load in debug) for self-consistency checks, a
clean closed body (`cube-0N.stl`, `make_apex_body` from `fossil_sign_bench_bodies`, a marching-cubes sphere) where
the sign must be exact; not `bunny.stl` (1.4 s to load in debug since issue #57, 19 s before; its tests are in
validation). A test that needs a large mesh, a fine grid or many cases goes to `src/tests/validation/`, often as a
full copy (`fossil_val_*`) of a regression test that keeps a smaller case. The device bench
`fossil_test_distance_bench` lives there too (built by `tests-nvidia`).

The budget is enforced: `fobis rule --ex makecoverage`, which CI runs, times every test after
`scripts/run_tests.sh` with `scripts/check_test_budgets.sh` and fails if one is over its budget (2 s, or its line in
`src/tests/budgets.txt`; times 2 on CI, whose runners are slower; a test over it is timed again and the faster run
counts). It always prints the five slowest tests. To check
a test of your own in that build:

```bash
fobis build --mode tests-gnu-debug --coverage && bash scripts/check_test_budgets.sh
```

An exception in `budgets.txt` needs a reason: a test that cannot get under 2 s with a smaller input belongs in
`src/tests/validation/`.

A test prints `Are all tests passed? T` and exits 0 on success; a `fossil_xfail_*` program must fail.

## Fortran Coding Style

- **Clarity over brevity**: `real :: gas_ideal_air` is better than `real :: gia`
- Single-character variable names only for loop counters
- Name all constants
- `implicit none` in every module and program
- Declare `intent` for all procedure arguments, ordered: pass arg → `inout` → `in` → `out` → optional
- Indent with two spaces (not tabs)
- No trailing whitespace; blank lines must contain no spaces
- Use `>, <, ==` instead of `.gt., .lt., .eq.`
- Avoid Windows-style CRLF line endings

### Recommended git whitespace settings

```ini
[color]
  ui = true
[color "diff"]
  whitespace = red reverse
[core]
  whitespace = fix,-indent-with-non-tab,trailing-space,cr-at-eol
```

## Commit style

Use [Conventional Commits](https://www.conventionalcommits.org/) so that `CHANGELOG.md` is generated automatically from the git log:

| Prefix | Purpose | Changelog section |
|--------|---------|-------------------|
| `feat:` | New feature or capability | New features |
| `fix:` | Bug fix | Bug fixes |
| `perf:` | Performance improvement | Performance |
| `refactor:` | Code restructuring | Refactoring |
| `docs:` | Documentation only | Documentation |
| `test:` | Tests | Testing |
| `build:` | Build system | Build system |
| `ci:` | CI/CD pipeline | CI/CD |
| `chore:` | Maintenance | Miscellaneous |

Append `!` for breaking changes (`feat!:`, `fix!:`). Reference issues with `#123` — they are auto-linked.

```
feat: add R32P kind parameter
fix: correct byte_size for character arrays (#42)
feat!: rename check_endian to init_endian
```

---

## Creating a release

Releases are fully automated via `scripts/release.sh` and GitHub Actions. The only steps needed are:

```bash
# Install git-cliff once
npx git-cliff@latest

# Then, to release:
scripts/release.sh patch   # v1.2.3 → v1.2.4
scripts/release.sh minor   # v1.2.3 → v1.3.0
scripts/release.sh major   # v1.2.3 → v2.0.0
scripts/release.sh v2.1.0  # explicit version
```

`release.sh` will ask for confirmation, then:

1. Regenerate `CHANGELOG.md` from the git log via [git-cliff](https://git-cliff.org/)
2. Update `VERSION`
3. Commit with `chore(release): vX.Y.Z`
4. Create an annotated git tag
5. Push commit + tag

Pushing the tag triggers the GitHub Actions release workflow, which automatically:
- Runs the regression suite and uploads coverage (run the validation suite before tagging)
- Builds this documentation site and deploys it to GitHub Pages
- Packages a versioned tarball
- Publishes a GitHub release with the changelog section as release notes
