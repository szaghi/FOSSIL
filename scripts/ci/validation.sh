#!/usr/bin/env bash
# The validation suite (issue #45) on one toolchain, run nightly by .github/workflows/validation.yml through
# fortran-toolchains' matrix workflow (in the gnu-16.2 image: FC, TOOLCHAIN, ... of its environment contract), and as is
# locally or in `docker run <image>`. The heavy programs of src/tests/validation/, built -O2 into exe-validation/,
# run -j $(nproc) by scripts/run_validation.sh (their logs: exe-validation/*.log, the workflow's artifact).
set -euo pipefail
echo "toolchain: ${TOOLCHAIN:-local}, FC=${FC:-gfortran}"
if [ -f src/third_party/.deps_config.ini ]; then fobis fetch; fi
fobis build --mode validation-gnu --fc "${FC:-gfortran}" # the only validation mode: gfortran
./scripts/run_validation.sh -j "$(nproc)"
