#!/usr/bin/env bash
# Fixed entrypoint for running OpticTrace.jl's test suite, so the exact
# command stays constant across runs (see .claude/settings.json's allowlist,
# which permits this script without a prompt). Full output goes to a temp
# log for follow-up grep'ing; only the tail is printed here.
cd "$(dirname "${BASH_SOURCE[0]}")"

log="$(mktemp -t optictrace-test).log"
julia --project=. -e 'using Pkg; Pkg.test()' > "$log" 2>&1
status=$?

echo "EXIT=$status"
echo "full log: $log"
tail -n 150 "$log"

exit $status
