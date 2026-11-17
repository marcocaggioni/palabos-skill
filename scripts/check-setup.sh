#!/usr/bin/env bash
# Check the machine is ready to run the Palabos Couette example.
# Exits 0 when ready, 1 otherwise. No Palabos, no GPU, no license
# needed to pass the Copilot/skill parts — only the solver parts fail.
set -u
fail=0
say() { printf '%-28s %s\n' "$1" "$2"; }

check() {  # name, command...
  local name="$1"; shift
  if command -v "$1" >/dev/null 2>&1; then say "$name" "OK ($1)"; else say "$name" "MISS"; fail=1; fi
}

check "copilot CLI" copilot
check "uv" uv
check "cmake" cmake
check "C++ compiler" g++
check "make" make
check "MPI" mpirun
check "python3" python3
python3 -c "import yaml, matplotlib" 2>/dev/null && say "python yaml+matplotlib" "OK" || { say "python yaml+matplotlib" "MISS (pip install pyyaml matplotlib)"; fail=1; }

if [ -n "${PALABOS_ROOT:-}" ] && [ -d "$PALABOS_ROOT" ]; then
  say "Palabos (PALABOS_ROOT)" "OK ($PALABOS_ROOT)"
elif [ -d "$HOME/palabos" ]; then
  say "Palabos (~/palabos)" "OK"
else
  say "Palabos checkout" "MISS (git clone https://gitlab.com/unigespc/palabos.git ~/palabos)"
  fail=1
fi

[ -f skills/palabos-simulation/SKILL.md ] && say "skill (skills/)" "OK" || { say "skill (skills/)" "MISS"; fail=1; }
[ -f .github/skills/palabos-simulation/SKILL.md ] && say "skill (.github/skills/)" "OK" || { say "skill (.github/skills/)" "MISS"; fail=1; }
[ -f examples/couette/case.yaml ] && say "Couette case.yaml" "OK" || { say "Couette case.yaml" "MISS"; fail=1; }

if [ "$fail" -eq 0 ]; then echo "READY"; else echo "NOT READY — fix the MISS lines above"; fi
exit "$fail"
