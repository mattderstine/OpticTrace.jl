# GitHub Copilot Instructions – Julia Project

## Language & Environment
- This project uses **Julia** version 1.12+ as the primary programming language
- Package manager: built-in `Pkg` with `Project.toml` and `Manifest.toml`
- Editor: VS Code with the Julia extension

## Code Style
- Follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- Use 4-space indentation
- Use `camelCase` for variables and functions, `PascalCase` for types and structs
- Prefer descriptive names over abbreviations
- Add docstrings to all public functions using Julia's `"""..."""` format
- Keep functions short and focused (single responsibility)

## Julia-Specific Conventions
- Prefer multiple dispatch over if/else type branching
- Avoid global mutable state; use function arguments and return values
- Write type-stable functions where possible — avoid `Any` unless necessary
- Use `@inbounds` and `@simd` only after profiling, not preemptively
- Prefer in-place functions (with `!` suffix) for performance-critical code
- Use broadcasting (`.`) instead of explicit loops where idiomatic

## Error Handling
- Use `try/catch` sparingly — prefer validation with informative error messages
- Throw descriptive exceptions: `throw(ArgumentError("message"))`
- Validate inputs at function boundaries

## Testing
- Tests live in `/test` with `Test` stdlib (`@test`, `@testset`)
- Every public function should have corresponding tests
- Use descriptive `@testset` names

## Performance
- Avoid type instability — use `@code_warntype` to check
- Prefer stack-allocated types (immutable structs) where possible
- Use `BenchmarkTools.jl` (`@benchmark`) for performance measurement
- Profile with `Profile` stdlib before optimizing

## Dependencies
- List key packages here, e.g.:
  - `DataFrames.jl` – tabular data
  - `Plots.jl` / `Makie.jl` – visualization
  - `Flux.jl` – machine learning
- Do not add new dependencies without updating `Project.toml`

## Project Structure
- `src/` – main source code, entry point is `src/MyProject.jl`
- `test/` – test files mirroring `src/` structure
- `scripts/` – one-off scripts and experiments (not part of the module)
- `docs/` – documentation (Documenter.jl)
- `data/` – raw and processed data (not committed to git if large)

## What to Avoid
- Do not use Python-style OOP patterns — use multiple dispatch instead
- Do not write type-unstable functions for hot paths
- Do not use `eval` or `@eval` unless absolutely necessary
- Do not ignore compiler warnings about type inference