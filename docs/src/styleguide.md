# Style guide
The following lists a few coding conventions for Trixi:

  * Modules, types, structs with `CamelCase`.
  * Functions, variables with lowercase `snake_case`.
  * Indentation with 2 spaces (*never* tabs!), line continuations indented with 4 spaces.
  * Maximum line length (strictly): **100**.
  * Functions that mutate their *input* are named with a trailing `!`.
  * Functions order their parameters [similar to Julia Base](https://docs.julialang.org/en/v1/manual/style-guide/#Write-functions-with-argument-ordering-similar-to-Julia-Base-1).
  * Prefer `for i in ...` to `for i = ...` for better semantic clarity and greater flexibility.
  * Executable code should only use ASCII characters.
  * Docstrings and comments can and should use Unicode characters where it helps understanding.
  * Multiline expressions should be explicitly grouped by parentheses and not
    rely on Julia's implicit line continuation syntax.
  * When naming multiple functions of a single or similar category, prefer to put the
    *general classification* first and the *specialization* second. Example: Use `flux_central`
    instead of `central_flux`. This helps when searching for available functions on the REPL
    (e.g., when trying to find all flux functions).

Based on that, and personal experience, a formatting tool with a few helpful
options is included in `utils/julia-format.jl`. Note, however, that this tool is
not yet optimal, as it re-indents too greedily.

This is a list of handy style guides that are mostly consistent with each
other and this guide, and which have been used as a basis:

  * [https://www.juliaopt.org/JuMP.jl/stable/style/](https://www.juliaopt.org/JuMP.jl/stable/style/)
  * [https://github.com/jrevels/YASGuide](https://github.com/jrevels/YASGuide)

