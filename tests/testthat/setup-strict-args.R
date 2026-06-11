# Surface R partial-argument-matching collisions as warnings during test
# runs. Every model method in this ecosystem has a `par` formal, so a
# stray kwarg like `p = ...` would otherwise silently bind to `par` and
# corrupt the call. See the rlang monorepo CLAUDE.md, "R partial-argument-
# matching footgun".
options(warnPartialMatchArgs = TRUE,
        warnPartialMatchAttr = TRUE,
        warnPartialMatchDollar = TRUE)
