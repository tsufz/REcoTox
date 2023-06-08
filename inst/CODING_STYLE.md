# Coding style conventions

Author: Tobias Schulze

Date: 06 June 2023

The following coding style conventions apply to the `MZquant` programming environment.

## Importing other packages

-   No full import, only selected import of functions

-   Use base-R wherever possible or write small helperfunctions to decrease dependencies

-   New R-version issues should be addressed (ensured by development version checks in CI)

## Semantic package versioning

The package follows **Semantic Versioning**:

The version number MAJOR.MINOR.PATCH, it is defined:

-   MAJOR version when package is not backward compatible (change of API)
-   MINOR version when functionality with backward compatibility was added
-   PATCH version when committing backward compatible bug fixes or light issues

## Semantic commits and version triggers

To ensure reproducible commits, standardized and semantic commits are required. A future use will be the automated generation of the changelog out of the commit messages.

Following the ideas developed for [Python Packages](https://py-pkgs.org/07-releasing-versioning#automatic-version-bumping "Python Packages"), this schema should be used:

```         
<type>(optional <scope>): short summary in present tense

(optional <body>: explains motivation for the change)

(option <footer>: note BREAKING CHANGES here, and issues to be closed)
```

<type> refers to the kind of change made and is usually one of:

-   *feat*: A new feature for the user.
-   *fix*: A bug fix for the user.
-   *docs*: Changes to documentation.
-   *style*: Changes not affecting the meaning of the code (white-space, formatting, missing semi-colons, etc).
-   *refactor*: A code change neither fixing a bug nor adds a feature.
-   *perf*: A code change improving the performance.
-   *test*: Changes in the test framework.
-   *build*: Changes in the build process or tools (including CI/CD)

<scope> is an optional keyword that provides context for where the change was made. It can be anything relevant to your package or development workflow (e.g., it could be the module or function name affected by the change).

***Some examples***

A <type> of *fix* with <scope> *modelling* triggers a **patch** version bump:

```         
$ git commit -m "fix(modelling): fix confusing error message"
```

A <type> of *feat* with <scope> *package* triggers a **minor** version bump:

```         
$ git commit -m "feat(package): add example data to package"
```

A <type> of *feat* with <scope> *system_settings* triggers a **BREAKING CHANGE** and **major** version bump:

```         
$ git commit -m "feat(system_setting): use yaml style settings replacing json style"
$
$ BREAKING CHANGE: json style setting files will not work anymore
```

## Helper-functions

-   At least the exported functions must be documented. Use `roclets` in front of each function.

## Function names

-   All function names must be written in `lowercase` and separated by an `underscore` (e.g. `plot_function`).

## General coding

-   External functions are called by `tibble::tibble()` to avoid problems with masked functions.
-   `Tidyverse` style coding is recommended.
-   Use semantic `variable` names and not abbreviation, register variables in `zzz.R` in `global_Variables`
-   Write important settings and results in the `.MZquant.env` for provenance.
