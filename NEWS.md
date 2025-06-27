# methyldeconv 1.1.0

* Added unified signature matrix getter functions for all deconvolution methods (EpiDISH, Houseman, MethylCC, MethylResolver, MethAtlas), with consistent data frame output and documentation.
* All deconvolution methods now check for overlap between input and reference CpGs, issuing warnings if overlap is low.
* Extended and improved vignettes: new vignettes for signatures, visualization, and contributing; improved getting started guide; new navbar dropdown in documentation.
* Clarified expected input structure for visualization functions in documentation.
* Enhanced test suite: more robust tests for CpG overlap, output structure, and method extensions.
* Improved and reorganized documentation: new manual pages, improved roxygen2 docs, and reorganized documentation files.
* CI/CD improvements: manual installation of preprocessCore for compatibility; simplified getting started vignette for CI stability.
* Increased CpG count in overlap tests to avoid singular matrix errors.
* Various bugfixes and clarifications in vignettes and documentation.

# methyldeconv 0.0.3

* rebranded methylDeconv to methyldeconv

# methylDeconv 0.0.2

* First version of package.

# methylDeconv 0.0.1

* Added a `NEWS.md` file to track changes to the package.
