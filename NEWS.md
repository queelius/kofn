# kofn 0.2.0

* Generalized masking: `loglik_masked()` and `rdata_masked()` extend masked
  cause-of-failure from series (k=1) to arbitrary k-out-of-m systems.

* Periodic inspection (Scheme 1): `loglik_scheme1()` and `fit_scheme1()`
  support general k via the critical-state density engine.

* Performance: eagerly precomputed critical/functioning states in
  `coherent_system()` and vectorized `f_sys_general()` give ~70x speedup
  for Scheme 1 likelihood evaluation.

* Exported `multistart_mle()` for users building custom likelihoods.

* New vignettes: symmetry-breaking comparison, periodic inspection,
  observation schemes, general k-out-of-m.

# kofn 0.1.0

* Initial release with exponential and Weibull parallel system MLE.
* EM algorithm for Weibull parallel systems.
* Coherent system infrastructure via minimal path/cut sets.
* Fisher information comparison across observation schemes.
