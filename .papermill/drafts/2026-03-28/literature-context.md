# Literature Context

## Key Prior Art and Positioning

### Identifiability Theory
- Nowik (1990): Formal identifiability conditions for coherent systems. Showed non-identifiability of individual parameters from system data. Our Proposition 1 provides the explicit permutation invariance proof for k-out-of-n.
- Tsiatis (1975): Classical non-identifiability in competing risks. Our work extends to general k.
- Cox (1959): Series system analysis foundations.

### Masked Cause-of-Failure (Series Systems)
- Miyakawa (1984): First masked data MLE for series systems.
- Lin/Usher/Guess (1993): Exact MLE with masked data.
- Usher (1996): Weibull extension.
- Craiu/Duchesne (2004): EM for masked competing risks.
Our generalized masking (Definition 2) extends this from k=1 to general k.

### k-out-of-n Specific
- Komarova (2017): Nonparametric identification from autopsy data. Our work is parametric but compares multiple mechanisms.
- Rodrigues/Pereira/Polpo (2019): Masked data for coherent systems. Closest prior work; we go further with the symmetry framing and quantitative comparison.
- Sarhan/El-Bassiouny (2003): Parallel system masked data (2 components). We handle general m and k.
- Dembinska/Jasinski (2021): Discrete component lifetimes in k-out-of-n. Different setting (discrete).

### System Signatures
- Samaniego (2007): System signatures for component inference.
- Bhattacharya/Samaniego (2010), Ng/Navarro/Balakrishnan (2017), Yang/Ng/Balakrishnan (2019): EM approaches using known/unknown signatures.

### Partial Autopsy Precedent
- Flehinger/Reiser/Yashchin (1998): Two-stage diagnosis model, closest to our partial autopsy. We generalize from their specific two-stage to arbitrary r inspections.
- Meilijson (1981): Established that autopsy (knowing which component failed) restores identifiability for series systems.

### Label Switching Connection
- Stephens (2000), Jasra/Holmes/Stephens (2005): Label switching in Bayesian mixtures. Same mathematical structure (likelihood invariant under permutation), different solution approach (post-processing vs. data augmentation).

### Composite Likelihood
- Varin/Reid/Firth (2011): Framework for our Scheme 1 analysis.

## Novelty Assessment
1. Symmetry-breaking framing: HIGH. No prior work frames k-out-of-n estimation explicitly as a symmetry problem.
2. Generalized masking for k>1: MODERATE-HIGH. Extension of well-known series framework.
3. Partial autopsy spectrum: HIGH. The r=1 phase transition result is new.
4. Quantitative comparison of 6 mechanisms: HIGH. No systematic comparison exists.
5. Delta insensitivity finding: HIGH. Unexpected and practically important.
