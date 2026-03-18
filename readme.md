# Guruswami-Sudan List Decoder Implementation in SageMath

This repository provides a complete implementation and experimental analysis of the Guruswami-Sudan (GS) list decoding algorithm for Reed-Solomon (RS) codes. The implementation is written in SageMath to leverage its robust algebraic capabilities over finite fields. 

Traditional unique decoding algorithms, such as Berlekamp-Massey, are strictly limited by the unique decoding radius of $(n-k)/2$. The Guruswami-Sudan algorithm effectively recovers messages from corrupted codewords even when the noise level exceeds this bound, theoretically approaching the Johnson bound.

## Environment Requirements
This project **must be run in a Linux environment**. 

### Prerequisites
* **OS:** Linux (Ubuntu, Debian, Fedora, etc. or WSL on Windows)
* **Software:** [SageMath](https://www.sagemath.org/) (Version 10.7 or compatible) 

## Features

* **Bivariate Polynomial Interpolation:** Constructs a curve $Q(x,y)$ passing through received points with a configurable multiplicity parameter $m$.
* **Root Finding & Factorization:** Extracts valid message polynomials by factoring the interpolated bivariate polynomial into the form $y-P(x)$.
* **Automated Experiment Suite:** Includes a complete testing environment (`GSExperimentSuite`) to simulate noisy channels, evaluate runtime, map decoding radii, and test adversarial signal collisions.

## Repository Structure
* **`RS.sage`**: The core source code containing the `GSDecoder` implementation and the `GSExperimentSuite`.
* **`experiment1_logs.txt`**: Generated log outputs from the execution of the first experiment.
* **`exp1_radius_comparison.png` / `exp1_radius.png`**: Visualizations of the success probability against the number of errors across varying multiplicities $m$.
* **`exp2_runtime.png`**: A chart illustrating the growth in runtime cost as the multiplicity parameter increases.
* **`exp3_listsize.png`**: A histogram showing the distribution of output list sizes under random noise models.

## Experimental Insights

This implementation was  tested over the finite field $GF(2^8)$ with block length $n=64$. Key findings from the research include:

1. **Decoding Radius vs. Multiplicity (Exp 1):** Increasing $m$ monotonically increases the error-correction capability past the unique decoding bound. For a code of $RS(64,16)$, $m=1$ corrects 27 errors, and $m=4$ corrects 31 errors (unique bound is 24). 
2. **Runtime Complexity (Exp 2):** The interpolation matrix growth introduces a steep computational cost. Runtime scales non-linearly, approaching $O(m^6)$ , making smaller multiplicities ($m=2$ or $m=3$) more practical for real-time applications.
3. **List Size Behavior (Exp 3 & 4):** Despite being a "list" decoder, random symmetric noise almost universally produces a list size of exactly 1.
4. **Rate Dependency (Exp 4):** The algorithm shows massive error-correction gains for low-rate codes ($R \le 0.25$) but offers negligible improvements over unique decoders for high-rate codes.
5. **Signal Collision Validation (Exp 5):** Under adversarial signal collision (randomly interleaving two codewords), the decoder successfully demonstrates true list decoding by recovering both constituent messages simultaneously with a 100% double capture rate.

## Usage

You can run the main simulation block by executing the Sage script in your terminal. The bottom of the `RS.sage` file contains the execution block for the experiment suite:

```python
sage RS.sage
```
## References
V. Guruswami, A. Rudra, and M. Sudan, Essential Coding Theory, Draft avail-
able at https://cse.buffalo.edu/faculty/atri/courses/coding-theory/book/
web-coding-book.pdf, 2025. (Chapters 12-13)

## Author

Nadav Cremisi

