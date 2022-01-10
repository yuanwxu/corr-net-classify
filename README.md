# CACONET

CACONET is a computational framework that can be used to distinguish
between "diseased" and "healthy" microbial correlation networks
inferred from relative abundance data. It can also be used to identify
potential signature interactions characteristic of the networks,
offering possible targets for further biological and clinical
research. CACONET consists of an inference component for
compositional-aware correlation inference, a classification component
for the classification of the correlation networks, and an explanation
component for extraction of signature interactions from the
classifier. Please refer to the [publication](https://doi.org/10.1093/bioinformatics/btab879).

## Data preparation

Currently, CACONET supports binary classification. Suppose samples
from diseased and healthy subjects were collected and sequenced for
microbiome profiling, resulting in OTU tables containing relative taxa
abundances, `X_0` and `X_1` for healthy and disease, respectively. An
example demonstrating data preparation using 16S rRNA amplicon
sequencing of CRC fecal microbiome, downloaded from the [The
Microbiome Quality Control project](https://www.hmpdacc.org/MBQC/),
can be found in `mbqc_baseline_data.R` in the `MBQC` directory.

## Inference

Correlation network inference is done separately for `X_0` and `X_1`,
using a hierarchical Bayesian method,
[BAnOCC](https://doi.org/10.1371/journal.pcbi.1005852). Users are
advised to refer to the original publication and relevant
[tutorials](https://github.com/biobakery/banocc) for setting suitable
parameter values for their own data. Example of using the above real
data to obtain posterior correlation networks corresponding to healthy
and diseased samples is provided in `banocc_mbqc.R`

## Classification

CACONET then performs graph-level classification using the combined
posterior correlation networks from the inference step, thereby
incorporating posterior uncertainty. The
[DGCNN](https://www.aaai.org/ocs/index.php/AAAI/AAAI18/paper/viewPaper/17146)
algorithm was used and implemented with the Python library
[StellarGraph](https://stellargraph.readthedocs.io/en/stable/).

## Explanation

A greedy algorithm was used to find the top `n` most important nodes
that best differentiate the correlation networks. We can then extract
these nodes from the result of the inference step and visually examine
their associations. The file `dgcnn_mbqc_nnode.py` contains modules
for classification and explanation. The file `process_output_mbqc.R`
contains several useful functions for visualization and diagnostics.



