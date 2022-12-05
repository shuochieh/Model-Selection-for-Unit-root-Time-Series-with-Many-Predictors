# Model-Selection-for-Unit-root-Time-Series-with-Many-Predictors

# Data

Data are provided in the public domain. ⚠️ Remember to supplement ⚠️

# Code
## Abstract
Codes for the simulation studies are in the simulation folder, whereas codes for real data analysis are in the real data folder.

We use R 4.0.2 to run these codes, along with packages <code>glmnet, Ohit, fGarch, RcppEigen, doParallel, ggplot2, ggfortify</code>.
For both the simulation and real data analysis, we run the codes on the university's computing clusters using multiple cores. One can adjust the number of cores to use in the codes by changing the value of <code>no_clusters</code> in the codes. ⚠️ remember to change Nslots to no_clusters in simulation codes ⚠️

## Reproducibility workflow

### Simulations
⚠️ Remember to include the file that aggregates the saved results ⚠️


### Real data analysis

- Draw Figure 1
- Clean data
- Draw Figure 2
- Rolling-window prediction
