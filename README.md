## [TimeFlow: a density-driven pseudotime method for flow cytometry data analysis ](https://)
**Authors**: Margarita Liarou, Thomas Matthes, and and Stéphane Marchand-Maillet.

Pre-print is available [**here**](https://).


We developed TimeFlow, a new pseudotime computation method for the analysis of multi-dimensional flow cytometry data. TimeFlow orders the cells within a sample from the least to the most differentiated along their maturation pathway. It tracks cell transitions over a graph following smooth changes in the cell population density. We applied TimeFlow on healthy human bone marrow samples to model the temporal dynamics of twenty surface protein markers for monocytes, neutrophils, erythrocytes and B cells.

![workflow](Figures/TimeFlow-Overview.png)

## Datasets
We have made available in the `Pre-processed-datasets` the dataset of P1-Monocytes. You can access all our new flow cytometry datasets, obtained from the bone marrow of three healthy patients, at [this Open Science Framework repository](https://osf.io/rs59m/)

All datasets have been pre-processed (see Supplementary Section S5) and stored in CSV format. The twenty first columns of each dataset correspond to the twenty CD markers and the last column contains the gating labels. There are fifteen datasets with linear trajectories: P1/2/3_Mono/Neu/Ery/Bcells.csv and three datasets with branching trajectories: P1/2/3-BM. The CD markers include: CD200, CD45, CD45RA, CD64, CD3, CD15, CD133, CD117, CD56, HLA.DR, CD19, CD33, CD34, CD371, CD7, CD16, CD123, CD36, CD38.
Note that the gating labels are only used for visualization and evaluation purposes and not during pseudotime computation. 

## Pseudotime analysis 
To reproduce the analysis for Monocytes presented in Results Section 3, as well as the Supplementary Figures S1-S5, please follow the Tutorials in the `P1-Monocytes-Analysis`. `P1-Monocytes-Analysis/Tutorial-1-Density-Estimation.ipynb` shows how to use a normalizing flow model to compute the probability density function of the observed single-cell data. `P1-Monocytes-Analysis/Tutorial-2-Pseudotime-Computation.ipynb` shows how to compute the cell pseudotime using TimeFlow and `P1-Monocytes-Analysis/Tutorial-3-Marker-Dynamics.ipynb` shows how to fit and visualize the evolution of CD markers along pseudotime for the linear monocytic trajectory. 

## Requirements

The Python/PyTorch requirements are the following:

- Python version: 3.11.3

Python Package versions:

- numpy==1.24.3
- pandas==2.2.2
- sklearn==1.4.2
- igraph==0.10.8
- torch==2.0.1+cpu
- seaborn==0.12.2
- matplotlib==3.7.1
- pygam==0.8.0

##  Citation
You can cite as it follows:

## Contact
TimeFlow is under continuous development and upgrades will follow soon. Please contact us at **margarita.liarou@unige.ch** for any question about the manuscript or support in using TimeFlow for your research. 
 
## License
TimeFlow: a density-driven pseudotime method for flow cytometry data analysis is licensed under the Creative Commons Zero v1.0 Universal License. More information can be found [**here**](https://github.com/MargaritaLiarou1/TimeFlow/blob/main/LICENSE).