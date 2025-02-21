## [TimeFlow: a density-driven pseudotime method for flow cytometry data analysis ](https://www.biorxiv.org/content/10.1101/2025.02.16.638508v1)
**Authors**: Margarita Liarou, Thomas Matthes, and Stéphane Marchand-Maillet.

Pre-print available [**here**](https://www.biorxiv.org/content/10.1101/2025.02.16.638508v1).


We developed TimeFlow, a new pseudotime computation method for the analysis of multi-dimensional flow cytometry data. TimeFlow orders the cells within a sample from the least to the most differentiated along their maturation pathway. It tracks cell transitions over a graph following smooth changes in the cell population density. We applied TimeFlow on healthy human bone marrow samples to model the temporal dynamics of twenty surface protein markers for monocytes, neutrophils, erythrocytes and B-cells.

![workflow](Figures/TimeFlow-Overview.png)

## Datasets
We have made available in the `Pre-processed-datasets` the dataset of P1-Monocytes. All new flow cytometry datasets (obtained from the bone marrow of three healthy patients) will be available upon manuscript acceptance.

All datasets have been pre-processed (see Supplementary Section S5) and stored in CSV format. The twenty first columns of each dataset correspond to the twenty CD markers and the last column contains the gating labels. There are fifteen datasets with linear trajectories: P1/2/3_Mono/Neu/Ery/Bcells.csv and three datasets with branching trajectories: P1/2/3-BM. The CD markers include: CD200, CD45, CD45RA, CD64, CD3, CD15, CD133, CD117, CD56, HLA.DR, CD19, CD33, CD34, CD371, CD7, CD16, CD123, CD36, CD38.
Note that the gating labels are only used for visualization and evaluation purposes and not during pseudotime computation. 

## Pseudotime analysis 
To reproduce the analysis for Monocytes presented in Results Section 3, as well as the Supplementary Figures S1-S5, please follow the Tutorials in the `P1-Monocytes-Analysis`. `P1-Monocytes-Analysis/Tutorial-1-Density-Estimation.ipynb` shows how to use a normalizing flow model to compute the probability density function of the observed single-cell data. In `P1-Monocytes-Analysis/Tutorial-2-Pseudotime-Computation.ipynb` we compute the cell pseudotime using TimeFlow and in `P1-Monocytes-Analysis/Tutorial-3-Marker-Dynamics.ipynb` we show how to fit and visualize the evolution of CD markers along pseudotime for the linear monocytic trajectory. 

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
The BibTeX for TimeFlow is the following:

```
@article {Liarou2025.02.16.638508,
	author = {Liarou, Margarita and Matthes, Thomas and Marchand-Maillet, St{\'e}phane},
	title = {TimeFlow: a density-driven pseudotime method for flow cytometry data analysis},
	elocation-id = {2025.02.16.638508},
	year = {2025},
	doi = {10.1101/2025.02.16.638508},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Pseudotime methods order cells undergoing differentiation from the least to most differentiated. We developed TimeFlow, a new method for computing pseudotime in multi dimensional flow cytometry datasets. TimeFlow tracks the differentiation path of each cell on a graph by following smooth changes in the cell population density. To compute the probability density function of the cells, it uses a normalizing flow model. We profiled bone marrow samples from three healthy patients using a 20-color antibody panel for flow cytometry and prepared datasets that ranged from 5,000 to 600,000 cells and included monocytes, neutrophils, erythrocytes and B-cells at various maturation stages. TimeFlow computed fine-grained pseudotime for all the datasets, and the cell orderings were consistent with prior knowledge of human hematopoiesis. Experiments showed its potential in generalizing across patients and unseen cell states. We compared our method to 11 other pseudotime methods using in-house and public datasets and found very good performance for both linear and branching trajectories. TimeFlow{\textquoteright}s pseudotemporal orderings are useful for modelling the dynamics of cell surface proteins along linear trajectories. The biologically meaningful results in branching trajectories suggest the possibility of future applications with automated cell lineage detection. Code is available at https://github.com/MargaritaLiarou1/TimeFlow and bone marrow data will be accessible upon acceptance.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2025/02/20/2025.02.16.638508},
	eprint = {https://www.biorxiv.org/content/early/2025/02/20/2025.02.16.638508.full.pdf},
	journal = {bioRxiv}
}

```

## Contact
Please contact us at *margarita.liarou@unige.ch* for any question about TimeFlow. 

## License
TimeFlow: a density-driven pseudotime method for flow cytometry data analysis is licensed under the Creative Commons Zero v1.0 Universal License. More information can be found [**here**](https://github.com/MargaritaLiarou1/TimeFlow/blob/main/LICENSE).