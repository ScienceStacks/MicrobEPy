<img src="microbepy_logo.png" alt="drawing" width="200"/>

# MicrobEPy (Pronounced "microbe pie")

Microbes such as bacteria and archaea are fundamental to life on earth, especially to human health. Phenotypes of microbes such as growth rate evolve rapidly, often as a result of interactions with other microbes in a community. 

Understanding the evolutionary dynamics of microbe communities requires answering two questions: (a) What mutations are acquired over time by each species in the community? (b) How do mutations affect phenotypes? MicrobeEPy is a python package for data driven discovery to answer these questions.

Some biology background is required to understand the data and the science questions addressed. The microbial community in this study consists of two organisms, a bacteria Desulfovibrio vulgaris Hildenborough (hereafter, DVH), and an archaea, Methanococcus maripaludis (hereafter, MMP). There is no known natural habitat in which these organisms co-exist. However, by growing both with lactate, each assists the metabolic cycle of the other. When DVH metabolizes lactate, hydrogen is produced. MMP uses hydrogen in its metabolism to produce methane. By consuming hydrogen, MMP accelerates the rate of DVH metabolic reactions. And so, a virtuous cycle develops  (see Stolyar et al., 2007).

The rate at which this virtuous cycle proceeds depends on the many factors. The ability of DVH to transport lactate across its cell membrane; the efficiency of DVH lactate metabolism; the efficiency of the MMP hydrogen transport; and the efficiency of MMP metabolism. A community contains millions even billions of cells of DVH and MMP. Over time, there are mutations in both DVH and MMP to increase metabolism. This efficiency is reflected in the rate at which cells grow (referred to as growth rate) and the total mass of the living cells produced (referred to as yield).

A study was conducted in which two DVH variants and two MMP variants were growth under two different conditions (shaken or unshaken) for a total of eight different experimental conditions. There were three replications of each experimental condition, and hence a total of 24 evolutionary lines. Cells were extracted from lines over time. An extraction, or community, may be cells with the same mutations (isolates) or an extraction may contain cells with a diversity of mutations.

Three kinds of data were collected in this study. Community data details how the cells being exampled were obtained, such as the evolutionary line, the time at which they were extracted, and if they are isolates. Mutation data (obtained from DNA sequencing and a subsequent analysis pipeline) specifies the position in the genome where DNA base pairs have changed. Culture data describes the phenotype of the cells when they are growth together, their growth rate and yield.

To answer the question “How are mutations acquired over time?”, we use the community and mutation data. Our primary technique is a correlation analysis that quantifies the frequency with which a pair of mutations occur in communities over time. We also analyze the co-occurrence of mutations in isolates; that is, for cells with the same genome, is it common for two mutations to occur together. For example, to increase growth rate DVH needs to be increase the transport of lactate and increase the efficiency of lactate metabolism. These are different mutations, and so it is interesting if they co-occur in the same isolates.

To answer the question “How do mutations affect growth rate and yield?”, we use community, mutation, and culture data. We pursued a number of analyses including simple least squares linear regression, forward regression, and tree regression. These approaches failed to provide insight largely because of high variability in the data (because evolution itself is a highly variable process). We subsequently used classification (e.g., tree classifiers) to characterize “low” and “high” values of phenotype. This ultimately led to substantive results, including the identification of some important mutations.

The python package MicrobEPy represents the accumulation of the insights we acquired from doing the foregoing analysis. The package takes as input data that are structured as a SQL database (e.g., using sqlite). MicrobEPy provides services for data access, analysis, and visualization. Data access can be done through a SQL query or with helper functions that provide simple mechanisms for filtering rows and columns. The DataProvider class provides two pandas DataFrames. The first contains information about the community (row index) and mutations (column names). The second contains information about the community (row index) and phenotype (columns).

MicrobEPy provides extensive capabilities for analysis of data about microbe communities. This includes: calculations of correlations between mutations and communities; predictive models using linear regression, forward regression, and tree regression; and procedures for searching for patterns of mutations.

Our work identified a number of specialized visualizations that seemed particularly effective at revealing insights into microbe communities. The visualizations provided by MicrobEPy include: a plot of  “phenotype space” that relates growth rate to yield, plots of mutation frequencies, heatmaps of mutation correlations, and visual assessments of the effectiveness of phenotype predictions.

## Usage Notes
1. The project can be used as a git submodule.
1. The "project directory" refers to the folder containing this repository.
