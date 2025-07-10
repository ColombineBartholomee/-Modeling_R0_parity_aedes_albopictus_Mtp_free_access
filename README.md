# -Modeling_R0_parity_aedes_albopictus_Mtp_free_access

R scripts developed in the frame of the V2MOC projet (Impact de la végétalisation des villes sur le risque vectoriel chez l’Homme et les animaux), as part the PhD thesis of Colombine Bartholomée (IRD, UMR MIVEGEC). It was funded by the RIVOC project: by the Occitanie region and the University of Montpellier. This folder contains the codes used to model the parturity of Aedes albopictus females based on environmental variables, as well as the calculation of the basic reproduction number (R₀) using field data. Two formulas are implemented for R₀ estimation: one from Garret-Jones in 1964 and used by Anne Lise Tran in 2005 and another from Poletti (2011).

This folder contains:

    01_Modele_Univ.R: performs univariate modeling of the probability of a female being parous.
    02_Univ_Viz.R: visualizes the results of the univariate modeling.
    3_Multi_GLMM.R: conducts multivariate analysis based on the results of the univariate modeling.
    04_Viz_Multi_GLMM.R: visualizes the results of the multivariate analysis.
    05_GAMMS_Predictions: realization of GAMMs for predictions of parity rate and abundance of females/trap
    06_Calcul_R0.R: presents the results of the two formulas used for R₀ estimation across months and environments, using two methods to compute confidence intervals: based on the standard error of the final R₀ or using the confidence intervals of each parameter involved in the R₀ calculation.
    07_Sensibility_Analysis_Sobol_Indices.R: provides an estimation of sobol indices for different parameters and different method used.
    08_Sensibility_Analysis_Variance_Fast99: provides an estimation of sensibility indices for different parameters and different method use and variance, suing the method fast 99.
    09_Visualise_R0 : allows to visualise plots of R0.
    10_Shiny_App_Exposition_Vizualise_Results: application shiny allowing to check all plots.
    functions.R: contains functions for visualizing univariate results.
    functions_R0.R: contains functions for calculating R₀ using the two formulas.
