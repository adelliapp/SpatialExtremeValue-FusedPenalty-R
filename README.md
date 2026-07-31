# Rainfall Extreme Value Modeling Research
Objectives
- Develop Generalized Pareto Distribution (GPD) models using fused lasso and fused ridge regularizations for modeling extreme rainfall.
- Predict return levels to estimate the future risks of extreme weather events.
- Model spatial dependency structures across different locations using a Gaussian copula.
- Compare model performance between fused lasso and fused ridge regularizations to identify the best model based on Takeuchi Information Criterion (TIC) values.

Methodology
- Approach: Extreme rainfall predictions in North Sumatra are analyzed using the GPD framework via the Peak Over Threshold (POT) method.
- Data Sources: Daily rainfall data from 2014 to 2024 collected from 5 BMKG meteorological stations in North Sumatra.
- Analysis Steps:Data preprocessing, Threshold determination, Exceedance data extraction, Parameter estimation using Maximum Likelihood Estimation (MLE) with regularization, Spatial dependency modeling using Gaussian copula, Model evaluation using TIC.
- Risk Assessment: Return level calculations for 10, 25, and 100-year return periods to illustrate extreme rainfall risks.

Result
- Regularized GPD Models: Parameter estimates for the Generalized Pareto Distribution integrated with fused lasso and fused ridge regularizations.
- Return Level Estimates: Calculated return levels representing the quantified risks of extreme rainfall events across specific return periods.
- Best-Performing Model: Identification of the optimal regularization approach determined by the lowest Takeuchi Information Criterion (TIC) value.
- Return Level Heatmap Visualization: Spatial mapping of extreme rainfall risks based on the best-performing GPD model.
- Spatial Dependency Structure: Modeled inter-location dependencies capturing the spatial connectivity of extreme weather events via Gaussian copula.


# SpatialExtremeValue-FusedPenalty-R
R code for spatial extreme value modeling using GPD with Fused Lasso and Fused Ridge regularization. Includes spatial dependence modeling by Gaussian Copula and return level estimation for 10, 25, and 100-year return periods and the visualization.
