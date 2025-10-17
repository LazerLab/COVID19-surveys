This repository contains the code for Quintana-Mathé and Uslu et al., *"Using Opt-In Non-Probability Surveys to Estimate COVID-19 Infection and Vaccination Rates"*, forthcoming in the Journal of Survey Statistics and Methodology (JSSAM).

The code structure is as follows:

- `infections.R` and `vaccinations.R` contain the main code to calculate national and state level vaccination and cumulative infection rates from the CSP, Axios-Ipsos, and the CDC data, generating Figures 1, 2, S2, S4, S5, S6 and Tables S3, S4, S7 and S8. These scripts also include the code to generate in-text numbers, for the third validation criteria, the data-driven confidence intervals, and for Figure S3 on the association of vaccination and cumulative infection rates with trust.  
- `functions_analysis.R` and `functions_read_data.R` contain helper functions for the rest of the code.
- `get_census.R` uses the census API to obtain the population totals for the calculations.

In COVID19-surveys/SM/:

- `vaccinations_unwt.R` and `infections_unwt.R` generate supplementary Figures S8 and S9 (same as Figure 1 and 2 but without using weights)
- `robustness_infections.R` generates Figure S1
- `sample_composition.R` generates Table S6
- `Facebook_recruitment.R` focuses on the comparison between the PureSpectrum and Facebook samples, and generates Figure 4 and the regressions of Tables S1 and S2.
- `downsample_trust_infections.R` and `downsample_trust_vax.R` includes the code for the simulations downsampling by trust levels, generating Figure 3. 

Data sources:

- CDC vaccination data: https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Age-and-Sex-Trends-in-the-Uni/5i5k-6cmh/
- Axios-Ipsos data: https://ropercenter.cornell.edu
- CSP/CHIP50 data: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ZUJU6H (https://doi.org/10.7910/DVN/ZUJU6H)
- NYT COVID-19 data: https://raw.githubusercontent.com/nytimes/covid-19-data/master/

Notes:

The CSP/CHIP50 dataset includes only the variables used in the main text analyses to prevent re-identification through individual demographic or response data. This causes some parts of the code in the Supplementary Material (COVID19-surveys/SM/) to not work properly as they need individual-level demographic data for reweighting, marginal effects etc.