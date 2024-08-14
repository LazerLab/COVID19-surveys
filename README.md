This repository contains the full code for the *Large, Non-Probability Surveys Can Produce Valid Estimates of COVID-19 Vaccination and Cumulative Infection Rates* paper.

<br>

The main scripts are:

`get_census_pop.R`  extracts population counts from the 2020 US Census, that are used to estimate the official percentages of vaccination and cumulative infection rate thorough the code.  

`functions_read_data.R` and `functions_analysis.R` contain functions used thorough the code. 

`vaccinations.R` : Generates Figure 1, the text numbers of the vaccination section, Extended data Figures 1 and 2, part of Extended data Figure 9,  and SI Tables 2 and 6.

`infections.R`: Generates Figure 2, the text numbers of the infections sections, Extended data Figures 3, 4, and 5, SI section 3,  SI Tables 1 and 5, part of Extended data Figure 9. 

<br>

The SI folder contains the code for the SI and some extended data figures:

`SI2_infections.R`: Extended data Figure 6, SI Section 2.

`vaccinations_unwt.R`: Extended data Figure 7.

`infections_unwt.R`: Extended data Figure 8.

`sample_composition.R`: SI Table 4.

`SI4_downsample_trust_infections.R` and `SI4_downsample_trust_vax.R`: SI Section 4, Extended data Figure 10.

`SI5_downsample_state_vax.R`: SI Section 5.

`SI6_Facebook_recruitment.R`: SI Section 6. 
