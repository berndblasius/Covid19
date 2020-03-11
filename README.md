# Covid19
Wang et al SEIR-model for 2019/20 Covid-19 outbreak


**Description**: A julia implementation of the SEIR-like model describing the Coronavirus Covid-19 outbreak 2019/20,
based on the publication: Wang et al. (2020) Evolving epidemiology and impact of non-pharmaceutical interventions on the outbreak of coronavirus disease 2019 in Wuhan, China. MedRxive https://doi.org/10.1101/2020.03.03.20030593


**Disclaimer**: 
*This code is not meant to be a qualatitive tool for predicting the current outbreak dynamics in different parts of the world*:
- The publication by Wang et al. was posted on a preprint server and has not yet been passed peer review.
- I did not contact these authors, but instead just implemented the code from their publication, as I understood it, into Julia.
- The parameter values from the paper have been fitted to the situation in Wuhan. There are many reasons to doubt whether these parameters can simply be transfered to the situation in Europe, or other regions in the world.

The motivations for this code is merely to give people a starting ground for running and exloring outbreak simulations.  

**Description**: The model is a variant of a standard 
Susceptible-Exposed-Infectious-Recovered (SEIR) model.
Additionally, the model distinguishes between, the reported number of infected (I) and the unreported number of infected (A). Individuals in the A-class can cause further infections, in the same way as the reported (or acertained) infected, but they are 'invisible' to mitigation and hospitalization. The model also includes a class of hospitalized (or quarantined) individuals (H). Only acertained infected can enter the H-class and individuals from the H-class cannot infect susceptibles. Additionally the model allows to consider inflows and outflows of people into the region.

**Some Results** Taking the parameter values that were fitted to Wuhan, it is possible to simulate the first weeks of the outbreak in Germany and (to some extend) of Italy. Playing around a bit (no serious parameter fitting yet) shows that for capturing the German outbreak requires either to increase the transmission rate from b=1.75 to r=2, or to increase the acertainment rate from r=0.19 to r=0.32. This means that if the transmission rate in Germany is similar to that in Wuhan, maybe the rate of detection of infected is larger.  
The calculated basic reproduction rate with R0=3.8 is much larger than the typical assumed values of R0=2..2.5. This may indicate that the model is either not parameterized well enough or that the model structure, or model simplifications (e.g. neglecting heterogeneity in contact networks), are not really adequate to capture the dynamics (at least in Europe).  
Corresponding to the inflated R0-value, the model also yields a too large fraction of totally infected individuals compared to values given by the experts.
 