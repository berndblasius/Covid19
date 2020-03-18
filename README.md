# Covid19 
## Modelling and data analysis tools to study the 2020 COVID-19 outbreak

Algorithms and code, written by me, to analyse outbreak data and so epidemiological simulations - in my favourite computer language [Julia](https://julialang.org).

**Data**: The folder `data` contains time series data taken from John Hopkins CSSE repository
https://github.com/CSSEGISandData/COVID-19  
Additionally, the folder contains a file with arrival dates of the virus in different countries (calculated by the function `arrival_time` in the src-folder).

**Results**: The folder `figures` contains some figures with analysis results.  
Further results of the analysis can be found at my [web page](https://www.staff.uni-oldenburg.de/bernd.blasius/project/corona). 



**Code**: The folder `src` contains algorithms for simulation of outbreak and analysis of data. 


##### read_data.jl  
Functions to read-in outbreak time series for a given country. Data are taken from the [John Hopkins CSSE repository](https://github.com/CSSEGISandData/COVID-19). I will try to update the data regularly. The functions include some minor preprocessing, such as interpolation for known irregularities in the John Hopkins data set.
The returned number of infected gives the prevalence, that is, the number of active cases, in contrast to the often given cumulative number of cases.  
Additionally a function to calculate arrival dates in countries.

##### corona_tsa.jl  
Some functions for time series analysis of outbreak data. The focus is on plotting current dynamics and trends of the outbreaks in different countries.
In particular, daily incidence, visualization of log-transformed data and analysis of local trends. 

##### Wang_SEIAHR.jl
A julia implementation of the SEIR-like model describing the Coronavirus COVID-19 outbreak 2019/20,
based on the publication: Wang et al. (2020) Evolving epidemiology and impact of non-pharmaceutical interventions on the outbreak of coronavirus disease 2019 in Wuhan, China. MedRxive https://doi.org/10.1101/2020.03.03.20030593


**Disclaimer**: 
*This code is not meant to be a quantitive tool for predicting the current outbreak dynamics in different parts of the world*:
- The publication by Wang et al. was posted on a preprint server and has not yet been passed peer review.
- I did not contact these authors, but instead just implemented the code from their publication, as I understood it, into Julia.
- The parameter values from the paper have been fitted to the situation in Wuhan. There are many reasons to doubt whether these parameters can simply be transfered to the situation in Europe, or other regions in the world.

The motivations for this code is merely to give people a starting ground for running and exploring outbreak simulations.  

**Description**: The model is a variant of a standard 
Susceptible-Exposed-Infectious-Recovered (SEIR) model.
Additionally, the model distinguishes between the reported (or acertained) number of infected (I) and the unreported number of infected (A). Individuals in the A-class can cause further infections, in the same way as the reported infected, but they are 'invisible' to mitigation and hospitalization. The model also includes a class of hospitalized (or quarantined) individuals (H). Only acertained infected can enter the H-class and individuals from the H-class cannot infect susceptibles. Additionally the model allows to consider inflows and outflows of people into the region.


 