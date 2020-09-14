# functions to load the data
# from John Hopkins CSSE repository
# https://github.com/CSSEGISandData/COVID-19
#
# the data should already be downloaded from the repository and be on 
# your local hard drive. Make sure to adapt the path variables below


using CSV, DataFrames, Dates

function smooth(x)
# make log-transform and smooth data
    lx = log10.(x)
    b = (1/2)*ones(2)     # flat (rectangular) window -> moving average
    lx = filtfilt(b,lx)
    10 .^ lx
end

# read-in covid-19 time series
# for infected, the funcion returns the prevalence (active cases)
# this is calculated by subtracting dead and recovered from cumulative cases
function covid_ts(country, starting_date)

  df_confirmed = CSV.read(csse_path*"time_series_covid19_confirmed_global.csv")
  df_death     = CSV.read(csse_path*"time_series_covid19_deaths_global.csv")
  df_recovered = CSV.read(csse_path*"time_series_covid19_recovered_global.csv")
  first_date = Date(2020,1,22)              # first entry in database: Jan 22, 2020
  countries = df_confirmed[:,2]
  countriesR = df_recovered[:,2]  # recovered still have a different file format
  

  inds = findall(x->x==country, countries)
  indsR = findall(x->x==country, countriesR)
  maxdays = size(df_confirmed)[2] - 4          # number of reported days
  days = Dates.value(starting_date - first_date) # days we are interested
  I = zeros(Int,maxdays-days)
  D = zeros(Int,maxdays-days)
  R = zeros(Int,maxdays-days)
  for i in inds
      I .+= collect(df_confirmed[i, 5 + days:end])  # sum over all Provinces/States in a country
      D .+= collect(df_death[i,     5 + days:end]) 
  end

  for i in indsR
      R .+= collect(df_recovered[i, 5 + days:end]) 
  end
  I = I .- R .- D
  1:length(I),I,R,D
end

covid_ts(country) = covid_ts(country, Date(2020,1,22))


# calculate arrival times: first records of confirmed cases > th for all countries
function arrival_times(th=0, pandemie_start=Date(2019,12,1))
  # th: threshold for arrival
  # pandemie_start: hypothetical starting date of the pandemie
  df_confirmed = CSV.read(csse_path*"time_series_covid19_confirmed_global.csv")
  df_death     = CSV.read(csse_path*"time_series_covid19_deaths_global.csv")
  df_recovered = CSV.read(csse_path*"time_series_covid19_recovered_global.csv")

  first_date = Date(2020,1,22)          # first entry in database: Jan 22, 2020
  countries = df_confirmed[:,2]
  maxdays = size(df_confirmed)[2] - 4          # number of reported days

  println("n countries: ", length(unique(countries)))

  firsts = zeros(Int,length(countries))
  for i = 1:length(countries)
      R = collect(df_confirmed[i,5:end])
      x = findfirst(x->x>th, R)
      firsts[i] = x == nothing ? maxdays : x 
  end

  # collect first records over Provinces/States of a country
  # probably can be done with two clever lines of code 
  # here the hard way..
  p = sortperm(countries)
  countries = countries[p]
  firsts = firsts[p]

  sep_countries = [countries[1]]
  sep_firsts    = [firsts[1]]
  c = countries[1]
  i = firsts[1]
  for j=2:length(countries)
      c1 = countries[j]
      i1 = firsts[j]
      if c1==c # same countrie
          if i1 < i
             sep_firsts[end] = i1 
             i = i1
          end
      else
         push!(sep_countries,c1)
         push!(sep_firsts,i1)
         c = c1
         i = i1
      end
  end
 
  arrivals = first_date .+ Day.(sep_firsts .- 1)
  time_diff = Dates.value.(arrivals .- pandemie_start)  # time difference in days to start of pandemie
  
  # return vector of countries, dates of arrival, and time_difference
  sep_countries, arrivals, time_diff
end


function read_all_data()
  df_confirmed = CSV.read(csse_path*"time_series_covid19_confirmed_global.csv")
  df_death     = CSV.read(csse_path*"time_series_covid19_deaths_global.csv")
  #df_recovered = CSV.read(csse_path*"time_series_covid19_recovered_global.csv")
  first_date = Date(2020,1,22)              # first entry in database: Jan 22, 2020
  countries = df_confirmed[:,2]

  unique_countries = unique(countries)
  ncountries = length(unique_countries)
  maxdays = size(df_confirmed)[2] - 4          # number of reported days
  I = zeros(Int,maxdays,ncountries)
  D = zeros(Int,maxdays,ncountries)
  #R = zeros(Int,maxdays,ncountries)

  j = 0
  for c in unique_countries
      j += 1
      inds = findall(x->x==c, countries)
      for i in inds
          I[:,j] .+= collect(df_confirmed[i, 5:end])  # sum over all Provinces/States in a country
          D[:,j] .+= collect(df_death[i,     5:end]) 
          #R[:,j] .+= collect(df_recovered[i, 5:end]) 
      end
  end
  1:length(I),I,D,unique_countries
  #1:length(I),I,D,R,unique_countries

end


# this is a path to the Johns Hopkins COVID19 repository
# it is assumed that you have a clone for this on your computer
csse_path = "../../COVID-19/csse_covid_19_data/csse_covid_19_time_series/"


#end # module
