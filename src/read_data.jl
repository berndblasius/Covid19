# two function to get the data
# data based on John Hopkins CSSE repository
# https://github.com/CSSEGISandData/COVID-19

module C

using CSV, DataFrames, Dates

# read-in covid-19 time series
# for infected, the funcion returns the prevalence (active cases)
# this is calculated by subtracting dead and recovered from cumulative cases
function covid_ts(country, starting_date)

  df_confirmed = CSV.read("../data/time_series_19-covid-Confirmed.csv")
  df_death     = CSV.read("../data/time_series_19-covid-Deaths.csv")
  df_recovered = CSV.read("../data/time_series_19-covid-Recovered.csv")
  first_date = Date(2020,1,22)              # first entry in database: Jan 22, 2020
  countries = df_confirmed[:,2]

  inds = findall(x->x==country, countries)
  maxdays = size(df_confirmed)[2] - 4          # number of reported days
  days = Dates.value(starting_date - first_date) # days we are interested
  I = zeros(Int,maxdays-days)
  D = zeros(Int,maxdays-days)
  R = zeros(Int,maxdays-days)
  for i in inds
      I .+= collect(df_confirmed[i, 5 + days:end])  # sum over all Provinces/States in a country
      D .+= collect(df_death[i,     5 + days:end]) 
      R .+= collect(df_recovered[i, 5 + days:end]) 
      
      # hack to interpolate wrong data on March 12
      ind_error = Dates.value(Date(2020,3,12)-starting_date) + 1  # time difference
      I[ind_error] = ceil(Int,sqrt(I[ind_error-1]*I[ind_error+1]))
      R[ind_error] = ceil(Int,sqrt(R[ind_error-1]*R[ind_error+1]))
      D[ind_error] = ceil(Int,sqrt(D[ind_error-1]*D[ind_error+1]))

  end
  I = I .- R .- D
  1:length(I),I,R,D
end

covid_ts(country) = covid_ts(country, Date(2020,1,22))


# calculate arrival times: first records of confirmed cases > th for all countries
function arrival_times(th=0, pandemie_start=Date(2019,12,1))
  # th: threshold for arrival
  # pandemie_start: hypothetical starting date of the pandemie
  df_confirmed = CSV.read("../data/time_series_19-covid-Confirmed.csv")
  df_death     = CSV.read("../data/time_series_19-covid-Deaths.csv")
  df_recovered = CSV.read("../data/time_series_19-covid-Recovered.csv")
  first_date = Date(2020,1,22)          # first entry in database: Jan 22, 2020
  countries = df_confirmed[:,2]

  firsts = zeros(Int,length(countries))
  for i = 1:length(countries)
      R = collect(df_confirmed[i,5:end])
      x = findfirst(x->x>th, R)
      firsts[i] = x == nothing ? length(names(df_confirmed))-4 : x 
  end

  # collect first records over Provinces/States of a country
  p = sortperm(countries)
  countries = countries[p]
  firsts = firsts[p]

  sep_countries = [countries[1]]
  sep_firsts    = [firsts[1]]
  c = countries[1]
  i = firsts[1]
  for i=2:length(countries)
      c1 = countries[i]
      i1 = firsts[i]
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
 
  # print into csv-file
  arrivals = first_date .+ Day.(sep_firsts .- 1)
  time_diff = Dates.value.(arrivals .- pandemie_start)  # time difference in days to start of pandemie
  df = DataFrame(country=sep_countries, arrival_time=arrivals, time_dist=time_diff)  # put into dataframe
  CSV.write("../arrival_times.csv", df) 

  # return vector of countries, dates of arrival
  sep_countries, sep_firsts, time_diff
end

arrival_times()

end
