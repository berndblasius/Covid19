# two function to get the data
# data based on John Hopkins CSSE repository
# https://github.com/CSSEGISandData/COVID-19


using CSV, Dates

# read-in covid-19 time series
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
      I .+= collect(df_confirmed[i, 5 + days:end]) 
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


