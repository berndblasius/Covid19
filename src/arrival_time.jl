# nlyse COVID-19 arrival time in countries
# this corresponds to Fig.2 in 
# Blasius B. (2020) Power-law distribution in the number of confirmed COVID-19 cases. Chaos
# https://doi.org/10.1063/5.0013031

module A

using CSV, DataFrames, Polynomials, PyPlot, StatsBase
using LsqFit

include("read_data.jl")
include("tools.jl")



function time_exponents()
# compute time dependent exponents
   t,C_ts,D_ts,countries = read_all_data()

   (days, ncountries) = size(C_ts)
   println("days ", days)
   println("total _countries ", ncountries)

   # calculate exponent from log-likelihood
   expsC = zeros(days)
   expsD = zeros(days)
   for t=1:days
       C = C_ts[t,:]
       C = C[findall(x->x>0, C)]  # remove zeros
       mu, st = estimate(C)
       expsC[t] = mu

       D = D_ts[t,:]
       D = D[findall(x->x>0, D)]  # remove zeros
       mu, st = estimate(D)
       expsD[t] = mu
   end
   expsC, expsD
end


function prevalence_overtime(fig, fig_x, fig_y1,fig_width_y)
   t,I_ts,D_ts,countries = read_all_data()

   (days, ncountries) = size(I_ts)
   println("days ", days)
   println("total _countries ", ncountries)
   
   first_date = Date(2020,1,22)  # first entry in database: Jan 22, 2020
  
   # make sure here to include only the data until the time span
   # that should be investigated
   last_date = first_date + Day(days-1)
   println(last_date)

   cols = ["k", "r", "b", "g", "m", "y"]
   times = [0,14,28,42,56,70]   # vector of time delays
   times = times .+ (days - 76)  
   dates = string.(last_date .- Day.(times))
   println(dates)

   for i = 1: length(times)
       t = times[i] # time difference
       c = cols[i]  # color

       I = I_ts[end .- t,:]
       ind = findall(x->x>0, I)
       I = sort(I[ind], rev=true)
       nI = length(I)
       println("countries, n ",nI)

       D = D_ts[end .- t,:]
       ind = findall(x->x>0, D)
       D = sort(D[ind], rev=true)
       nD = length(D)

       # log -bin Histogram
       ex1I,ey1I = log2histogram(I)
       x1I = log.(ex1I)
       y1I = log.(ey1I)
       ex1D,ey1D = log2histogram(D)
       x1D = log.(ex1D)
       y1D = log.(ey1D)

   ax = PyPlot.axes([fig_x, fig_y1, 0.8, fig_width_y])
       loglog(ex1I,ey1I,"-", color=c)
       scatter(ex1I,ey1I,s=50, color=c, label="distribution on "*dates[i])
       myfont = Dict("fontname"=>"Sans","fontsize"=>20)
       xticks(fontsize=14)
       yticks(fontsize=14)
       xlim(1,1e6)
       xlabel("cases, n",fontsize=16)
       ylabel(L"frequency, $P(n)$",fontsize=16, labelpad=10)
       legend(loc="lower left",fancybox="true", fontsize=12)
       ax1 = ax.transAxes
       text(-0.15, 0.95, "a)", fontsize=18, fontweight="bold", transform = ax1)

   ax = PyPlot.axes([0.67, 0.82, 0.26, 0.15]) # Inset: cumulative distribution 
       loglog(I,1:nI,"-", color=c)
       ylim(1,200)
       xlim(1,5e5)
       xticks(fontsize=12)
       yticks(fontsize=12)
       myfont = Dict("fontname"=>"Sans","fontsize"=>10)
       xlabel("cases, n", fontsize=14)
       ylabel(L"$N \cdot C(n)$", fontsize=14)
   end

end



function arrival_histogram(fig, fig_x, fig_y2,fig_width_y)

  c,d,t = arrival_times()
  ind = sortperm(t)  
  t = t[ind]; c=c[ind]; d = d[ind]
  # first entry Jan 22 2020
  # at this day 6 countries already invaded, China, Japan, South Korea, Thailand, Taiwan, US)
  t = t .- t[1]

  h = StatsBase.fit(Histogram, t,  nbins=100, closed=:left)
  h_weight = h.weights
  h_bins = h.edges[1][2:end]

  println(size(h_bins))
  println(size(h_weight))
  @. model(x, p) = p[1]*exp(x*p[2])
  p0 = [0.5, 0.5]
  expfit = curve_fit(model, collect(h_bins)[1:57], h_weight[1:57], p0)
  println("exp fit: ",coef(expfit))
  cf = coef(expfit)

  alphaC, alphaD =  time_exponents()
  ndays = length(alphaC)

  h_dates = Date(2020,1,21) + Day.(h_bins)

  ind_march22 = Dates.value(Date(2020,3,22) - Date(2020,1,22)) + 1

  show_days= Dates.value(Date(2020,4,05) - Date(2020,1,22)) + 1

   majorformatter = matplotlib.dates.DateFormatter("%d.%m")
   majorlocator = matplotlib.dates.DayLocator(interval=7)

   ax = PyPlot.axes([fig_x, fig_y2, 0.8, fig_width_y])
       bar(h_dates,h_weight)
       plot(h_dates[1:61], cf[1]*exp.(cf[2]*h_bins[1:61]),"k", linewidth=2,
       linestyle="dashed")
       plot([Date(2020,3,22),Date(2020,3,22)], [0,15],"-r", linewidth=2)
       ax.xaxis_date()
       ax.xaxis.set_major_formatter(majorformatter)
       ax.xaxis.set_major_locator(majorlocator)
       xticks(Date(2020,1,26) .+ Day.(0:7:10*7),fontsize=14)
       yticks(fontsize=14)
       xlim(Date(2020,1,21),Date(2020,4,05))
       ylim(0,15)
       ylabel("countries",fontsize=16, labelpad=10)
       xlabel("arrival date", fontsize=16)
       for label in ax.get_xticklabels()
           label.set_ha("right")
           label.set_rotation(30)
       end
       ax1 = ax.transAxes
       text(-0.15, 0.95, "b)", fontsize=18, fontweight="bold", transform = ax1)

   ax1 = PyPlot.axes([0.25, 0.3, 0.28, 0.15]) # Inset: cumulative distribution 
       plot(1:ndays,alphaC, linewidth=2, color="black")
       # plot the exponent of deaths from day 26, the first day with at least
       # 5 countries with a confirmed death case
       plot(26:ndays,alphaD[26:ndays], linewidth=2, color="blue")
       plot([ind_march22, ind_march22],[1,2],"r-")
       xticks(fontsize=12)
       yticks(fontsize=12)
       xlim(0,show_days)
       ylim(1,2)
       ax1.yaxis.grid()
       xlabel("time (days)", fontsize=14)
       ylabel(L"exponent $\mu$", fontsize=14, labelpad=10)
end


function make_plot()
    fig_sep_y   = 0.1
    fig_width_y = 0.4
    fig_width_y = 0.39
    fig_x = 0.15
    fig_y2=0.1
    fig_y1=fig_y2 + fig_width_y + fig_sep_y

    fig_Iy = 0.25
    fig_Iwidth_y = 0.1
    fig_Iwidth_y = 0.12
    fig_Iwidth_x = 0.22
    fig_Iwidth_x = 0.26
    
    fig=figure("prevalence over time1", figsize=(8,12))
    clf()
    prevalence_overtime(fig, fig_x, fig_y1,fig_width_y)
    arrival_histogram(fig, fig_x, fig_y2,fig_width_y)

   #savefig("../figures/distribution_time.png")
end

make_plot()

end  #module
