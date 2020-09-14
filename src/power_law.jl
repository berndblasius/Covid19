# Source code for Fig.1
# power-law distribution in the number of reported COVID-19 cases
# in countries worldwide and in counties in the US
# this corresponds to Fig.1 in 
# Blasius B. (2020) Power-law distribution in the number of confirmed COVID-19 cases. Chaos
# https://doi.org/10.1063/5.0013031

module A

using CSV, DataFrames, Polynomials, PyPlot, StatsBase

include("read_data.jl")
include("tools.jl")


function distribution_countries()
   println("Read countries")
   t,C_ts,D_ts,countries = read_all_data()
   (days, ncountries) = size(C_ts)
   println("days ", days)
   println("total _countries ", ncountries)

   NNN=61 # this is March 22
   C = C_ts[NNN,:]  
   D = D_ts[NNN,:] 
   
   ind = findall(x->x>0, C)
   C = sort(C[ind], rev=true)
   nC = length(C)
   nmaxC = maximum(C)
   println("countries with case, n ",nC)
   println("max case: ",nmaxC)

   ind = findall(x->x>0, D)
   D = sort(D[ind], rev=true)
   nD = length(D)
   nmaxD = maximum(D)
   println("countries with death: ",nD)
   println("maxdeath: ", nmaxD)

   # Histogram with bins 2,4,8,16,..
   ex1,ey1 = log2histogram(C)
   x1 = log.(ex1)
   y1 = log.(ey1)
   p = Polynomials.fit(x1,y1,1)
   muC=-coeffs(p)[2]
   println("linear fit ",p)
   ey1fit = exp.(p.(x1))

   alphaC = 1 + length(C) / sum(log.(C)) # calculate exponent by log-likelihood
   println("alphaC ", alphaC)
   alphaC, stC = estimate(C)
   println("alphaC ", alphaC, " pm ", stC)

   # #########
   ex1D,ey1D = log2histogram(D)
   x1D = log.(ex1D)
   y1D = log.(ey1D)
   p = Polynomials.fit(x1D,y1D,1)
   muD=-coeffs(p)[2]
   println("linear fit ",p)
   ey1fitD = exp.(p.(x1D))

   alphaD = 1 + length(D) / sum(log.(D)) # calculate exponent by log-likelihood
   println("alphaD ", alphaD)
   alphaD, stD = estimate(D)
   println("alphaD ", alphaD, " pm ", stD)

   ax = PyPlot.axes([fig_x1, fig_y1, fig_width_x, fig_width_y])
       scatter(ex1,ey1,s=50,label="confirmed cases")
       plot(ex1,ey1,"-",linewidth=2)
       loglog(ex1,ey1fit,"k", label=L"regression line, $\mu=1.18$") 
       ylabel(L"frequency, $P_C(n)$",fontsize=18,labelpad=14)
       title("countries", fontsize=16)
       xlim(0.8,2e5)
       ylim(1e-7,1)
       xticks(fontsize=16)
       yticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],fontsize=16)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       ax1 = ax.transAxes
       text(-0.14, 0.95, "a)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix1, fig_y1+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,5.5, length=50)
       loglog(cc, cumulpdf.(cc,alphaC,nmaxC), color="black", linewidth=2)
       scatter(C,(1:nC)/nC,s=20)
       myfont = Dict("fontname"=>"Sans","fontsize"=>14)
       xticks(fontsize=14)
       yticks(fontsize=14)
       xlim(1,1e5)
       ylim(5e-3,1)
       xlabel("cases, n",fontsize=16)
       ylabel("C(n)",fontsize=16)

   ax = PyPlot.axes([fig_x1, fig_y2, fig_width_x, fig_width_y])
       scatter(ex1D,ey1D,s=50,label="confirmed deaths")
       loglog(ex1D,ey1D,"-",linewidth=2)
       loglog(ex1D,ey1fitD,"k", label=L"regression line, $\mu=1.35$") 
       xticks(fontsize=16)
       yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],fontsize=16)
       ylim(1e-6,1)
       xlim(0.8,1e4)
       xlabel("cases, n",fontsize=18)
       ylabel(L"frequency, $P_D(n)$",fontsize=18,labelpad=14)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       ax1 = ax.transAxes
       text(-0.14, 0.95, "b)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix1, fig_y2+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,5, length=50)
       loglog(cc, cumulpdf.(cc,muD,1.5e4), color="black", linewidth=2)
       scatter(D,(1:nD)/nD,s=20)
       myfont = Dict("fontname"=>"Sans","fontsize"=>14)
       xticks(fontsize=14)
       yticks(fontsize=14)
       xlim(1,1e4)
       ylim(1e-2,1)
       xlabel("cases, n",fontsize=16)
       ylabel("C(n)",fontsize=16)
end

function distribution_counties()
   println(); println()
   println("Read US counties")
   # make sure to have access to the data:
   # this means, clone into the Johns Hopkins Covid19 repository
   # and replace the folder variable below with the correct folder to the data
   folder = "../../COVID-19/csse_covid_19_data/csse_covid_19_daily_reports/"
   df = CSV.read(folder*"03-31-2020.csv")

   countries = df.Country_Region
   C = df.Confirmed
   D = df.Deaths

   ind = findall(x->x=="US", countries)
   C = C[ind]
   D = D[ind]

   ind = findall(x->x>0, C)
   C = sort(C[ind], rev=true)
   nC = length(C)
   println("no of entries: ", nC)
   println("max ", maximum(C), " min ", minimum(C))
   println(C[1:3])
   n_ones = length(findall(x->x==1, C))
   println("counties with one case: ", n_ones)

   # Histogram with bins 2,4,8,16,..
   ex1,ey1 = log2histogram(C)
   x1 = log.(ex1)
   y1 = log.(ey1)
   p = Polynomials.fit(x1,y1,1)
   muC=-coeffs(p)[2]
   println("linear fit log2 ",p)
   ey1fit_log2 = exp.(p.(x1))

   alphaC = 1 + length(C) / sum(log.(C)) # calculate exponent by log-likelihood
   println("alphaC ", alphaC)
   alphaC, stC = estimate(C)
   println("alphaC ", alphaC, " pm ", stC)


   ind = findall(x->x>0, D)
   D = sort(D[ind], rev=true)
   nD = length(D)
   println("no of entries: ", nD)
   println("max ", maximum(D), " min ", minimum(D))
   println(D[1:3])
   n_ones = length(findall(x->x==1, D))
   println("counties with one case: ", n_ones)

   # Histogram with bins 2,4,8,16,..
   ex1D,ey1D = log2histogram(D)
   x1D = log.(ex1D)
   y1D = log.(ey1D)
   p = Polynomials.fit(x1D,y1D,1)
   muD=-coeffs(p)[2]
   println("linear fit log2 deaths ",p)
   ey1fit_log2D = exp.(p.(x1D))

   alphaD = 1 + length(D) / sum(log.(D)) # calculate exponent by log-likelihood
   println("alphaD ", alphaD)
   alphaD, stD = estimate(D)
   println("alphaD ", alphaD, " pm ", stD)


   ax = PyPlot.axes([fig_x2, fig_y1, fig_width_x, fig_width_y])
       scatter(ex1,ey1,s=50,label="confirmed cases")
       plot(ex1,ey1,"-",linewidth=2)
       loglog(ex1,ey1fit_log2,"k", label=L"regression line, $\mu=1.58$") 
       title("US counties", fontsize=16)
       xlim(0.8,1e5)
       ylim(1e-8,1)
       xticks(fontsize=16)
       yticks([1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],fontsize=16)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       ax1 = ax.transAxes
       text(-0.14, 0.95, "c)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix2, fig_y1+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,5.0, length=50)
       loglog(cc, cumulpdf.(cc,muC,7e4), color="black", linewidth=2)
       scatter(C,(1:nC)/nC,s=20)
       loglog(C,(1:nC)/nC,"-")
       myfont = Dict("fontname"=>"Sans","fontsize"=>14)
       xticks(fontsize=14)
       yticks(fontsize=14)
       xlim(1,8e4)
       ylim(3e-4,1)
       xlabel("cases, n",fontsize=16)
       ylabel("C(n)",fontsize=16)

   ax = PyPlot.axes([fig_x2, fig_y2, fig_width_x, fig_width_y])
       scatter(ex1D,ey1D,s=50,label="confirmed deaths")
       loglog(ex1D,ey1D,"-",linewidth=2)
       loglog(ex1D,ey1fit_log2D,"k", label=L"regression line, $\mu=1.83$") 
       xticks(fontsize=16)
       yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],fontsize=16)
       xlim(0.8,1.5e3)
       ylim(1e-6,1)
       xlabel("cases, n",fontsize=18)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       ax1 = ax.transAxes
       text(-0.14, 0.95, "d)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix2, fig_y2+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,3, length=50)
       loglog(cc, cumulpdf.(cc,muD,3e3), color="black", linewidth=2)
       scatter(D,(1:nD)/nD,s=20)
       xticks([1,10,100,1000],fontsize=14)
       yticks(fontsize=14)
       xlim(1,1.5e3)
       ylim(1e-3,1)
       xlabel("cases, n",fontsize=16)
       ylabel("C(n)",fontsize=16)
end


fig_Iy = 0.25
fig_Ix1 = 0.37
fig_Ix2 = 0.84
fig_Iwidth_y = 0.12
fig_Iwidth_x = 0.12

fig_sep_y   = 0.06
fig_width_y = 0.4
fig_width_x = 0.4
fig_x1 = 0.1
fig_x2 = 0.57
fig_y2=0.08
fig_y1=fig_y2 + fig_width_y + fig_sep_y

figure("distribution", figsize=(12,8))
clf()
distribution_countries()
distribution_counties()
#savefig("../figures/fig1.pdf")


end # module
