# Simple metapopulation model for coupled countries with SIR dynamics
# The aim is to see whether we find similar power-law distribution in the
# prevalence, as in the COVID-19 data
#
# see my paper https://arxiv.org/abs/2004.00940

module C

using Dates, PyPlot, Random, StatsBase, DifferentialEquations, Polynomials
using LsqFit


function SIRD(du,u,p,t)  # standard SIR model + dead class
  S,I,R,D = u
  b,g,m = p
  N = S+I+R+D
  du[1] = dS = -b*S*I/N       # susceptible
  du[2] = dI = b*S*I/N  - g*I # infected
  du[3] = dR = (1-m)*g*I      # recovered
  du[4] = dD = m*g*I          # dead
end

function SIRDmeta(du,u,p,t)  # SIRD metapopulation model (uncoupled countries)
  b,g,m,N = p
  j = 1
  for i=1:N
      S,I,R,D = u[j:j+3]
      Nc = S+I+R+D
      du[j] = dS = -b*S*I/Nc       # susceptible
      du[j+1] = dI = b*S*I/Nc - g*I # infected
      du[j+2] = dR = (1-m)*g*I      # recovered
      du[j+3] = dD = m*g*I          # dead
      j = j+4
  end
end

cumulpdf(x,alpha,xmax) =  x^(1-alpha) - xmax^(1-alpha)

function loghistogram(ts, nb=20)
   his = fit(Histogram, ts, nbins=nb)
   x1 = his.edges[1][2:end]
   ind = findall(x->x != 0, his.weights)  # remove bins with entry 0
   w = his.weights[ind]
   w = w ./ float(sum(w))
   x1 = x1[ind]
   y1 = log.(w) .- x1
   ex1 = exp.(x1)
   ey1 = exp.(y1)
   x1,y1,ex1,ey1
end

function nlys_discrete_SIR()
    Ti = 6.0         # infectious period
    gamma = 1. / Ti  # recovery rate
    m = 0.01         # case fatality rate
    b = 0.4          # contact rate

    Nc  = 50 * 1e6  # population size 

    I0 = 1
    S0 = Nc - I0
    R0 = 0
    D0 = 0

    ndays = 60    # runtime
    N = 200       # max number of countries
    p = 0.0005    # invasion probability
    firsts = [0]

    alphaI = zeros(ndays)
    alphaD = zeros(ndays)
    N1 = 1 
    u0 = [S0;I0;R0;D0]   # initial values
    tspan = (0.0,1.0)    # run the SIR-model for bursts of 1 day
    for i=1:ndays
       par = (b,gamma,m,N1)
       prob = ODEProblem(SIRDmeta,u0,tspan,par)
       sol = solve(prob)
       u0 = sol[:,end]

       for j=N1+1:N  # infections of other countries
           if minimum(rand(N1)) < p
               #push!(u0,1)
               push!(u0,S0)
               push!(u0,I0)
               push!(u0,R0)
               push!(u0,D0)
               N1 += 1
               push!(firsts,i)  # invasion at time i
           end
       end
       # calculate time-dependent critical exponent by log-likelihood
       tmp = reshape(u0,(4,:)) # vector of infected
       Iv = tmp[2,:] # vector of infected
       Rv = tmp[3,:] # vector of dead
       Dv = tmp[4,:] # vector of dead
       Iv = Iv .+ Rv .+ Dv
       ind = findall(x->x>=1,Dv)
       Dv = Dv[ind]
       alphaI[i] = 1 .+ length(Iv) ./ sum(log.(Iv)) # calculate exponent by log-likelihood
       alphaD[i] = 1 .+ length(Dv) ./ sum(log.(Dv)) # calculate exponent by log-likelihood
    end

    u0 = reshape(u0,(4,:))
    S = u0[1,:]
    I = u0[2,:]
    R = u0[3,:]
    D = u0[4,:]
    C = I .+ R .+ D   # total number of cases

    ind = findall(x->x>0, C)
    C = sort(C[ind], rev=true)
    nC = length(C)
   
    ind = findall(x->x>=1, R)
    R = sort(R[ind], rev=true)
    nR = length(R)

    ind = findall(x->x>=1, D)
    D = sort(D[ind], rev=true)
    nD = length(D)


    #  *********** Arrival time histogram  ***********************

    println(N1, " countries infected")
    h = fit(Histogram, firsts,  nbins=60, closed=:left)
    h_weight = h.weights
    h_bins = h.edges[1][1:end-1]

    @. model(x, p) = p[1]*exp(x*p[2])
    p0 = [0.5, 0.5]
    expfit = curve_fit(model, collect(h_bins), h_weight, p0)
    println(coef(expfit))
    cf = coef(expfit)

   figure("meta arrival", figsize=(12,8))
   clf()
   ax = PyPlot.axes([0.1, 0.1, 0.8, 0.8])
       bar(h_bins,h_weight)
       plot(h_bins, cf[1]*exp.(cf[2]*h_bins),"k", linewidth=2,
       linestyle="dashed")
       ylabel("countries",fontsize=20, labelpad=10)
       xticks(fontsize=20)
       yticks(fontsize=20)
       xlim(0,61)
       legend(bbox_to_anchor=(0.03, 0.2), loc="upper left", fontsize=16, fancybox="true")
       xlabel("arrival time (days)",fontsize=20)
   ax = PyPlot.axes([0.22, 0.55, 0.28, 0.3]) # Inset: cumulative distribution 
       plot(1:ndays,alphaI, linewidth=2, color="black")
       plot(1:ndays,alphaD, linewidth=2, color="blue")
       xticks(fontsize=16)
       yticks(fontsize=16)
       xlim(0,60)
       ylim(1,2)
       ax.yaxis.grid()
       xlabel("time (days)", fontsize=16)
       ylabel(L"exponent, $\mu$", fontsize=16, labelpad=10)
   savefig("../paper_figures/meta_arrival.png")


   # ************ distribution   ******************************** 
   x1,y1,ex1,ey1 = loghistogram(log.(C), 25)
   pol = polyfit(x1,y1,1)
   muC=-coeffs(pol)[2]
   println("linear fit muC:",muC)
   pol = polyfit(x1,y1,1)
   ey1fit = exp.(polyval(pol,x1))
   alpha = 1 + length(C) / sum(log.(C)) # calculate exponent by log-likelihood
   println("alpha ", alpha)
   
   x1R,y1R,ex1R,ey1R = loghistogram(log.(R), 15)
   pol = polyfit(x1R,y1R,1)
   muR=-coeffs(pol)[2]
   println("linear fit muR: ",muR)
   pol = polyfit(x1R,y1R,1)
   ey1fitR = exp.(polyval(pol,x1R))
   alphaR = 1 + length(R) / sum(log.(R)) # calculate exponent by log-likelihood
   println("alphaR ", alphaR)
   
   x1D,y1D,ex1D,ey1D = loghistogram(log.(D), 15)
   pol = polyfit(x1D,y1D,1)
   muD=-coeffs(pol)[2]
   println("linear fit muD: ",muD)
   pol = polyfit(x1D,y1D,1)
   ey1fitD = exp.(polyval(pol,x1D))
   alphaD = 1 + length(D) / sum(log.(D)) # calculate exponent by log-likelihood
   println("alphaD ", alphaD)
   
   fig_sep_y   = 0.06
   fig_width_y = 0.4
   fig_Iy = 0.25
   fig_Ix = 0.68
   fig_Ix = 0.66
   fig_Iwidth_y = 0.1
   fig_Iwidth_y = 0.12
   fig_Iwidth_x = 0.22
   fig_Iwidth_x = 0.26

   fig_x = 0.18
   fig_y2=0.08
   fig_y1=fig_y2 + fig_width_y + fig_sep_y

   figure("meta2", figsize=(8,12))
   clf()
   ax = PyPlot.axes([fig_x, fig_y1, 0.8, fig_width_y])
       scatter(ex1,ey1,s=50,label="total detected cases")
       loglog(ex1,ey1fit,"k", label=L"regression line, $\mu=1.22$") 
       xlim(1,5e6)
       ylabel(L"frequency, $P_C(n)$",fontsize=18,labelpad=10)
       xticks(fontsize=16)
       yticks(fontsize=16)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       #text(0.02, 0.85, "a", fontsize=18, fontweight="bold", transform = ax1)
       text(-0.2, 0.95, "a)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix, fig_y1+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,7, length=50)
       loglog(cc, N1*cumulpdf.(cc,muC,1e7), color="black", linewidth=2)
       scatter(C,1:N1,s=20)
       xticks(fontsize=16)
       yticks(fontsize=16)
       xlim(1,1e7)
       ylim(1,200)
       xlabel("cases, n",fontsize=16)
       ylabel("C(n)",fontsize=16)

   ax = PyPlot.axes([fig_x, fig_y2, 0.8, fig_width_y])
       scatter(ex1D,ey1D,s=50,label="deaths")
       loglog(ex1D,ey1fitD,"k", label=L"regression line, $\mu=1.23$") 
       xticks(fontsize=16)
       yticks(fontsize=16)
       xlim(1,5e6)
       xlabel("cases, n",fontsize=18)
       ylabel(L"frequency, $P_D(n)$",fontsize=18,labelpad=10)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       ax1 = ax.transAxes
       text(-0.2, 0.95, "b)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix, fig_y2+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,5.5, length=50)
       loglog(cc, nD*cumulpdf.(cc,muD,1e4), color="black", linewidth=2)
       scatter(D,1:nD,s=20)
       myfont = Dict("fontname"=>"Sans","fontsize"=>14)
       xticks(fontsize=16)
       yticks(fontsize=16)
       xlim(1,1e5)
       ylim(1,100)
       xlabel("cases, n",fontsize=16)
       ylabel("C(n)",fontsize=16)
       #savefig("../paper_figures/meta.png")



    # ***************  Lorenz curve  ************************
    figure("Lorenz curve")
    clf()
    cc = (1:nC) / nC
    wealth = cumsum(C) / maximum(cumsum(C))
    plot(cc, wealth)
    xlabel("fraction of countries")
    ylabel("fraction of cases")


    B = 1 - sum(wealth) / nC
    G = 1 - 2*B
    println("Gini coefficient ", G)

end



function run_SIR()
# run a simple SIR model

  tmax = 120.0  # all time units in days

  Ti = 6.0        # infectious period
  gamma = 1. / Ti  # recovery rate
  m = 0.01         # case fatality rate
  b = 0.4          # contact rate
  
  println("initial growth rate ", b - gamma)
  println("R0 ", b/gamma)
  println("T12 ",log(2)/(b - gamma))
  println("b ",b)

  N  = 50 * 1e6  # population size 
  I0 = 1
  S0 = N - I0
  R0 = 0
  D0 = 0

  u0 = [S0;I0;R0;D0]   # initial condistion
  tspan = (0.0,tmax)
  p = (b,gamma,m)
  prob = ODEProblem(SIRD,u0,tspan,p)
  sol = solve(prob)
  t = 0:0.1:tmax  # use interpolation of results
  res = sol(t)  
  S = res[1,:]
  I = res[2,:]
  R = res[3,:]
  D = res[4,:]

  C = I .+ R .+ D   # total number of cases

  figure("SIR", figsize=(12,8))
  clf()
  ax = gca()
  plot(t,C,"k-", label="total cases", linewidth=2)
  plot(t,S,"g-", label="susceptible", linewidth=2)
  plot([],[],"r-",label="infected", linewidth=2)
  plot(t,R,"b-", label="recovered", linewidth=2)
  plot([],[],"m-",label="deaths")
  xticks(fontsize=16)
  yticks(fontsize=16)
  ylim(0,6e7)
  xlim(0,120)
  ylabel("total cases, recovered and susceptible",labelpad=15,fontsize=18)
  xlabel("time (days)",fontsize=18,labelpad=10)
  legend(loc="lower left",fancybox="true", fontsize=16)
  ax1 = ax.twinx()
  plot(t,I,"r-",label="infected", linewidth=2)
  plot(t,D,"m-",label="deaths", linewidth=2)
  yticks(fontsize=16)
  ylim(0,1.2e7)
  xlim(0,120)
  ylabel("infected and deaths",fontsize=18,labelpad=15)
  PyPlot.tight_layout()
  #savefig("../paper_figures/SIRD.png")

end


function generate_rnd()
# test loglog histograms with artificially generated random numbers
# sampled from a power-law distribution
    n = 100       # sample size
    mu = 1.25     # critical exponent
    C =  rand(n) .^ (-1/(mu-1))   # power-law distributed numbers
    C = sort(C, rev=true)
    nC = length(C)

    x1,y1,ex1,ey1 = loghistogram(log.(C), 20)

    pol = polyfit(x1,y1,1)
    muC=-coeffs(pol)[2]
    println("linear fit muC:",muC)
    ey1fit = exp.(polyval(pol,x1))
    alpha = 1 + length(C) / sum(log.(C)) # calculate exponent by log-likelihood
    println("alpha ", alpha)

    #######      # 2nd sample with different critical exponent
    mu = 2.0
    D =  rand(n) .^ (-1/(mu-1))
    D = sort(D, rev=true)
    nD = length(D)

    x1D,y1D,ex1D,ey1D = loghistogram(log.(D), 20)

    pol = polyfit(x1D,y1D,1)
    muD=-coeffs(pol)[2]
    println("linear fit muD:",muD)
    ey1fitD = exp.(polyval(pol,x1D))
    alphaD = 1 + length(D) / sum(log.(D)) # calculate exponent by log-likelihood
    println("alphaD ", alphaD)

   fig_Iy = 0.25
   fig_Ix = 0.63
   fig_Iwidth_y = 0.12
   fig_Iwidth_x = 0.26
   fig_sep_y   = 0.06
   fig_width_y = 0.4
   fig_x = 0.15
   fig_y2=0.08
   fig_y1=fig_y2 + fig_width_y + fig_sep_y

   figure("rnd", figsize=(8,12))
   clf()
   ax = PyPlot.axes([fig_x, fig_y1, 0.8, fig_width_y])
       scatter(ex1,ey1,s=50,label=L"random numbers, $\mu=1.25$")
       loglog(ex1,ey1fit,"k", label=L"regression line, $\mu=1.15$") 
       xlim(1,2e6)
       ylabel("frequency, P(n)",fontsize=16,labelpad=10)
       xticks(fontsize=14)
       yticks(fontsize=14)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=12)
       ax1 = ax.transAxes
       text(-0.15, 0.95, "a)", fontsize=18, fontweight="bold", transform = ax1)
       PyPlot.tight_layout()
   ax = PyPlot.axes([fig_Ix, fig_y1+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,7, length=50)
       loglog(cc, nC*cumulpdf.(cc,muC,1e6), color="black", linewidth=2)
       scatter(C,1:nC,s=30)
       xticks(fontsize=14)
       yticks(fontsize=14)
       xlim(1,1e7)
       ylim(1,200)
       xlabel("cases, n",fontsize=14)
       ylabel("C(n)",fontsize=14)
       PyPlot.tight_layout()

   ax = PyPlot.axes([fig_x, fig_y2, 0.8, fig_width_y])
       scatter(ex1D,ey1D,s=50,label=L"random numbers, $\mu=2.0$")
       loglog(ex1D,ey1fitD,"k", label=L"regression line, $\mu=1.63$") 
       xticks(fontsize=14)
       yticks(fontsize=14)
       xlim(1,2e3)
       ylim(5e-6,5e-1)
       xlabel("cases, n",fontsize=16)
       ylabel("frequency, P(n)",fontsize=16,labelpad=10)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=12)
       ax1 = ax.transAxes
       text(-0.15, 0.95, "b)", fontsize=18, fontweight="bold", transform = ax1)
       PyPlot.tight_layout()
   ax = PyPlot.axes([fig_Ix, fig_y2+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,3.5, length=50)
       loglog(cc, nD*cumulpdf.(cc,muD,1e3), color="black", linewidth=2)
       scatter(D,1:nD,s=30)
       myfont = Dict("fontname"=>"Sans","fontsize"=>14)
       xticks(fontsize=14)
       yticks(fontsize=14)
       xlim(1,5e3)
       ylim(1,100)
       xlabel("cases, n",fontsize=14)
       ylabel("C(n)",fontsize=14)
       PyPlot.tight_layout()
       #savefig("../paper_figures/rnd.png")

end

Random.seed!(0)

nlys_discrete_SIR()
#run_SIR()
#generate_rnd()


end   # module

