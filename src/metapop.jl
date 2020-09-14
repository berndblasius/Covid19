# Poower-law distributions of infections in a simple metapopulation model
# this corresponds to Fig.3 in 
# Blasius B. (2020) Power-law distribution in the number of confirmed COVID-19 cases. Chaos
# https://doi.org/10.1063/5.0013031


module C

using Dates, LsqFit, PyPlot, Random, StatsBase, DifferentialEquations, Polynomials

include("tools.jl")


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

    ndays = 75    # total runtime

    ndays = 80    # total runtime
    nB = 62   # runtime until states are saved
    nB1 = 80   # runtime until states are saved
    nB2 = 65   # runtime until states are saved
    nB3 = 50   # runtime until states are saved
    nB4 = 35   # runtime until states are saved
    nB5 = 20   # runtime until states are saved
    nB6 = 5   # runtime until states are saved


    N = 200       # max number of countries
    p = 0.0005    # invasion probability
    p = 0.0006    # invasion probability
    firsts = [0]

    # local growth rate
    r = b - gamma
    C_max = exp(r*(ndays+1))
    println("Cmax: ",C_max)

    alphaI = zeros(ndays)
    alphaD = zeros(ndays)
    C = zeros(0)
    D = zeros(0)
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
       Rv = tmp[3,:] # vector of recovered
       Dv = tmp[4,:] # vector of dead
       Cv = Iv .+ Rv .+ Dv   # all cases
       if i == nB
          C = Cv
          D = Dv
       end
       if i == nB1
           println("here")
          C1 = Cv
          D1 = Dv
       elseif i == nB2
          C2 = Cv
          D2 = Dv
       elseif i == nB3
          C3 = Cv
          D3 = Dv
       elseif i == nB4
          C4 = Cv
          D4 = Dv
       elseif i == nB5
          C5 = Cv
          D5 = Dv
       elseif i == nB6
          C6 = Cv
          D6 = Dv
       end
       ind = findall(x->x>=1,Dv)
       Dv = Dv[ind]
       ind = findall(x->x>=1,Cv)
       Cv = Cv[ind]
       if length(Cv) >= 5
          mu,st = estimate(Cv)
          alphaI[i] = mu
       end
       if length(Dv) >= 5
          mu,st = estimate(Dv)
          alphaD[i] = mu
       end
    end


   # ************ distribution   ******************************** 
   ind = findall(x->x>=0.5, C)
   C = sort(C[ind], rev=true)
   nC = length(C)
   println("maximum C: ",maximum(C))
   ex1,ey1 = log2histogram(C)
   x1 = log.(ex1)
   y1 = log.(ey1)
   pol = Polynomials.fit(x1,y1,1)
   muC=-coeffs(pol)[2]
   println("linear fit muC:",muC)
   pol = Polynomials.fit(x1,y1,1)
   ey1fit = exp.(pol.(x1))
   #alpha = 1 + length(C) / sum(log.(C)) # calculate exponent by log-likelihood
   #println("alpha ", alphaI[62])
   
   
   println("length D ", length( findall(x->x>=1,D )))
   ind = findall(x->x>=0.5, D)
   D = sort(D[ind], rev=true)
   nD = length(D)
   ex1D,ey1D = log2histogram(D)
   x1D = log.(ex1D)
   y1D = log.(ey1D)
   pol = Polynomials.fit(x1D,y1D,1)
   muD=-coeffs(pol)[2]
   println("linear fit muD: ",muD)
   pol = Polynomials.fit(x1D,y1D,1)
   ey1fitD = exp.(pol.(x1D))

    C1 = sort(C1[findall(x->x>0.5, C1)],rev=true)
    exC1,eyC1 = log2histogram(C1)
    nC1 = length(C1)
 
    C2 = sort(C2[findall(x->x>0.5, C2)],rev=true)
    exC2,eyC2 = log2histogram(C2)
    nC2 = length(C2)

    C3 = sort(C3[findall(x->x>0.5, C3)],rev=true)
    exC3,eyC3 = log2histogram(C3)
    nC3 = length(C3)
    
    C4 = sort(C4[findall(x->x>0.5, C4)],rev=true)
    exC4,eyC4 = log2histogram(C4)
    nC4 = length(C4)

    C5 = sort(C5[findall(x->x>0.5, C5)],rev=true)
    exC5,eyC5 = log2histogram(C5)
    nC5 = length(C5)
    
    #  *********** Arrival time histogram  ***********************
    println(N1, " countries infected")
    h = StatsBase.fit(Histogram, firsts,  nbins=100, closed=:left)
    h_weight = h.weights
    h_bins = h.edges[1][1:end-1]

    @. model(x, p) = p[1]*exp(x*p[2])
    p0 = [0.5, 0.5]
    expfit = curve_fit(model, collect(h_bins[1:62]), h_weight[1:62], p0)
    println("Spreading rate: ",coef(expfit))
    cf = coef(expfit)


   ax = PyPlot.axes([fig_x2, fig_y1, fig_width_x, fig_width_y])
       loglog(exC1,eyC1,"-", color="k")
       scatter(exC1,eyC1,s=50, color="k", label="arrival time: $nB1 days")

       loglog(exC2,eyC2,"-", color="r")
       scatter(exC2,eyC2,s=50, color="r", label="arrival time: $nB2 days")

       loglog(exC3,eyC3,"-", color="b")
       scatter(exC3,eyC3,s=50, color="b", label="arrival time: $nB3 days")
       loglog(exC4,eyC4,"-", color="g")
       scatter(exC4,eyC4,s=50, color="g", label="arrival time: $nB4 days")
       loglog(exC5,eyC5,"-", color="m")
       scatter(exC5,eyC5,s=50, color="m", label="arrival time: $nB5 days")
       legend(loc="lower left",fancybox="true", fontsize=12)
       xticks(fontsize=14)
       yticks([1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],fontsize=14)
       xlim(0.8,1e8)
       ylim(5e-9,1)
       xlabel("cases, n",fontsize=16)
       ax1 = ax.transAxes
       text(-0.14, 0.95, "c)", fontsize=20, fontweight="bold", transform = ax1)

   ax = PyPlot.axes([fig_Ix3, fig_y1+fig_Iy, fig_Iwidth_x, fig_Iwidth_y]) # Inset: cumulative distribution 
       loglog(C1,(1:nC1),"-", color="k")
       loglog(C2,(1:nC2),"-", color="r")
       loglog(C3,(1:nC3),"-", color="b")
       loglog(C4,(1:nC4),"-", color="g")
       loglog(C5,(1:nC5),"-", color="m")
       ylim(1,200)
       xlim(1,1e8)
       xticks(fontsize=12)
       yticks(fontsize=12)
       xlabel("cases, n", fontsize=14)
       ylabel(L"$N \cdot C(n)$", fontsize=14, labelpad=0)


   ax = PyPlot.axes([fig_x2, fig_y2, fig_width_x, fig_width_y])
       bar(h_bins,h_weight)
       plot(h_bins[1:63], cf[1]*exp.(cf[2]*h_bins[1:63]),"k", linewidth=2, linestyle="dashed")
       plot([62,62],[0,14],"-r",linewidth=2)
       ylabel("countries",fontsize=16, labelpad=5)
       xticks(fontsize=14)
       yticks(0:2:12,fontsize=14)
       xlim(0,ndays)
       ylim(0,14)
       xlabel("arrival time (days)",fontsize=16)
       ax1 = ax.transAxes
       text(-0.14, 0.95, "d)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix2, fig_y2+fig_Iy, fig_Iwidth_x1, fig_Iwidth_y1]) # Inset: cumulative distribution 
       ind = findfirst(x->x>0, alphaI)
       plot(ind:ndays,alphaI[ind:ndays], linewidth=2, color="black")
       ind = findfirst(x->x>0, alphaD)
       plot(ind:ndays,alphaD[ind:ndays], linewidth=2, color="blue")
       plot([62,62],[1,2],"r-")
       xticks([0,20,40,60,80],fontsize=12)
       yticks([1,1.2,1.4,1.6,1.8,2.0],fontsize=12)
       xlim(10,ndays)
       ylim(1,1.8)
       ax.yaxis.grid()
       xlabel("time (days)", fontsize=14)
       ylabel(L"exponent, $\mu$", fontsize=14, labelpad=5)

   # **************************************************************************
   # plot distribution
   ax = PyPlot.axes([fig_x1, fig_y1, fig_width_x, fig_width_y])
       loglog(ex1,ey1,"-",linewidth=2)
       scatter(ex1,ey1,s=50,label="total detected cases")
       loglog(ex1,ey1fit,"k", label=L"regression line, $\mu=1.18$") 
       xlim(0.8,1e7)
       xlabel("cases, n",fontsize=16)
       ylabel(L"frequency, $P_C(n)$",fontsize=16,labelpad=10)
       xticks(fontsize=14)
       yticks([1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],fontsize=14)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       ax1 = ax.transAxes
       text(-0.2, 0.95, "a)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix1, fig_y1+fig_Iy, fig_Iwidth_x1, fig_Iwidth_y1]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,7, length=50)
       loglog(cc, cumulpdf.(cc,muC,5e6), color="black", linewidth=2)
       scatter(C,(1:nC)./nC,s=20)
       xticks(fontsize=12)
       yticks(fontsize=12)
       xlim(1,1e7)
       ylim(5e-3,1)
       xlabel("cases, n",fontsize=14)
       ylabel("C(n)",fontsize=14)

   ax = PyPlot.axes([fig_x1, fig_y2, fig_width_x, fig_width_y])
       loglog(ex1D,ey1D,"-",linewidth=2)
       scatter(ex1D,ey1D,s=50,label="deaths")
       loglog(ex1D,ey1fitD,"k", label=L"regression line, $\mu=1.31$") 
       xticks(fontsize=14)
       yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],fontsize=14)
       xlim(0.8,5e4)
       xlabel("cases, n",fontsize=16)
       ylabel(L"frequency, $P_D(n)$",fontsize=16,labelpad=10)
       legend(bbox_to_anchor=(0.01, 0.25), loc="upper left", fancybox="true", fontsize=14)
       ax1 = ax.transAxes
       text(-0.2, 0.95, "b)", fontsize=20, fontweight="bold", transform = ax1)
   ax = PyPlot.axes([fig_Ix1, fig_y2+fig_Iy, fig_Iwidth_x1, fig_Iwidth_y1]) # Inset: cumulative distribution 
       cc = 10 .^ range(0,4.2, length=50)
       loglog(cc, cumulpdf.(cc,muD,5e4), color="black", linewidth=2)
       scatter(D,(1:nD)/nD,s=20)
       xticks(fontsize=12)
       yticks(fontsize=12)
       xlim(1,6e4)
       ylim(1e-2,1)
       xlabel("cases, n",fontsize=14)
       ylabel("C(n)",fontsize=14)



end

Random.seed!(12)

fig_Iy = 0.25
fig_Ix1 = 0.34
fig_Ix3 = 0.83
fig_Ix2 = 0.63
fig_Iwidth_y  = 0.12
fig_Iwidth_x  = 0.13
fig_Iwidth_y1 = 0.13
fig_Iwidth_x1 = 0.15

fig_sep_y   = 0.06
fig_sep_y   = 0.08
fig_width_y = 0.4
fig_width_x = 0.4
fig_x1 = 0.1
fig_x2 = 0.57
fig_y2=0.08
fig_y1=fig_y2 + fig_width_y + fig_sep_y

figure("distribution_meta1", figsize=(12,8))
clf()
nlys_discrete_SIR()
#savefig("../figures/meta_all.png")


end # module



