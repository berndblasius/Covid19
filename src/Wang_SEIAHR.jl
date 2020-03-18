# SEIR-like model for Corona outbreak 2019/20
#
# based on Wang et al. Evolving epidemiology and impact of non-pharmaceutical
# interventions on the outbreak of coronavirus disease 2019 in Wuhan, China.
# MedRxive (2020) doi.org/10.1101/2020.03.03.2003059

module C

using Dates, PyPlot, DifferentialEquations, Polynomials

include("read_data.jl")


function rhs(du,u,p,t)
  S,E,I,A,H,R = u
  b,alpha,n,De,Dq,Di,Dh,r = p
  N = S+E+I+A+H+R
  du[1] = dS = -b*S*(I+alpha*A)/N + n - n*S/(N-I-H)   # susceptible
  du[2] = dE = b*S*(I+alpha*A)/N - E/De - n*E/(N-I-H) # latent
  du[3] = dI = r*E/De - I/Dq - I/Di                   # reported infectious
  du[4] = dA = (1-r)*E/De - A/Di - n*A/(N-I-H)        # unreported infectious
  du[5] = dH = I/Dq - H/Dh                            # hospitalized
  du[6] = dR = (I+A)/Di + H/Dh - n*R/(N-I-H)          # recovered
end


function nlys(scenario)
  tmax = 100.0  # all time units in days

  # The Wang et al. study changes parameters withing 4 periods
  # period 1: Jan 1-10   # before Chunyun (Chinese New Year)
  # period 2: Jan 11-22  # Chunyun period massive population movement
                         # hospitals starting to be crowded
  # period 3: Jan 23 - Feb 1 # implementation of social distancing
  # period 4: Feb 2-18   # improvement in medical resources

  b =  1.75  # transmission rate [1.75, 1.75, 0.58, 0.15]
  alpha = 1  # transmission rate ratio of unacertained ov acertained cases
  r =  0.19  # acertainment rate [0.19, 0.19, 0.22, 0.17]
  De =  5.2  # latent period
  Di =  2.3  # infectious period
  Dq =  10   # time from onset of illness to hospitalization [10, 7, 5, 2]
  Dh =  30   # hospitalization period
  n  = 0     # daily inbound and outbound size [500000, 800000, 0 , 0]
 
  if scenario == :wuhan
      N  = 10 * 1e6  # population size Wuhan
      E0 = 346
      I0 = 80
      A0 = 80
      H0 = 27
  elseif scenario == :de  
      N  = 82 * 1e6  # population size Germany
      E0 = 500   # this is just put in by hand to make a smooth start
      I0 = 16
      A0 = 16    # intially equal number of reported and unreported 
      H0 = 0     # assume initially there are no hospitalized 

      r = 0.22 # slightly increase the acertainment rate 
      #Di =  4.2  # infectious period
  elseif scenario == :it
      N  = 60 * 1e6  # population size Italy
      E0 = 1000
      I0 = 20
      A0 = 200
      H0 = 0
      b = 2.0 # to fit the data 
      #r = 0.32 # or we increase the acertainment rate a lot
  elseif scenario == :us
      N = 327 * 1e6  # population size US
      E0 = 346
      I0 = 80
      A0 = 80
      H0 = 27
  elseif scenario == :OL   # try the model on a city-scale
      N  = 170000.0 # population size of Oldenburg (my city in Germany)
      E0 = 50
      I0 = 5
      A0 = 5
      H0 = 0
      #b = 2.0 # to fit the data in the first 10 days in Germany,
      r = 0.32 # or we increase the acertainment rate a lot
  end
  R0 = 0
  S0 = N - E0 - I0 - A0 - H0 - R0


  u0 = [S0;E0;I0;A0;H0;R0]
  tspan = (0.0,tmax)
  p = (b,alpha,n,De,Dq,Di,Dh,r)
  prob = ODEProblem(rhs,u0,tspan,p)
  sol = solve(prob)
  t = 0:0.1:tmax  # use interpolation of results
  res = sol(t)  
  S = res[1,:]
  E = res[2,:]
  I = res[3,:]
  A = res[4,:]
  H = res[5,:]
  R = res[6,:]

  # effective reproductive number  (formula from Wang et al)
  R0 = Di .* b ./ (A .+ I) .* (alpha .* A .+ Dq .* I ./ (Di .+ Dq))  

  t_de,I_de,R_de,D_de = covid_ts("Germany", Date(2020,2,24))
  t_it,I_it,R_it,D_it = covid_ts("Italy",   Date(2020,2,21))


  ## Plots ****************************************************

  if scenario == :OL  # this scenario calculates the epidemic course for a small city
      figure("ts_OL")
      clf()
      subplot(3,1,1)
      plot(t,I)
      plot(t,I .+ A, "r")
      ylabel("infected")
      xlim(0,50)
      subplot(3,1,2)
      plot(t, log10.(I))
      plot(t, log10.(I .+ A),"r")
      ylabel("log10(infected)")
      xlim(0,50)
      subplot(3,1,3)
      plot(t,R)  
      ylabel("recovered")
      xlabel("time (days)")
  end

  if scenario == :de    # simulate scenario for Germany
      figure("ts_ge")
      clf()

      subplot(3,2,1)  # show the whole outbreak
      #plot(t_it .+ 6 ,I_it,"g*")  # time delay Italy cases
      plot(t_de,I_de,"r*")        # German cases
      plot(t,I)                   # model for acertained infected
      plot(t,A .+ I .+ H, "g")    # model of real (including non-acerted) infected
      xlim(0,60)
      ylabel("infected")

      subplot(3,2,2)   # show only the first weeks
      #plot(t_it .+ 6 ,I_it,"g*")  # time delay Italy cases
      plot(t_de,I_de,"r*")
      plot(t,I)
      xlim(0,30)     
      ylim(0,10000)
      ylabel("infected, inital period")

      subplot(3,2,3)  # whole outbreak on a half-logarithmic scale
      #plot(t_it .+ 6,log10.(I_it),"g*")
      plot(t_de,log10.(I_de),"r*")
      plot(t, log10.(I))
      plot(t, log10.(A), "g")
      ylabel("log10(infected)")
      xlim(0,60)

      subplot(3,2,4)  # first weeks on a half-logarithmic scale
      #plot(t_it .+ 6,log10.(I_it),"g*")
      plot(t_de,log10.(I_de),"r*")
      plot(t, log10.(I))
      xlabel("time (days)")
      ylabel("log10(infected)")
      xlim(0,30)
      ylim(0,4.5)

      subplot(3,2,5)  # other model stages
      plot(t,E)
      plot(t,R)
      plot(t_de,R_de,"r*")
      plot(t,H)
      plot(t,S)
      ylabel("Other states")
      xlim(0,60)
      xlabel("time (days)")

      subplot(3,2,6)  # other model stages
      plot(t,E)
      plot(t,R,"r")
      plot(t_de,R_de,"r*")
      plot(t,H)
      #plot(t,S)
      ylabel("Other states")
      xlim(0,30)
      ylim(0,100)
      xlabel("time (days)")


      #subplot(3,2,6)
      #plot(t,R0)  # basic reproductive number according to the formula of Wang et al.
      #xlim(0,60)
      #ylabel("R0")
      #xlabel("time (days)")

      #savefig("outbreak_Germany.png")



  end


  if scenario == :it   # analysis for Italy (it seems the is a slight reduction
                       # in growth rates after first week)
      figure("ts_it")
      clf()

      subplot(2,2,1)
      #plot(1:14,I_Italy,"b*")
      plot(t_it,I_it,"r*")
      plot(t,I)
      plot(t,A, "g")
      xlim(0,80)
      #ylim(0,2000)
      #plot(t,Rt)
      ylabel("infected")
      #xlim(0,10)

      subplot(2,2,2)
      #plot(1:14,I_Italy,"b*")
      plot(t_it,I_it,"r*")
      plot(t,I)
      xlim(0,20)
      ylim(0,15000)
      ylabel("infected, inital period")
      subplot(2,2,3)
      #plot(1:14,log10.(I_Italy),"b*")
      plot(t_it,log10.(I_it),"r*")
      plot(t, log10.(I))
      plot(t, log10.(A), "g")
      xlabel("time (days)")
      ylabel("log10(infected)")
      #xlim(0,50)

      subplot(2,2,4)
      plot(t,E)
      plot(t,R)
      plot(t,H)
      plot(t,S)
      xlim(0,80)
  end

end


#scenario = :wuhan # Wuhan
#scenario = :de    # Germany
#scenario = :it    # Italy
#scenario = :us    # US
#scenario = :OL    # Oldenburg (my city in Germany)
nlys(:de)


end   # module

