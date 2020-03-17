# Time series analysis of 2020 worldwide coronavirus outbreak
#
# data based on John Hopkins CSSE repository
# https://github.com/CSSEGISandData/COVID-19

module C

using  Dates, PyPlot, DifferentialEquations, Polynomials

include("read_data.jl")


function local_fit(x,y,di)
# local slope of graph y(x) from regression in running window
# window width: 2*di+1
  len = length(x)
  slope = zeros(len)
  for i=1:len
    i1 = max(1,i-di)
    i2 = min(len,i+di)
    p = polyfit(x[i1:i2],y[i1:i2],1)
    slope[i] = coeffs(p)[2]
  end
  slope
end

# doubling time when slope is taken from log10-transformed data
doubling_time(slope) = log(2.)/log(10.)/slope



function nlys()
   
   t_de,I_de,R_de,D_de = covid_ts("Germany",       Date(2020,2,23))
   t_it,I_it,R_it,D_it = covid_ts("Italy",         Date(2020,2,20))
   t_es,I_es,R_es,D_es = covid_ts("Spain",         Date(2020,2,24))
   t_fr,I_fr,R_fr,D_fr = covid_ts("France",        Date(2020,2,25))
   t_uk,I_uk,R_uk,D_uk = covid_ts("United Kingdom",Date(2020,2,27))
   t_at,I_at,R_at,D_at = covid_ts("Austria",       Date(2020,2,29))
   t_nl,I_nl,R_nl,D_nl = covid_ts("Netherlands",   Date(2020,2,23))
   t_ch,I_ch,R_ch,D_ch = covid_ts("Switzerland",   Date(2020,2,25))
   t_cn,I_cn,R_cn,D_cn = covid_ts("China",         Date(2020,1,22))
   t_jp,I_jp,R_jp,D_jp = covid_ts("Japan",         Date(2020,1,24))
   t_kr,I_kr,R_kr,D_kr = covid_ts("Korea, South",  Date(2020,2,17))
   t_us,I_us,R_us,D_us = covid_ts("US",            Date(2020,2,22))
   t_ir,I_ir,R_ir,D_ir = covid_ts("Iran",          Date(2020,2,19))
   t_il,I_il,R_il,D_il = covid_ts("Israel",        Date(2020,2,27))


  lI_de = log10.(I_de)
  lI_it = log10.(I_it)
  lI_fr = log10.(I_fr)
  lI_es = log10.(I_es)
  lI_jp = log10.(I_jp)
  lI_cn = log10.(I_cn)
  lI_kr = log10.(I_kr)
  lI_ir = log10.(I_ir)
  lI_us = log10.(I_us)
  lI_ch = log10.(I_ch)
  lI_il = log10.(I_il)
  lI_at = log10.(I_at)
  lI_uk = log10.(I_uk)
  lI_nl = log10.(I_nl)

  slope_de = local_fit(t_de,lI_de,2)
  slope_it = local_fit(t_it,lI_it,2)
  slope_es = local_fit(t_es,lI_es,2)
  slope_fr = local_fit(t_fr,lI_fr,2)
  slope_jp = local_fit(t_jp,lI_jp,2)
  slope_cn = local_fit(t_cn,lI_cn,2)
  slope_kr = local_fit(t_kr,lI_kr,2)
  slope_ir = local_fit(t_ir,lI_ir,2)
  show(lI_us)
  println(length(lI_us))
  show(t_us)
  println()
  
  slope_us = local_fit(t_us,lI_us,2)
  slope_ch = local_fit(t_ch,lI_ch,2)
  slope_at = local_fit(t_at,lI_at,2)
  slope_il = local_fit(t_il,lI_il,2)
  slope_uk = local_fit(t_uk,lI_uk,2)
  slope_nl = local_fit(t_nl,lI_nl,2)

  p_tot_de = polyfit(t_de,lI_de,1)  # fit full slope Germany
  coef = coeffs(p_tot_de)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de= polyval(p_tot_de, t_de)

  p_tot_it1 = polyfit(t_it[4:9],lI_it[4:9],1)  # fit slope 1st phase Italy
  coef = coeffs(p_tot_it1)[2]
  println("slope Italy Phase 1 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it1 = polyval(p_tot_it1, t_it)

  p_tot_it = polyfit(t_it[10:end],lI_it[10:end],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it)[2]
  println("slope Italy Phase 2 ",coef, " T12 ", doubling_time(coef))
  len = length(lI_it)
  lI_fit_it = polyval(p_tot_it, t_it)

  p_tot_es = polyfit(t_es,lI_es,1)  # fit full slope Spain
  lI_fit_es= polyval(p_tot_es, t_es)

  p_tot_fr = polyfit(t_fr,lI_fr,1)  # fit full slope France
  lI_fit_fr= polyval(p_tot_fr, t_fr)

  p_tot_cn = polyfit(t_cn,lI_cn,1)  # fit full slope China
  lI_fit_cn= polyval(p_tot_cn, t_cn)

  p_tot_kr = polyfit(t_kr,lI_kr,1)  # fit full slope South Korean
  lI_fit_kr= polyval(p_tot_kr, t_kr)

  p_tot_jp = polyfit(t_jp,lI_jp,1)  # fit full slope Japan
  lI_fit_jp= polyval(p_tot_jp, t_jp)

  p_tot_ir = polyfit(t_ir,lI_ir,1)  # fit full slope Iran
  lI_fit_ir= polyval(p_tot_ir, t_ir)

  p_tot_us = polyfit(t_us,lI_us,1)  # fit full slope USA
  lI_fit_us= polyval(p_tot_us, t_us)

  p_tot_ch = polyfit(t_ch,lI_ch,1)  # fit full slope Switzerland
  lI_fit_ch= polyval(p_tot_ch, t_ch)

  p_tot_at = polyfit(t_at,lI_at,1)  # fit full slope Austria
  lI_fit_at= polyval(p_tot_at, t_at)

  p_tot_il = polyfit(t_il,lI_il,1)  # fit full slope Israel
  lI_fit_il= polyval(p_tot_il, t_il)

  p_tot_uk = polyfit(t_uk,lI_uk,1)  # fit full slope UK
  lI_fit_uk= polyval(p_tot_uk, t_uk)


# ********************************************
  figure("Outbreak Europe1")
  clf()
  subplot(3,4,1)         # Germany
      plot(t_de,I_de,"r*")        
      ylabel("infected")
      title("Germany")
  subplot(3,4,5)
      plot(t_de,lI_de,"r*")        
      plot(t_de,lI_fit_de,"k")
      ylabel("log10 infected")
  ax = subplot(3,4,9)
      plot(t_de,slope_de,"r*-")
      #ax.yaxis.grid(which="major", color="r", linestyle="-", linewidth=2)
      ax.yaxis.grid()
      xlabel("time (days)")
      ylabel("slope")

  subplot(3,4,2)          # Italy
      plot(t_it,I_it,"r*")        
      #ylabel("infected Italy")
      title("Italy")
  subplot(3,4,6)
      plot(t_it,log10.(I_it),"r*")        
      plot(t_it,lI_fit_it,"k")
      plot(t_it,lI_fit_it1,"k")
      #ylabel("log10 infected Italy")
  ax=subplot(3,4,10)
      plot(t_it,slope_it,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")

  subplot(3,4,3)              # Spain
      plot(t_es,I_es,"r*")      
      #ylabel("infected Spain")
      title("Spain")
  subplot(3,4,7)
      plot(t_es,lI_es,"r*")        
      plot(t_es,lI_fit_es,"k")
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_es,slope_es,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")

  subplot(3,4,4)              # France
      plot(t_fr,I_fr,"r*")      
      #ylabel("infected")
      title("France")
  subplot(3,4,8)
      plot(t_fr,lI_fr,"r*")        
      plot(t_fr,lI_fit_fr,"k")
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_fr,slope_fr,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")

  #savefig("outbreak_europe.png")

  
# ********************************************
  figure("Outbreak Asia1")
  clf()
  subplot(3,4,1)         # China
      plot(t_cn,I_cn,"r*")        
      ylabel("infected")
      title("China")
  subplot(3,4,5)
      plot(t_cn,lI_cn,"r*")        
      plot(t_cn,lI_fit_cn,"k")
      ylabel("log10 infected")
  ax = subplot(3,4,9)
      plot(t_cn,slope_cn,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      ylabel("slope")

  subplot(3,4,2)          # Japan
      plot(t_jp,I_jp,"r*")        
      #ylabel("infected Italy")
      title("Japan")
  subplot(3,4,6)
      plot(t_jp,log10.(I_jp),"r*")        
      plot(t_jp,lI_fit_jp,"k")
      #ylabel("log10 infected Italy")
  ax=subplot(3,4,10)
      plot(t_jp,slope_jp,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")

  subplot(3,4,3)              # South Korean
      plot(t_kr,I_kr,"r*")      
      #ylabel("infected Spain")
      title("South Korea")
  subplot(3,4,7)
      plot(t_kr,lI_kr,"r*")        
      plot(t_kr,lI_fit_kr,"k")
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_kr,slope_kr,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")

  subplot(3,4,4)              # Iran
      plot(t_ir,I_ir,"r*")      
      #ylabel("infected")
      title("Iran")
  subplot(3,4,8)
      plot(t_ir,lI_ir,"r*")        
      plot(t_ir,lI_fit_ir,"k")
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_ir,slope_ir,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")

  #savefig("outbreak_europe.png")

# ********************************************
  figure("Outbreak World1")
  clf()
  subplot(3,4,1)         # USA
      plot(t_us,I_us,"r*")        
      ylabel("infected")
      title("USA")
  subplot(3,4,5)
      plot(t_us,lI_us,"r*")        
      plot(t_us,lI_fit_us,"k")
      ylabel("log10 infected")
  ax = subplot(3,4,9)
      plot(t_us,slope_us,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      ylabel("slope")

  subplot(3,4,2)         # Israel
      plot(t_il,I_il,"r*")        
      ylabel("infected")
      title("Israel")
  subplot(3,4,6)
      plot(t_il,lI_il,"r*")        
      plot(t_il,lI_fit_il,"k")
      ylabel("log10 infected")
  ax = subplot(3,4,10)
      plot(t_il,slope_il,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      ylabel("slope")

  subplot(3,4,3)              # UK
      plot(t_uk,I_uk,"r*")      
      #ylabel("infected Spain")
      title("UK")
  subplot(3,4,7)
      plot(t_uk,lI_uk,"r*")        
      plot(t_uk,lI_fit_uk,"k")
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_uk,slope_uk,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")





# ********************************************
  figure("Outbreak Small Europe1")
  clf()

  subplot(3,4,3)              # Switzerland
      plot(t_ch,I_ch,"r*")      
      #ylabel("infected Spain")
      title("Switzerland")
  subplot(3,4,7)
      plot(t_ch,lI_ch,"r*")        
      plot(t_ch,lI_fit_ch,"k")
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_ch,slope_ch,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")

  subplot(3,4,4)              # Austria
      plot(t_at,I_at,"r*")      
      #ylabel("infected")
      title("Austria")
  subplot(3,4,8)
      plot(t_at,lI_at,"r*")        
      plot(t_at,lI_fit_at,"k")
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_at,slope_at,"r*-")
      ax.yaxis.grid()
      xlabel("time (days)")
      #ylabel("slope")


   figure("ts")
   clf()
   subplot(2,1,1)
   plot(t_us,I_us,"*-")
   
   subplot(2,1,2)
   semilogy(t_us,I_us,"*-")
   semilogy(t_us,R_us)
   semilogy(t_us,D_us)

end


nlys()

end   # module

