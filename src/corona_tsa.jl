# Time series analysis of 2020 worldwide coronavirus outbreak
#
# data based on John Hopkins CSSE repository
# https://github.com/CSSEGISandData/COVID-19

module C

using  Dates, DSP, PyPlot, DifferentialEquations, Polynomials

include("read_data.jl")

function smooth_log(x)
# make log-transform and smooth data
    lx = log10.(x)
    b = (1/2)*ones(2)     # flat (rectangular) window -> moving average
    lx = filtfilt(b,lx)
    10 .^ lx, lx
end

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
  slope, log(10.) .* slope
end

# doubling time when slope is taken from log10-transformed data
doubling_time(slope) = log(2.)/log(10.)/slope



function nlys()

   today = " March 19"
   
   t_de,I_de,R_de,D_de = covid_ts("Germany",       Date(2020,2,23))
   t_it,I_it,R_it,D_it = covid_ts("Italy",         Date(2020,2,20))
   t_es,I_es,R_es,D_es = covid_ts("Spain",         Date(2020,2,25))
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
   t_no,I_no,R_no,D_no = covid_ts("Norway",        Date(2020,2,28))
   t_ca,I_ca,R_ca,D_ca = covid_ts("Canada",        Date(2020,2,27))


  I_de, lI_de = smooth_log(I_de)
  I_it, lI_it = smooth_log(I_it)
  I_fr, lI_fr = smooth_log(I_fr)
  I_es, lI_es = smooth_log(I_es)
  I_jp, lI_jp = smooth_log(I_jp)
  I_cn, lI_cn = smooth_log(I_cn)
  I_kr, lI_kr = smooth_log(I_kr)
  I_ir, lI_ir = smooth_log(I_ir)
  I_us, lI_us = smooth_log(I_us)
  I_ch, lI_ch = smooth_log(I_ch)
  I_il, lI_il = smooth_log(I_il)
  I_at, lI_at = smooth_log(I_at)
  I_uk, lI_uk = smooth_log(I_uk)
  I_nl, lI_nl = smooth_log(I_nl)
  I_no, lI_no = smooth_log(I_no)
  I_ca, lI_ca = smooth_log(I_ca)

  di = 2 # length of running window: 2*di+1
  slope_de,growth_de = local_fit(t_de,lI_de,di)
  slope_it,growth_it = local_fit(t_it,lI_it,di)
  slope_es,growth_es = local_fit(t_es,lI_es,di)
  slope_fr,growth_fr = local_fit(t_fr,lI_fr,di)
  slope_jp,growth_jp = local_fit(t_jp,lI_jp,di)
  slope_cn,growth_cn = local_fit(t_cn,lI_cn,di)
  slope_kr,growth_kr = local_fit(t_kr,lI_kr,di)
  slope_ir,growth_ir = local_fit(t_ir,lI_ir,di)
  slope_us,growth_us = local_fit(t_us,lI_us,di)
  slope_ch,growth_ch = local_fit(t_ch,lI_ch,di)
  slope_at,growth_at = local_fit(t_at,lI_at,di)
  slope_il,growth_il = local_fit(t_il,lI_il,di)
  slope_uk,growth_uk = local_fit(t_uk,lI_uk,di)
  slope_nl,growth_nl = local_fit(t_nl,lI_nl,di)
  slope_no,growth_no = local_fit(t_no,lI_no,di)
  slope_ca,growth_ca = local_fit(t_ca,lI_ca,di)

  p_tot_de1 = polyfit(t_de[4:12],lI_de[4:12],1)  # fit full slope Germany
  coef = coeffs(p_tot_de1)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de1= polyval(p_tot_de1, t_de)

  p_tot_de2 = polyfit(t_de[12:end],lI_de[12:end],1)  # fit full slope Germany
  coef = coeffs(p_tot_de2)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de2= polyval(p_tot_de2, t_de)

  p_tot_it1 = polyfit(t_it[4:9],lI_it[4:9],1)  # fit slope 1st phase Italy
  coef = coeffs(p_tot_it1)[2]
  println("slope Italy Phase 1 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it1 = polyval(p_tot_it1, t_it)

  p_tot_it = polyfit(t_it[10:end],lI_it[10:end],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it)[2]
  println("slope Italy Phase 2 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it = polyval(p_tot_it, t_it)

  p_tot_es1 = polyfit(t_es[6:18],lI_es[6:18],1)  # fit full slope Spain
  lI_fit_es1= polyval(p_tot_es1, t_es)
  p_tot_es2 = polyfit(t_es[19:24],lI_es[19:24],1)  # fit full slope Spain
  lI_fit_es2= polyval(p_tot_es2, t_es)

  p_tot_fr1 = polyfit(t_fr[3:14],lI_fr[3:14],1)  # fit full slope France
  lI_fit_fr1= polyval(p_tot_fr1, t_fr)
  p_tot_fr2 = polyfit(t_fr[14:end],lI_fr[14:end],1)  # fit full slope France
  lI_fit_fr2= polyval(p_tot_fr2, t_fr)

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

  p_tot_ch = polyfit(t_ch[11:22],lI_ch[11:22],1)  # fit full slope Switzerland
  lI_fit_ch= polyval(p_tot_ch, t_ch)

  p_tot_at = polyfit(t_at,lI_at,1)  # fit full slope Austria
  lI_fit_at= polyval(p_tot_at, t_at)

  p_tot_il = polyfit(t_il,lI_il,1)  # fit full slope Israel
  lI_fit_il= polyval(p_tot_il, t_il)

  p_tot_uk = polyfit(t_uk,lI_uk,1)  # fit full slope UK
  lI_fit_uk= polyval(p_tot_uk, t_uk)

  p_tot_no = polyfit(t_no,lI_no,1)  # fit full slope Norway
  lI_fit_no= polyval(p_tot_no, t_no)

  p_tot_nl = polyfit(t_nl[15:23],lI_nl[15:23],1)  # fit full slope Netherlands
  lI_fit_nl= polyval(p_tot_nl, t_nl)

  p_tot_ca = polyfit(t_ca,lI_ca,1)  # fit full slope Norway
  lI_fit_ca= polyval(p_tot_ca, t_ca)

  #majorformatter = matplotlib.dates.DateFormatter("%d.%m.%Y")
  #majorlocator = matplotlib.dates.DayLocator(interval=1)

# ********************************************
  fig = figure("Outbreak Europe",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Germany
      plot(t_de,I_de,"r*")        
      ylabel("active cases",fontsize=12)
     
      title("Germany,"* today)
  ax=subplot(3,4,5)
     # p1 = plot_date(Date(2020,2,17),y,linestyle="-",marker="None",label="test")
      semilogy(t_de,I_de,"r*")        
      semilogy(t_de,10 .^ lI_fit_de1,"k")
      semilogy(t_de,10 .^ lI_fit_de2,"k")
      ax.yaxis.grid()
      ylim(5,5e4)
      ylabel("active cases",labelpad=20,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_de,growth_de,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      #ax.yaxis.grid(which="major", color="r", linestyle="-", linewidth=2)
      ax.yaxis.grid()
     # annotate("slope 0.107",xy=[20,2],
     # xycoords="axes fraction"
     # )
      #ax.xaxis.set_major_formatter(majorformatter)
      #fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      xlabel("time (days)",fontsize=12)
      ylabel("growth rate",labelpad=20,fontsize=12)

  subplot(3,4,2)          # Italy
      plot(t_it,I_it,"r*")        
      #ylabel("infected Italy")
      title("Italy,"*today)
  ax=subplot(3,4,6)
      semilogy(t_it,I_it,"r*")        
      semilogy(t_it,10 .^ lI_fit_it1,"k")
      semilogy(t_it,10 .^ lI_fit_it,"k")
      ylim(5,5e4)
      ax.yaxis.grid()
  ax=subplot(3,4,10)
      plot(t_it,growth_it,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("slope")

  subplot(3,4,3)              # Spain
      plot(t_es,I_es,"r*")      
      title("Spain,"*today)
  ax=subplot(3,4,7)
      semilogy(t_es,I_es,"r*")        
      semilogy(t_es,10 .^ lI_fit_es1,"k")
      semilogy(t_es,10 .^ lI_fit_es2,"k")
      ax.yaxis.grid()
      ylim(5,5e4)
   ax= subplot(3,4,11)
      plot(t_es,growth_es,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,4)              # France
      plot(t_fr,I_fr,"r*")      
      title("France,"*today)
  ax=subplot(3,4,8)
      semilogy(t_fr,I_fr,"r*")        
      semilogy(t_fr,10 .^ lI_fit_fr1,"k")
      semilogy(t_fr,10 .^ lI_fit_fr2,"k")
      ax.yaxis.grid()
      ylim(5,5e4)
  ax = subplot(3,4,12)
      plot(t_fr,growth_fr,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  PyPlot.tight_layout()
  savefig("../figures/outbreak_europe.png")

  
# ********************************************
  fig=figure("Outbreak Asia",figsize=(12,8))
  clf()
  subplot(3,4,1)         # China
      plot(t_cn,I_cn,"r*")        
      ylabel("active cases",fontsize=12)
      title("China,"*today)
  ax=subplot(3,4,5)
      semilogy(t_cn,I_cn,"r*")        
      #plot(t_cn,lI_fit_cn,"k")
      ylim(1e2,1e5)
      ax.yaxis.grid()
      ylabel("active cases", labelpad=20,fontsize=12)
  ax = subplot(3,4,9)
      plot(t_cn,growth_cn,"r*-")
      ylim(-0.1,0.5)
      #yticks([-0.05,0,0.05,0.1,0.15,0.2])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=10)

  subplot(3,4,2)          # Japan
      plot(t_jp,I_jp,"r*")        
      #ylabel("infected Italy")
      title("Japan,"*today)
  ax=subplot(3,4,6)
      semilogy(t_jp,I_jp,"r*")        
      semilogy(t_jp,10 .^ lI_fit_jp,"k")
      ax.yaxis.grid()
      ylim(5,2e4)
      #ylabel("log10 infected Italy")
  ax=subplot(3,4,10)
      plot(t_jp,growth_jp,"r*-")
      ylim(-0.1,0.5)
      #yticks([-0.05,0,0.05,0.1,0.15,0.2])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,3)              # South Korean
      plot(t_kr,I_kr,"r*")      
      #ylabel("infected Spain")
      title("South Korea,"*today)
  ax=subplot(3,4,7)
      semilogy(t_kr,I_kr,"r*")        
      ax.yaxis.grid()
      ylim(5,2e4)
      #plot(t_kr,lI_fit_kr,"k")
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_kr,growth_kr,"r*-")
      ylim(-0.1,0.5)
      #yticks([-0.05,0,0.05,0.1,0.15,0.2])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,4)              # Iran
      plot(t_ir,I_ir,"r*")      
      #ylabel("infected")
      title("Iran,"*today)
  ax=subplot(3,4,8)
      semilogy(t_ir,I_ir,"r*")        
      semilogy(t_ir,10 .^ lI_fit_ir,"k")
      ax.yaxis.grid()
      ylim(5,2e4)
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_ir,growth_ir,"r*-")
      ylim(-0.1,0.5)
      #yticks([-0.05,0,0.05,0.1,0.15,0.2])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")
  PyPlot.tight_layout()
  savefig("../figures/outbreak_Asia.png")

# ********************************************
  fig=figure("Outbreak World",figsize=(12,8))
  clf()
  subplot(3,4,1)         # USA
      plot(t_us,I_us,"r*")        
      ylabel("active cases",fontsize=12)
      title("USA,"*today)
  ax=subplot(3,4,5)
      semilogy(t_us,I_us,"r*")        
      semilogy(t_us,10 .^ lI_fit_us,"k")
      ax.yaxis.grid()
      ylim(5,2e4)
      ylabel("active cases",fontsize=12,labelpad=12)
  ax = subplot(3,4,9)
      plot(t_us,growth_us,"r*-")
      ax.yaxis.grid()
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=12)

  subplot(3,4,2)              # Canada
      plot(t_ca,I_ca,"r*")      
      #ylabel("infected")
      title("Canada,"* today)
  ax=subplot(3,4,6)
      semilogy(t_ca,I_ca,"r*")        
      semilogy(t_ca,10 .^ lI_fit_ca,"k")
      ax.yaxis.grid()
      ylim(5,2e4)
      #ylabel("log10 infected")
  ax = subplot(3,4,10)
      plot(t_ca,growth_ca,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,3)              # UK
      plot(t_uk,I_uk,"r*")      
      #ylabel("infected Spain")
      title("UK,"*today)
  ax=subplot(3,4,7)
      semilogy(t_uk,I_uk,"r*")        
      semilogy(t_uk,10 .^ lI_fit_uk,"k")
      ax.yaxis.grid()
      ylim(5,2e4)
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_uk,growth_uk,"r*-")
      ax.yaxis.grid()
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")
 

  subplot(3,4,4)         # Israel
      plot(t_il,I_il,"r*")        
      title("Israel,"*today)
  ax=subplot(3,4,8)
      semilogy(t_il,I_il,"r*")        
      semilogy(t_il,10 .^ lI_fit_il,"k")
      ax.yaxis.grid()
      ylim(1,1e3)
  ax = subplot(3,4,12)
      plot(t_il,growth_il,"r*-")
      ax.yaxis.grid()
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      xlabel("time (days)",fontsize=12)
      
  PyPlot.tight_layout()
  savefig("../figures/outbreak_World.png")





# ********************************************
  fig=figure("Outbreak Small Europe",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Netherlands
      plot(t_nl,I_nl,"r*")        
      ylabel("active cases",fontsize=12)
     
      title("Netherlands,"* today)
  ax=subplot(3,4,5)
      semilogy(t_nl,I_nl,"r*")        
      semilogy(t_nl,10 .^ lI_fit_nl,"k")
      ax.yaxis.grid()
      ylim(5,1e4)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_nl,growth_nl,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      #ax.yaxis.grid(which="major", color="r", linestyle="-", linewidth=2)
      ax.yaxis.grid()
     # annotate("growth 0.107",xy=[20,2],
     # xycoords="axes fraction"
     # )
      #ax.xaxis.set_major_formatter(majorformatter)
      #fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=20)

  subplot(3,4,2)          # Norway
      plot(t_no,I_no,"r*")        
      #ylabel("infected Italy")
      title("Norway,"*today)
  ax=subplot(3,4,6)
      semilogy(t_no,I_no,"r*")        
      semilogy(t_no,10 .^ lI_fit_no,"k")
      ax.yaxis.grid()
      ylim(5,1e4)
  ax=subplot(3,4,10)
      plot(t_no,growth_no,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,3)              # Switzerland
      plot(t_ch,I_ch,"r*")      
      #ylabel("infected Spain")
      title("Switzerland,"* today)
  ax=subplot(3,4,7)
      semilogy(t_ch,I_ch,"r*")        
      semilogy(t_ch,10 .^ lI_fit_ch,"k")
      ax.yaxis.grid()
      ylim(5,1e4)
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_ch,growth_ch,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,4)              # Austria
      plot(t_at,I_at,"r*")      
      #ylabel("infected")
      title("Austria,"* today)
  ax=subplot(3,4,8)
      semilogy(t_at,I_at,"r*")        
      semilogy(t_at,10 .^ lI_fit_at,"k")
      ax.yaxis.grid()
      ylim(5,1e4)
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_at,growth_at,"r*-")
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")
  PyPlot.tight_layout()
  savefig("../figures/outbreak_europe1.png")


# ********************************************

   lt_it = log10.(t_it)
   p_tot_it = polyfit(lt_it,lI_it,1)  # fit full slope Germany
   println("loglog It ", coeffs(p_tot_it)[2])
   I_fit_it = 10 .^ polyval(p_tot_it, lt_it)

   lt_de = log10.(t_de)
   p_tot_de = polyfit(lt_de,lI_de,1)  # fit full slope Germany
   println("loglog It ", coeffs(p_tot_de)[2])
   I_fit_de = 10 .^ polyval(p_tot_de, lt_de)
   
   figure("ts")
   clf()
   subplot(2,1,1)
   loglog(t_it,I_it,"*-")
   loglog(t_it,I_fit_it,"-")
   
   subplot(2,1,2)
   loglog(t_de,I_de,"*-")
   loglog(t_de,I_fit_de,"-")
   #semilogy(t_es,I_es,"r*-")

end

# **********************************************************
# **********************************************************
# **********************************************************


function nlys_ItDe()
  today = " March 18"
  t_it,I_it,R_it,D_it = covid_ts("Italy",Date(2020,2,21))
  I_it, lI_it = smooth_log(I_it)
  di = 2
  slope_it,growth_it = local_fit(t_it,lI_it,di)
  
  p_tot_it1 = polyfit(t_it[3:8],lI_it[3:8],1)  # fit slope 1st phase Italy
  coef = coeffs(p_tot_it1)[2]
  println("slope Italy Phase 1 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it1 = polyval(p_tot_it1, t_it)

  p_tot_it2 = polyfit(t_it[9:22],lI_it[9:22],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it2)[2]
  println("slope Italy Phase 2 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it2 = polyval(p_tot_it2, t_it)

  p_tot_it3 = polyfit(t_it[23:end],lI_it[23:end],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it3)[2]
  println("slope Italy Phase 3 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it3 = polyval(p_tot_it3, t_it)

  # ********************
  t_de,I_de,R_de,D_de = covid_ts("Germany",  Date(2020,2,24))
  t_de = t_de .+ 3 # bring to the same time axis
  I_de, lI_de = smooth_log(I_de)
  di = 2
  slope_de, growth_de = local_fit(t_de,lI_de,di)

  p_tot_de1 = polyfit(t_de[4:11],lI_de[4:11],1)  # fit full slope Germany
  coef = coeffs(p_tot_de1)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de1= polyval(p_tot_de1, t_de)

  p_tot_de2 = polyfit(t_de[12:end],lI_de[12:end],1)  # fit full slope Germany
  coef = coeffs(p_tot_de2)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de2= polyval(p_tot_de2, t_de)


   figure("Germany_Italy", figsize=(12,6))
   clf()
   subplot(1,2,1)
      semilogy(t_it[1:12],10 .^ lI_fit_it1[1:12],"k")
      semilogy(t_it[5:25],10 .^ lI_fit_it2[5:25],"k")
      semilogy(t_it[19:28],10 .^ lI_fit_it3[19:28],"k")
      semilogy(t_it,I_it,"b*", label="Italy", markersize=8)        
      semilogy(t_de[1:15],10 .^ lI_fit_de1[1:15],"k")
      semilogy(t_de[8:end],10 .^ lI_fit_de2[8:end],"k")
      semilogy(t_de,I_de,"r*",label="Germany",markersize=8)        
      #plot([10,10],[8,5e4],"grey",linewidth=1)
      #ylim(0.5,5)
      ylim(8,5e4)
      ylabel("active cases",fontsize=14)
      xlabel("time (Feb 21 - March 19)", fontsize=14)
      legend(loc="upper left",fancybox="true")
   ax = subplot(1,2,2)
      plot(t_it,growth_it,"b*-",label="Italy",markersize=8)
      plot(t_de,growth_de,"r*-",label="Germany",markersize=8)
      #plot([10,10],[0,0.6],"grey",linewidth=1)
      #ylim(0,0.22)
      #yticks([0,0.05,0.1,0.15,0.2])
      ylim(0,0.6)
      #yticks([0,0.05,0.1,0.15,0.2])
      ax.yaxis.grid()
      #xlabel("time (days)",fontsize=14)
      xlabel("time (Feb 21 - March 18)", fontsize=14)
      ylabel("growth rate",fontsize=14)
      legend(loc="upper right",fancybox="true")
  PyPlot.tight_layout()
  savefig("../figures/ItalyGermany_March19.png")


end

#country = "Italy"; date = Date(2020,2,20)

#nlys_country(country, date)
nlys_ItDe()

nlys()

end   # module

