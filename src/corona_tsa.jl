# Time series analysis of 2020 worldwide coronavirus outbreak
#
# At first data were taken from Johns Hopkins CSSE repository
# https://github.com/CSSEGISandData/COVID-19
#
# but since this data source for non-understandable reasons on March 23 totally
# changed their data format, without reporting confirmed cases any more
# I changed to datahub 
# https://github.com/datasets/covid-19/tree/master/data.

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
    #p = polyfit(x[i1:i2],y[i1:i2],1)
    p = fit(x[i1:i2],y[i1:i2],1)
    cp = coeffs(p)
    if length(cp) == 1
        slope[i] = 0.0
    else
      slope[i] = coeffs(p)[2]
    end
  end
  slope, log(10.) .* slope
end

# doubling time when slope is taken from log10-transformed data
doubling_time(slope) = log(2.)/log(10.)/slope



# **********************************************************
function nlys_Asia()
   di = 3

   t_cn,I_cn,R_cn,D_cn = covid_ts("China",         Date(2020,1,22))
   #(D_cn, C_cn, R_cn, dates_cn) =  read_datahub("China",Date(2020,1,22))
   #println(R_cn)
   #t_cn = 1:length(C_cn)
   #I_cn = C_cn #  .- R_cn .- D_cn
   I_cn, lI_cn = smooth_log(I_cn)
   #println(t_cn)
   #println(lI_cn)
   slope_cn,growth_cn = local_fit(t_cn,lI_cn,di)

   di = 2
   t_jp,I_jp,R_jp,D_jp = covid_ts("Japan",         Date(2020,1,24))
   #(D_jp, C_jp, R_jp, dates_jp) =  read_datahub("Japan",Date(2020,1,24))
   #t_jp = 1:length(C_jp)
   #I_jp = C_jp .- R_jp .- D_jp
   I_jp, lI_jp = smooth_log(I_jp)
   slope_jp,growth_jp = local_fit(t_jp,lI_jp,di)

   t_kr,I_kr,R_kr,D_kr = covid_ts("Korea, South",  Date(2020,2,17))
   #(D_kr, C_kr, R_kr, dates_kr) =  read_datahub("Korea, South",Date(2020,2,17))
   #t_kr = 1:length(C_kr)
   #I_kr = C_kr .- R_kr .- D_kr
   I_kr, lI_kr = smooth_log(I_kr)
   slope_kr,growth_kr = local_fit(t_kr,lI_kr,di)

   t_ir,I_ir,R_ir,D_ir = covid_ts("Iran",          Date(2020,2,22))
   #(D_ir, C_ir, R_ir, dates_ir) =  read_datahub("Iran",Date(2020,2,19))
   #t_ir = 1:length(C_ir)
   #I_ir = C_ir .- R_ir .- D_ir
   I_ir, lI_ir = smooth_log(I_ir)
   slope_ir,growth_ir = local_fit(t_ir,lI_ir,di)


  fig=figure("Outbreak Asia",figsize=(12,8))
  clf()
  subplot(3,4,1)         # China
      plot(t_cn,I_cn,"r*")        
      #plot(dates_nl,I_nl,"r*")        
      ylabel("active cases",fontsize=12)
     
      title("China,"* today)
  ax1=subplot(3,4,5)
      #semilogy(dates_nl,I_nl,"r*", label=country_it, markersize=8)        
      semilogy(t_cn,I_cn,"r*")        
      #semilogy(t_nl,10 .^ lI_fit_nl,"k")
      #ax1.xaxis_date()
      #ax1.xaxis.set_major_formatter(majorformatter)
      #ax1.xaxis.set_major_locator(majorlocator)
      #fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      #xlim(Date(2020,2,21),Date(2020,3,26))
      #ax1.xaxis.grid()
      ax1.yaxis.grid()
      ylim(5,5e5)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_cn,growth_cn,"r*-")
      ylim(-0.1,0.5)
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

  subplot(3,4,2)          # Japan
      plot(t_jp,I_jp,"r*")        
      #ylabel("infected Italy")
      title("Japan,"*today)
  ax=subplot(3,4,6)
      semilogy(t_jp,I_jp,"r*")        
      #semilogy(t_no,10 .^ lI_fit_no,"k")
      ax.yaxis.grid()
      ylim(5,1e5)
  ax=subplot(3,4,10)
      plot(t_jp,growth_jp,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,3)              # Korea, South
      plot(t_kr,I_kr,"r*")      
      #ylabel("infected Spain")
      title("Korea, South,"* today)
  ax=subplot(3,4,7)
      semilogy(t_kr,I_kr,"r*")        
     # semilogy(t_ch,10 .^ lI_fit_ch,"k")
      ax.yaxis.grid()
      ylim(5,1e5)
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_kr,growth_kr,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,4)              # Iran
      plot(t_ir,I_ir,"r*")      
      #ylabel("infected")
      title("Iran,"* today)
  ax=subplot(3,4,8)
      semilogy(t_ir,I_ir,"r*")        
     # semilogy(t_at,10 .^ lI_fit_at,"k")
      ax.yaxis.grid()
      ylim(5,1e5)
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_ir,growth_ir,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")
  PyPlot.tight_layout()
  savefig(folder*"outbreak_Asia.png")

end


# **********************************************************
function nlys_Europe()
   di = 2

   t_de,I_de,R_de,D_de = covid_ts("Germany",       Date(2020,2,23))
   #(D_de, C_de, R_de, dates_de) =  read_datahub("Germany",Date(2020,2,24))
   #t_de = 1:length(C_de)
   #I_de = C_de #.- R_de .- D_de
   I_de, lI_de = smooth_log(I_de)
   slope_de,growth_de = local_fit(t_de,lI_de,di)

   t_it,I_it,R_it,D_it = covid_ts("Italy",         Date(2020,2,21))
   #(D_it, C_it, R_it, dates_it) =  read_datahub("Italy",Date(2020,2,21))
   #t_it = 1:length(C_it)
   #I_it = C_it #.- R_it .- D_it
   I_it, lI_it = smooth_log(I_it)
   slope_it,growth_it = local_fit(t_it,lI_it,di)

   t_fr,I_fr,R_fr,D_fr = covid_ts("France",        Date(2020,2,27))
   #(D_fr, C_fr, R_fr, dates_fr) =  read_datahub("France",Date(2020,2,25))
   #t_fr = 1:length(C_fr)
   #I_fr = C_fr #.- R_fr .- D_fr
   I_fr, lI_fr = smooth_log(I_fr)
   slope_fr,growth_fr = local_fit(t_fr,lI_fr,di)

   t_es,I_es,R_es,D_es = covid_ts("Spain",         Date(2020,2,25))
   #(D_es, C_es, R_es, dates_es) =  read_datahub("Spain",Date(2020,2,25))
   #t_es = 1:length(C_es)
   #I_es = C_es #.- R_es .- D_es
   I_es, lI_es = smooth_log(I_es)
   slope_es,growth_es = local_fit(t_es,lI_es,di)


  fig=figure("Outbreak Europe",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Germany
      plot(t_de,I_de,"r*")        
      #plot(dates_nl,I_nl,"r*")        
      ylabel("active cases",fontsize=12)
     
      title("Germany,"* today)
  ax1=subplot(3,4,5)
      #semilogy(dates_nl,I_nl,"r*", label=country_it, markersize=8)        
      semilogy(t_de,I_de,"r*")        
      #semilogy(t_nl,10 .^ lI_fit_nl,"k")
      #ax1.xaxis_date()
      #ax1.xaxis.set_major_formatter(majorformatter)
      #ax1.xaxis.set_major_locator(majorlocator)
      #fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      #xlim(Date(2020,2,21),Date(2020,3,26))
      #ax1.xaxis.grid()
      ax1.yaxis.grid()
      ylim(5,2e5)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_de,growth_de,"r*-")
      ylim(-0.1,0.5)
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

  subplot(3,4,2)          # Italy
      plot(t_it,I_it,"r*")        
      #ylabel("infected Italy")
      title("Italy,"*today)
  ax=subplot(3,4,6)
      semilogy(t_it,I_it,"r*")        
      #semilogy(t_no,10 .^ lI_fit_no,"k")
      ax.yaxis.grid()
      ylim(5,2e5)
  ax=subplot(3,4,10)
      plot(t_it,growth_it,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,3)              # France
      plot(t_fr,I_fr,"r*")      
      #ylabel("infected Spain")
      title("France,"* today)
  ax=subplot(3,4,7)
      semilogy(t_fr,I_fr,"r*")        
     # semilogy(t_ch,10 .^ lI_fit_ch,"k")
      ax.yaxis.grid()
      ylim(5,2e5)
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_fr,growth_fr,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,4)              # Spain
      plot(t_es,I_es,"r*")      
      #ylabel("infected")
      title("Spain,"* today)
  ax=subplot(3,4,8)
      semilogy(t_es,I_es,"r*")        
     # semilogy(t_at,10 .^ lI_fit_at,"k")
      ax.yaxis.grid()
      ylim(5,2e5)
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_es,growth_es,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")
  PyPlot.tight_layout()
  savefig(folder*"outbreak_Europe.png")

end

# **********************************************************
# **********************************************************
# **********************************************************
function nlys_World()
   di = 2

   t_us,I_us,R_us,D_us = covid_ts("US",            Date(2020,2,22))
   I_us, lI_us = smooth_log(I_us)
   slope_us,growth_us = local_fit(t_us,lI_us,di)

   t_ca,I_ca,R_ca,D_ca = covid_ts("Canada",        Date(2020,2,27))
   I_ca, lI_ca = smooth_log(I_ca)
   slope_ca,growth_ca = local_fit(t_ca,lI_ca,di)

   t_uk,I_uk,R_uk,D_uk = covid_ts("United Kingdom",Date(2020,2,27))
   I_uk, lI_uk = smooth_log(I_uk)
   slope_uk,growth_uk = local_fit(t_uk,lI_uk,di)

   t_il,I_il,R_il,D_il = covid_ts("Israel",        Date(2020,2,27))
   I_il, lI_il = smooth_log(I_il)
   slope_il,growth_il = local_fit(t_il,lI_il,di)

  fig=figure("Outbreak World",figsize=(12,8))
  clf()
  subplot(3,4,1)         # US
      plot(t_us,I_us,"r*")        
      #plot(dates_nl,I_nl,"r*")        
      ylabel("active cases",fontsize=12)
     
      title("US,"* today)
  ax1=subplot(3,4,5)
      semilogy(t_us,I_us,"r*")        
      ax1.yaxis.grid()
      ylim(5,5e6)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_us,growth_us,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=20)

  subplot(3,4,2)          # Canada
      plot(t_ca,I_ca,"r*")        
      title("Canada,"*today)
  ax=subplot(3,4,6)
      semilogy(t_ca,I_ca,"r*")        
      ax.yaxis.grid()
      ylim(5,5e4)
  ax=subplot(3,4,10)
      plot(t_ca,growth_ca,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,3)              # United Kingdom 
      plot(t_uk,I_uk,"r*")      
      title("United Kingdom,"* today)
  ax=subplot(3,4,7)
      semilogy(t_uk,I_uk,"r*")        
      ax.yaxis.grid()
      ylim(5,5e5)
  ax= subplot(3,4,11)
      plot(t_uk,growth_uk,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,4)              # Israel
      plot(t_il,I_il,"r*")      
      title("Israel,"* today)
  ax=subplot(3,4,8)
      semilogy(t_il,I_il,"r*")        
      ax.yaxis.grid()
      ylim(5,5e4)
  ax = subplot(3,4,12)
      plot(t_il,growth_il,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
  PyPlot.tight_layout()
  savefig(folder*"outbreak_World.png")

end

# **********************************************************
# **********************************************************
# **********************************************************
function nlys_SouthAmerica()
   di = 2

   t_br,I_br,R_br,D_br = covid_ts("Brazil",        Date(2020,2,22))
   I_br, lI_br = smooth_log(I_br)
   slope_br,growth_br = local_fit(t_br,lI_br,di)

   t_cl,I_cl,R_cl,D_cl = covid_ts("Chile",        Date(2020,2,27))
   I_cl, lI_cl = smooth_log(I_cl)
   slope_cl,growth_cl = local_fit(t_cl,lI_cl,di)

   t_ar,I_ar,R_ar,D_ar = covid_ts("Argentina",Date(2020,2,27))
   I_ar, lI_ar = smooth_log(I_ar)
   slope_ar,growth_ar = local_fit(t_ar,lI_ar,di)

   t_ec,I_ec,R_ec,D_ec = covid_ts("Ecuador",        Date(2020,2,27))
   I_ec, lI_ec = smooth_log(I_ec)
   slope_ec,growth_ec = local_fit(t_ec,lI_ec,di)

  fig=figure("Outbreak South America",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Brazil
      plot(t_br,I_br,"r*")        
      ylabel("active cases",fontsize=12)
      title("Brazil,"* today)
  ax1=subplot(3,4,5)
      semilogy(t_br,I_br,"r*")        
      ax1.yaxis.grid()
      ylim(5,5e5)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_br,growth_br,"r*-")
      ylim(0,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=20)

  subplot(3,4,2)          # Chile
      plot(t_cl,I_cl,"r*")        
      title("Chile,"*today)
  ax=subplot(3,4,6)
      semilogy(t_cl,I_cl,"r*")        
      ax.yaxis.grid()
      ylim(5,5e4)
  ax=subplot(3,4,10)
      plot(t_cl,growth_cl,"r*-")
      ylim(0,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,3)              # Argentina 
      plot(t_ar,I_ar,"r*")      
      title("Argentina,"* today)
  ax=subplot(3,4,7)
      semilogy(t_ar,I_ar,"r*")        
      ax.yaxis.grid()
      ylim(5,1e5)
  ax= subplot(3,4,11)
      plot(t_ar,growth_ar,"r*-")
      ylim(0,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,4)              # Ecuador
      plot(t_ec,I_ec,"r*")      
      title("Ecuador,"* today)
  ax=subplot(3,4,8)
      semilogy(t_ec,I_ec,"r*")        
      ax.yaxis.grid()
      ylim(5,5e4)
  ax = subplot(3,4,12)
      plot(t_ec,growth_ec,"r*-")
      ylim(0,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
  PyPlot.tight_layout()
  savefig(folder*"outbreak_SouthAmerica.png")

end



# **********************************************************
# **********************************************************
# **********************************************************
function nlys_Africa()
   di = 2

   t_eg,I_eg,R_eg,D_eg = covid_ts("Egypt",        Date(2020,3,8))
   I_eg, lI_eg = smooth_log(I_eg)
   slope_eg,growth_eg = local_fit(t_eg,lI_eg,di)

   t_sa,I_sa,R_sa,D_sa = covid_ts("South Africa",        Date(2020,3,8))
   I_sa, lI_sa = smooth_log(I_sa)
   slope_sa,growth_sa = local_fit(t_sa,lI_sa,di)

   t_al,I_al,R_al,D_al = covid_ts("Algeria",Date(2020,3,8))
   I_al, lI_al = smooth_log(I_al)
   slope_al,growth_al = local_fit(t_al,lI_al,di)

   t_ca,I_ca,R_ca,D_ca = covid_ts("Cameroon",   Date(2020,3,12))
   I_ca, lI_ca = smooth_log(I_ca)
   slope_ca,growth_ca = local_fit(t_ca,lI_ca,di)

  fig=figure("Africa",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Egypt
      plot(t_eg,I_eg,"r*")        
      ylabel("active cases",fontsize=12)
      title("Egypt,"* today)
  ax1=subplot(3,4,5)
      semilogy(t_eg,I_eg,"r*")        
      ax1.yaxis.grid()
      ylim(5,5e5)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_eg,growth_eg,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=20)

  subplot(3,4,2)          # South Africa
      plot(t_sa,I_sa,"r*")        
      title("South Africa,"*today)
  ax=subplot(3,4,6)
      semilogy(t_sa,I_sa,"r*")        
      ax.yaxis.grid()
      ylim(5,5e4)
  ax=subplot(3,4,10)
      plot(t_sa,growth_sa,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,3)              # Algeria
      plot(t_al,I_al,"r*")      
      title("Algeria,"* today)
  ax=subplot(3,4,7)
      semilogy(t_al,I_al,"r*")        
      ax.yaxis.grid()
      ylim(5,1e5)
  ax= subplot(3,4,11)
      plot(t_al,growth_al,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,4)              # Cameroon
      plot(t_ca,I_ca,"r*")      
      title("Cameroon,"* today)
  ax=subplot(3,4,8)
      semilogy(t_ca,I_ca,"r*")        
      ax.yaxis.grid()
      ylim(5,1e4)
  ax = subplot(3,4,12)
      plot(t_ca,growth_ca,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
  PyPlot.tight_layout()
  savefig(folder*"outbreak_Africa.png")

end


# **********************************************************
# **********************************************************
# **********************************************************
function nlys_Russia()
   di = 2

   t_ru,I_ru,R_ru,D_ru = covid_ts("Russia",        Date(2020,2,29))
   I_ru, lI_ru = smooth_log(I_ru)
   slope_ru,growth_ru = local_fit(t_ru,lI_ru,di)

   t_cz,I_cz,R_cz,D_cz = covid_ts("Czechia",        Date(2020,3,12))
   I_cz, lI_cz = smooth_log(I_cz)
   slope_cz,growth_cz = local_fit(t_cz,lI_cz,di)

   t_hu,I_hu,R_hu,D_hu = covid_ts("Hungary",        Date(2020,3,5))
   I_hu, lI_hu = smooth_log(I_hu)
   slope_hu,growth_hu = local_fit(t_hu,lI_hu,di)

   t_pl,I_pl,R_pl,D_pl = covid_ts("Poland",Date(2020,2,27))
   I_pl, lI_pl = smooth_log(I_pl)
   slope_pl,growth_pl = local_fit(t_pl,lI_pl,di)


  fig=figure("Outbreak Russia",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Russia
      plot(t_ru,I_ru,"r*")        
      ylabel("active cases",fontsize=12)
      title("Russia,"* today)
  ax1=subplot(3,4,5)
      semilogy(t_ru,I_ru,"r*")        
      ax1.yaxis.grid()
      ylim(5,5e5)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_ru,growth_ru,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=20)

  subplot(3,4,2)          # Czechia 
      plot(t_cz,I_cz,"r*")        
      title("Czechia,"*today)
  ax=subplot(3,4,6)
      semilogy(t_cz,I_cz,"r*")        
      ax.yaxis.grid()
      ylim(5,5e4)
  ax=subplot(3,4,10)
      plot(t_cz,growth_cz,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,3)              # Hungary 
      plot(t_hu,I_hu,"r*")      
      title("Hungary,"* today)
  ax=subplot(3,4,7)
      semilogy(t_hu,I_hu,"r*")        
      ax.yaxis.grid()
      ylim(5,1e5)
  ax= subplot(3,4,11)
      plot(t_hu,growth_hu,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,4)              # Poland
      plot(t_pl,I_pl,"r*")      
      title("Poland,"* today)
  ax=subplot(3,4,8)
      semilogy(t_pl,I_pl,"r*")        
      ax.yaxis.grid()
      ylim(5,1e5)
  ax = subplot(3,4,12)
      plot(t_pl,growth_pl,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
  PyPlot.tight_layout()
  savefig(folder*"outbreak_Russia.png")

end





# **********************************************************
# **********************************************************
# **********************************************************
function nlys_Australia()
   di = 2

   t_au,I_au,R_au,D_au = covid_ts("Australia",        Date(2020,2,29))
   I_au, lI_au = smooth_log(I_au)
   slope_au,growth_au = local_fit(t_au,lI_au,di)

   t_nz,I_nz,R_nz,D_nz = covid_ts("New Zealand",        Date(2020,3,12))
   I_nz, lI_nz = smooth_log(I_nz)
   slope_nz,growth_nz = local_fit(t_nz,lI_nz,di)

   t_ph,I_ph,R_ph,D_ph = covid_ts("Philippines",        Date(2020,3,5))
   I_ph, lI_ph = smooth_log(I_ph)
   slope_ph,growth_ph = local_fit(t_ph,lI_ph,di)

   t_ml,I_ml,R_ml,D_ml = covid_ts("Malaysia",Date(2020,2,27))
   I_ml, lI_ml = smooth_log(I_ml)
   slope_ml,growth_ml = local_fit(t_ml,lI_ml,di)


  fig=figure("Outbreak Australia",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Australia
      plot(t_au,I_au,"r*")        
      ylabel("active cases",fontsize=12)
      title("Australia,"* today)
  ax1=subplot(3,4,5)
      semilogy(t_au,I_au,"r*")        
      ax1.yaxis.grid()
      ylim(5,5e5)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_au,growth_au,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      ylabel("growth",fontsize=12,labelpad=20)

  subplot(3,4,2)          # New Zealand
      plot(t_nz,I_nz,"r*")        
      title("New Zealand,"*today)
  ax=subplot(3,4,6)
      semilogy(t_nz,I_nz,"r*")        
      ax.yaxis.grid()
      ylim(5,5e4)
  ax=subplot(3,4,10)
      plot(t_nz,growth_nz,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,3)              # Philippines 
      plot(t_ph,I_ph,"r*")      
      title("Philippines,"* today)
  ax=subplot(3,4,7)
      semilogy(t_ph,I_ph,"r*")        
      ax.yaxis.grid()
      ylim(5,1e5)
  ax= subplot(3,4,11)
      plot(t_ph,growth_ph,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)

  subplot(3,4,4)              # Malaysia
      plot(t_ml,I_ml,"r*")      
      title("Malaysia,"* today)
  ax=subplot(3,4,8)
      semilogy(t_ml,I_ml,"r*")        
      ax.yaxis.grid()
      ylim(5,1e4)
  ax = subplot(3,4,12)
      plot(t_ml,growth_ml,"r*-")
      ylim(-0.1,0.5)
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
  PyPlot.tight_layout()
  savefig(folder*"outbreak_Australia.png")

end

# **********************************************************
# **********************************************************
# **********************************************************
# **********************************************************
function nlys_smallEurope()
   di = 2

   t_nl,I_nl,R_nl,D_nl = covid_ts("Netherlands",   Date(2020,2,23))
   #(D_nl, C_nl, R_nl, dates_nl) =  read_datahub("Netherlands",Date(2020,2,23))
   #t_nl = 1:length(C_nl)
   #I_nl = C_nl .- R_nl .- D_nl
   I_nl, lI_nl = smooth_log(I_nl)
   slope_nl,growth_nl = local_fit(t_nl,lI_nl,di)

   t_no,I_no,R_no,D_no = covid_ts("Norway",        Date(2020,2,28))
   #(D_no, C_no, R_no, dates_no) =  read_datahub("Norway",Date(2020,2,28))
   #t_no = 1:length(C_no)
   #I_no = C_no .- R_no .- D_no
   I_no, lI_no = smooth_log(I_no)
   slope_no,growth_no = local_fit(t_no,lI_no,di)

   t_ch,I_ch,R_ch,D_ch = covid_ts("Switzerland",   Date(2020,2,27))
   #(D_ch, C_ch, R_ch, dates_ch) =  read_datahub("Switzerland",Date(2020,2,25))
   #t_ch = 1:length(C_ch)
   #I_ch = C_ch .- R_ch .- D_ch
   I_ch, lI_ch = smooth_log(I_ch)
   slope_ch,growth_ch = local_fit(t_ch,lI_ch,di)

   t_at,I_at,R_at,D_at = covid_ts("Austria",       Date(2020,2,29))
   #(D_at, C_at, R_at, dates_at) =  read_datahub("Austria",Date(2020,2,29))
   #t_at = 1:length(C_at)
   #I_at = C_at .- R_at .- D_at
   I_at, lI_at = smooth_log(I_at)
   slope_at,growth_at = local_fit(t_at,lI_at,di)

   t_dk,I_dk,R_dk,D_dk = covid_ts("Denmark",        Date(2020,2,27))
   #(D_dk, C_dk, R_dk, dates_dk) =  read_datahub("Denmark",Date(2020,2,27))
   #t_dk = 1:length(C_dk)
   #I_dk = C_dk .- R_dk .- D_dk
   I_dk, lI_dk = smooth_log(I_dk)
   slope_dk,growth_dk = local_fit(t_dk,lI_dk,di)

   t_se,I_se,R_se,D_se = covid_ts("Sweden",        Date(2020,2,27))
   #(D_dk, C_dk, R_dk, dates_dk) =  read_datahub("Denmark",Date(2020,2,27))
   #t_dk = 1:length(C_dk)
   #I_dk = C_dk .- R_dk .- D_dk
   I_se, lI_se = smooth_log(I_se)
   slope_se,growth_se = local_fit(t_se,lI_se,di)

   t_tr,I_tr,R_tr,D_tr = covid_ts("Turkey",        Date(2020,2,27))
   #(D_dk, C_dk, R_dk, dates_dk) =  read_datahub("Denmark",Date(2020,2,27))
   #t_dk = 1:length(C_dk)
   #I_dk = C_dk .- R_dk .- D_dk
   I_tr, lI_tr = smooth_log(I_tr)
   slope_tr,growth_tr = local_fit(t_tr,lI_tr,di)

   t_fi,I_fi,R_fi,D_fi = covid_ts("Finland",        Date(2020,2,27))
   #(D_dk, C_dk, R_dk, dates_dk) =  read_datahub("Denmark",Date(2020,2,27))
   #t_dk = 1:length(C_dk)
   #I_dk = C_dk .- R_dk .- D_dk
   I_fi, lI_fi = smooth_log(I_fi)
   slope_fi,growth_fi = local_fit(t_fi,lI_fi,di)


# ********************************************
  majorformatter = matplotlib.dates.DateFormatter("%d.%m.%Y")
  #majorformatter = matplotlib.dates.DateFormatter("%d.%m")
  majorlocator = matplotlib.dates.DayLocator(interval=7)
  fig=figure("Outbreak Small Europe",figsize=(12,8))
  clf()
  subplot(3,4,1)         # Netherlands
      plot(t_nl,I_nl,"r*")        
      #plot(dates_nl,I_nl,"r*")        
      ylabel("active cases",fontsize=12)
     
      title("Netherlands,"* today)
  ax1=subplot(3,4,5)
      #semilogy(dates_nl,I_nl,"r*", label=country_it, markersize=8)        
      semilogy(t_nl,I_nl,"r*")        
      #semilogy(t_nl,10 .^ lI_fit_nl,"k")
      #ax1.xaxis_date()
      #ax1.xaxis.set_major_formatter(majorformatter)
      #ax1.xaxis.set_major_locator(majorlocator)
      #fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      #xlim(Date(2020,2,21),Date(2020,3,26))
      #ax1.xaxis.grid()
      ax1.yaxis.grid()
      ylim(5,5e4)
      ylabel("active cases",labelpad=16,fontsize=12)
   ax = subplot(3,4,9)
      plot(t_nl,growth_nl,"r*-")
      ylim(-0.1,0.5)
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

  subplot(3,4,2)          # Turkey
      plot(t_tr,I_tr,"r*")        
      #ylabel("infected Italy")
      title("Turkey,"*today)
  ax=subplot(3,4,6)
      semilogy(t_tr,I_tr,"r*")        
      #semilogy(t_dk,10 .^ lI_fit_dk,"k")
      ax.yaxis.grid()
      ylim(5,1e5)
  ax=subplot(3,4,10)
      plot(t_tr,growth_tr,"r*-")
      ylim(-0.1,0.5)
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
     # semilogy(t_ch,10 .^ lI_fit_ch,"k")
      ax.yaxis.grid()
      ylim(5,5e4)
      #ylabel("log10 infected Spain")
  ax= subplot(3,4,11)
      plot(t_ch,growth_ch,"r*-")
      ylim(-0.1,0.5)
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
     # semilogy(t_at,10 .^ lI_fit_at,"k")
      ax.yaxis.grid()
      ylim(5,5e4)
      #ylabel("log10 infected")
  ax = subplot(3,4,12)
      plot(t_at,growth_at,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")
  PyPlot.tight_layout()
  savefig(folder*"outbreak_europe1.png")

# ********************************************
  fig=figure("Outbreak Small Europe2",figsize=(12,8))
  clf()

  subplot(3,4,1)          # Sweden
      plot(t_se,I_se,"r*")        
      #ylabel("infected Italy")
      title("Sweden,"*today)
  ax=subplot(3,4,5)
      semilogy(t_se,I_se,"r*")        
      #semilogy(t_dk,10 .^ lI_fit_dk,"k")
      ax.yaxis.grid()
      ylim(5,1e5)
  ax=subplot(3,4,9)
      plot(t_se,growth_se,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,2)          # Denmark
      plot(t_dk,I_dk,"r*")        
      #ylabel("infected Italy")
      title("Denmark,"*today)
  ax=subplot(3,4,6)
      semilogy(t_dk,I_dk,"r*")        
      #semilogy(t_dk,10 .^ lI_fit_dk,"k")
      ax.yaxis.grid()
      ylim(5,1e4)
  ax=subplot(3,4,10)
      plot(t_dk,growth_dk,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")


  subplot(3,4,3)          # Norway
      plot(t_no,I_no,"r*")        
      #ylabel("infected Italy")
      title("Norway,"*today)
  ax=subplot(3,4,7)
      semilogy(t_no,I_no,"r*")        
      #semilogy(t_no,10 .^ lI_fit_no,"k")
      ax.yaxis.grid()
      ylim(5,1e4)
  ax=subplot(3,4,11)
      plot(t_no,growth_no,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  subplot(3,4,4)          # Finland
      plot(t_fi,I_fi,"r*")        
      #ylabel("infected Italy")
      title("Finland,"*today)
  ax=subplot(3,4,8)
      semilogy(t_fi,I_fi,"r*")        
      #semilogy(t_no,10 .^ lI_fit_no,"k")
      ax.yaxis.grid()
      ylim(5,1e4)
  ax=subplot(3,4,12)
      plot(t_fi,growth_fi,"r*-")
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2,0.25])
      ax.yaxis.grid()
      xlabel("time (days)",fontsize=12)
      #ylabel("growth")

  PyPlot.tight_layout()
  savefig(folder*"outbreak_europe2.png")
end



# **********************************************************
# **********************************************************
# **********************************************************
function nlys_ItDe()

  t_it,I_it,R_it,D_it = covid_ts("Italy",Date(2020,2,21))
  #(D_it, C_it, R_it, dates_it) =  read_datahub("Italy",Date(2020,2,21))
  #(D_it, I_it, dates_it) =  read_datahub("Italy",Date(2020,2,21))
  #t_it = 1:length(I_it)
  #I_it = C_it .- R_it .- D_it
  dates_it = Date(2020,2,21) .+ Day.(t_it .- 1)
  I_it, lI_it = smooth_log(I_it)
  #lI_it = log10.(I_it)
  di = 2
  slope_it,growth_it = local_fit(t_it,lI_it,di)

  
  p_tot_it1 = fit(t_it[3:8],lI_it[3:8],1)  # fit slope 1st phase Italy
  coef = coeffs(p_tot_it1)[2]
  println("slope Italy Phase 1 ",coef, " T12 ", doubling_time(coef))
  #lI_fit_it1 = polyval(p_tot_it1, t_it)
  lI_fit_it1 = p_tot_it1.(t_it)

  p_tot_it2 = fit(t_it[9:22],lI_it[9:22],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it2)[2]
  println("slope Italy Phase 2 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it2 = p_tot_it2.(t_it)

  p_tot_it3 = fit(t_it[23:35],lI_it[23:35],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it3)[2]
  println("slope Italy Phase 3 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it3 = p_tot_it3.(t_it)

  # ********************
  t_de,I_de,R_de,D_de = covid_ts("Germany",  Date(2020,2,24))
  #(D_de, C_de, R_de, dates_de) =  read_datahub("Germany",Date(2020,2,24))
  #(D_de, I_de, dates_de) =  read_datahub("Germany",Date(2020,2,24))
  #t_de = 1:length(I_de)
  #I_de = C_de .- R_de .- D_de
  dates_de = Date(2020,2,24) .+ Day.(t_de .- 1 )
  I_de, lI_de = smooth_log(I_de)
  #lI_de = log10.(I_de)
  di = 2
  slope_de,growth_de = local_fit(t_de,lI_de,di)
  t_de = t_de .+ 3 # bring to the same time axis
  
  p_tot_de1 = fit(t_de[4:11],lI_de[4:11],1)  # fit full slope Germany
  coef = coeffs(p_tot_de1)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de1= p_tot_de1.(t_de)

  p_tot_de2 = fit(t_de[12:25],lI_de[12:25],1)  # fit full slope Germany
  coef = coeffs(p_tot_de2)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de2= p_tot_de2.(t_de)

  p_tot_de3 = fit(t_de[26:34],lI_de[26:34],1)  # fit full slope Germany
  coef = coeffs(p_tot_de3)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de3= p_tot_de3.(t_de)

  p_tot_de4 = fit(t_de[35:end],lI_de[35:end],1)  # fit full slope Germany
  coef = coeffs(p_tot_de4)[2]
  println("slope Germany ",coef, " T12 ", doubling_time(coef))
  lI_fit_de4= p_tot_de4.(t_de)

   majorformatter = matplotlib.dates.DateFormatter("%d.%m.%Y")
   #majorformatter = matplotlib.dates.DateFormatter("%d.%m")
   majorlocator = matplotlib.dates.DayLocator(interval=7)
   lang = "en"
   if lang == "en"
       country_de = "Germany"
       country_it = "Italy"
   else
       country_de = "Deutschland"
       country_it = "Italien"
   end
   fig=figure("Germany_Italy", figsize=(12,6))
   clf()
   subplot(1,2,1)
      #semilogy(t_it[1:12],10 .^ lI_fit_it1[1:12],"k")
      #semilogy(t_it[5:25],10 .^ lI_fit_it2[5:25],"k")
      #semilogy(t_it[19:28],10 .^ lI_fit_it3[19:28],"k")
      #semilogy(t_it,I_it,"b*", label="Italy", markersize=8)        
      #semilogy(t_de[1:15],10 .^ lI_fit_de1[1:15],"k")
      #semilogy(t_de[8:end],10 .^ lI_fit_de2[8:end],"k")
      #semilogy(t_de,I_de,"r*",label="Germany",markersize=8)        
      semilogy(dates_it[1:12],10 .^ lI_fit_it1[1:12],"k")
      semilogy(dates_it[5:25],10 .^ lI_fit_it2[5:25],"k")
      #semilogy(dates_it[19:37],10 .^ lI_fit_it3[19:37],"k")
      semilogy(dates_it,I_it,"b*", label=country_it, markersize=8)        
      #semilogy(dates_it,I_it,"b*", label="Italy", markersize=8)        
      semilogy(dates_de[1:15],10 .^ lI_fit_de1[1:15],"k")
      semilogy(dates_de[8:25],10 .^ lI_fit_de2[8:25],"k")
      semilogy(dates_de[23:35],10 .^ lI_fit_de3[23:35],"k")
      #semilogy(dates_de[32:end],10 .^ lI_fit_de4[32:end],"k")
      #semilogy(dates_de,I_de,"r*",label="Germany",markersize=8)        
      semilogy(dates_de,I_de,"r*",label=country_de,markersize=8)        
      ax1 = gca()
      ax1.xaxis_date()
      ax1.xaxis.set_major_formatter(majorformatter)
      ax1.xaxis.set_major_locator(majorlocator)
      fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      xlim(Date(2020,2,21),today_date)
      ylim(8,2e5)
      ax1.xaxis.grid()
      ax1.yaxis.grid()
      if lang == "en"
       # ylabel("active cases",fontsize=14)
        xlabel("time",fontsize=14)
        #xlabel("time (Feb 21 - March 19)", fontsize=14)
      else
        #xlabel("Zeit", fontsize=14)
        ylabel("aktive Fälle",fontsize=14)
      end
      #xlabel("time (Feb 21 - March 19)", fontsize=14)
      legend(loc="upper left",fancybox="true")
   subplot(1,2,2)
      #plot(t_it,growth_it,"b*-",label="Italy",markersize=8)
      #plot(t_de,growth_de,"r*-",label="Germany",markersize=8)
        plot(dates_it,growth_it,"b*-",label=country_it,markersize=8)
        plot(dates_de,growth_de,"r*-",label=country_de,markersize=8)
      #plot([10,10],[0,0.6],"grey",linewidth=1)
      #ylim(0,0.22)
      #yticks([0,0.05,0.1,0.15,0.2])
      ax1 = gca()
      ax1.xaxis_date()
      ax1.xaxis.set_major_formatter(majorformatter)
      ax1.xaxis.set_major_locator(majorlocator)
      fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      ylim(-0.1,0.5)
      xlim(Date(2020,2,21),today_date)
      #yticks([0,0.05,0.1,0.15,0.2])
      ax1.xaxis.grid()
      ax1.yaxis.grid()
      if lang == "en"
        #xlabel("time",fontsize=14)
        ylabel("growth rate",fontsize=14)
        #xlabel("time (Feb 21 - March 19)", fontsize=14)
      else
        #xlabel("Zeit", fontsize=14)
        ylabel("Wachstumsrate",fontsize=14)
      end
      legend(loc="upper right",fancybox="true")
  PyPlot.tight_layout()
  savefig(folder*"ItalyGermany_March.png")
  #savefig("../figures/ItalienDeutschland_March.png")
end


# **********************************************************
function nlys_ItUS()
  #(D_it, C_it, R_it, dates_it) =  read_datahub("Italy",Date(2020,2,22))
  #(D_it, I_it, dates_it) =  read_datahub("Italy",Date(2020,2,22))
  #t_it = 1:length(I_it)
  t_it,I_it,R_it,D_it = covid_ts("Italy",Date(2020,2,21))
  #I_it = C_it .- R_it .- D_it
  dates_it = Date(2020,2,21) .+ Day.(t_it .- 1)
  I_it, lI_it = smooth_log(I_it)
  #lI_it = log10.(I_it)
  di = 2
  slope_it,growth_it = local_fit(t_it,lI_it,di)
  
  p_tot_it1 = fit(t_it[2:7],lI_it[2:7],1)  # fit slope 1st phase Italy
  coef = coeffs(p_tot_it1)[2]
  println("slope Italy Phase 1 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it1 = p_tot_it1.(t_it)

  p_tot_it2 = fit(t_it[8:21],lI_it[8:21],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it2)[2]
  println("slope Italy Phase 2 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it2 = p_tot_it2.(t_it)

  p_tot_it3 = fit(t_it[22:end],lI_it[22:end],1)  # fit slope 2nd phase Italy
  coef = coeffs(p_tot_it3)[2]
  println("slope Italy Phase 3 ",coef, " T12 ", doubling_time(coef))
  lI_fit_it3 = p_tot_it3.(t_it)


  # ********************  US
  t_us,I_us,R_us,D_us = covid_ts("US",            Date(2020,2,22))
  #(D_us, C_us, R_us, dates_us) =  read_datahub("US",Date(2020,2,22))
  #(D_us, I_us, dates_us) =  read_datahub("US",Date(2020,2,22))
  #t_us = 1:length(I_us)
  #I_us = C_us # for the US we just take the confirmed cases
  dates_us = Date(2020,2,22) .+ Day.(t_us .- 1 )
  I_us, lI_us = smooth_log(I_us)
  #lI_us = log10.(I_us)
  p_tot_us = fit(t_us[11:26],lI_us[11:26],1)  # fit full slope USA
  lI_fit_us= p_tot_us.(t_us)
  coef = coeffs(p_tot_us)[2]
  slope_us, growth_us = local_fit(t_us,lI_us,di)
  println("slope US ",coef, " T12 ", doubling_time(coef))


   majorformatter = matplotlib.dates.DateFormatter("%d.%m.%Y")
   #majorformatter = matplotlib.dates.DateFormatter("%d.%m")
   majorlocator = matplotlib.dates.DayLocator(interval=7)

   fig = figure("ItalyUS", figsize=(12,6))
   clf()
   subplot(1,2,1)
      #semilogy(t_it[1:12],10 .^ lI_fit_it1[1:12],"k")
      #semilogy(t_it[5:25],10 .^ lI_fit_it2[5:25],"k")
      #semilogy(t_it[19:28],10 .^ lI_fit_it3[19:28],"k")
      #semilogy(t_it,I_it,"b*", label="Italy", markersize=8)        
      #semilogy(t_us,I_us,"r*",label="US",markersize=8)        
      #semilogy(t_us[8:27],10 .^ lI_fit_us[8:27],"k")
      semilogy(dates_it[1:11],10 .^ lI_fit_it1[1:11],"k")
      semilogy(dates_it[4:24],10 .^ lI_fit_it2[4:24],"k")
      semilogy(dates_it[18:31],10 .^ lI_fit_it3[18:31],"k")
      semilogy(dates_it,I_it,"b*", label="Italy", markersize=8)        
      semilogy(dates_us,I_us,"r*",label="US",markersize=8)        
      semilogy(dates_us[8:30],10 .^ lI_fit_us[8:30],"k")
      ax1 = gca()
      ax1.xaxis_date()
      ax1.xaxis.set_major_formatter(majorformatter)
      ax1.xaxis.set_major_locator(majorlocator)
      ax1.xaxis.grid()
      ax1.yaxis.grid()
      fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      xlim(Date(2020,2,21),today_date)
      ylim(8,1e6)
      ylabel("active cases",fontsize=14)
      #xlabel("time (Feb 21 - March 19)", fontsize=14)
      legend(loc="upper left",fancybox="true")
   subplot(1,2,2)
      #plot(t_it,growth_it,"b*-",label="Italy",markersize=8)
      #plot(t_us,growth_us,"r*-",label="US",markersize=8)
      plot(dates_it,growth_it,"b*-",label="Italy",markersize=8)
      plot(dates_us,growth_us,"r*-",label="US",markersize=8)
      #plot([10,10],[0,0.6],"grey",linewidth=1)
      #ylim(0,0.22)
      #yticks([0,0.05,0.1,0.15,0.2])
      ax1 = gca()
      ax1.xaxis_date()
      ax1.xaxis.set_major_formatter(majorformatter)
      ax1.xaxis.set_major_locator(majorlocator)
      ax1.xaxis.grid()
      ax1.yaxis.grid()
      fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
      xlim(Date(2020,2,21),today_date)
      ylim(-0.1,0.5)
      #yticks([0,0.05,0.1,0.15,0.2])
      #xlabel("time (days)",fontsize=14)
      #xlabel("time (Feb 21 - March 19)", fontsize=14)
      ylabel("growth rate",fontsize=14)
      legend(loc="upper right",fancybox="true")
  PyPlot.tight_layout()
  savefig(folder*"/ItalyUS_March.png")


end


# ********************************************
function corona_scaling() 

   t_it,I_it,R_it,D_it = covid_ts("Italy",Date(2020,2,20))
   I_it, lI_it = smooth_log(I_it)
   lt_it = log10.(t_it)
   p_tot_it = fit(lt_it,lI_it,1)  # fit full slope Germany
   println("loglog It ", coeffs(p_tot_it)[2])
   I_fit_it = 10 .^ p_tot_it.(lt_it)

   #lt_de = log10.(t_de)
   #p_tot_de = fit(lt_de,lI_de,1)  # fit full slope Germany
   #println("loglog It ", coeffs(p_tot_de)[2])
   #I_fit_de = 10 .^ polyval(p_tot_de, lt_de)
   
   figure("scaling")
   clf()
   #subplot(2,1,1)
   title("Power-law scaling of COVID-19 outbreak in Italy")
   loglog(t_it,I_it,"*-", markersize=8)
   loglog(t_it,I_fit_it,"-k", label="regression line, slope=2.75")
   ax = gca()
   ax.yaxis.grid()
   ax.xaxis.grid()
   xlabel("time (days)",fontsize=12)
   ylabel("active cases (Italy)",fontsize=12)
   legend(loc="upper left",fancybox="true")
   savefig(folder*"Italy_scaling.png")

end



#nlys_country(country, date)
today = " August 19"
today_date = Date(2020,8,20)
folder = "../figures/"

nlys_Asia()
nlys_Europe()
nlys_World()
nlys_SouthAmerica()
nlys_Africa()
nlys_Russia()
nlys_Australia()
nlys_smallEurope()
nlys_ItDe()
#nlys_ItUS()

#corona_scaling()

# Russia, Poland
# Indonesia, Phillippines, Malaysia
# India, Pakistan

end   # module

