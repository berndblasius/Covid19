using Optim

#cumulpdf(x,alpha,xmax) =  x^(1-alpha) - xmax^(1-alpha)

cumulpdf(x,a,xmax) = (x^(1-a)-xmax^(1-a)) / (1-xmax^(1-a))


function loghistogram(ts, nb=20)
# histogram of log-transformed values   
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

function log2histogram1(ts)
# histogram with logarithmic binning
# bin positions at powers of two, 2,4,8,16,...
    maxt = maximum(ts)    
    maxind = ceil(Int,log2(maxt))
    bins = 2.0 .^ (1:maxind)   # bins
    nbins = length(bins)
    w = zeros(nbins)        # weights
    for i=1:length(ts)
       t = ts[i]
       if t == 1
           w[1] += 1
       else
           ind = ceil(Int,log2(t)) 
           w[ind] +=1
       end
    end
    ind = findall(x->x != 0, w)  # remove bins with entry 0
    w = w[ind]
    bins = bins[ind]
    distr = w ./ bins  # get distribution by normalization
    distr = distr / sum(distr)

    bins, distr
end

function log2histogram(ts)
# histogram with logarithmic binning
# bin positions at powers of two, 1,2,4,8,16,...
    maxt = maximum(ts)    
    maxind = ceil(Int,log2(maxt))
    bins = 2.0 .^ (0:maxind)   # bins
    nbins = length(bins)
    weight = zeros(nbins)        # weights
    widths = zeros(nbins)        # widths of bins
    for t in ts
       ind = ceil(Int,log2(t)) + 1
       weight[ind] +=1
    end
    widths[1] = 0.5
    widths[2:end] = bins[1:end-1] 
    ind = findall(x->x != 0, weight)  # remove bins with entry 0
    weight = weight[ind]
    widths = widths[ind]
    bins   = bins[ind]
    distr = weight ./ widths  # get distribution by normalization
    distr = distr / sum(distr)

    bins, distr
end


# log likelihood
like(mu,g,nmax) = log((1-mu)/(nmax^(1-mu)-1)) - mu*g

like_inf(mu,g) = log(mu-1) - mu*g

function standard_err(M,mu,nmax) 
# standard error of exponent
   mu1 = mu-1.0    
   denom = 1-nmax^mu1
   1.0 / sqrt(1.0/(mu1*mu1) - nmax^mu1 * log(log(nmax))/(denom*denom)) /
   sqrt(M)
end

function KS_test(C,n,mu,nmax)
    # Kolmogorov-Smirnov test
   diff = abs.(cumulpdf.(C,mu,nmax) .- (1:n)/n)
   maximum(diff)
end


function estimate(C)
   nmax = maximum(C)
   n = length(C)
   estimate(C,n,nmax)
end

function estimate(C,n,nmax)
   # estimate fitting values for truncated power law
   # critical exponent, standard deviation Kolmogorov-Smirnov distance
   g = sum(log.(C)) / n   # log(geometric mean)
   #g = sum(log.(2*C)) / n   # try correction from Clauset
   # find maximum using Brent's method
   mu = optimize(x->-like(x,g,nmax), 1.01, 3.0).minimizer
  # mu = optimize(x->-like(x,g,nmax), 1.01, 3.0, GoldenSection()).minimizer
   st = standard_err(n,mu,nmax) 
   #d = KS_test(C,n,mu,nmax)
   #mu, st, d
   mu, st
end

function randPL(n,mu,nmax)
# generate sample of random numbers taken from a truncated power law  
   mu1 = 1.0 - mu
   r = (1.0 .- rand(n) .* (1.0 - nmax^mu1)) .^ (1.0/mu1)
   sort(r, rev=true)
end 
