from math import *
import numpy as N

def fit(x, y, sig=None):
    #Given a set of data points x, y,
    #with individual standard deviations sig,
    #fit them to a straight line y = a + bx by minimizing chi-sq.
    #Returned are a, b and their respective probable uncertainties siga and sigb,
    #the chi-square chi2, and the scatter sigdat.
    #If mwt=0 on input, then the standard deviations are assumed to be unavailable:
    #the normalization of chi2 is to unit standard deviation on all points.
    if sig is None or sig == 1.0:
	mwt=False
    else:
	mwt=True
    sx=0.0
    sy=0.0
    st2=0.0
    b=0.0
    ndata = len(x)
    if mwt:
        #Accumulate sums ...
        ss=0.0
        for i in range(ndata):
            #...with weights
            wt=1.0/sig[i]**2
            ss += wt
            sx += x[i]*wt
            sy += y[i]*wt
    else:
        for i in range(ndata):
            #...or without weights.
            sx += x[i]
            sy += y[i]
        ss=ndata
    sxoss=sx/ss
    if mwt:
        for i in range(ndata):
            t=(x[i]-sxoss)/sig[i]
            st2 += t*t
            b += t*y[i]/sig[i]
    else:
        for i in range(ndata):
            t=x[i]-sxoss
            st2 += t*t
            b += t*y[i]
    b /= st2
    #Solve for a, b, siga, and sigb.
    a=(sy-sx*b)/ss
    siga=sqrt((1.0+sx*sx/(ss*st2))/ss)
    sigb=sqrt(1.0/st2)
    #Calculate chi2.
    chi2=0.0
    if mwt:
        for i in range(ndata):
            chi2 += ((y[i] - a - b*x[i]) / sig[i])**2
    else:
        #For unweighted data evaluate typical sig using chi2,
        #and adjust the standard deviations.
        for i in range(ndata):
            chi2 += (y[i] - a - b*x[i])**2
        sigdat=sqrt(chi2/(ndata-2))
        siga *= sigdat
        sigb *= sigdat
    return a, b, siga, sigb, chi2


def fit_slope(x, y, sig, mwt, a=0.0, siga=0.0):
    #Given a set of data points x, y,
    #with individual standard deviations sig,
    #fit them to a straight line y = a + bx by minimizing chi-sq,
    #where the intercept of the line is fixed.
    #Returned are b and its probable uncertainty sigb,
    #the chi-square chi2, and the scatter sigdat.
    #If mwt=0 on input, then the standard deviations are assumed to be unavailable:
    #the normalization of chi2 is to unit standard deviation on all points.
    sx=0.0
    sy=0.0
    sxx=0.0
    sxy=0.0
    b=0.0
    ndata = len(x)
    if mwt:
        #Accumulate sums ...
        ss = 0.0
        for i in range(ndata):
            #...with weights
            wt=1.0/sig[i]**2
            ss += wt
            sx += x[i]*wt
            sy += y[i]*wt
            sxx += x[i]*x[i]*wt
            sxy += x[i]*y[i]*wt
    else:
        for i in range(ndata):
            #...or without weights.
            sx += x[i]
            sy += y[i]
            sxx += x[i]*x[i]
            sxy += x[i]*y[i]
        ss=ndata
    #Solve for b and sigb.
    b = (sxy - a*sx) / sxx
    sigb = sqrt(1.0 / sxx + (siga * sx / sxx)**2)
    chi2=0.0
    #Calculate chi2.
    if mwt:
        for i in range(ndata):
            chi2 += ((y[i] - a - b*x[i]) / sig[i])**2
    else:
        #For unweighted data evaluate typical sig using chi2,
        #and adjust the standard deviations.
        for i in range(ndata):
            chi2 += (y[i] - a - b*x[i])**2
        sigdat=sqrt(chi2/(ndata-1))
        sigb *= sigdat
    return b, sigb, chi2


def fit_intercept(x, y, sig, mwt, b=0.0, sigb=0.0):
    #Given a set of data points x, y,
    #with individual standard deviations sig,
    #fit them to a straight line y = a + bx by minimizing chi-sq,
    #where the intercept of the line is fixed.
    #Returned are a and its probable uncertainty siga,
    #the chi-square chi2, and the scatter sigdat.
    #If mwt=0 on input, then the standard deviations are assumed to be unavailable:
    #the normalization of chi2 is to unit standard deviation on all points.
    sx=0.0
    sy=0.0
    a=0.0
    ndata = len(x)
    if mwt:
        #Accumulate sums ...
        ss = 0.0
        for i in range(ndata):
            #...with weights
            wt=1.0/sig[i]**2
            ss += wt
            sx += x[i]*wt
            sy += y[i]*wt
    else:
        for i in range(ndata):
            #...or without weights.
            sx += x[i]
            sy += y[i]
        ss=ndata
    #Solve for a and siga.
    a = (sy - b*sx) / ss
    siga = sqrt(1.0 / ss + (sigb * sx / ss)**2)
    chi2=0.0
    #Calculate chi2.
    if mwt:
        for i in range(ndata):
            chi2 += ((y[i] - a - b*x[i]) / sig[i])**2
    else:
        #For unweighted data evaluate typical sig using chi2,
        #and adjust the standard deviations.
        for i in range(ndata):
            chi2 += (y[i] - a - b*x[i])**2
        sigdat=sqrt(chi2/(ndata-1))
        siga *= sigdat
    return a, siga, chi2


def fit_i(x, y, sig, int_scat=0.0):
    #Given a set of data points x, y,
    #with individual measurement error standard deviations sig,
    #fit them to a straight line relation y = a + bx
    #with intrinsic scatter int_scat by minimizing chi-sq.
    #Returned are a, b and their respective probable uncertainties siga and sigb,
    #the chi-square chi2, and the scatter sigdat.
    #If mwt=0 on input, then the standard deviations are assumed to be unavailable:
    #the normalization of chi2 is to unit standard deviation on all points.
    sigtot = N.sqrt(sig**2 + int_scat**2)
    return fit(x, y, sigtot, 1)


def fit_slope_i(x, y, sig, int_scat=0.0, a=0.0, siga=0.0):
    #Given a set of data points x, y,
    #with individual standard deviations sig,
    #fit them to a straight line y = a + bx
    #with intrinsic scatter int_scat by minimizing chi-sq.
    #where the intercept of the line is fixed.
    #Returned are b and its probable uncertainty sigb,
    #the chi-square chi2, and the scatter sigdat.
    #If mwt=0 on input, then the standard deviations are assumed to be unavailable:
    #the normalization of chi2 is to unit standard deviation on all points.
    sigtot = N.sqrt(sig**2 + int_scat**2)
    return fit_slope(x, y, sigtot, 1, a, siga)


def fit_intercept_i(x, y, sig, int_scat, b=0.0, sigb=0.0):
    #Given a set of data points x, y,
    #with individual standard deviations sig,
    #fit them to a straight line y = a + bx by minimizing chi-sq,
    #where the intercept of the line is fixed.
    #Returned are a and its probable uncertainty siga,
    #the chi-square chi2, and the scatter sigdat.
    #If mwt=0 on input, then the standard deviations are assumed to be unavailable:
    #the normalization of chi2 is to unit standard deviation on all points.
    sigtot = N.sqrt(sig**2 + int_scat**2)
    return fit_intercept(x, y, sigtot, 1, b, sigb)


def ttest(n1, ave1, var1, n2, ave2, var2):
    df = n1+n2-2
    svar = ((n1-1)*var1+(n2-1)*var2)/df
    t = (ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2))
    prob = betai(0.5*df,0.5,df/(df+t**2))
    return t, prob


def tutest(n1, ave1, var1, n2, ave2, var2):
    t = (ave1-ave2)/sqrt(var1/n1+var2/n2)
    df = (var1/n1+var2/n2)**2/((var1/n1)**2/(n1-1)+(var2/n2)**2/(n2-1))
    prob = betai(0.5*df,0.5,df/(df+t**2))
    return t, prob

def ftest(d1, d2):
    m1 = mean(d1)
    v1 = variance(d1, m=m1)
    m2 = mean(d2)
    v2 = variance(d2, m=m2)
    if v1 > v2: 
	vr = v1/v2
	df1 = len(d1)-1
	df2 = len(d2)-1
    else:
	vr = v2/v1
	df1 = len(d2)-1
	df2 = len(d1)-1
    prob = 2.0*betai(0.5*df2, 0.5*df1, df2/(df2+df1*vr))
    if prob > 1.0: prob = 2.0 - prob
    return prob
		  
def betai(a, b, x):
    if (x < 0.0 or x > 1.0):
        print "Bad x in routine betai"
    if (x == 0.0 or x == 1.0):
        bt=0.0
    else:
        bt = exp(gammln(a+b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.0-x))
    if (x < (a+1.0)/(a+b+2.0)):
        return bt * betacf(a,b,x)/a
    else:
        return 1.0 - bt * betacf(b,a,1.0-x)/b


def betacf(a, b, x):
    qab=a+b
    qap=a+1.0
    qam=a-1.0
    c=1.0
    d=1.0-qab*x/qap
    FPMIN =  1.0e-30
    EPS = 3.0e-7
    if (abs(d) < FPMIN): d=FPMIN
    d=1.0/d
    h=d
    for m in range(1, 101):
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.0+aa*d
        if (abs(d) < FPMIN): d=FPMIN
        c=1.0+aa/c
        if (abs(c) < FPMIN): c=FPMIN
        d=1.0/d
        h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.0+aa*d
        if (abs(d) < FPMIN): d=FPMIN
        c=1.0+aa/c
        if (abs(c) < FPMIN): c=FPMIN
        d=1.0/d
        delta=d*c
        h *= delta
        if (abs(delta-1.0) < EPS):  break
    if (m > 100): print "a or b too big, or MAXIT too small in betacf"
    return h


def gammln(xx):
    cof = [76.18009172947146,-86.50532032941677,
           24.01409824083091,-1.231739572450155,
           0.1208650973866179e-2,-0.5395239384953e-5]
    y=x=xx
    tmp=x+5.5
    tmp -= (x+0.5)*log(tmp)
    ser=1.000000000190015
    for j in range(6):
        y = y+1
        ser += cof[j]/y
    return -tmp+log(2.5066282746310005*ser/x)


def gammp(a, x):
    if (x < 0.0 or a <= 0.0):
        print "Invalid arguments in routine gammp"
    if (x < (a+1.0)):
        gamser, gln = gser(a,x)
        return gamser
    else:
        gammcf, gln = gcf(a,x)
        return 1.0-gammcf


def gser(a, x):
    gln=gammln(a);
    if (x <= 0.0):
        if (x < 0.0):
            print "x less than 0 in routine gser"
        gamser=0.0
        return
    else:
        ap=a
        delta=sum=1.0/a
        for n in range(1, 101):
            ap +=1
            delta *= x/ap
            sum += delta
            if (abs(delta) < abs(sum)*3.0e-7):
                gamser=sum*exp(-x+a*log(x)-(gln))
                break
        if n>100: print "a too large, ITMAX too small in routine gser"
    return gamser, gln


def gcf(a, x):
    gln=gammln(a)
    b=x+1.0-a
    FPMIN = 1.0e-30
    c=1.0/FPMIN
    d=1.0/b
    h=d
    for i in range(1, 101):
        an = -i*(i-a)
        b += 2.0
        d=an*d+b
        if (abs(d) < FPMIN): d=FPMIN
        c=b+an/c
        if (abs(c) < FPMIN): c=FPMIN
        d=1.0/d
        delta=d*c
        h *= delta
        if (abs(delta-1.0) < 3.0e-7): break
    if (i > 100): print "a too large, ITMAX too small in gcf"
    gammcf=exp(-x+a*log(x)-(gln))*h;
    return gammcf, gln


def erff(x):
    e = gammp(0.5,x*x)
    if x < 0.0:  e = -e
    return e 

def erffc(x):
    return 1.0 - erff(x)

def mean(x, s=None):
    # calculates the weighted mean of x with errors s
    if s is None:
        return N.sum(x) / len(x)
    else:
        w = 1.0/s**2
        sumw = N.sum(w)
        sumx = N.sum(w*x)
        return sumx/sumw

def variance(x, s=None, m=None):
    # calculates the weighted variance of x with errors s and mean m
    if m is None:  m = mean(x, s)
    if s is None:
        sumdx2 = N.sum((x - m)**2)
        return sumdx2/len(x)
    else:
        w = 1.0/s**2
        sumw = N.sum(w)
        sumdx2 = N.sum(w * (x - m)**2)
        return sumdx2/sumw

def stderr(x, s, m):
    # calculates the standard error of x with errors s about mean m
    w = 1.0/s**2
    sumw = N.sum(w)
    sumw2 = N.sum(w*w)
    sumdx2 = N.sum(w * (x - m)**2)
    return sqrt(sumdx2/sumw2)
