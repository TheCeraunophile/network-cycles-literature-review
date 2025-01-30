import csv
from sympy import *
import networkx as nx
import matplotlib.pyplot as plt
from math import *
from random import *
import csv
from functools import reduce

def plfit(x, *varargin):
    vec     = []
    sample  = []
    xminx   = []
    limit   = []
    finite  = False
    nosmall = False
    nowarn  = False

    # parse command-line parameters trap for bad input
    i=0
    while i<len(varargin):
        argok = 1
        if type(varargin[i])==str:
            if varargin[i]=='range':
                Range = varargin[i+1]
                if Range[1]>Range[0]:
                    argok=0
                    vec=[]
                try:
                    vec=[X*float(Range[2])+Range[0] for X in range(int((Range[1]-Range[0])/Range[2]))]


                except:
                    argok=0
                    vec=[]


                if Range[0]>=Range[1]:
                    argok=0
                    vec=[]
                    i-=1

                i+=1


            elif varargin[i]== 'sample':
                sample  = varargin[i+1]
                i = i + 1
            elif varargin[i]==  'limit':
                limit   = varargin[i+1]
                i = i + 1
            elif varargin[i]==  'xmin':
                xminx   = varargin[i+1]
                i = i + 1
            elif varargin[i]==  'finite':       finite  = True
            elif varargin[i]==  'nowarn':       nowarn  = True
            elif varargin[i]==  'nosmall':      nosmall = True
            else: argok=0


        if not argok:
            print('(PLFIT) Ignoring invalid argument #',i+1)

        i = i+1

    if vec!=[] and (type(vec)!=list or min(vec)<=1):
        print('(PLFIT) Error: ''range'' argument must contain a vector or minimum <= 1. using default.\n')

        vec = []

    if sample!=[] and sample<2:
        print('(PLFIT) Error: ''sample'' argument must be a positive integer > 1. using default.\n')
        sample = []

    if limit!=[] and limit<min(x):
        print('(PLFIT) Error: ''limit'' argument must be a positive value >= 1. using default.\n')
        limit = []

    if xminx!=[] and xminx>=max(x):
        print('(PLFIT) Error: ''xmin'' argument must be a positive value < max(x). using default behavior.\n')
        xminx = []



    # select method (discrete or continuous) for fitting
    if     reduce(lambda X,Y:X==True and floor(Y)==float(Y),x,True): f_dattype = 'INTS'
    elif reduce(lambda X,Y:X==True and (type(Y)==int or type(Y)==float or type(Y)==int),x,True):    f_dattype = 'REAL'
    else:                 f_dattype = 'UNKN'

    if f_dattype=='INTS' and min(x) > 1000 and len(x)>100:
        f_dattype = 'REAL'


    # estimate xmin and alpha, accordingly

    if f_dattype== 'REAL':
        xmins = unique(x)
        xmins.sort()
        xmins = xmins[0:-1]
        if xminx!=[]:

            xmins = [min([X for X in xmins if X>=xminx])]


        if limit!=[]:
            xmins=[X for X in xmins if X<=limit]
            if xmins==[]: xmins = [min(x)]

        if sample!=[]:
            step = float(len(xmins))/(sample-1)
            index_curr=0
            new_xmins=[]
            for i in range (0,sample):
                if round(index_curr)==len(xmins): index_curr-=1
                new_xmins.append(xmins[int(round(index_curr))])
                index_curr+=step
            xmins = unique(new_xmins)
            xmins.sort()



        dat   = []
        z     = sorted(x)

        for xm in range(0,len(xmins)):
            xmin = xmins[xm]
            z    = [X for X in z if X>=xmin]

            n    = len(z)
            # estimate alpha using direct MLE

            a    = float(n) / sum([log(float(X)/xmin) for X in z])
            if nosmall:
                if (a-1)/sqrt(n) > 0.1 and dat!=[]:
                    xm = len(xmins)+1
                    break


            # compute KS statistic
            #cx   = map(lambda X:float(X)/n,range(0,n))
            cf   = [1-pow((float(xmin)/X),a) for X in z]
            dat.append( max( [abs(cf[X]-float(X)/n) for X in range(0,n)]))
        D     = min(dat)
        xmin  = xmins[dat.index(D)]
        z     = [X for X in x if X>=xmin]
        z.sort()
        n     = len(z)
        alpha = 1 + n / sum([log(float(X)/xmin) for X in z])
        if finite: alpha = alpha*float(n-1)/n+1./n  # finite-size correction
        if n < 50 and not finite and not nowarn:
            print('(PLFIT) Warning: finite-size bias may be present.\n')

        L = n*log((alpha-1)/xmin) - alpha*sum([log(float(X)/xmin) for X in z])
    elif f_dattype== 'INTS':

        x=list(map(int,x))
        if vec==[]:
            for X in range(150,351):
                vec.append(X/100.)    # covers range of most practical
                                    # scaling parameters
        zvec = list(map(zeta, vec))

        xmins = unique(x)
        xmins.sort()
        xmins = xmins[0:-1]
        if xminx!=[]:
            xmins = [min([X for X in xmins if X>=xminx])]

        if limit!=[]:
            limit = round(limit)
            xmins=[X for X in xmins if X<=limit]
            if xmins==[]: xmins = [min(x)]

        if sample!=[]:
            step = float(len(xmins))/(sample-1)
            index_curr=0
            new_xmins=[]
            for i in range (0,sample):
                if round(index_curr)==len(xmins): index_curr-=1
                new_xmins.append(xmins[int(round(index_curr))])
                index_curr+=step
            xmins = unique(new_xmins)
            xmins.sort()

        if xmins==[]:
            print('(PLFIT) Error: x must contain at least two unique values.\n')
            alpha = 'Not a Number'
            xmin = x[0]
            D = 'Not a Number'
            return [alpha,xmin,D]

        xmax   = max(x)

        z      = x
        z.sort()
        datA=[]
        datB=[]

        for xm in range(0,len(xmins)):
            xmin = xmins[xm]
            z    = [X for X in z if X>=xmin]
            n    = len(z)
            # estimate alpha via direct maximization of likelihood function

            # force iterative calculation
            L       = []
            slogz   = sum(map(log,z))
            xminvec = list(map(float,list(range(1,xmin))))
            for k in range(0,len(vec)):
                L.append(-vec[k]*float(slogz) - float(n)*log(float(zvec[k]) - sum([pow(float(X),-vec[k]) for X in xminvec])))


            I = L.index(max(L))
            # compute KS statistic
            fit = reduce(lambda X,Y: X+[Y+X[-1]],\
                         ([pow(X,-vec[I])/(float(zvec[I])-sum([pow(X,-vec[I]) for X in list(map(float,list(range(1,xmin))))])) for X in range(xmin,xmax+1)]),[0])[1:]
            cdi=[]
            for XM in range(xmin,xmax+1):
                cdi.append(len([X for X in z if floor(X)<=XM])/float(n))

            datA.append(max( [abs(fit[X] - cdi[X]) for X in range(0,xmax-xmin+1)]))
            datB.append(vec[I])
        # select the index for the minimum value of D
        I = datA.index(min(datA))
        xmin  = xmins[I]
        z     = [X for X in x if X>=xmin]
        n     = len(z)
        alpha = datB[I]
        if finite: alpha = alpha*(n-1.)/n+1./n  # finite-size correction
        if n < 50 and not finite and not nowarn:
            print('(PLFIT) Warning: finite-size bias may be present.\n')

        L     = -alpha*sum(map(log,z)) - n*log(zvec[vec.index(max([X for X in vec if X<=alpha]))] - \
                                              sum([pow(X,-alpha) for X in range(1,xmin)]))
    else:
        print('(PLFIT) Error: x must contain only reals or only integers.\n')
        alpha = []
        xmin  = []
        L     = []

    return [alpha,xmin,L]


def zeta(s):
    """
    Riemann zeta function, real argument
    """
    if not isinstance(s, (float, int)):
        try:
            s = float(s)
        except (ValueError, TypeError):
            try:
                s = complex(s)
                if not s.imag:
                    return complex(zeta(s.real))
            except (ValueError, TypeError):
                pass
            raise NotImplementedError
    if s == 1:
        raise ValueError("zeta(1) pole")
    if s >= 27:
        return 1.0 + 2.0**(-s) + 3.0**(-s)
    n = int(s)
    if n == s:
        if n >= 0:
            return _zeta_int[n]
        if not (n % 2):
            return 0.0
    if s <= 0.0:
        return 0
    if s <= 2.0:
        if s <= 1.0:
            return _polyval(_zeta_0,s)/(s-1)
        return _polyval(_zeta_1,s)/(s-1)
    z = _polyval(_zeta_P,s) / _polyval(_zeta_Q,s)
    return 1.0 + 2.0**(-s) + 3.0**(-s) + 4.0**(-s)*z


def plplot(x,xmin,alpha):
        # select method (discrete or continuous) for fitting
    if     reduce(lambda X,Y:X==True and floor(Y)==float(Y),x,True): f_dattype = 'INTS'
    elif reduce(lambda X,Y:X==True and (type(Y)==int or type(Y)==float or type(Y)==int),x,True):    f_dattype = 'REAL'
    else:                 f_dattype = 'UNKN'

    if f_dattype=='INTS' and min(x) > 1000 and len(x)>100:
        f_dattype = 'REAL'
    plt.close()
    plt.ion()
    h=[[],[]]
    # estimate xmin and alpha, accordingly
    if f_dattype== 'REAL':
        n = len(x)
        c1 = sorted(x)
        c2 = [X/float(n) for X in range(n,0,-1)]
        q = sorted([X for X in x if X>=xmin])
        cf = [pow(float(X)/xmin,1.-alpha) for X in q]
        cf = [X*float(c2[c1.index(q[0])]) for X in cf]

        h[0]=plt.loglog(c1, c2, 'bo',markersize=8,markerfacecolor=[1,1,1],markeredgecolor=[0,0,1])
        h[1]=plt.loglog(q, cf, 'k--',linewidth=2)

        xr1 = pow(10,floor(log(min(x),10)))
        xr2 = pow(10,ceil(log(min(x),10)))
        yr1 = pow(10,floor(log(1./n,10)))
        yr2 = 1


        plt.axhspan(ymin=yr1,ymax=yr2,xmin=xr1,xmax=xr2)
        plt.ylabel('Pr(X >= x)',fontsize=16);
        plt.xlabel('x',fontsize=16)
        plt.draw()

    elif f_dattype== 'INTS':
        n = len(x)
        q = sorted(unique(x))
        c=[]
        for Q in q:
            c.append(len([X for X in x if floor(X)==Q])/float(n))
        c1 = q+[q[-1]+1]
        c2 = [1.-Z for Z in reduce(lambda X,Y: X+[Y+X[-1]],c,[0])]
        c2 = [X for X in c2 if float(X)>=pow(10,-10.)]
        c1 = c1[0:len(c2)]
        print("alpha:", alpha, type(alpha))
        print("xmin:", xmin, type(xmin))
        print("q:", q, type(q), "Contents:", q)
        cf = [pow(X,-alpha)/(float(zeta(alpha)) - sum([pow(Y,-alpha) for Y in range(1,xmin)])) for X in range(xmin,q[-1]+1)]
        cf1 = list(range(xmin,q[-1]+2))
        cf2 = [1.-Z for Z in reduce(lambda X,Y: X+[Y+X[-1]],cf,[0])]
        cf2 = [X*float(c2[c1.index(xmin)]) for X in cf2]



        h[0]=plt.loglog(c1, c2, 'bo',markersize=8,markerfacecolor=[1,1,1],markeredgecolor=[0,0,1])
#         print c1
#         print c2
        h[1]=plt.loglog(cf1, cf2, 'k--',linewidth=2)
#         print cf1
#         print cf2

        xr1 = pow(10,floor(log(min(x),10)))
        xr2 = pow(10,ceil(log(min(x),10)))
        yr1 = pow(10,floor(log(1./n,10)))
        yr2 = 1


        #plt.axhspan(ymin=yr1,ymax=yr2,xmin=xr1,xmax=xr2)
        #plt.ylabel('Pr(X >= x)',fontsize=16);
        #plt.xlabel('x',fontsize=16)
        #plt.draw()
    return h


def unique(seq):
    # not order preserving
    result_set = {}
    for item in seq:
        result_set[item] = None
    return list(result_set.keys())


def _polyval(coeffs, x):
    p = coeffs[0]
    for c in coeffs[1:]:
        p = c + x*p
    return p

_zeta_int = [\
-0.5,
0.0,
1.6449340668482264365,1.2020569031595942854,1.0823232337111381915,
1.0369277551433699263,1.0173430619844491397,1.0083492773819228268,
1.0040773561979443394,1.0020083928260822144,1.0009945751278180853,
1.0004941886041194646,1.0002460865533080483,1.0001227133475784891,
1.0000612481350587048,1.0000305882363070205,1.0000152822594086519,
1.0000076371976378998,1.0000038172932649998,1.0000019082127165539,
1.0000009539620338728,1.0000004769329867878,1.0000002384505027277,
1.0000001192199259653,1.0000000596081890513,1.0000000298035035147,
1.0000000149015548284]

_zeta_P = [-3.50000000087575873, -0.701274355654678147,
-0.0672313458590012612, -0.00398731457954257841,
-0.000160948723019303141, -4.67633010038383371e-6,
-1.02078104417700585e-7, -1.68030037095896287e-9,
-1.85231868742346722e-11][::-1]

_zeta_Q = [1.00000000000000000, -0.936552848762465319,
-0.0588835413263763741, -0.00441498861482948666,
-0.000143416758067432622, -5.10691659585090782e-6,
-9.58813053268913799e-8, -1.72963791443181972e-9,
-1.83527919681474132e-11][::-1]

_zeta_1 = [3.03768838606128127e-10, -1.21924525236601262e-8,
2.01201845887608893e-7, -1.53917240683468381e-6,
-5.09890411005967954e-7, 0.000122464707271619326,
-0.000905721539353130232, -0.00239315326074843037,
0.084239750013159168, 0.418938517907442414, 0.500000001921884009]

_zeta_0 = [-3.46092485016748794e-10, -6.42610089468292485e-9,
1.76409071536679773e-7, -1.47141263991560698e-6, -6.38880222546167613e-7,
0.000122641099800668209, -0.000905894913516772796, -0.00239303348507992713,
0.0842396947501199816, 0.418938533204660256, 0.500000000000000052]


def Average(lst):
    return sum(lst) / len(lst)

def degree_sequence(g):
    degree_sequence=sorted(list(dict(nx.degree(g)).values()),reverse=True) # degree sequence
    #print "degree_sequence", degree_sequence
    return degree_sequence

def degree_distribution(g,name):
    distinct_degrees=[]
    degreeSeqForPL=[]
    degree_count=[]
    N=g.number_of_nodes()
    deg_seq=degree_sequence(g)
    for ds in deg_seq:
        if ds not in distinct_degrees:
            distinct_degrees.append(ds)
    print("distDeg=",distinct_degrees)
    for i in distinct_degrees:
        degree_count.append(deg_seq.count(i))
    print("degCount=",degree_count)
    yl=[]
############################################## Emadi #############################################
    #print("distDeg",distinct_degrees)
    #print ("[",)
    temp = len(distinct_degrees)
    for ii in range(0, temp):
        bardia = distinct_degrees[ii]
        for iii in range(0, degree_count[ii]):
            #print (distinct_degrees[ii],",",)
            degreeSeqForPL.append(distinct_degrees[ii])
    #print ("]")
    #print ("degCount",degree_count)
    #print ("degree_distribution",degreeSeqForPL)

############################################## Emadi #############################################
    for i in range(max(degree_count)+2):
        yl.append(i)
    #print(yl)
    width = 0.5
    fig, ax = plt.subplots(figsize=(15,8), dpi = 1000)
    rects1 = ax.bar(distinct_degrees, degree_count, width, data=None, color='black')
    rects2 = ax.plot(distinct_degrees, degree_count, width, data=None, color='gray')
    ax.set_title("Stepwise Delta=7", weight='bold', size=20)
    ax.set_ylabel("count", weight='bold', size=17)
    ax.set_xlabel("degrees", weight='bold', size=13)
    ax.set_xticks(distinct_degrees)
    ax.set_yticks(degree_count)
    #ax.set_xlim(0,45)
    #ax.set_ylim(0,200)
    for tick in ax.get_xticklabels():
        tick.set_rotation(-90)
    plt.xticks(weight='bold', size=13)
    plt.yticks(weight='bold', size=15)
    #plt.savefig("degree_Dist_"+name+".png")


def main(graph):
    deglist = degree_sequence(graph)
    print("degree_sequence= ", deglist)
    # degree_distribution(g,"BA.2_network_relaxed")

    average = Average(degree_sequence(graph))
    print("Average Degree =", round(average, 2))
    xxx = deglist
    [alpha, xmin, L] = plfit(xxx)

    h = plplot(xxx,xmin,alpha)
    # plt.title("Power Law Stepwise Delta=4", family= 'times new roman', weight='bold', size=20)
    plt.xlabel("log(k)", weight='bold', size=20)
    plt.ylabel("log(P(k))", weight='bold', size=20)
    plt.xticks(weight='bold', size=18)
    plt.yticks(weight='bold', size=18)
    plt.savefig("power_law_Stepwise_Delta_10.png", dpi = 1000, bbox_inches='tight')
    plt.show()
    print("Alpha=",alpha)
    print("Xmin=",xmin)
    print("L=",L)
