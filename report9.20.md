




$$
R_{SFR}\propto \begin{cases} (1+z)^{3.44} & z\leq z_{peak} \\ (1+z_{peak})^{3.44} & z\geq z_{peak} \\\end{cases} \quad \tag{1}
$$

~~~python
    # define the star formation rate segment function:
def RSFR(z):
    zpeak=1
    if z<=zpeak:
        Rsfr=(1+z)**(3.44)
        return Rsfr
    elif z>=zpeak:
        Rsfr=(1+zpeak)**(3.44)
        return Rsfr
~~~

$$
\Theta(\epsilon,z)=\frac{\hat\Gamma(\alpha+2,\epsilon^{\beta}10^{0.15\beta z})}{\Gamma(\alpha+2)} \quad \tag{2}
$$





~~~python
    # define the complete gamma function:
def comGammaFunc(v):
    gamma=quad(lambda t:t**(v-1)*np.e**(-t),0,float("inf"))
    return gamma[0]
    #define the incomplete gamma function:
def incomGammaFunc(v,z):
    sgamma=quad(lambda u:u**(v-1)*np.e**(-u),0,z)[0]
    bgamma=quad(lambda u:u**(v-1)*np.e**(-u),z,float('inf'))[0]
    gamma=comGammaFunc(v)-bgamma
    return gamma,bgamma,sgamma
    #and define the Seitafunction:
def SeitaFunc(eps,z,alpha,beta):
    Seita1=incomGammaFunc(alpha+2,eps**beta*10**(0.15*beta*z))[0]
    Seita2=comGammaFunc(alpha+2)
    Seita=Seita1/Seita2
    return Seita
~~~



<center>$R_{GRB}(z)=\rho_0 R_{SFR}(z)\Theta(\epsilon,z)  \quad \tag{3}$

</center>

~~~python
    # define the grb rate function:
def RGRB(z,eps,alpha,beta,rho):
    A=1/(33.30270146296203)
    RGRB=A*rho*RSFR(z)*SeitaFunc(eps,z,alpha,beta)
    return RGRB
~~~

In the code above, A is the normalized constant,

<center>$A=(\int_{-\infin}^{+\infin}R_{GRB}(z)\ dz)^{-1}  \quad \tag{4}$

</center>





<center>$\Phi(L_{\gamma})=\frac{A_{L}}{\sqrt{2\pi}\sigma_{L_{\gamma}}}exp[-\frac{(logL_{\gamma}-logL_{c})^2}{2\sigma^2_{L_{\gamma}}}] \quad \tag{5}$</center>

Where $A_{L}​$ is normalized constant,

<center>$A_{L}=(\int_{0}^{\infin}\frac{1}{\sqrt{2\pi}\sigma_{L_{\gamma}}}exp[-\frac{(logL_{\gamma}-logL_{c})^2}{2\sigma^2_{L_{\gamma}}}]dL_{\gamma})^{-1} \quad \tag{6}$</center>

~~~python
# Here during the defination, normalized constant A_{L} is included:
def Luminosityfunction(L_gamma):
    L_critical=10**(49.69)      #unit is erg
    sigma_L=0.4
    A_L=1/(1.7235434382660358e+50)
    luminosityfunc=A_L*np.exp(-(np.log10(L_gamma)-np.log10(
        L_critical))**2/(2*sigma_L**2))/(np.sqrt(2*np.pi)*sigma_L)
    return luminosityfunc
~~~



<center>$\Psi(\theta_{j})=\frac{A_{\theta}}{\sqrt{2\pi}\sigma_{\theta}}exp[{-\frac{(log\theta_{j}-log\theta_{c})^2}{2\sigma^2_{\theta}}}] \quad \tag{7}$</center>

~~~python
# Define the angle distribution as log-normal distribution:
def thetalogdistri(theta_jet):
    theta_critical=10**(-1.27)
    sigema_theta=0.6
    A_theta=1/0.32112249370542306
    Psi=A_theta*np.exp(-(np.log10(theta_jet)-np.log10(theta_critical))**2/
               (2*sigema_theta**2))/(np.sqrt(2*np.pi)*sigema_theta)
    return Psi
~~~



<center>$P(z,L_{\gamma},\theta_{j})=\frac{L_{\gamma}}{4\pi D^2_{L}(z)(1-cos\theta_{j})k} \quad \tag{8}$</center>

<center>$E_{p}\times(1+z)=E'_{P}=200keV(L/10^{52})^{1/2}/C \quad \tag{9}$</center>

~~~python
# Define peak flux P:
def P(z,L_gamma,theta_jet):
    L=L_gamma/(1-np.cos(theta_jet/180*np.pi))
    C=random.uniform(0.1,1)
    ep=200*(L/10**52)**0.5/C/(1+z)
    P=L/(4*np.pi*dl(z)**2*nk(ep,z,-1.1,-2.2,15,150))
    return P
~~~

<center>$k=\frac{\int_{1/(1+z)}^{10^4/(1+z)}EN(E)dE}{\int_{E_{min}}^{E_{max}}N(E)dE} \quad \tag{10}$</center>

~~~python
# calculate the k-corrention in photons.s-1.cm-2:
def nk(Epeak,Z,alpha,beita,bandmin,bandmax):
    a1=quad(lambda E:E*NE(E,Epeak,alpha,beita),1/(1+Z),10**4/(1+Z))
    a2=quad(lambda E:NE(E,Epeak,alpha,beita),bandmin,bandmax)
    k=a1[0]/a2[0]
    # return k
    return k*1.6*10**(-9)          #transform kev to erg
~~~

<center>$\eta_t(P)=\begin{cases} P^2 & P <\ 0.45 \quad photons \ cm^{-2}\ s^{-1} \\ 0.67(1.0-0.4/P)^{0.52}, & P\geq \ 0.45 \quad photons \ cm^{-2}\ s^{-1} \end{cases} \quad \tag{11}$</center>

~~~python
# BAT trigger probability:
def eta_t(P):
    if P<0.45:
        eta_t=P**2
        return eta_t
    elif P>=0.45:
        eta_t=0.67*(1.0-0.4/P)**0.52
        return eta_t
~~~



<center>$\eta_{z}(P)=0.26+0.032e^{1.61logP} \quad \tag{12}$</center>

~~~python
# weak dependence of probability on the observed peak flux:
def eta_z(P):
    eta_z=0.26+0.032*np.e**(1.61*np.log10(P))
    return eta_z
~~~



<center>$\eta_{a}(\theta_{j})=\frac{\Omega}{4\pi}(1-cos\theta_{j}) \quad \tag{13}$</center>

~~~python
def eta_a(theta_jet):
    eta_a=1.4*(1-np.cos(theta_jet))/(4*np.pi)      #where 1.4 sr is instrument solid angle
    return eta_a
~~~



**Monte Carlo mimic & selection effect:**

~~~python
#set the parameters:
eps=0.4
alpha=-1.1
beta=-2.2
rho=1
samplenumber=170

#---------------------------------------------------------------------------------------
zmocksample=[]
lmocksample=[]
thetamocksample=[]
lmin=46
lmax=52
P_mock=[]
i=0
while i<samplenumber:
    a=random.uniform(10**(-2.5),10**(0.))
    b=random.uniform(0,2.1)
    theta=a/np.pi*180
    if b<=thetalogdistri(a):
        c=random.uniform(0,1)
        if eta_a(a)>c:
            while True:
                a=random.uniform(0.01, 10)
                b=random.uniform(0,2.5)
                z=a
                if b<=RGRB(a,eps,alpha,beta,rho):
                    break
                else:
                    continue
            while True:
                a=random.uniform(10.**lmin,10.**lmax)
                b=random.uniform(0,6*10.**(-51))
                l=a
                if b<=Luminosityfunction(a):
                    break
                else:
                    continue
            p=P(z,l,theta)
            if p<80:  #annexed condition
                a=random.uniform(0,1)
                if 0<a<eta_t(p):
                    a=random.uniform(0,1)
                    if 0<a<eta_z(p):
                        P_mock=np.append(P_mock,p)
                        zmocksample=np.append(zmocksample,z)
                        thetamocksample=np.append(thetamocksample,theta)
                        lmocksample=np.append(lmocksample,l/(1-np.cos(theta/180*np.pi)))
                        i=i+1
~~~

Get the probability distribution of mock sample:

~~~python
zmock= zmocksample
lmock=lmocksample
thetamock=thetamocksample
pmock=P_mock

proz=[]
prol=[]
protheta=[]
prothetaz=[]
prop=[]


#---------------------------------------------------------------------------------------
# for z:
for i in range(20):
    counter=0
    for j in range(170):
        if 10**(i*0.05)-1<=zmock[j]<=10**((i+1)*0.05)-1:
            counter=counter+1
    proz=np.append(proz,counter)
proz=proz/np.sum(proz)

# --------------------------------------------------------------------------------------
# for luminosity:
for i in range(20):
    counter=0
    for j in range(170):
        if 10.**(49+i*0.3)<=lmock[j]<=10.**(49+(i+1)*0.3):
            counter=counter+1
    prol=np.append(prol,counter)
prol=prol/np.sum(prol)

# --------------------------------------------------------------------------------------
# for cumulate p:
num=40
for i in range(num):
    counter=0
    for j in range(samplenumber):
        if pmock[j]>=10**(-1+3/num*i):
            counter=counter+1
    prop=np.append(prop,counter)
prop=prop/prop[0]

# --------------------------------------------------------------------------------------
# randomly pick 77 grbs subsample from mock sample with 170 grbs:
i=0
subsamplez=[]
subsampletheta=[]
min=-2.5
max=0
while i<77:
    a=random.randint(0,169)
    subsamplez=np.append(subsamplez,zmock[a])
    subsampletheta=np.append(subsampletheta,thetamock[a])
    i=i+1

for i in range(10):
    counter=0
    for j in range(77):
        if 10.**(-2.5+0.25*i)<=subsampletheta[j]/180*np.pi<=10.**(-2.5+0.25*(i+1)):
            counter=counter+1
    protheta=np.append(protheta,counter)
protheta=protheta/np.sum(protheta)

for i in range(10):
    counter=0
    for j in range(77):
        if 10**(i*0.1)-1<=subsamplez[j]<=10**((i+1)*0.1)-1:
            counter=counter+1
    prothetaz=np.append(prothetaz,counter)
prothetaz=prothetaz/np.sum(prothetaz)

~~~

![Z-pdf](/Users/dingding/Desktop/calculate/9.18/z_pdf.eps)

![mockz](/Users/dingding/Desktop/calculate/9.18/mockz.eps)



![L-pdf](/Users/dingding/Desktop/calculate/9.18/L-pdf.eps)

![mockL](/Users/dingding/Desktop/calculate/9.18/mockl.eps)

![theta-pdf](/Users/dingding/Desktop/calculate/9.18/theta_pdf.eps)



![mocktheta](/Users/dingding/Desktop/calculate/9.18/mocktheta.eps)

![thetapdf2](/Users/dingding/Desktop/calculate/9.18/conditiontheta_pdf.eps)

![mocktheta2](/Users/dingding/Desktop/calculate/9.18/conditionmocktheta.eps)





**Result **

![](/Users/dingding/Desktop/calculate/9.17//same/z9.17.eps)

![](/Users/dingding/Desktop/calculate/9.17//same/luminosity9.17.eps)

![](/Users/dingding/Desktop/calculate/9.17//same/z-ldistribution.eps)

![](/Users/dingding/Desktop/calculate/9.17//same/theta9.17.eps)

![](/Users/dingding/Desktop/calculate/9.17//same/zdistributioninsubsample9.17.eps)

![](/Users/dingding/Desktop/calculate/9.17//same/z-thetadistribution9.17.eps)





