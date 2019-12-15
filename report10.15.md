

Z- distribution
$$
R_{SFR}\propto \begin{cases} (1+z)^{3.44} & z\leq z_{peak} \\ (1+z_{peak})^{3.44} & z\geq z_{peak} \\\end{cases} \quad \tag{1}
$$

$$
\Theta(\epsilon,z)=\frac{\hat\Gamma(\alpha+2,\epsilon^{\beta}10^{0.15\beta z})}{\Gamma(\alpha+2)} \quad \tag{2}
$$

$$
R_{GRB}(z)=\rho_0 R_{SFR}(z)\Theta(\epsilon,z)  \quad \tag{3}
$$

$$
A=(\int_{-\infin}^{+\infin}R_{GRB}(z)\ dz)^{-1}  \quad \tag{4}
$$


Result (Pdf & mock sample Pdf):

![ZPDF](/Users/dingding%201/Desktop/calculate/10.11/z_pdf.eps)

![zmock](/Users/dingding%201/Desktop/calculate/10.11/mockzpdf1.eps)





$L_{\gamma}$-distribution:
$$
\Phi(L_{\gamma})=\frac{A_{L}}{\sqrt{2\pi}\sigma_{L_{\gamma}}}exp[-\frac{(logL_{\gamma}-logL_{c})^2}{2\sigma^2_{L_{\gamma}}}] \quad \tag{5}
$$

$$
A_{L}=(\int_{10^{46}}^{10^{52}}\frac{1}{\sqrt{2\pi}\sigma_{L_{\gamma}}}exp[-\frac{(logL_{\gamma}-logL_{c})^2}{2\sigma^2_{L_{\gamma}}}]dL_{\gamma})^{-1} \quad \tag{6}
$$


Result (Pdf & mock sample pdf)

![pdfl](/Users/dingding 1/Desktop/calculate/9.18/L-pdf.eps)

![mockl](/Users/dingding 1/Desktop/calculate/10.11/mockLpdf1.eps)

\theta_{j}-distribution:
$$
\Psi(\theta_{j})=\frac{A_{\theta}}{\sqrt{2\pi}\sigma_{\theta}}exp[{-\frac{(log\theta_{j}-log\theta_{c})^2}{2\sigma^2_{\theta}}}] \quad \tag{7}
$$
Result(Pdf & mock sample Pdf):

![thetapdf](/Users/dingding 1/Desktop/calculate/10.11/thetapdf.eps)

![mocktheta](/Users/dingding 1/Desktop/calculate/10.11/mocktheta.eps)


$$
P(z,L_{\gamma},\theta_{j})=\frac{L_{\gamma}}{4\pi D^2_{L}(z)(1-cos\theta_{j})k} \quad \tag{8}
$$

$$
E_{p}\times(1+z)=E'_{P}=200keV(L/10^{52})^{1/2}/C \quad \tag{9}
$$

$$
k=\frac{\int_{1/(1+z)}^{10^4/(1+z)}EN(E)dE}{\int_{E_{min}}^{E_{max}}N(E)dE} \quad \tag{10}
$$

$$
\eta_t(P)=\begin{cases} P^2 & P <\ 0.45 \quad photons \ cm^{-2}\ s^{-1} \\ 0.67(1.0-0.4/P)^{0.52}, & P\geq \ 0.45 \quad photons \ cm^{-2}\ s^{-1} \end{cases} \quad \tag{11}
$$

