# Model definitions and manuscript correspondence

## GWTC-4 BP2P

The official posterior provides a joint vector for each realization j:

`alpha_1, alpha_2, m_break, delta_m1, delta_m2, lambda_0, lambda_1,
m1_low, m2_low, m_high, mu_1, mu_2, sigma_1, sigma_2, beta_q,
rate_0, kappa, mu_chi, sigma_chi, mu_tilt, sigma_tilt, zeta_tilt`.

Select complete rows, not independently sampled marginal distributions. The
historical loader requires finite mapped parameters and finite `variance < 1`;
3,962 posterior rows pass for the recorded input. Select 500 without replacement
with NumPy `default_rng(20260918)`, then assign uint64 child seeds from that same
root generator. Hyperparameters stay fixed within each catalog. Source variables
are redrawn conditionally on that vector.

The mass model (manuscript A1–A2) is

\[
p(m_1|\Lambda)\propto[\lambda_0 B(m_1)+\lambda_1 G_1(m_1)
+(1-\lambda_0-\lambda_1)G_2(m_1)]S(m_1|m_{1,\min},\delta m_1).
\]

The continuous broken power law is proportional to
`(m1/m_break)**(-alpha_1)` below the break and
`(m1/m_break)**(-alpha_2)` above. Each component is normalized on
`[m1_low,m_high]` before the common taper; the final mixture is normalized again.
The posterior used here fixes `m_high=300` solar masses. It does not imply every
source reaches that mass or that this is a detector-frame mass limit.

\[
p(q|m_1,\Lambda)\propto q^{\beta_q}S(m_1q|m_{2,\min},\delta m_2),
\quad m_2=q m_1,\quad m_{2,\min}/m_1\leq q\leq1.
\]

The Planck taper is zero below the lower bound, one above lower bound + width,
and `1/[1+exp(width/x+width/(x-width))]` between them. The current implementation
uses 20,000 primary-mass points and 1,024 logarithmic primary-mass bins with 2,048
conditional q points, selecting the nearest primary-mass bin in log space.

Spin magnitudes are iid truncated Gaussians on [0,1]. A shared Bernoulli mixture
assignment per binary selects two truncated-Gaussian tilt cosines on [-1,1]
with probability `zeta_tilt`, or two uniform cosines otherwise. Spin azimuths
are uniform on [0,2pi]. The Gaussian parameters need not correspond to a peak
exactly aligned at cos(tilt)=1; they are taken from the posterior.

For low z (A10–A11):

\[
R_{\mathrm{low},j}(z)=R_{0,j}(1+z)^{\kappa_j},\quad
w_j(z)=\frac{(1+z)^{\kappa_j}}{1+z}\frac{dV_c}{dz}.
\]

The Poisson mean is `R0_j * years * integral(w_j dz)`. Here R0 has units
Gpc^-3 source-year^-1, duration is in observer years, and the volume is full sky.
Low-z sampling uses a 20,000-point linear-z grid, trapezoidal integration, and
inverse-CDF interpolation. Planck18 supplies dVc/dz with the per-steradian
Astropy result multiplied by 4pi sr and converted to Gpc^3.

## Manuscript high-z extension (A12–A20)

Set gamma=3.2, z_turn=1.9, kappa_high=5.6, z_join=1.5, z_max=30:

\[
g(z)=\frac{(1+z)^{3.2}}{1+[(1+z)/2.9]^{8.8}},\quad
h(z)=\begin{cases}(1+z)^{3.2},&z\leq1.5,\\
2.5^{3.2}g(z)/g(1.5),&z>1.5.\end{cases}
\]

The **reference** h is continuous at 1.5 and asymptotes to (1+z)^-5.6.
The two branches need not have equal derivatives. z_turn is a turnover parameter,
not the exact rate maximum. The old whole-range Madau–Dickinson resampling version
and the later 100-catalog fresh-realization experiment are not the manuscript preset.

Define the shape integrals, with R0 factored out:

\[
I_{\rm low}=\int_0^{1.5}\frac{h(z)}{1+z}\frac{dV_c}{dz}dz,
\quad I_{\rm high}=\int_{1.5}^{30}\frac{h(z)}{1+z}\frac{dV_c}{dz}dz.
\]

Their units are Gpc^3; retaining R0 instead would give rate integrals and the
same ratio. With Planck18 they are 1822.0072111311952 and 4102.560973818097,
giving ratio 2.2516710958959476. For each retained low catalog:

\[
N_{{\rm high},j}|N_{{\rm old},j}\sim
\mathrm{Poisson}(N_{{\rm old},j} I_{\rm high}/I_{\rm low}).
\]

This explicitly conditions the high count on the **realized** low count. It is
not the independent Poisson count from R0_j times an absolute high-z volume.
It produces additional count covariance between low and high components.

Sample high-z redshifts from `h(z)/(1+z) * dVc/dz`, normalized only on (1.5,30],
using 30,000 points uniform in log(1+z). Trapezoid integration is with respect
to z on that nonuniform grid, not an integral with respect to log(1+z).
All high-z sources use the matching low catalog's mass and spin posterior draw.
The extension root seed is 20260920; it supplies 500 uint64 child seeds in catalog
order. Each worker first draws its Poisson count, then m1, q, spins, redshift,
time, sky position and orientation, retaining the historical RNG ordering.

**Reference versus actual low-z evolution:** retained rows use kappa_j, not
fixed gamma=3.2. Continuity is a statement about h(z); it is not a demonstrated
per-realization continuity of `R0_j (1+z)^kappa_j` into the appended distribution.
The paper's Appendix A makes this reference-rate distinction. A future model
matched separately to every kappa_j would be a scientific change and would not
reproduce the published ensemble prescription.

## Stored fields and frame convention

| Columns | Definition / units |
| --- | --- |
| z | Cosmological redshift |
| m1_Msun, m2_Msun | Source-frame component masses, solar masses |
| q | m2/m1 |
| a_1, a_2 | Dimensionless spin magnitudes |
| cos_tilt_1, cos_tilt_2 | Spin tilts relative to orbital angular momentum |
| phi_1_rad, phi_2_rad | Spin azimuths, radians; not orbital/coalescence phase |
| t_c_detector_yr | Uniform in [0,Tcat], observer-frame years |
| ra_rad | Uniform in [0,2pi], radians |
| sin_dec, cos_iota | Uniform in [-1,1] |
| psi_rad | Uniform in [0,pi], radians |

Derive detector masses as (1+z)*source mass, distance with
`Planck18.luminosity_distance(z)`, and the waveform's aligned spin projections
as `chi1L=a_1*cos_tilt_1`, `chi2L=a_2*cos_tilt_2` (A8). No orbital reference
phase, distance, detector masses, source ID or selection weight is stored in the
historical 15-column CSV. In combined files the first N_old rows are original;
subsequent rows are high-z. Preserve that order if component identity is needed.

Planck18 uses H0=67.66 km/s/Mpc, Om0=0.30966, Ob0=0.04897, Tcmb0=2.7255 K,
Neff=3.046 and neutrino masses [0,0,0.06] eV. Use Astropy's named cosmology to
retain its radiation and neutrino treatment.

This is a phenomenological high-z extrapolation, not a GWTC-4 measurement at
high redshift. There is no mass/spin-redshift evolution or detector selection
in this population generator.
