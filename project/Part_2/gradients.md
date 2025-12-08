# -----------------------------------------------------

# 1. Objective Function

# -----------------------------------------------------

The waveform is parameterized by **phases**:

[
s_k = \frac{1}{\sqrt{N}} e^{j x_k}, \qquad k=1,\dots,N \equiv nPulse.
]

Define the aperiodic autocorrelation using FFT-IFFT:

[
c = \mathrm{IFFT}\big( S \odot S^{!*} \big),
]

where (S = \mathrm{FFT}(s, nfft)).
The objective is:

[
J(x) = \log!\left(\sum_{n} | w_n, c_n |^{p} \right)
]

Define the interior quantity:

[
J_{\text{SL}} = \sum_n |w_n c_n|^p.
]

Then

[
J = \log(J_{\text{SL}}).
]

So the gradient is:

[
\nabla_x J = \frac{1}{J_{\text{SL}}} \nabla_x J_{\text{SL}}.
]

The whole problem reduces to differentiating (J_{\text{SL}}).

---

# -----------------------------------------------------

# 2. Useful rewriting of the ACF

# -----------------------------------------------------

The autocorrelation is:

[
c = \mathrm{IFFT}( S S^{!*} ).
]

Using convolution theorem:

[
c = s \star s^{!*},
]

i.e., the autocorrelation is a **circular convolution**:

[
c_n = \sum_{m} s_m , s^{!*}_{(m-n)}.
]

This is important because the derivative of a convolution is another convolution.

---

# -----------------------------------------------------

# 3. Derivative of the objective wrt autocorrelation values

# -----------------------------------------------------

Consider:

[
J_{SL} = \sum_n |w_n c_n|^p = \sum_n w_n^p |c_n|^p.
]

The Wirtinger derivative w.r.t. (c_n^*) is:

[
\frac{\partial |c_n|^p}{\partial c_n^*}
= \frac{p}{2} |c_n|^{p-2} c_n.
]

Thus

[
\frac{\partial J_{SL}}{\partial c_n^*}
= \frac{p}{2} w_n^p |c_n|^{p-2} c_n.
]

Define the intermediate variable used in MATLAB:

[
t_n \equiv |c_n|^{p-2} c_n,
]
[
w_n^p = w_n \quad \text{(in your mask design, } w_n = w_n^p\text{)}.
]

In MATLAB:

```matlab
temp = (abs(ccorr).^(p-2)).*ccorr;
tempsl = fft(wsl.*temp);   % = FFT(w_n t_n)
```

---

# -----------------------------------------------------

# 4. Chain rule: derivative of c wrt s

# -----------------------------------------------------

The circular autocorrelation is:

[
c = s \star s^{!*}.
]

The derivative of a correlation w.r.t. the signal is:

[
\frac{\partial c}{\partial s^*} = s \star \delta
]
[
\frac{\partial c^*}{\partial s^*} = s \star \delta
]

The differential of (J_{SL}) is:

[
dJ_{SL}
= \sum_n \frac{\partial J_{SL}}{\partial c_n^*} , dc_n^*.
]

Using convolution identity and FFT/IFFT duality:

[
dc^* = \mathrm{IFFT}( S \odot dS^* ).
]

After substitution:

[
\nabla_{s^*} J_{SL}
= p , s \star (w \odot t).
]

In FFT form:

[
\nabla_{s^*} J_{SL}
= p , \mathrm{IFFT}\big( S \cdot \mathrm{FFT}( w \odot t ) \big).
]

This matches MATLAB:

```matlab
c = ifft(fft(s,nfft) .* tempsl);   % c = IFFT(S * FFT(w t))
g_s = p*c(1:N);                    % gradient wrt complex s
```

BUT MATLAB uses:

```matlab
g = p*imag(conj(s).*c(1:N));
```

Why?

Because the derivative must be w.r.t. **real phase parameters** (x_k), not the complex amplitudes (s_k).

---

# -----------------------------------------------------

# 5. Chain rule: gradient wrt phases x

# -----------------------------------------------------

Each signal sample is:

[
s_k = \frac{1}{\sqrt{N}} e^{j x_k}.
]

Derivative:

[
\frac{\partial s_k}{\partial x_k}
= j s_k.
]
[
\frac{\partial s_k^*}{\partial x_k}
= -j s_k^*.
]

The total derivative:

[
\frac{\partial J_{SL}}{\partial x_k}
= 2\operatorname{Im}\left( s_k^* \frac{\partial J_{SL}}{\partial s_k^*} \right).
]

Thus:

[
\nabla_x J_{SL}
= 2,\Im!\left( s^* \odot \nabla_{s^*} J_{SL} \right).
]

Substitute the previously computed gradient:

[
\nabla_x J_{SL}
= 2,\Im!\left(
s^* \odot
\big[ p, c(1:N) \big]
\right).
]

This exactly matches MATLAB:

```matlab
g = p*imag(conj(s).*c(1:N));
VJ = 2*g/Jsl;
```

Finally:

[
\nabla_x J
= \frac{1}{J_{SL}} \nabla_x J_{SL}
= \frac{2p}{J_{SL}}
\Im!\left( s^* \odot c(1:N) \right).
]

---

# -----------------------------------------------------

# ✔ Final Gradient (closed-form)

# -----------------------------------------------------

[
\boxed{
\nabla_x J
= \frac{2p}{\sum_n |w_n c_n|^p },
\Im!\left(
s^* \odot \mathrm{IFFT}!\left(
S \odot \mathrm{FFT}( w \odot |c|^{p-2} c )
\right)_{1:N}
\right)
}
]

This is **exactly** your MATLAB code:

```matlab
temp   = (abs(ccorr).^(p-2)).*ccorr;
tempsl = fft(wsl.*temp);
c      = ifft(fft(s,nfft) .* tempsl);
g      = p*imag(conj(s).*c(1:nPulse));
VJ     = 2*g/Jsl;
```
