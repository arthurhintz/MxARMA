# output arma

    Code
      mxarma.fit(y1, ar = c(1, 2), ma = c(1, 2))$model
    Output
             Estimate Std. Error z value Pr(>|z|)
      alpha    0.3751     0.3475  1.0794   0.2804
      phi1     0.9823     0.7607  1.2912   0.1966
      phi2    -0.1457     0.6100  0.2389   0.8112
      theta1   0.3133     0.7610  0.4116   0.6806
      theta2  -0.0880     0.3764  0.2337   0.8153

---

    Code
      mxarma.fit(y2, ar = c(1), ma = c(1, 2))$model
    Output
             Estimate Std. Error z value Pr(>|z|)
      alpha    0.5287     0.0416 12.6977    0e+00
      phi1     0.5735     0.0352 16.3112    0e+00
      theta1   0.7368     0.0413 17.8214    0e+00
      theta2   0.1319     0.0392  3.3624    8e-04

---

    Code
      mxarma.fit(y3, ar = c(1, 2), ma = c(1))$model
    Output
             Estimate Std. Error z value Pr(>|z|)
      alpha    0.5483     0.0462 11.8755   0.0000
      phi1     0.6492     0.0479 13.5555   0.0000
      phi2     0.1456     0.0466  3.1251   0.0018
      theta1   0.6844     0.0375 18.2741   0.0000

---

    Code
      mxarma.fit(y4, ar = c(1), ma = c(1))$model
    Output
             Estimate Std. Error z value Pr(>|z|)
      alpha    0.5551     0.0300 18.5184        0
      phi1     0.5426     0.0242 22.4126        0
      theta1   0.7441     0.0195 38.1780        0

# output ar

    Code
      mxarma.fit(y1, ar = c(1, 2))$model
    Warning <simpleWarning>
      no non-missing arguments to max; returning -Inf
    Output
            Estimate Std. Error z value Pr(>|z|)
      alpha   0.5714     0.0423 13.5022        0
      phi1    0.5938     0.0263 22.5508        0
      phi2    0.1668     0.0263  6.3412        0

---

    Code
      mxarma.fit(y2, ar = c(1))$model
    Warning <simpleWarning>
      no non-missing arguments to max; returning -Inf
    Output
            Estimate Std. Error z value Pr(>|z|)
      alpha   0.5022     0.0244 20.5889        0
      phi1    0.5983     0.0210 28.5382        0

---

    Code
      mxarma.fit(y3, ar = c(1, 2, 3))$model
    Warning <simpleWarning>
      no non-missing arguments to max; returning -Inf
    Output
            Estimate Std. Error z value Pr(>|z|)
      alpha   0.5366     0.0719  7.4671    0e+00
      phi1    0.6324     0.0277 22.8382    0e+00
      phi2    0.1666     0.0325  5.1293    0e+00
      phi3    0.0981     0.0276  3.5617    4e-04

# output ma

    Code
      mxarma.fit(y1, ar = c(1, 2))$model
    Warning <simpleWarning>
      no non-missing arguments to max; returning -Inf
    Output
            Estimate Std. Error z value Pr(>|z|)
      alpha   0.2880     0.0152 18.9917        0
      phi1    0.6013     0.0260 23.1119        0
      phi2   -0.2079     0.0260  7.9868        0

---

    Code
      mxarma.fit(y2, ar = c(1))$model
    Warning <simpleWarning>
      no non-missing arguments to max; returning -Inf
    Output
            Estimate Std. Error z value Pr(>|z|)
      alpha   0.3261     0.0149 21.9367        0
      phi1    0.4147     0.0221 18.7320        0

---

    Code
      mxarma.fit(y3, ar = c(1, 2, 3))$model
    Warning <simpleWarning>
      no non-missing arguments to max; returning -Inf
    Output
            Estimate Std. Error z value Pr(>|z|)
      alpha   0.2770     0.0175 15.8701    0.000
      phi1    0.6288     0.0277 22.7143    0.000
      phi2   -0.1668     0.0321  5.1946    0.000
      phi3    0.0414     0.0277  1.4945    0.135

