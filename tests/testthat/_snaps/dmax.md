# pdf_maxwell works correctly

    Code
      dmax(1, c(2, 5))
    Output
      [1] 0.29479494 0.02465028

---

    Code
      dmax(2, 1)
    Output
      [1] 0.07962814

---

    Code
      dmax(5, c(1.5, 0.9, 0.8))
    Output
      [1] 1.723877e-05 9.536039e-16 3.976352e-20

---

    Code
      dmax(10, c(1, 0.912))
    Output
      [1] 1.639681e-53 1.408517e-64

---

    Code
      dmax(20, c(3.6789))
    Output
      [1] 1.183744e-15

---

    Code
      dmax(50, c(10, 20))
    Output
      [1] 1.215535e-13 3.545640e-04

