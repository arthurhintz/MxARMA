# cdf_maxwell works correctly

    Code
      pmax(1, c(2, 5))
    Output
      [1] 0.111999714 0.008386613

---

    Code
      pmax(2, c(1))
    Output
      [1] 0.9829499

---

    Code
      pmax(5, c(1.5, 0.9, 0.8))
    Output
      [1] 0.9999968 1.0000000 1.0000000

---

    Code
      pmax(10, c(1, 0.912))
    Output
      [1] 1 1

---

    Code
      pmax(20, c(3.6789))
    Output
      [1] 1

---

    Code
      pmax(50, c(10, 20))
    Output
      [1] 1.0000000 0.9988199

