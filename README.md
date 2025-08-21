# PC5215_Numerical-Recipes
NUS Course PC5215 - Numerical Recipes course materials repository, course instructor: Prof. Wang Jian-Sheng.

**Course Website:** https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/cz5101.html

**INSTRUCTOR:** Prof. [Wang Jian-Sheng](https://www.physics.nus.edu.sg/~phywjs)

**TEXTBOOKS**: "Introducing Python", B. Lubanovic; ``[Numerical Recipes](http://www.numerical.recipes/) in C (1992), the Art of Scientific Computing,'' W. H. Press, etc.

**MODULE DESCRIPTION**: Covers computational techniques for the solution of problems arising in physics, with an emphasis on molecular simulation and modelling. Topics will be from the text, "Numerical Recipes", Press et al, supplemented with example problems in materials and condensed matter physics. The textbook is in C, but this course will use Python for programming. This course insures that graduate students intend to do research in computational physics will have sufficient background in computational methods and programming experience. Assessment: 4 labs assignments 40%, midterm 15%, final 45%.

**COURSE OUTLINE (lecture slides)**

week 0, 4-8 Aug, no class

[week 1 ](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec1.pptx), 11,14 Aug, introduction, Python

week 2, 18,21 Aug, Chap 1, Python continued, floating point representation, error, accuracy, etc

[week 3](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec2.pptx), 25,28 Aug, Chap 2, LU decomposition, [lab 1](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/lab1.pdf)

[week 4](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec3.pptx), 1,4 Sep, Chap 3, interpolation and extrapolation

[week 5](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec4.pptx), 8,11 Sep, Chap 4, integration, gaussian quadrature,

[week 6](http://www.physics.nus.edu.sg/~phywjs/BeijingWorkshop.html), 15,18 Sep, Chap 7, Monte Carlo integration

week 7, recess week, no class

[week 8](http://www.physics.nus.edu.sg/~phywjs/BeijingWorkshop.html), 29 Sep,2 Oct, Monte Carlo method, ([midterm test](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/midterm-2024-ans.pdf) at Thursday lecture time)

week 9, 6,9 Oct, Monte Carlo, continued

[week 10](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec7.pptx), 13,16 Oct, Chap 10, optimization, conjugate gradient

[week 11](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec9.pptx), 23 Oct, Chap 15, modeling of data (least square) (Monday is Deepavali)

[week 12](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec10.pptx), 27,30 Oct, Chap 16, ordinary differential equations, molecular dynamics

[week 13](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec12.pptx), 3,6 Nov, ODE, continued

[week 14](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec11.pptx), 10,13 Nov, Chap 17, boundary value problem, revision (last week)

\---

[Review Problems](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/prob-review.pdf) (with partial answers)

[2008](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final08-answers.pdf), [2010](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final10-answers.pdf), [2011](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final11-answers.pdf), [2012](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final12-answers.pdf), [2013](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final13-answers.pdf), [2014](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final14-answers.pdf), [2015](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final15-answers.pdf), [2016](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final16-answers.pdf), [2017](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final17-answers.pdf), [2018](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final18-answers.pdf), [2019](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final19-answers.pdf), [2020](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final20-answers.pdf), [2022](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final22-answers.pdf), [2023](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final23-answers.pdf), [2024](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/final24-answers.pdf). final exam and answers.

Some sample codes: [stability.py](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/stability.py), [ludcmp.py](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/ludcmp.py), [polint.py](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/polint.py), [trapzd.py](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/trapzd.py), [mnist.py](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/mnist.py), [conjugategradient.py](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/CG.py), [nfit.dat](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/nfit.dat), [testieee.c](https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/testieee.c).

Some interesting sites related to the course: [LAPACK](http://www.netlib.org/lapack/), an easy to read article on [conjugate gradient method](http://www-2.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf), fast Fourier transform [FFTW](http://www.fftw.org/), [symplectic integrator](https://arxiv.org/abs/1607.03882).

We recommend to install Anaconda ([www.anaconda.com](https://www.anaconda.com/download/success)) with Jupyter notebook for Python programming.
