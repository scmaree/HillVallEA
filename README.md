# HillVallEA (GECCO 2018)
S.C. Maree, T. Alderliesten, D. Thierens, P.A.N. Bosman


The Hill-Valley Evolutionary Algorithm (HillVallEA) is a real-valued evolutionary algorithm specifically aimed for multi-modal optimization. It has been described initially in 


>**Real-Valued Evolutionary Multi-Modal Optimization driven by Hill-Valley Clustering In Proceedings of the Genetic and Evolutionary Computation Conference**
> S.C. Maree, T. Alderliesten, D. Thierens, and P.A.N. Bosman. 
> GECCO-2018, ACM Press, New York, New York, 2018.

and benchmarked in the following technical reports

> **Benchmarking the Hill-Valley Evolutionary Algorithm for the GECCO 2018 Competition on Niching**
>S.C. Maree, T. Alderliesten, D. Thierens, P.A.N. Bosman
> arXiv preprint [arXiv:1807.00188](https://arxiv.org/abs/1807.00188), 2018.

and

> **Benchmarking HillVallEA for the GECCO 2019 Competition on Multimodal Optimization**
>S.C. Maree, T. Alderliesten, P.A.N. Bosman
> arXiv preprint [arXiv:1907.10988](https://arxiv.org/abs/1907.10988), 2019
.


---

# Getting started
Start by making a clone of the repository,

``` git clone https://github.com/SCMaree/HillVallEA```



Two example scripts have been provided to demonstrate the use of HillVallEA. Call `make` to build using your favorite compiler. This builds the two provided example scripts, `example_simple` and `example_cec2013_benchmark`. 


The script `example_simple` uses HillVallEA to solve the Six Hump Camel Back problem and `example_simple.cpp` can function as a guideline in order to implement your own problem. 


The script `example_cec2013_benchmark` runs HillVallEA on the problems of the [CEC2013 niching benchmark](https://github.com/mikeagn/CEC2013/)
reproduces the obtained peak ratio and static f1 averaged over a number of runs as stated in the above mentioned technical report.

To clean up after compilation, call `make clean`. 

