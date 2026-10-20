LotSizingProblemAlgorithms

This repository contains the implementation of the efficient math-programming-based heuristics (matheuristics) developed for the Integrated Lot Sizing and Scheduling Problem on Parallel Machines with Sequence-Dependent and Non-Triangular Setup Costs and Times (ILSSP-NT).

The algorithms featured in this codebase include hybrid Relax-and-Fix (RF), Fix-and-Optimize (FO), Path-Relinking (PR), and Kernel Search (KS) approaches specifically tailored to handle shortcut items and non-triangular setup behaviors.

📄 Reference / Citation
If you use this code or build upon this research, please cite the following paper:

Carvalho, D. M., & Nascimento, M. C. V. (2022). Hybrid matheuristics to solve the integrated lot sizing and scheduling problem on parallel machines with sequence-dependent and non-triangular setup. European Journal of Operational Research, 296(2), 158–173. DOI: 10.1016/j.ejor.2021.03.050

BibTeX Entry:

Snippet de código
@article{CarvalhoNascimento2022,
  title = {Hybrid matheuristics to solve the integrated lot sizing and scheduling problem on parallel machines with sequence-dependent and non-triangular setup},
  author = {Carvalho, Desiree M. and Nascimento, Mari{\'a} C. V.},
  journal = {European Journal of Operational Research},
  volume = {296},
  number = {2},
  pages = {158--173},
  year = {2022},
  publisher = {Elsevier},
  doi = {10.1016/j.ejor.2021.03.050}
}
⚙️ Prerequisites & Dependencies
To compile and execute these algorithms, you will need:

A C++ compiler supporting modern standards (e.g., C++11 or higher).

IBM ILOG CPLEX Optimizer (Tested with version 12.10).

Note: A valid CPLEX license is required to run the optimization subproblems. Students and academic researchers can typically obtain a free license through the IBM Academic Initiative.

🛠️ Repository Structure
src/ — Contains the C++ source files for the matheuristic frameworks (RF, FO, PR, KS).

instances/ — Benchmark instance files adapted from the literature for single and parallel machine configurations.

LICENSE — MIT License (permits free use, modification, and distribution with attribution).
