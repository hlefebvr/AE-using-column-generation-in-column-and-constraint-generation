# Using Column Generation in Column-and-Constraint Generation for Adjustable Robust Optimization

![GitHub](https://img.shields.io/github/license/hlefebvr/AE-using-column-generation-in-column-and-constraint-generation)
![GitHub issues](https://img.shields.io/github/issues-raw/hlefebvr/AE-using-column-generation-in-column-and-constraint-generation)
![Repo status](https://www.repostatus.org/badges/latest/active.svg)

This repository contains the code and data used for the paper 
[Lefebvre, H., Schmidt, M., and Thürauf, J. (2023) "Using Column Generation in Column-and-Constraint Generation for Adjustable Robust Optimization"](https://optimization-online.org/XXXXXXXXX).

## How to cite

If you are using this repository for your research, please cite our preprint.

```latex
@misc{Lefebvre2023,
  title={Using Column Generation in Column-and-Constraint Generation for Adjustable Robust Optimization},
  author={Henri Lefebvre and Martin Schmidt and Johannes Thürauf},
  year={2023},
  url = {https://optimization-online.org/XXXXXXXXX}
}
```

## How to use

### Dependencies

Compiling requires the following dependencies*:

- [CMake](https://cmake.org/) version 3.22.1.
- [GCC](https://gcc.gnu.org/) version 11.4.0 (for `__has_include` preprocessor).
- [Gurobi](https://www.gurobi.com/) version 10.0.1.
- [idol](https://github.com/hlefebvr/idol) version 0.3.5-alpha. **Idol will be downloaded automatically, you do not need to install it.**

Running requires the following dependencies*:

- [Gurobi](https://www.gurobi.com/) version 10.0.1.
- [idol](https://github.com/hlefebvr/idol) version 0.3.5-alpha. **Idol will be downloaded automatically, you do not need to install it.**

_* reported versions are not minimal required versions but are ones which have been tested._

### Compiling

#### Finding Gurobi

*(Recommended)* First, be sure to have your `GUROBI_HOME` environment variable configured according to [the official Gurobi guidelines](https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer-).

*(Alternatively)* You may use the `GUROBI_DIR` CMake variable, please refer to the [idol documentation page](https://hlefebvr.github.io/idol/installation/options.html).

#### Create Makefile with CMake

By running the following command, CMake will generate an appropriate Makefile. Additionally, it will download the idol
library. Thus, an internet connection is required for this step.

```shell
mkdir build
cd build
cmake ..
```

#### Compile with make

From the build directory, run the following.

```shell
make
```

### Running

Running our code is done as follows:

```shell
./build/FLP/FLP_solve <PATH_TO_INSTANCE> <STD_PHASE_TIME_LIMIT> <GAMMA> <USE_HEURISTIC>
./build/JSP/JSP_solve <PATH_TO_INSTANCE> <STD_PHASE_TIME_LIMIT> <GAMMA> <USE_HEURISTIC>
```

Arguments are as follows:

- `PATH_TO_INSTANCE`: the absolute path to an instance file.
- `STD_PHASE_TIME_LIMIT`: the time limit for the first phase of the algorithm (i.e., the time limit after which we switch to the column generation approach).
- `GAMMA`: the uncertainty budget $\Gamma$.
- `USE_HEURISRIC`: `true` if the "integer master" heuristic should be used, `false` otherwise.
