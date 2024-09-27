# ppi-prediction-from-rna-seq

Construct a model to predict protein-protein interactions using RNA-seq data by identifying co-expressed gene groups and integrating with proteomics data.

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Abstract

There is an overwhelming amount of public RNA sequencing data available online for free. Since RNA-seq is generally much cheaper than mass spectrometry, if we are interested in predicting the behavior of proteins, we may wish to do so directly from RNA-seq data.

The purpose of this project is to construct a model and pipeline that identifies co-expressed groups of genes under particular conditions at the mRNA level, and integrates this with proteomics data to learn protein-protein interactions from mRNA co-expression and other relevant features.

Participants are encouraged to use any public sequencing data they find suitable, along with any quantitative or computational approaches. You may find it helpful to use [GTEx](https://gtexportal.org/home/) data, as well as data from [STRING](https://string-db.org/) and [BioGRID](https://thebiogrid.org/).

## Installation

Provide instructions on how to install and set up the project, such as installing dependencies and preparing the environment.

```bash
# Example command to install dependencies (Python)
pip install project-dependencies

# Example command to install dependencies (R)
install.packages("project-dependencies")
```

## Quick Start

Provide a basic usage example or minimal code snippet that demonstrates how to use the project.

```python
# Example usage (Python)
import my_project

demo = my_project.example_function()
print(demo)
```
```r
# Example usage (R)
library(my_project)

demo <- example_function()
print(demo)
```

## Usage

Add detailed information and examples on how to use the project, covering its major features and functions.

```python
# More usage examples (Python)
import my_project

demo = my_project.advanced_function(parameter1='value1')
print(demo)
```
```r
# More usage examples (R)
library(demoProject)

demo <- advanced_function(parameter1 = "value1")
print(demo)
```

## Contribute

Contributions are welcome! If you'd like to contribute, please open an issue or submit a pull request. See the [contribution guidelines](CONTRIBUTING.md) for more information.

## Support

If you have any issues or need help, please open an [issue](https://github.com/hackbio-ca/ppi-prediction-from-rna-seq/issues) or contact the project maintainers.

## License

This project is licensed under the [MIT License](LICENSE).
