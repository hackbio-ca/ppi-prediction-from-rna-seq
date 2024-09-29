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

## Contribution Instructions for Hackathon Group Members

1. Clone the repository. Either use SSH or HTTPS:

```bash
git clone git@github.com:hackbio-ca/ppi-prediction-from-rna-seq.git
git clone https://github.com/hackbio-ca/ppi-prediction-from-rna-seq.git
```

2. Making a branch and implementing changes:

```bash
git checkout -b your-name
```
Use `git status` to check you're on the correct branch. You do not want to implement any features on the `main` branch.
```
git add filename
git commit -m "Changes to filename"
git push origin your-name
```

3. Make a pull request. Navigate to your branch on the github repository page, and click **Compare & pull request**. Make sure you compare `main` and the `your-name` branch.

4. The pull request will be reviewed and hopefully accepted. At this point, github will tell you that the `your-name` branch can be safely deleted.

5. Navigate back to your local `main` branch and update, then cleanup.

```bash
git checkout main
git pull origin main
git branch -d your-name
```

6. Upon implementing a new feature, simply repeat from step 2.

## Sources

| **Publication** | Cell line |
| :--- | :---: |
| Johnson *et al.*, (2021) [^1] | HEK293T |
| Johnson *et al.*, (2021) [^1]| Jurkat |
| Johnson *et al.*, (2021) [^1] | HUVEC |
| Huttlin *et al.*, (2021) [^2] | HEK293T |
| Huttlin *et al.*, (2021) [^2] | HCT116 |
| Göös *et al.*, (2022) [^3] | HEK293 |
| Khoroshkin *et al.*, (2024) [^4] | K562 |
| Banks *et al.*, (2014) [^5] | HEK293 |

[^1]: [Johnson, K. L. *et al.* Revealing protein-protein interactions at the transcriptome scale by sequencing. *Molecular Cell* **81,** 3877 (2021).](https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00574-8)
[^2]: [Huttlin, E. L. *et al.* Dual proteome-scale networks reveal cell-specific remodeling of the human interactome. *Cell* **184**, 11 (2021).](https://doi.org/10.1016/j.cell.2021.04.011)
[^3]: [Göös, H. *et al.* Human transcription factor protein interaction networks. *Nature Communications* **13**, (2022).](https://www.nature.com/articles/s41467-022-28341-5)
[^4]: [Khoroshkin, M. *et al.* Systematic identification of post-transcriptional regulatory modules. *Nature Communications* **15**, (2024).](https://www.nature.com/articles/s41467-024-52215-7)
[^5]: [Banks, C. *et al.* Controlling for Gene Expression Changes in Transcription Factor Protein Networks. *Molecular & Cellular Proteomics* **13**, (2014).](https://www.mcponline.org/article/S1535-9476(20)33081-4/fulltext)
