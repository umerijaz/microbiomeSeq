# microbiomeSeq
An R package for microbial community analysis in an environmental context

This package is developed to enhance the available statistical analysis procedures in R by providing more 
analysis produre and visualisation of results for microbial communities data obtained from 16S rRNA. See full details of the usage and dependencies at
[microbiomeSeq](http://userweb.eng.gla.ac.uk/umer.ijaz/projects/microbiomeSeq_Tutorial.html) tutorial.

#### **Main Features**

* [Alpha and Beta diversity](https://github.com/umerijaz/microbiomeSeq/wiki/Alpha-and-Beta-Diversity)
* [Differential abundance analysis](https://github.com/umerijaz/microbiomeSeq/wiki/Differential-Abundance).
* [Effects of environmental variables on community structure](https://github.com/umerijaz/microbiomeSeq/wiki/Community-Environment-relationship).
* [Co-occurence pattern analysis in community data](https://github.com/umerijaz/microbiomeSeq/wiki/Co-ocurrence-Pattern-Analysis).

![alt text](https://github.com/umerijaz/microbiomeSeq/blob/master/disclaimer.png)

#### **Installation** 

Install the package with its associated dependencies and load it for usage in R.
```
library(devtools)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)
```

##### **Data format/requirement**
The data is required to be a phyloseq object (of `phyloseq-class`) comprising taxa abundance information, taxonomy assignment, 
sample data which is a combination of measured environmental variables together with any categorical variables present in 
the samples. If the phylogenetic tree is available, it can also be part but not so relevant for most of the functionality 
implemented here so far.
We choose to use this format since we can have enormous options for 
manipulating the data as we progress with the analysis and
visualisations. Details of format and comprehensive manipulations 
of phyloseq objects are available at [phyloseq page](https://github.com/joey711/phyloseq).
