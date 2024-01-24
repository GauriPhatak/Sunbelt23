# Sunbelt23
REpository for Sunbelt presentation code

# Notes on the Simulation code:

# Observe all ergm terms for known community netowrks.
 Look at which terms are common to communities vs overall networks.
 using ERGM to generate networks with covariates.
 Short introductoin on ERGM models.

 ergm package recognizes missing tie information and codes it appropriately.

 can use ergm simulate function to simulate the model.

 can use estimates from existing community detection based networks as parameters to generate other networks.
 Compare the comm det methods on the generated networks.
# This will help validate if the ergm model used to fit the CD network captures the relationships within and between communities.

 Datasets to be used:

 Wastewater dataset

 DBLP dataset: Emailed the researcher who owns it

 Weibo?  dataset : https://www.kdd.org/kdd-cup/view/kdd-cup-2012-track-1

 Socio patterns datasets: http://www.sociopatterns.org/datasets/

 Protein sequencing data? https://string-db.org/cgi/input?sessionId=bhrcNwbbYsYP&input_page_show_search=on

 Biological datasets: http://dp.univr.it/~laudanna/LCTST/downloads/index.html

 Urbanity global network: https://figshare.com/articles/dataset/Global_Urban_Network_Dataset/22124219

 Largish network with ground truth and attributes: https://github.com/he-tiantian/Attributed-Graph-Data

 Papers with code: https://paperswithcode.com/datasets?mod=graphs

 Pajek datasets: http://vlado.fmf.uni-lj.si/pub/networks/data/default.htm also : http://vladowiki.fmf.uni-lj.si/doku.php?id=pajek:data:urls:index

# Look at clique based community detection methods.look at triads and higher relationship between nodes within and between communities


 Covariates assigned to each node. Type of ovariates: continous. Maybe from mixture of distributions.
 Multinomial distributions.
 How to relate covariates with the cluster assignments.
 Cov assignment based on cluster assignment. Differ the level of dependence between complete to no dependence.

https://stats.stackexchange.com/questions/397526/generate-value-of-variables-for-given-correlation-coefficient
https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables/15040#15040
https://stats.stackexchange.com/questions/13382/how-to-define-a-distribution-such-that-draws-from-it-correlate-with-a-draw-from/13384#13384

