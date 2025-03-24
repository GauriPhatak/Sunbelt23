# Community detection and Vector auto regression models research.

Network systems are ubiquitous across various domains, and analyzing these systems is crucial
for understanding complex relationships among multiple variables. Social networks, formed through
relationships between individuals, offer valuable insights into structural patterns and key nodes within
a system. Social network analysis (SNA) enables the identification of information flow, disease propa-
gation, and other dynamic processes, thereby supporting data-driven decision-making based on network
structures. Community detection within social networks identifies latent groupings that can inform
targeted interventions. Incorporating community structures into public policy design enhances the ef-
fectiveness of group-specific actions. As a well-established area of study in social network analysis,
community detection often incorporates node attributes in conjunction with network structure to im-
prove the identification of communities. However, missing data — whether pertaining to relationships
between nodes or node attributes — frequently poses a challenge in real-world networks. This research
investigates the use of Affiliation Graph Models (AGM) for addressing the issue of missing data in com-
munity detection. Specifically, we explore methods for simultaneously imputing missing node attribute
values and detecting community structures in networks where some node attributes are absent. Our
approach aims to predict missing information by leveraging the community affiliations of the nodes,
enhancing both the accuracy of community detection and the completeness of network data.

The motivating application of this project is to investigate the mechanisms of pathogen spread
across the state of Oregon. To achieve this, we utilize spatio-temporal data derived from COVID-
19 measurements in wastewater samples collected at various locations across the state. Community
detection is employed as a method to group nodes based on their similarity, while change-point detection
and Granger causality networks — developed using penalty-based vector auto regression (VAR) models
— are utilized as approaches to find important locations and assess if these locations change over time.
In the first part of our research, we extend community detection methods based on the AGM to identify
clusters of similar locations. This analysis facilitates a better understanding of the health status of one
location in relation to another, enabling more effective allocation of resources. In the second part, we
identify locations with the greatest influence on surrounding areas and the broader state. This allows us
to pinpoint regions where pathogens are introduced early, providing insights into the dynamics of disease
spread and opportunities for targeted interventions.


## Dataset repositories:
 DBLP dataset: Emailed the researcher who owns it

 Weibo?  dataset : https://www.kdd.org/kdd-cup/view/kdd-cup-2012-track-1

 Socio patterns datasets: http://www.sociopatterns.org/datasets/

 Protein sequencing data? https://string-db.org/cgi/input?sessionId=bhrcNwbbYsYP&input_page_show_search=on

 Biological datasets: http://dp.univr.it/~laudanna/LCTST/downloads/index.html

 Urbanity global network: https://figshare.com/articles/dataset/Global_Urban_Network_Dataset/22124219

 Largish network with ground truth and attributes: https://github.com/he-tiantian/Attributed-Graph-Data

 Papers with code: https://paperswithcode.com/datasets?mod=graphs

 Pajek datasets: http://vlado.fmf.uni-lj.si/pub/networks/data/default.htm also : http://vladowiki.fmf.uni-lj.si/doku.php?id=pajek:data:urls:index
