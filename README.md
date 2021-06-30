# MCC
## Construct distance-based multivariate control charts

<br>
These codes produce distance-based multivariate control charts (<i>d<sub>t</sub></i> and <i>d<sup>b</sup><sub>t</sub></i>) from a data frame of abundance data (e.g., fish count, benthic cover) and a data frame of factors.<br>
<br>
The rows of the two data frames should match, with each row representing a sample. The abundance data should only contain numerical values (counts, density, cover etc.).<br>
For standard MCCs, factor data frame should contain "year" and "site" columns, with “year” representing survey year (or any appropriate survey frequencies) and “site” representing individual survey sites.<br>
For regional MCCs, an additional column “region” is required.<br>
If the data frames only contain data from a single site or region, the "site/region" column is still required.<br>
<br>
Three examples are included along with example datasets.
<li>mcc_example_1.R</li>
Standard MCCs for annual fish survey program with 15 survey sites.  MCCs are constructed for each site with respective confidence bounds.
<li>mcc_example_2.R</li>
Standard MCCs for annual fish survey program with 15 survey sites separated into three regions.  MCCs are constructed for each region with respective confidence bounds.
<li>Mcc_example_3.R</li>
MCCs after fitting covariates, with surveys from randomized sites.  MCCs are separately constructed for three regions, but shared confidence bounds for the three regions are used for plotting.<br>
<br>
References<br>
Standard MCCs:<br>
Anderson MJ and Thompson AA. 2004. Multivariate control charts for ecological and environmental monitoring. Ecological Applications 14:1921–1935; DOI 10.1890/03-5379<br>
MCCs after fitting covariates:<br>
Fukunaga A and Kosaki RK. 2017. Use of multivariate control charts to assess the status of reef fish assemblages in the Northwestern Hawaiian Islands. PeerJ 5:e3651; DOI 10.7717/peerj.3651

