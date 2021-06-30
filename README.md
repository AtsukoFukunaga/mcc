# mcc
## Construct distance-based multivariate control charts

These codes produce distance-based multivariate control charts (dt and dbt) from a data frame of abundance data (e.g. fish count, benthic cover) and a data frame of factors.<br>
<br>
The rows of the two data frames should match, with each row representing a sample. The abundance data should only contain numerical values (counts, density, cover etc.) and the factor data frame should contain "year" and "site" columns.  If the data frames only contain data from a single site, the "site" column is still required.
