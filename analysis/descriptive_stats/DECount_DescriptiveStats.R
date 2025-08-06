library(tidyverse)
library(stringr)
astro <- read.csv("tableS22.csv")
micro <- read.csv("tableS18.csv")
oligo<- read.csv("tables26.csv")

#genes tested
astro %>%
   pull(gene)%>%unique()%>%length()

micro %>%
   pull(gene)%>%unique()%>%length()

oligo %>%
   pull(gene)%>%unique()%>%length()

#DE in any condition
astro %>%
   filter(lfsr<0.05)%>%
   pull(gene)%>%unique()%>%length()

micro %>%
   filter(lfsr<0.05)%>%
   pull(gene)%>%unique()%>%length()

oligo %>%
   filter(lfsr<0.05)%>%
   pull(gene)%>%unique()%>%length()
#DE in >2 condition
astro %>%
   group_by(gene)%>%
   filter(lfsr<0.05)%>%
   summarize(n_condition=n())%>%
   filter(n_condition>=2)%>%
   pull(gene)%>%length()

micro %>%
   group_by(gene)%>%
   filter(lfsr<0.05)%>%
   summarize(n_condition=n())%>%
   filter(n_condition>=2)%>%
   pull(gene)%>%length()

oligo %>%
   group_by(gene)%>%
   filter(lfsr<0.05)%>%
   summarize(n_condition=n())%>%
   filter(n_condition>=2)%>%
   pull(gene)%>%length()

#DE in 2 or more studies
astro %>%
   mutate(study = str_extract(comparison, "^[^_]+")) %>%
   filter(lfsr < 0.05) %>%
   group_by(gene) %>%
   summarize(n_study = n_distinct(study)) %>%  
   filter(n_study > 2)%>%
   pull(gene)%>%length()  

micro %>%
   mutate(study = str_extract(comparison, "^[^_]+")) %>%
   filter(lfsr < 0.05) %>%
   group_by(gene) %>%
   summarize(n_study = n_distinct(study)) %>%  
   filter(n_study > 2)%>%
   pull(gene)%>%length()  

oligo %>%
   mutate(study = str_extract(comparison, "^[^_]+")) %>%
   filter(lfsr < 0.05) %>%
   group_by(gene) %>%
   summarize(n_study = n_distinct(study)) %>%  
   filter(n_study > 2)%>%
   pull(gene)%>%length()  
