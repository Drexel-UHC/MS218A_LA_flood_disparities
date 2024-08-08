
{
  library(tidyverse)
  library(dfoptim)
  library(optimx)
  library(epitools)
  library(lme4)
  library(lmerTest)
  library(kableExtra)
  library(Hmisc)
  library(broom.mixed)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggthemes)
  library(rgeos)
  library(RColorBrewer)
  library(lubridate)
  library(scales)
  library(rbenchmark)
  remotes::install_github("bbolker/broom.mixed")
}

###Codebook-type notes
#L1 = level 1 = city
#L3 = level 3 = neighborhood


###read in analysis data

#all flood events (main analysis)
all <- read_csv("../Data/Analysis data/analysis_ge1_event_2024_03_13.csv")

#repeat/multiple flood events (sensitivity analysis)
all_rpt <- read_csv("../Data/Analysis data/analysis_ge2_event_2024_03_13.csv")

#Check N
all %>%
  ungroup() %>%
  summarise(SALID1 = n_distinct(SALID1),
            SALID2 = n_distinct(SALID2))

all_rpt %>%
  summarise(SALID1 = n_distinct(SALID1),
            SALID2 = n_distinct(SALID2))

#Check nhood (L3) area size for main text
all %>%
  summarise(area = median(BECADAREAL3),
            popL3 = median(CNSPOPL3),
            N = n())

sqrt(0.36)


###Table 1 summary stats
#L3 level vars

{
  ######Country specific L3-level vars
  #L3s
  (ta1 <-all %>%
     group_by(country) %>%
     distinct(id) %>%
     tally() %>%
     dplyr::rename(L3s = n)
  )
  
  #L3 pop density 
  (tabpd3 <- all %>%
      mutate(contvar = BECPOPDENSL3/1000) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(pdl3 = cont_all) %>%
      dplyr::select(country, L3s, pdl3)
  )
  
  (ta2 <- full_join(ta1, tabpd3))
  
  #L3 education
  (tabedu3 <- all %>%
      mutate(contvar = CNSMINPR_L3) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(edul3 = cont_all) %>%
      dplyr::select(country, L3s, edul3)
  )
  
  (ta3 <- full_join(ta2, tabedu3))
  
  #L3 int. density
  (tabint3 <- all %>%
      mutate(contvar = BECADINTDENSL3) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 0),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 0),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(intl3 = cont_all) %>%
      dplyr::select(country, L3s, intl3)
  )
  
  (ta4 <- full_join(ta3, tabint3))
  
  #L3 green
  (tabgreen3 <- all %>%
      mutate(contvar = BECMEDNDVINW2017L3) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 2),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 2),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(greenl3 = cont_all) %>%
      dplyr::select(country, L3s, greenl3)
  )
  
  (ta5 <- full_join(ta4, tabgreen3))
  
  #L3 distance
  (tabdist3 <- all %>%
      mutate(contvar = BECNBHDCENTL3) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(distl3 = cont_all) %>%
      dplyr::select(country, L3s, distl3)
  )

  (ta6 <- full_join(ta5, tabdist3))
  
  #L3 coastal
  (tabcst3 <- all %>%
      mutate(contvar = BECCOASTL3) %>%
      group_by(country) %>%
      mutate(l3_n_flood = sum(BECCOASTL3),
             l3_n_total = n()) %>%
      filter(row_number() == 1) %>%
      summarise(cont_md = l3_n_flood,
                cont_iqr = format(round(100*(l3_n_flood/l3_n_total), 1), nsmall = 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, " %)", sep = "")) %>%
      dplyr::rename(cstl3 = cont_all) %>%
      dplyr::select(country, cstl3)
  )
  
  (ta7 <- full_join(ta6, tabcst3))
  
  #L3 altitude
  (tabalt3 <- all %>%
      mutate(contvar = BECELEVATIONMEDIANL3) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 0),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 0),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(altl3 = cont_all) %>%
      dplyr::select(country, L3s, altl3)
  )
  
  (ta8 <- full_join(ta7, tabalt3))
  
  #L3 slope
  (tabslp3 <- all %>%
      mutate(contvar = BECSLOPEMEDIANL3) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(slpl3 = cont_all) %>%
      dplyr::select(country, L3s, slpl3)
  )
  
  (ta9_L3 <- full_join(ta8, tabslp3))
  
  #####L3 vars only, but total
  
  #L3s
  (ta1 <-all %>%
      distinct(id) %>%
      tally() %>%
      dplyr::rename(L3s = n)
  )
  
  #L3 pop density 
  (tabpd3 <- all %>%
      mutate(contvar = BECPOPDENSL3/1000) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(pdl3 = cont_all) %>%
      dplyr::select(country, L3s, pdl3)
  )
  
  (ta2 <- full_join(ta1, tabpd3))
  
  #L3 education
  (tabedu3 <- all %>%
      mutate(contvar = CNSMINPR_L3) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(edul3 = cont_all) %>%
      dplyr::select(country, L3s, edul3)
  )
  
  (ta3 <- full_join(ta2, tabedu3))
  
  #L3 int. density
  (tabint3 <- all %>%
      mutate(contvar = BECADINTDENSL3) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 0),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 0),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(intl3 = cont_all) %>%
      dplyr::select(country, L3s, intl3)
  )
  
  (ta4 <- full_join(ta3, tabint3))
  
  #L3 green
  (tabgreen3 <- all %>%
      mutate(contvar = BECMEDNDVINW2017L3) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 2),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 2),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(greenl3 = cont_all) %>%
      dplyr::select(country, L3s, greenl3)
  )
  
  (ta5 <- full_join(ta4, tabgreen3))
  
  #L3 distance
  (tabdist3 <- all %>%
      mutate(contvar = BECNBHDCENTL3) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(distl3 = cont_all) %>%
      dplyr::select(country, L3s, distl3)
  )
  
  (ta6 <- full_join(ta5, tabdist3))
  
  #L3 coastal
  (tabcst3 <- all %>%
      mutate(contvar = BECCOASTL3) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      mutate(l3_n_flood = sum(BECCOASTL3),
             l3_n_total = n()) %>%
      filter(row_number() == 1) %>%
      summarise(cont_md = l3_n_flood,
                cont_iqr = format(round(100*(l3_n_flood/l3_n_total), 1), nsmall = 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, " %)", sep = "")) %>%
      dplyr::rename(cstl3 = cont_all) %>%
      dplyr::select(country, cstl3)
  )
  
  (ta7 <- full_join(ta6, tabcst3))
  
  #L3 altitude
  (tabalt3 <- all %>%
      mutate(contvar = BECELEVATIONMEDIANL3) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 0),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 0),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(altl3 = cont_all) %>%
      dplyr::select(country, L3s, altl3)
  )
  
  (ta8 <- full_join(ta7, tabalt3))
  
  #L3 slope
  (tabslp3 <- all %>%
      mutate(contvar = BECSLOPEMEDIANL3) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L3s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(slpl3 = cont_all) %>%
      dplyr::select(country, L3s, slpl3)
  )
  
  (ta8_L3_tot <- full_join(ta8, tabslp3) %>%
      select(country, everything()))
  
  (ta9 <- t(full_join(ta8_L3_tot, ta9_L3)))
  kable(ta9) %>%
  kable_styling() %>%
    cat(., file = "tab1_neighborhood_2024-03-13.html")
}


##City level vars (L1) for Table 1

{
  ##By country
  #N of L1s
  (ta1 <-all %>%
     group_by(country) %>%
     distinct(SALID1) %>%
     tally() %>%
     dplyr::rename(L1s = n)
  )
  
  #L1 pop 
  (tabp1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECTPOP2010L1UX/1000) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(popl1 = cont_all) %>%
      dplyr::select(country, L1s, popl1)
  )
  
  (ta2 <- full_join(ta1, tabp1))
  
  #L1 GDP 
  (tabgdp1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = SECGDPGPPC/1000) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(gdpl1 = cont_all) %>%
      dplyr::select(country, L1s, gdpl1)
  )
  
  (ta3 <- full_join(ta2, tabgdp1))
  
  
  #L1 pop dens
  (tabpd1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECPOPDENS2010L1UX/1000) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(pdl1 = cont_all) %>%
      dplyr::select(country, L1s, pdl1)
  )
  
  (ta5 <- full_join(ta3, tabpd1))
  
  #L1 education
  (tabedu1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = CNSMINPR_L1AD) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(edu1 = cont_all) %>%
      dplyr::select(country, L1s, edu1)
  )
  
  (ta6 <- full_join(ta5, tabedu1))
  
  #L1 int density
  (tabint1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECADINTDENSL1UX) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(int1 = cont_all) %>%
      dplyr::select(country, L1s, int1)
  )
  
  (ta7 <- full_join(ta6, tabint1))
  
  #L1 green
  (tabgreen1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECMEDNDVINW2010L1UX) %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 2),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 2),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(green1 = cont_all) %>%
      dplyr::select(country, L1s, green1)
  )
  
  (ta8 <- full_join(ta7, tabgreen1))


##Total
  
  #N of L1s
  (ta1 <-all %>%
     mutate(country = "Total") %>%
     group_by(country) %>%
     distinct(SALID1) %>%
     tally() %>%
     dplyr::rename(L1s = n)
  )
  
  #L1 pop 
  (tabp1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECTPOP2010L1UX/1000) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 0),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 0),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(popl1 = cont_all) %>%
      dplyr::select(country, L1s, popl1)
  )
  
  (ta2 <- full_join(ta1, tabp1))
  
  #L1 GDP 
  (tabgdp1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = SECGDPGPPC/1000) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(gdpl1 = cont_all) %>%
      dplyr::select(country, L1s, gdpl1)
  )
  
  (ta3 <- full_join(ta2, tabgdp1))
  
    #L1 pop dens
  (tabpd1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECPOPDENS2010L1UX/1000) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(pdl1 = cont_all) %>%
      dplyr::select(country, L1s, pdl1)
  )
  
  (ta5 <- full_join(ta3, tabpd1))
  
  #L1 education
  (tabedu1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = CNSMINPR_L1AD) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(edu1 = cont_all) %>%
      dplyr::select(country, L1s, edu1)
  )
  
  (ta6 <- full_join(ta5, tabedu1))
  
  #L1 int density
  (tabint1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECADINTDENSL1UX) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 1),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 1),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(int1 = cont_all) %>%
      dplyr::select(country, L1s, int1)
  )
  
  (ta7 <- full_join(ta6, tabint1))
  
  #L1 green
  (tabgreen1 <- all %>%
      group_by(SALID1) %>%
      filter(row_number()==1) %>%
      mutate(contvar = BECMEDNDVINW2010L1UX) %>%
      mutate(country = "Total") %>%
      group_by(country) %>%
      summarise(cont_md = round(mean(contvar, na.rm = TRUE), 2),
                cont_iqr = round(sd(contvar, na.rm = TRUE), 2),
                L1s = n()) %>%
      mutate(cont_all = paste(cont_md, " (", cont_iqr, ")", sep = "")) %>%
      dplyr::rename(green1 = cont_all) %>%
      dplyr::select(country, L1s, green1)
  )
  
  (ta8_total <- full_join(ta7, tabgreen1))
  (ta9 <- t(full_join(ta8_total, ta8)))
  kable(ta9) %>%
    kable_styling() %>%
    cat(., file = "tab1_city_2024-03-13.html")
}


#####Table 2: Exposure stats
{
(t2_tot <- all %>%
  ungroup() %>%
  dplyr::mutate(tot_L3s = n(),
                tot_pop_L3s = sum(CNSPOPL3)) %>%
  dplyr::group_by(any_floods_y1n0, tot_L3s, tot_pop_L3s) %>%
  dplyr::summarise(
    L3_N = n(),
    pop_flood_cat = sum(CNSPOPL3)
  ) %>%
    mutate(country = "Overall") %>%
  dplyr::mutate(flood_L3_perc = round(100*(L3_N/tot_L3s), 1),
                flood_pop_perc = round(100*(pop_flood_cat/tot_pop_L3s), 1)) %>%
    mutate(tot_pop_L3s = round(tot_pop_L3s/1000000, 1),
           pop_flood_cat = round(pop_flood_cat/1000000, 1)) %>%
    filter(any_floods_y1n0 == 1) %>%
    ungroup() %>%
    select(country, tot_L3s, tot_pop_L3s, pop_flood_cat, flood_pop_perc)) 
  
  kable(t2_tot) %>%
    kable_styling() %>%
    cat(., file = "tab2_overall_2023-11-13.html")
  

(t2_educ <- all %>%
  group_by(educ_quint) %>%
  dplyr::mutate(tot_L3s = n(),
                tot_pop_L3s = sum(CNSPOPL3)) %>%
  dplyr::group_by(educ_quint, any_floods_y1n0, tot_L3s, tot_pop_L3s) %>%
  dplyr::summarise(
    L3_N = n(),
    pop_flood_cat = sum(CNSPOPL3)
  ) %>%
  dplyr::mutate(flood_L3_perc = round(100*(L3_N/tot_L3s), 1),
                flood_pop_perc = round(100*(pop_flood_cat/tot_pop_L3s), 1)) %>%
  filter(any_floods_y1n0 == 1) %>%
    ungroup() %>%
    mutate(tot_pop_L3s = round(tot_pop_L3s/1000000, 1),
           pop_flood_cat = round(pop_flood_cat/1000000, 1)) %>%
    select(-any_floods_y1n0, -L3_N, -flood_L3_perc))
  
  kable(t2_educ) %>%
    kable_styling() %>%
    cat(., file = "tab2_educ_2023-11-13.html")

  (t2_country <- all %>%
      group_by(country) %>%
      dplyr::mutate(tot_L3s = n(),
                    tot_pop_L3s = sum(CNSPOPL3)) %>%
      dplyr::group_by(country, any_floods_y1n0, tot_L3s, tot_pop_L3s) %>%
      dplyr::summarise(
        L3_N = n(),
        pop_flood_cat = sum(CNSPOPL3)
      ) %>%
      dplyr::mutate(flood_L3_perc = round(100*(L3_N/tot_L3s), 1),
                    flood_pop_perc = round(100*(pop_flood_cat/tot_pop_L3s), 1)) %>%
      filter(any_floods_y1n0 == 1) %>%
      ungroup() %>%
      mutate(tot_pop_L3s = round(tot_pop_L3s/1000000, 1),
             pop_flood_cat = round(pop_flood_cat/1000000, 1)) %>%
      select(-any_floods_y1n0, -L3_N, -flood_L3_perc))

  kable(t2_country) %>%
    kable_styling() %>%
    cat(., file = "tab2_country_2024-03-13.html")  

}


##Figure 2: Bar plot of education quintiles

{

###Pop weighting, instead of count of neighborhoods
all_plot <- all %>%
  filter(!is.na(educ_quint)) %>%
  group_by(educ_quint) %>%
  mutate(quint_totpop = sum(CNSPOPL3)) %>%
  group_by(educ_quint, flood_n_cat, quint_totpop) %>%
  summarise(quint_pop_flood = sum(CNSPOPL3)) %>%
  filter(flood_n_cat > 0) %>%
  mutate(quint_pop_flood_perc = round(100*quint_pop_flood/quint_totpop, 1)) %>%
  dplyr::select(educ_quint, quint_totpop, quint_pop_flood, quint_pop_flood_perc, flood_n_cat)

all_plot$flood_n_cat <- factor(all_plot$flood_n_cat, levels = c(1,2), ordered = TRUE)

ggplot(all_plot, aes(x=educ_quint, y = quint_pop_flood_perc, fill = as.factor(flood_n_cat))) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("% of residents in flooded neighborhoods") +
  xlab("Neighborhood educational attainment (1=lowest)") +
  scale_fill_manual(name = "Observed floods", labels = c("Single event", "Multiple events"), values = c("#1f83b4", "#ba43b4")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.8,.85),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

#ggsave("./Figures/bar_multiple_final.png", plot = last_plot(), width = 7.5, height = 6.5, unit = "in")

}  

######Within-city disparities 
###Table 3 - Slope Index of Inequality
###create SVI quintiles by city, standardized and flipped to range (0 high, 1 low)

{
  all_std <- all %>%
    mutate(svi = CNSMINPR_L3) %>%
    group_by(SALID1) %>% 
    mutate(svi_quint=as.numeric(cut(svi, breaks=quantile(svi, probs=seq(0, 1, by=1/5),
                                                         na.rm=T), include.lowest=T)),
           svi_quint_f=factor(svi_quint, levels=c(1:5),
                              labels=c("Low", "Low-Medium", "Medium", "Medium-High", "High")),
           svi_std = case_when(
             svi_quint_f == "Low" ~ 1,
             svi_quint_f == "Low-Medium" ~ 0.75,
             svi_quint_f == "Medium" ~ 0.5,
             svi_quint_f == "Medium-High" ~ 0.25,
             svi_quint_f == "High" ~ 0
           ))
  
  View(select(all_std, SALID1, id, contains("svi")))
  
  
  ###FOR ANY FLOODING
  city_slopes_models<-all_std %>% 
    group_by(SALID1, country) %>% 
    mutate(city_flood = sum(any_floods_y1n0)) %>% ##limit to cities that have experienced flooding
    filter(city_flood > 0) %>%
    group_modify(~{
      sii <- glm(any_floods_y1n0 ~ svi_std, data = .x, family = poisson) #city-specific model is here
      tibble(
        type=rep(c("adj_std"), times=1),
        outcome=c(rep("sii", times=1)),
        models=list(sii))
    })
  
  head(city_slopes_models)
  
  #exponentiate
  rii_sii<-city_slopes_models %>% 
    filter(type=="adj_std") %>% 
    group_by(SALID1, country) %>% 
    group_modify(~{
      cbind(.x %>% filter(outcome=="sii") %>% 
              pull(models) %>% `[[`(1) %>% tidy %>% 
              filter(term=="svi_std") %>% 
              mutate(sii=exp(estimate),
                     sii_lci=exp(estimate-(1.96*std.error)),
                     sii_uci=exp(estimate+(1.96*std.error))) %>% 
              select(sii, sii_lci, sii_uci))
    }) %>% 
    mutate(SALID1=factor(SALID1))
  
  table2_sii <- rii_sii %>% 
    mutate(sii=round(sii, digits=2),
           sii_ci=paste0(round(sii_lci, digits=2), ",",
                         round(sii_uci, digits=2))) %>% 
    select(country, SALID1, sii, sii_ci) 
  
  ##Overall results
  (rii_sii_overall <- rii_sii %>%
      ungroup() %>%
      summarise(sii_median = quantile(sii, probs = 0.5),
                sii_p10 = quantile(sii, probs = 0.25),
                sii_p90 = quantile(sii, probs = 0.75))
  )
  
  ##Country results
  (rii_sii_quant <- rii_sii %>%
      group_by(country) %>%
      summarise(sii_median = quantile(sii, probs = 0.5),
                sii_p10 = quantile(sii, probs = 0.25),
                sii_p90 = quantile(sii, probs = 0.75))
  )
  
  ##By city size
  
  city_info <- all %>% #limit to L1/city-level vars
    select(contains("1")) %>%
    group_by(SALID1) %>%
    mutate(SALID1 = as.character(SALID1)) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    mutate(qpopl1 = ntile(BECTPOP2010L1UX, 5))
  
  rii_all <- left_join(rii_sii, city_info) %>%
    mutate(estimate = log(sii))
  
  (rii_sii_quant <- rii_all %>%
      group_by(qpopl1) %>%
      summarise(N = n(),
                sii_median = quantile(sii, probs = 0.5),
                sii_p10 = quantile(sii, probs = 0.25),
                sii_p90 = quantile(sii, probs = 0.75))
  )
  
  
  ##proportion of cities with sii > 1
  #Overall
  (rii_sii_overall <- rii_all %>%
      ungroup() %>%
      mutate(gt1 = if_else(sii > 1, 1, 0)) %>%
      summarise(mean = mean(gt1))
  )
  #by country
  (rii_sii_country <- rii_all %>%
      ungroup() %>%
      mutate(gt1 = if_else(sii > 1, 1, 0)) %>%
      group_by(country) %>%
      summarise(mean = mean(gt1))
  )
  #by city pop size  
  (rii_sii_country <- rii_all %>%
      ungroup() %>%
      mutate(gt1 = if_else(sii > 1, 1, 0)) %>%
      group_by(qpopl1) %>%
      summarise(mean = mean(gt1),
                N = n())
  )
  
}

####Modeling - Table 4
#Univariable - Table 4

{
  
###RUN BEFORE MODELING - bins climate zones into first-levels groups
all <- all %>% 
  mutate(climate_lab = case_when(
    climate_L1 %in% c(1:4) ~ "Tropical",
    climate_L1 %in% c(5:8) ~ "Arid",
    climate_L1 %in% c(9:17) ~ "Temperate")) %>%
  group_by(SALID1) 

all$climate_lab <- relevel(as.factor(all$climate_lab), ref = "Temperate")

m1 <- glmer(any_floods_y1n0 ~ zpopdens + country + (1|SALID1), data = all, family = binomial)
tidy(m1,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
  
m2 <- glmer(any_floods_y1n0 ~ zeduc + country + (1|SALID1), data = all, family = binomial)
tidy(m2,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m3 <- glmer(any_floods_y1n0 ~ zintdens + country + (1|SALID1), data = all, family = binomial)
tidy(m3,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m4 <- glmer(any_floods_y1n0 ~ zgreen + country + (1|SALID1), data = all, family = binomial)
tidy(m4,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m5 <- glmer(any_floods_y1n0 ~ zdistance + country + (1|SALID1), data = all, family = binomial)
tidy(m5,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m6 <- glmer(any_floods_y1n0 ~ BECCOASTL3 + country + (1|SALID1), data = all, family = binomial)
tidy(m6,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m7 <- glmer(any_floods_y1n0 ~ zelev + country + (1|SALID1), data = all, family = binomial)
tidy(m7,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m8 <- glmer(any_floods_y1n0 ~ zslope + country + (1|SALID1), data = all, family = binomial)
tidy(m8,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m9 <- glmer(any_floods_y1n0 ~ zpop_L1 + country + (1|SALID1), data = all, family = binomial)
tidy(m9,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m10 <- glmer(any_floods_y1n0 ~ zpopdens_L1 + country + (1|SALID1), data = all, family = binomial)
tidy(m10,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m11 <- glmer(any_floods_y1n0 ~ zeduc_L1 + country + (1|SALID1), data = all, family = binomial)
tidy(m11,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m12 <- glmer(any_floods_y1n0 ~ zintdens_L1 + country + (1|SALID1), data = all, family = binomial)
tidy(m12,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m13 <- glmer(any_floods_y1n0 ~ zgreen_L1 + country + (1|SALID1), data = all, family = binomial)
tidy(m13,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m14 <- glmer(any_floods_y1n0 ~ zgdp_L1 + country + (1|SALID1), data = all, family = binomial)
tidy(m14,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m15 <- glmer(any_floods_y1n0 ~ climate_lab + country + (1|SALID1), data = all, family = binomial)
tidy(m15,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

}


##multivariable - Table 4
{

m_all <- glmer(any_floods_y1n0 ~ zpopdens+zeduc+zintdens+zgreen+zdistance+BECCOASTL3+zelev+ zslope+
                 zpop_L1+zpopdens_L1+zeduc_L1+zintdens_L1+zgreen_L1+zgdp_L1+ climate_lab + country + (1|SALID1), 
               data = all, family = binomial)
m_all_output <- tidy(m_all,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

m_all_output_cln <- m_all_output %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(estimate = round(estimate, digits = 2),
         conf.low = round(conf.low, digits = 2),
         conf.high = round(conf.high, digits = 2)) %>%
  mutate(tab = paste0(conf.low, ", ", conf.high))

write_csv(m_all_output_cln, "./multiv_output_2024-03-13.csv")

}

#######################
###SUPPLEMENTAL MATERIAL
#######################

####Supp Table S2: Table of info for each city

coords_raw <- read_csv("../Data/L1coords.csv") %>%
  mutate(SALID1 = as.numeric(SALID1))

stab <- left_join(all, coords_raw) %>%
  dplyr::select(id, SALID1, ISO2, CNSPOPL3, L1Name, tot_events) %>%
  mutate(country_name = case_when(
    ISO2 == "AR" ~ "Argentina",
    ISO2 == "BR" ~ "Brazil",
    ISO2 == "CL" ~ "Chile",
    ISO2 == "CO" ~ "Colombia",
    ISO2 == "CR" ~ "Costa Rica",
    ISO2 == "GT" ~ "Guatemala",
    ISO2 == "MX" ~ "Mexico",
    ISO2 == "PA" ~ "Panama"
  ))

tab <- stab %>% 
  group_by(country_name, L1Name) %>%
  summarise(N = n(),
            L1pop = sum(CNSPOPL3),
            L3pop_med = round(median(CNSPOPL3),0),
            L1_max = max(tot_events))

write_csv(tab, "./Suppl_tab_2023-11-17.csv")

#####Supp Table S3: Sensitivity analyses
#Model associations but for cities with 2+ floods only

{
  ###RUN BEFORE MODELING
  all_rpt <- all_rpt %>% 
    mutate(climate_lab = case_when(
      climate_L1 %in% c(1:4) ~ "Tropical",
      climate_L1 %in% c(5:8) ~ "Arid",
      climate_L1 %in% c(9:17) ~ "Temperate")) %>%
    group_by(SALID1) 
  
  all_rpt$climate_lab <- relevel(as.factor(all_rpt$climate_lab), ref = "Temperate")
  
  m_all_rpt <- glmer(any_floods_y1n0 ~ zpopdens+zeduc+zintdens+zgreen+zdistance+BECCOASTL3+zelev+ zslope+
                       zpop_L1+zpopdens_L1+zeduc_L1+zintdens_L1+zgreen_L1+zgdp_L1+ climate_lab + country + (1|SALID1), 
                     data = all_rpt, family = binomial)
  m_all_rpt_output <- tidy(m_all_rpt,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
  
  m_all_rpt_output_cln <- m_all_rpt_output %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(estimate = round(estimate, digits = 2),
           conf.low = round(conf.low, digits = 2),
           conf.high = round(conf.high, digits = 2)) %>%
    mutate(tab = paste0(conf.low, ", ", conf.high))
  
  write_csv(m_all_rpt_output_cln, "./multiv_output_2+_2024-03-13.csv")
  
}


####Supp Figure S1: Correlation figure

{
library(scales)  # NEED FOR TRANSPARENCY IN PLOTS (alpha() FUNCTION)          #
library(gridExtra)  # for plotting side-by-side plots                         # 
library(tidyverse)
library(GGally)

# CONTEXTUAL VARIABLES, RE-SCALED BY 1 SD, CENTERED AROUND MEAN.
context_vars <- c("zpopdens","zeduc","zintdens","zgreen","zdistance","BECCOASTL3","zelev","zslope","zpop_L1","zpopdens_L1","zeduc_L1","zintdens_L1","zgreen_L1","zgdp_L1")
context_var_labels <- c("Neighb. pop. density", "Neighb. educ.", "Neighb. int. density", "Neighb. green", "Neighb. distance", "Neighb. coastal", "Neighb. elevation", "Neighb. slope", "City pop. size", "City pop. density" , "City educ.", "City int. density", "City green", "City GDP")


cbsa_data <- all[,context_vars]
theme_text_sizes <- function(){
  theme(title        = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = 'black', size = 11),
        axis.title.x = element_text(color = 'black', size = 11),
        axis.text.y  = element_text(color = 'black', size = 11),
        axis.text.x  = element_text(color = 'black', size = 11),
        strip.text   = element_text(color = "black", size = 11),
        legend.text  = element_text(color = "black", size = 11))
}
cbsa_corr_plot <- 
  ggpairs(data = cbsa_data,
          columns = context_vars,
          columnLabels = context_var_labels,
          upper = list(continuous = 
                         wrap("cor", size=3, color="black", stars=F)),
          lower = "blank",
          diag = "blank") +
  theme_bw() +
  theme_text_sizes() +
  theme(
    # HORIZONTAL FACET TEXT FOR Y AXIS
    strip.text.y = element_text(angle=0),
    # Vertical FACET TEXT FOR X AXIS
    strip.text.x = element_text(angle=90),
    # LIGHTER GRAY BACKGROUND COLOR FOR FACET PANEL LABELS
    strip.background = element_blank(),
    axis.line = element_blank())

cbsa_corr_plot
pdf("./Figures/corr_plot_2024-01-31.pdf", height=10, width=10)
cbsa_corr_plot
dev.off()

ggsave("./Figures/corr_plot.png", plot = cbsa_corr_plot, height=10, width=10)
}

