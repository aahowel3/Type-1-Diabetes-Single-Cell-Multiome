#use include = False to comment out an entier code chunk 
#fitting logistic regression with donor as random variable
library(lme4)
library(multcomp)
library(tidyr)
library(dplyr)
library(ggplot2)

#data load in from npod1_reactome_pathwaysignifigance.RMD
#scroll down for 7/2/25 updates

df_beta_kegg4 = df_beta_kegg3 %>%
  mutate(across(starts_with('KEGG'), ~ifelse( .x > 0, 1, 0)))
  
  group_mean<- aggregate(x= df_beta_kegg4$KEGG.ANTIGEN.PROCESSING.AND.PRESENTATION,
                         # Specify group indicator
                         by = list(df_beta_kegg4$group),      
                         # Specify function (i.e. mean)
                         FUN = lite::fitBernoulli)
  
  group_mean = unlist(group_mean$x[,"mle"])
  
  df_beta_kegg4 = df_beta_kegg4 |> 
    mutate(Value = case_when(group== "Aab" ~ group_mean[1],
                             group=="ND" ~ group_mean[2],
                             group=="T1Dearly" ~ group_mean[3],
                             group == "T1Dlate"~ group_mean[4]))
  
  df_beta_kegg4$Value = round(1-(round(df_beta_kegg4$Value, digits = 4)), digits = 4)
  
  df_beta_kegg3$Value = df_beta_kegg4$Value
  df_beta_kegg5 = df_beta_kegg3[df_beta_kegg3$group != "ND",]
  
  ggplot(df_beta_kegg5, aes(x = group, y = df_beta_kegg5$KEGG.ANTIGEN.PROCESSING.AND.PRESENTATION, fill= group)) + 
    labs(y="UCell Pathway Enrichment Score") + 
    geom_violin(trim=FALSE) +
    #y value here changes where label is placed otherwise plots a label per point 
    geom_text(aes(label = Value, y=0.12)) +
    ggtitle("KEGG ANTIGEN PROCESSING AND PRESENTATION") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_discrete(name = "Disease State", labels = c("Multiple AAB+", "Early Stage T1D", "Late Stage T1D"))
  
  
  
###same figure for a continuous plot 
  ###selected apthway jak-stat signaling  
  
  group_mean<- aggregate(x= df_beta_kegg4$KEGG.JAK.STAT.SIGNALING.PATHWAY,
                         # Specify group indicator
                         by = list(df_beta_kegg4$group),      
                         # Specify function (i.e. mean)
                         FUN = lite::fitBernoulli)
  
 
   group_mean = unlist(group_mean$x[,"mle"])
  
  df_beta_kegg4 = df_beta_kegg4 |> 
    mutate(Value = case_when(group== "Aab" ~ group_mean[1],
                             group=="ND" ~ group_mean[2],
                             group=="T1Dearly" ~ group_mean[3],
                             group == "T1Dlate"~ group_mean[4]))
  
  df_beta_kegg4$Value = round(1-(round(df_beta_kegg4$Value, digits = 4)), digits = 4)
  
  df_beta_kegg3$Value = df_beta_kegg4$Value
  df_beta_kegg5 = df_beta_kegg3[df_beta_kegg3$group != "ND",]
  
  ggplot(df_beta_kegg5, aes(x = group, y = df_beta_kegg5$KEGG.JAK.STAT.SIGNALING.PATHWAY, fill= group)) + 
    labs(y="UCell Pathway Enrichment Score") + 
    geom_violin(trim=FALSE) +
    geom_text(aes(label = Value, y=0.15)) +
    ggtitle("KEGG JAK-STAT SIGNALING PATHWAY") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_discrete(name = "Disease State", labels = c("Multiple AAB+", "Early Stage T1D", "Late Stage T1D"))
  
  
#yes really you needed to run all the inputs again for the reactome pathways
  #file loadings are in npod1_reactome_pathwaysignifigance.Rmd
  #just run the whole first code chunk
  
  
  df_beta_kegg4 = df_beta_kegg3 %>%
    mutate(across(-c(group, donor), ~ifelse( .x > 0, 1, 0)))
  
  group_mean<- aggregate(x= df_beta_kegg4$REACTOME.INTERFERON.SIGNALING,
                         # Specify group indicator
                         by = list(df_beta_kegg4$group),      
                         # Specify function (i.e. mean)
                         FUN = lite::fitBernoulli)
  
  
  group_mean = unlist(group_mean$x[,"mle"])
  
  df_beta_kegg4 = df_beta_kegg4 |> 
    mutate(Value = case_when(group== "Aab" ~ group_mean[1],
                             group=="ND" ~ group_mean[2],
                             group=="T1Dearly" ~ group_mean[3],
                             group == "T1Dlate"~ group_mean[4]))
  
  df_beta_kegg4$Value = round(1-(round(df_beta_kegg4$Value, digits = 4)), digits = 4)
  
  df_beta_kegg3$Value = df_beta_kegg4$Value
  df_beta_kegg5 = df_beta_kegg3[df_beta_kegg3$group != "ND",]
  
  ggplot(df_beta_kegg5, aes(x = group, y = df_beta_kegg5$REACTOME.INTERFERON.SIGNALING, fill= group)) + 
    labs(y="UCell Pathway Enrichment Score") + 
    geom_violin(trim=FALSE) +
    geom_text(aes(label = Value, y=0.17)) +
    ggtitle("REACTOME INTERFERON SIGNALING") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_discrete(name = "Disease State", labels = c("Multiple AAB+", "Early Stage T1D", "Late Stage T1D"))
  
  
  
######  
#####  
#####  
  df_beta_kegg5 = df_beta_kegg2[df_beta_kegg2$group != "ND",]
#plot updates 7-2-25 - reviewer wants to spread per donor  
#starting from  df_beta_kegg2 generated  from npod1_reactome_pathwaysignifigance.RMD
  # Step 1: Identify pathway columns (exclude metadata)
  pathway_cols <- df_beta_kegg5 %>%
    select(-donor, -group) %>%
    select(where(is.numeric)) %>%
    colnames()
  
  # Step 2: Calculate proportion of zeros per donor per pathway
  # Store all results in a donor-level table
  donor_zero_props <- df_beta_kegg5 %>%
    select(donor, group) %>%
    distinct()
  
  for (pathway in pathway_cols) {
    prop_df <- df_beta_kegg5 %>%
      group_by(donor) %>%
      summarise(
        !!paste0(pathway, "_zero_prop") :=
          sum(.data[[pathway]] == 0, na.rm = TRUE) /
          sum(!is.na(.data[[pathway]])),
        .groups = "drop"
      )
    
    donor_zero_props <- left_join(donor_zero_props, prop_df, by = "donor")
  }
  
  # Step 3: Join proportion columns back to original df (preserves all rows)
  df_beta_kegg5 <- left_join(df_beta_kegg5, donor_zero_props, by = c("donor", "group"))


  ######
  ####
  ####actual plitting 
  
  plot_kegg_violin <- function(df, pathway, zero_prop, title_text = pathway) {
    library(dplyr)
    library(ggplot2)
    
    # Clean column names just in case
    names(df) <- trimws(names(df))
    
    # Set group order
    df$group <- factor(df$group, levels = c("Aab", "T1Dearly", "T1Dlate"))
    
    # Order donors by group
    donor_order <- df %>%
      distinct(donor, group) %>%
      arrange(group, donor) %>%
      pull(donor)
    df$donor <- factor(df$donor, levels = donor_order)
    
    # Calculate label position and extract zero proportion per donor
    label_df <- df %>%
      group_by(donor) %>%
      summarise(
        zero_prop = first(.data[[zero_prop]]),
        y_pos = max(.data[[pathway]], na.rm = TRUE) + 0.02,
        .groups = "drop"
      ) %>%
      mutate(
        label_text = sprintf("%.2f", zero_prop),
        donor = factor(donor, levels = donor_order)
      )
    
    # Plot
    ggplot(df, aes(x = donor, y = .data[[pathway]], fill = group)) +
      geom_violin(trim = FALSE) +
      geom_text(
        data = label_df,
        aes(x = as.numeric(donor) + 0.2, y = y_pos, label = label_text),
        inherit.aes = FALSE, size = 3, angle = 90, hjust = 0
      ) +
      labs(y = "UCell Pathway Enrichment Score") +
      ggtitle(title_text) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.ticks.x = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      scale_fill_discrete(
        name = "Disease State",
        labels = c("Multiple AAB+", "Early Stage T1D", "Late Stage T1D")
      )
  }
  
  plot_kegg_violin(
    df = df_beta_kegg5,
    pathway = "KEGG.ANTIGEN.PROCESSING.AND.PRESENTATION",
    zero_prop = "KEGG.ANTIGEN.PROCESSING.AND.PRESENTATION"
  )
  
 
