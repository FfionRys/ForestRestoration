##packages##
library(ggplot2)
library(tidyverse)
library(treemap)
library(cowplot)
library(dplyr)
library(forcats)
library(patchwork)
library(scales)
#load data
df <- read_csv("Updated data.csv")
#give papers ID code
#Paper ID code
df<-df%>%
  group_by(Title) %>%
  mutate(StudyID = cur_group_id()) %>%
  ungroup()
#create df of papers using info needed for year and location plots
papers <- df %>%
  select(StudyID, Title, Authors, Year, Location) %>%
  distinct()
#year
papers_by_year <- papers %>%
  count(Year) %>%
  arrange(Year)
integer_breaks <- function(...) {
  function(x) {
    seq(floor(x[1]), ceiling(x[2]), by = 1)
  }
}
year_plot <- ggplot(papers_by_year, aes(x = as.factor(Year), y = n)) +
  geom_col(fill = "#1f78b4", width = 0.7) +
  labs(
    title = "Publications per Year",
    x = "Year",
    y = "No. of Papers"
  ) +
  scale_y_continuous(
    breaks = integer_breaks(),
    labels = number_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_cowplot()
#location
papers_by_location <- papers %>%
  count(Location) %>%
  arrange(desc(n))
location_plot <- ggplot(papers_by_location, aes(x = reorder(Location, n), y = n)) +
  geom_col(fill = "#1f78b4", width = 0.7) +
  coord_flip() +
  labs(
    title = "Locations of Meta-Analyses",
    x = "Location",
    y = "No. of Papers"
  ) +
  theme_cowplot()
#combine plots
plot_grid(year_plot, location_plot, labels=c("(a)","(b)"))
#create count data for variables
count_data <- df %>%
  group_by(Focus, variable) %>%
  summarise(paper_count = n(), .groups = "drop")
#faceted plot for variables
ggplot(count_data, aes(x = reorder(variable, -paper_count), y = paper_count)) +
  geom_bar(stat = "identity", fill = "#1f78b4") +
  facet_wrap(~ Focus, scales = "free") +  
  coord_flip() +
  labs(title = "Number of Papers per Variable by Focus",
       x = "Variable",
       y = "Number of Papers") +
  theme_cowplot()
#frequency table for biodiversity variables
freq_table <- biodiversity_df %>%
  group_by(variable) %>%
  summarise(Frequency = n())
freq_table$label <- paste0(freq_table$variable, "\n(", freq_table$Frequency, ")")
##grouped by richness and abundance##
freq_table$Category <- ifelse(grepl("abundance", freq_table$label, ignore.case = TRUE),
                              "Abundance",
                              "Richness")
#create biodiversity treemap
biodiv_treemap <- treemap(freq_table,
                          index = c("Category", "label"),  
                          vSize = "Frequency",            
                          vColor = "Frequency",              
                          type = "value",                    
                          title = "Biodiversity Variables",
                          palette = "YlGn",
                          fontsize.labels = c(0, 12),        
                          fontcolor.labels = c("transparent", "white"),
                          fontface.labels = c(1, 1),
                          border.col = "black",
                          align.labels = list(c("center", "center"), c("center", "center"))
)

#change labels in soil and carbon to full names
variable_names <- c(
  TN  = "Total Nitrogen",
  AP  = "Available Phosphorus",
  SMC = "Soil Moisture Content",
  TP  = "Total Phosphorus",
  SOC = "Soil Organic Carbon",
  AN  = "Available Nitrogen",
  BD  = "Bulk Density",
  TC  = "Total Carbon",
  ABGC= "Aboveground Carbon",
  SIC = "Soil Inorganic Carbon"
)

#create combined frequency table for soil + carbon
combined_freq <- df %>%
  filter(tolower(Focus) %in% c("carbon", "soil properties")) %>%
  group_by(variable) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  arrange(desc(Frequency)) %>%
  mutate(
    full_name = ifelse(is.na(variable_names[variable]),
                       variable,
                       variable_names[variable]),
    label = paste0(full_name, "\n(", Frequency, ")")
  )

#plot treemap for soil + carbon
soil_carbon_treemap <- treemap(
  combined_freq,
  index = "label",
  vSize = "Frequency",
  vColor = "Frequency",
  type  = "value",
  palette = "Greens",
  title = "Carbon and Soil Properties Variables",
  fontsize.labels = 12,
  fontcolor.labels = "white",
  fontface.labels = 1,
  border.col = "black"
)

#effect size analysis#
#change some column names and filter out NA sample sizes
df <- df %>%
  mutate(
    EffectSize = as.numeric(`Effect size %`),
    Lower = as.numeric(`Lower Limit`),
    Upper = as.numeric(`Upper Limit`),
    SampleSize = as.numeric(`sample size`),
    Control = tolower(str_trim(`Control`))
  ) %>%
  filter(!is.na(SampleSize))
#correct confidence intervals omitting StudyID 5 + 23
df_subset1 <- subset(df, StudyID %in% c(5, 23))
df_subset2 <- subset(df, !StudyID %in% c(5, 23))

df2 <- df_subset2 %>%
  mutate(
    Lower_CI = if_else(Lower < 0, EffectSize + Lower, EffectSize - Lower),
    Upper_CI = if_else(Upper > 0, EffectSize + Upper, EffectSize - Upper)
  )

df_subset1$Lower_CI <- df_subset1$Lower
df_subset1$Upper_CI <- df_subset1$Upper

df <- rbind(df_subset1, df2)
#global biodiversity dataset
df_global_bio <- df %>%
  filter(tolower(Focus) == "biodiversity" & tolower(Location) == "global")%>%
  group_by(variable, Control) %>%
  filter(SampleSize == max(SampleSize, na.rm = TRUE)) %>%
  ungroup()
#create forest plot of global biodiversity
global_bio_plot <- ggplot(df_global_bio,
                          aes(y = reorder(variable, EffectSize, median), x = EffectSize)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Upper_CI, label = SampleSize),   
            hjust = -0.3,  
            size = 3) +
  labs(title = NULL, y = NULL, x = NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed")
#china biodiversity dataframe
df_china_bio <- df %>%
  filter(tolower(Focus) == "biodiversity" & tolower(Location) == "china")%>%
  group_by(variable, Control) %>%
  filter(SampleSize == max(SampleSize, na.rm = TRUE)) %>%
  ungroup()
#create forest plot of china biodiversity
china_bio_plot <- ggplot(df_china_bio,
                         aes(y = reorder(variable, EffectSize, median), x = EffectSize)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Upper_CI, label = SampleSize),
            hjust = -0.3,  
            size = 3) +
  labs(title = NULL,
       y = NULL, x = NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed")
#combine bio plots
combined_plot <- plot_grid(global_bio_plot, china_bio_plot, labels=c("(a) Global","(b) China"),  label_size = 12, label_y = 1.02)
#adjust details
final_plot <- ggdraw() +
  draw_plot(combined_plot, y = 0, height = 0.93) +
  draw_label("Percent Change", x = 0.6, y = 0.03, hjust = 0.5, size = 12 )

##soil + carbon analysis##
#change label names again
variable_names2 <- c(
  TN = "Total N",
  AP = "Available P",
  SMC = "Soil Moisture Content",
  TP = "Total P",
  SOC = "Soil Organic C",
  AN = "Available N",
  BD ="Bulk Density",
  TC = "Total C",
  ABGC = "Aboveground C",
  SIC = "Soil Inorganic C"
)
#global dataframe
df_global_soil_carb <- df %>%
  filter(
    tolower(Focus) %in% c("carbon", "soil properties"),
    tolower(Location) == "global",
  ) %>%
  group_by(variable, Control) %>%
  filter(SampleSize == max(SampleSize, na.rm = TRUE)) %>%
  ungroup()
# Reorder variables 
df_global_soil_carb$variable <- fct_reorder(df_global_soil_carb$variable, df_global_soil_carb$EffectSize, .fun = median)
#use variable_names2 
df_global_soil_carb <- df_global_soil_carb %>%
  mutate(variable = recode(variable, !!!variable_names2))
#plot
global_nutrients_plot <- ggplot(df_global_soil_carb,
                                aes(y = variable, x = EffectSize)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Upper_CI, label = SampleSize),
            hjust = -0.3,  
            size = 3) +
  labs(title = NULL,
       y = NULL, x = "Percentage change") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Control, scales = "fixed")
#china dataframe
df_china_soil_carb <- df %>%
  filter(
    tolower(Focus) %in% c("carbon", "soil properties"),
    tolower(Location) == "china",
  ) %>%
  group_by(variable, Control) %>%
  filter(SampleSize == max(SampleSize, na.rm = TRUE)) %>%
  ungroup()
#reorder variables
df_china_soil_carb$variable <- fct_reorder(df_china_soil_carb$variable, df_china_soil_carb$EffectSize, .fun = median)
#give variable_names2
df_china_soil_carb <- df_china_soil_carb %>%
  mutate(variable = recode(variable, !!!variable_names2))
#plot
china_nutrients_plot <- ggplot(df_china_soil_carb,
                               aes(y = variable, x = EffectSize)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Upper_CI, label = SampleSize),
            hjust = -0.3, 
            size = 3) +
  labs(title = "China Soil Properties, Nutrients, and Carbon",
       y = NULL, x = "Percentage change") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Control, scales = "fixed")
#Add regions to combine plots
df_global_soil_carb <- df_global_soil_carb %>%
  mutate(Region = "Global")
df_china_soil_carb <- df_china_soil_carb %>%
  mutate(Region = "China")
#create combined df
df_combined <- bind_rows(df_global_soil_carb, df_china_soil_carb)
#plot
combined_soil_plot <- ggplot(df_combined,
                             aes(y = variable, x = EffectSize, color = Region)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI),
                 height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_text(aes(x = Upper_CI, label = SampleSize),
            hjust = -0.3, size = 3,
            position = position_dodge(width = 0.5),
            show.legend = FALSE) +  
  labs(title = NULL,
       y = NULL, x = "Percentage change") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Control, scales = "fixed")

#brazil and tropics biodiversity df#
df_tropic_brazil_bio <- df %>%
  filter(
    tolower(Focus) == "biodiversity" &
      tolower(Location) %in% c("tropics", "brazil")
  ) %>%
  group_by(variable, Control) %>%
  filter(SampleSize == max(SampleSize, na.rm = TRUE)) %>%
  ungroup()
#reorder variables
df_tropic_brazil_bio <- df_tropic_brazil_bio %>%
  mutate(variable = factor(variable, levels = c(
    "Plant Richness", "Plant Abundance",
    "Bird Richness", 
    "Mammal Richness",
    "Reptile Richness",
    "Invertebrate Richness", "Invertebrate Abundance",
    "Overall Richness", "Overall Abundance"
  )))
#plot
trop_braz_bio_plot <- ggplot(df_tropic_brazil_bio,
                             aes(y = variable, x = EffectSize)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Upper_CI, label = SampleSize), hjust = -0.3, size = 3) +
  labs(title = NULL, y = NULL, x = "Percent Change") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed")

#tropics and Brazil soil/carbon df#
df_tropic_brazil_soil <- df %>%
  filter(
    tolower(Focus) == "soil properties" &
      tolower(Location) %in% c("tropics", "brazil")
  ) %>%
  group_by(variable, Control) %>%
  filter(SampleSize == max(SampleSize, na.rm = TRUE)) %>%
  ungroup()
#rename variables
variable_names3 <- c(
  TN = "Total N",
  AP = "Available P",
  ABGC = "Aboveground C"
)

df_tropic_brazil_soil <- df_tropic_brazil_soil %>%
  mutate(variable = recode(variable, !!!variable_names3))
#plot
trop_braz_soil_plot <- ggplot(df_tropic_brazil_soil,
                              aes(y = variable, x = EffectSize)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Upper_CI, label = SampleSize),
            hjust = -0.3,  
            size = 3) +
  labs(title = NULL,
       y = NULL, x = "Percent Change") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Control, scales = "fixed")
#combine plots
combined_plot2 <- plot_grid(trop_braz_bio_plot, trop_braz_soil_plot, labels=c("(a) Biodiversity","(b) Soil and Carbon"),  label_size = 12, label_y = 1.03)
final_plot2 <- ggdraw() +
  draw_plot(combined_plot2, y = 0, height = 0.93)


##view data used in plots##
df_china_bio %>%
  select(variable, EffectSize, Lower_CI, Upper_CI, StudyID) %>%
  arrange(variable)
df_global_bio %>%
  select(variable, EffectSize, Lower_CI, Upper_CI, StudyID) %>%
  arrange(variable)
df_global_soil_carb %>%
  select(variable, EffectSize, Lower_CI, Upper_CI) %>%
  arrange(variable)
df_china_soil_carb %>%
  select(variable, EffectSize, Lower_CI, Upper_CI) %>%
  arrange(variable)

##read dataframe for confounding variables##
dfnew <- read_csv("paper_info_new.csv")
#give Study ID
dfnew<-dfnew%>%
  group_by(Title) %>%
  mutate(StudyID = cur_group_id()) %>%
  ungroup()
#omit hydrology and combine carbon + soil
df_plot <- dfnew %>%
  filter(tolower(Focus) != "hydrology") %>%
  mutate(Focus = case_when(
    tolower(Focus) %in% c("carbon", "soil properties") ~ "Carbon & Soil Properties",
    TRUE ~ Focus
  ))
#frequency table
freq_table2 <- df_plot %>%
  group_by(Variable, Focus) %>%
  summarise(Frequency = n(), .groups = "drop")
#plot
ggplot(freq_table2, aes(x = reorder(Variable, -Frequency), y = Frequency, fill = Focus)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Focus, scales = "free_x") +
  labs(
    title = NULL,
    x = "Confounding Variable",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#view variables and correlations
dfnew %>%
  filter(!is.na(Correlation)) %>%
  select(Variable, Correlation, Control, StudyID) %>%
  arrange(Variable) %>%
  print(n = Inf)


#citations#
citation()
version$version.string
citation("ggplot2")
citation("tidyverse")
citation("forcats")
citation("scales")
citation("patchwork")
citation("treemap")
citation("cowplot")
citation("dplyr")
