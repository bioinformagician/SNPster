data_path = "C:/Users/FrederikTolberg/OneDrive - Evosep/Skrivebord/vacation genomics project/gwas_catalog_v1.0-associations_e113_r2024-11-20.tsv"

library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

gwas_data = read.csv(data_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(gwas_data)
View(gwas_data)

gwas_data_filtered <- gwas_data %>%
  filter(!is.na(OR.or.BETA)) %>%
  mutate(STRONGEST.SNP.RISK.ALLELE = str_extract(STRONGEST.SNP.RISK.ALLELE, "(?<=\\-).*")) %>%
  filter(STRONGEST.SNP.RISK.ALLELE != "?") %>%
  filter(
    !is.na(RISK.ALLELE.FREQUENCY),
    str_detect(X95..CI..TEXT., "(?i)\\[.*\\]( unit increase| unit decrease)$") |
      str_detect(X95..CI..TEXT., "(?i)^\\[.*\\]$")
  ) %>%
  mutate(
    or = !str_detect(X95..CI..TEXT., "(?i)unit"),
    increase_or_decrease = case_when(
      str_detect(X95..CI..TEXT., "(?i)unit increase") ~ "increase",
      str_detect(X95..CI..TEXT., "(?i)unit decrease") ~ "decrease",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(REPLICATION.SAMPLE.SIZE),
    !str_detect(X95..CI..TEXT., "(?i)NR"),
    str_detect(SNPS, "(?i)rs"),
    OR.or.BETA < 100,
    nchar(STRONGEST.SNP.RISK.ALLELE) == 1
  ) %>%
  mutate(RISK.ALLELE.FREQUENCY = as.numeric(RISK.ALLELE.FREQUENCY)) %>% # Convert to numeric
  select(
    OR.or.BETA, STRONGEST.SNP.RISK.ALLELE, X95..CI..TEXT., DISEASE.TRAIT,
    SNPS, or, PUBMEDID, RISK.ALLELE.FREQUENCY, increase_or_decrease,MAPPED_GENE
  ) %>%
  mutate(
    population_score = case_when(
      # If 'or == TRUE', apply log transformation
      or ~ RISK.ALLELE.FREQUENCY * 2 * log(OR.or.BETA),
      # If 'or == FALSE' and 'increase_or_decrease == "decrease"', negate the value
      !or & increase_or_decrease == "decrease" ~ -RISK.ALLELE.FREQUENCY * 2 * OR.or.BETA,
      # If 'or == FALSE' and 'increase_or_decrease == "increase"', no negation
      !or & increase_or_decrease == "increase" ~ RISK.ALLELE.FREQUENCY * 2 * OR.or.BETA,
      # Default to NA or 0 if no condition is met
      TRUE ~ NA_real_
    )
  )
colnames(gwas_data_filtered) = c("odds_ratio", "risk_allele", "confidence_interval", "trait", "rsid", "or", "pmid", "risk_allele_frequency", "trait_direction","gene", "population_score")

all_or_traits = unique(gwas_data_filtered$trait[gwas_data_filtered$or==TRUE])

# Trait-Category Mapping
trait_category_mapping <- list(
  "Neurological Disorders" = c(
    # Original
    "Schizophrenia", "Restless legs syndrome", "Multiple sclerosis", 
    "Parkinson's disease", "Alzheimer's disease", 
    "Alzheimer's disease (late onset)", "Alzheimer's disease in APOE e4- carriers", 
    "Generalized epilepsy", "Essential tremor", "Myopia (pathological)", 
    "Idiopathic inflammatory myopathy", "Bipolar disorder", 
    "Neuromyelitis optica", "Neuromyelitis optica (AQP4-IgG-positive)", 
    "Creutzfeldt-Jakob disease (sporadic)", "Major depressive disorder", 
    "Frontotemporal dementia", "Frontotemporal dementia with GRN mutation",
    "Alzheimer's disease (APOE e4 interaction)",
    
    # Newly added
    "Cluster headache",
    "Decreased low contrast letter acuity in multiple sclerosis",
    "Hippocampal sclerosis of aging",
    "Oligoclonal band status in multiple sclerosis",
    "Opioid dependence",
    "Parkinson's disease (familial, age at onset)"
  ),
  
  "Autoimmune and Inflammatory Diseases" = c(
    # Original
    "Ulcerative colitis", "Psoriasis", "Sjögren's syndrome", 
    "Rheumatoid arthritis", "Systemic lupus erythematosus", 
    "Behcet's disease", "Takayasu arteritis", "Atopic march", 
    "Dermatomyositis or juvenile dermatomyositis", "Polymyositis", 
    "Vitiligo (early onset)", "Vitiligo (late onset)", 
    "Primary biliary cholangitis", "Primary sclerosing cholangitis", 
    "Allergic rhinitis", "Childhood steroid-sensitive nephrotic syndrome",
    "Adult asthma", "Adult onset Still's disease",
    "Allergic disease (asthma, hay fever or eczema)",
    "Arthritis (juvenile idiopathic)","Asthma",
    
    # Newly added
    "Birdshot chorioretinopathy",
    "Crohn's disease",
    "Eczema",
    "Hydrolysed wheat protein allergy",
    "Membranous nephropathy",
    "Peanut allergy",
    "Periodontitis",
    "Poor prognosis in Crohn's disease",
    "Refractory celiac disease type II",
    "Rheumatoid arthritis (rheumatoid factor and/or anti-cyclic citrullinated peptide seropositive)",
    "Thyroid-associated orbitopathy in graves' disease",
    "Vogt-Koyanagi-Harada syndrome"
  ),
  
  "Cancer" = c(
    # Original
    "Bladder cancer", "Breast cancer", "Esophageal cancer", 
    "Prostate cancer", "Lung cancer", "Colorectal cancer", 
    "Diffuse large B cell lymphoma", "Hodgkin's lymphoma", 
    "Pancreatic cancer", "Gallbladder cancer", 
    "Non-small cell lung cancer (survival)", 
    "Breast Cancer in BRCA1 mutation carriers", 
    "Breast cancer (estrogen-receptor negative, progesterone-receptor negative, and human epidermal growth factor-receptor negative)", 
    "Endometrial endometrioid carcinoma", "Non-melanoma skin cancer", 
    "B-cell acute lymphoblastic leukaemia (ETV6-RUNX1 positive)", 
    "B-cell acute lymphoblastic leukaemia (high-hyperdiploidy)", 
    "Neuroblastoma or malignant cutaneous melanoma", 
    "Colon cancer", "Metastasis at diagnosis in osteosarcoma",
    "B-cell acute lymphoblastic leukaemia",
    
    # Newly added
    "Breast cancer (estrogen-receptor negative)",
    "Cardia gastric cancer",
    "Chronic lymphocytic leukemia",
    "EGFR mutation-positive lung adenocarcinoma",
    "Endometrial cancer",
    "Gastric cancer",
    "Lung adenocarcinoma",
    "Melanoma",
    "Meningioma",
    "Multiple myeloma",
    "Myeloproliferative neoplasms",
    "Nodular sclerosis Hodgkin lymphoma",
    "Non-cardia gastric cancer",
    "Prostate cancer (survival)",
    "Survival in colorectal cancer (distant metastatic)",
    "Survival in colorectal cancer (non-distant metastatic)",
    "Urinary bladder cancer",
    "Uterine fibroids"
  ),
  
  "Metabolic Disorders" = c(
    # Original
    "Type 1 diabetes", "Type 2 diabetes", "Metabolic syndrome", 
    "Polycystic ovary syndrome", "Obesity (early onset extreme)", 
    "Gout", "Gout vs. Hyperuricemia", "Hyperuricemia", 
    "Proliferative diabetic retinopathy", 
    "Diabetes in response to antihypertensive drug treatment (treatment strategy interaction)",
    
    # Newly added
    "Renal overload gout",
    "Renal underexcretion gout"
  ),
  
  "CVD and Chronic Conditions" = c(
    # Original
    "Renal function and chronic kidney disease",
    "Chronic obstructive pulmonary disease", 
    "Chronic obstructive pulmonary disease (moderate to severe)", 
    "Chronic obstructive pulmonary disease (severe)", 
    "Coronary artery disease", "Coronary artery disease in type 1 diabetes", 
    "Hypertension", "Abdominal aortic aneurysm", 
    "Venous thromboembolism", 
    "Venous thromboembolism adjusted for sickle cell variant rs77121243-T", 
    "Aortic valve stenosis", "Retinal detachment or retinal break",
    "Acute-on-chronic liver failure in hepatitis B",  # <-- originally placed here in your old snippet, though it fits Infectious
    "Atrial fibrillation",  # Already in the old snippet
    "Bicuspid aortic valve without aortic dilation",
    
    # Newly added
    "Brugada syndrome",
    "Chronic kidney disease (chronic kidney disease vs normal or mildly reduced eGFR) in type 1 diabetes",
    "Chronic kidney disease (severe chronic kidney disease vs normal kidney function) in type 1 diabetes",
    "Chronic venous disease",
    "Coronary heart disease",
    "Idiopathic pulmonary fibrosis",
    "Intracerebral hemorrhage",
    "Ischemic stroke",
    "Kidney disease (early and late stages) in type 1 diabetes",
    "Kidney disease (early stage) in type 1 diabetes",
    "Kidney disease (end stage renal disease vs non-end stage renal disease) in type 1 diabetes",
    "Kidney disease (end stage renal disease vs normoalbuminuria) in type 1 diabetes",
    "Kidney disease (late stage) in type 1 diabetes",
    "Peripheral artery disease",
    "Resistance to antihypertensive treatment in hypertension",
    "Small vessel stroke"
  ),
  
  "Infectious Diseases" = c(
    # Original
    "Acute respiratory distress syndrome in sepsis",  # newly added, but shown in your snippet
    "Tuberculosis", "Chronic hepatitis B infection", 
    "Chronic hepatitis C infection", "Hepatitis B", 
    "Hepatocellular carcinoma in post hepatitis C eradication by interferon therapy", 
    "Liver injury in anti-retroviral drug treated HIV", 
    "Liver injury in anti-tuberculosis drug treatment", 
    "Liver injury in combined anti-retroviral and anti-tuberculosis drug-treated HIV with tuberculosis", 
    "Leishmaniasis (visceral)", "Buruli ulcer", 
    "Buruli ulcer (age at onset)", "Recurrence of mild malaria attacks", 
    "Recurrence of malaria infection (mild or asymptomatic)",
    
    # Newly added
    "Acute-on-chronic liver failure in hepatitis B" 
    # If you prefer it here (which is more logical given it's a complication of HBV),
    # you can *remove* it from Cardiovascular if that was an earlier misplacement.
  ),
  
  "Developmental and Structural Conditions" = c(
    # Original
    "Achilles tendinopathy", "Anterior cruciate ligament rupture", 
    "Sagittal craniosynostosis", "Temporomandibular joint disorder", 
    "Osteoarthritis", "Osteoarthritis (hand, severe)", 
    "Myopia (pathological)", "Adolescent idiopathic scoliosis" 
    # (already listed in old snippet but also newly mentioned)
  ),
  
  "Other Conditions" = c(
    # Original
    "Alcohol consumption (heavy vs. light/non-drinkers)", 
    "Acne (severe)", "End-stage renal disease in Type 1 diabetics", 
    "Longevity (85 years and older)", "Longevity (90 years and older)", 
    "Longevity (age >90th survival percentile)", 
    "Longevity (age >99th survival percentile)", 
    "Pre-treatment pain in head and neck squamous cell carcinoma", 
    "Response to platinum-based neoadjuvant chemotherapy in cervical cancer", 
    "Tooth agenesis", "Tooth agenesis (mandibular second premolars)", 
    "Tooth agenesis (maxillary lateral incisors)", 
    "Sasang constitutional medicine type (So-Eum)",
    "Birdshot chorioretinopathy",  # NOTE: was also placed in Autoimmune. If you'd rather keep it in Autoimmune, remove from here.
    
    # Newly added
    "Bortezomib-induced peripheral neuropathy in multiple myeloma",
    "Chronic central serous retinopathy",
    "Early spontaneous preterm birth",
    "Exfoliation glaucoma or exfoliation syndrome",
    "Gallstone disease",
    "Gene methylation in lung tissue",
    "Glaucoma (primary open-angle)",
    "Irritable bowel syndrome",
    "Macular telangiectasia type 2",
    "Post-term birth",
    "Response to Pazopanib in cancer (hepatotoxicity)",
    "Trastuzumab-induced cardiotoxicity in cancer",
    "Uterine fibroids" # If you prefer fibroids under "Cancers and Tumors," remove from here.
  )
)


# Mapping Traits

# Mapping Traits
trait_to_category <- sapply(all_or_traits, function(trait) {
  category <- names(trait_category_mapping)[sapply(trait_category_mapping, function(traits) trait %in% traits)]
  if (length(category) == 0) {
    "Other"
  } else {
    category[1]
  }
})

# Combine into a data frame
mapped_traits <- data.frame(
  trait = all_or_traits,
  Category = unname(trait_to_category)
)

View(mapped_traits)
# Ensure column names match before merging
colnames(mapped_traits) <- c("trait", "Category")

colnames(gwas_data_filtered)
colnames(mapped_traits)
# Perform a left join
gwas_data_filtered <- merge(gwas_data_filtered, mapped_traits, by = "trait", all.x = TRUE)
unique(gwas_data_filtered[gwas_data_filtered$or == F, ]$trait)

unique(gwas_data_filtered[gwas_data_filtered$Category=='Other' & gwas_data_filtered$or==T, ]$trait)

category_counts <- gwas_data_filtered %>%
  filter(or == TRUE) %>%
  count(Category) %>%
  arrange(desc(n))
category_counts

write.csv(gwas_data_filtered, "C:/Users/FrederikTolberg/OneDrive - Evosep/Skrivebord/vacation genomics project/gwas_data_filtered.csv")


personal_data = fread("C:/Users/FrederikTolberg/OneDrive - Evosep/Skrivebord/vacation genomics project/genome_Frederik_FangelTolberg_v5_Full_20241115141826/genome_Frederik_FangelTolberg_v5_Full_20241117223640.txt")


#rs12913832  blue eyes
snp_vector <- c("i3003626", "rs1815739", "rs12913832", "rs1042522", "rs53576")
vec <- c("rs1805007", "rs1805008")
vec <- c("rs3791686")


personal_data[personal_data$`# rsid` %in% vec,]

vec %in% personal_data$`# rsid`

gwas_data$SNP_ID_CURRENT = paste0("rs", gwas_data$SNP_ID_CURRENT)
snp_vector %in% gwas_data$SNP_ID_CURRENT

personal_data[personal_data$`# rsid` %in% "rs53576",]

colnames(personal_data) = c("rsid", "chromosome", "position", "genotype")
#merge
merged_data = merge(gwas_data_filtered, personal_data, by = "rsid", all = FALSE)

colnames(merged_data) = c("rsid", "odds_ratio", "risk_allele", "confidence_interval", "trait", "or", "pmid", "population_score", "chromosome", "position", "dunno", "personal_genotype")


counts_for_each_trait <- merged_data %>%
  group_by(trait) %>%
  summarise(count = n())


View(counts_for_each_trait)

merged_data_pruned <- merged_data %>%
  rowwise() %>%
  mutate(
    risk_status = case_when(
      # Check if the genotype is homozygous for the risk allele
      str_count(personal_genotype, risk_allele) == 2 ~ 2,
      # Check if the genotype contains the risk allele
      str_detect(personal_genotype, risk_allele) ~ 1,
      # Otherwise, the risk allele is not present
      TRUE ~ 0
    )
  ) %>%
  ungroup()

counts_long <- merged_data_pruned %>%
  # Filter rows with valid odds_ratio and short trait length
  filter(
    !is.na(odds_ratio),
    !is.na(trait),
  ) %>%
  # Create a risk_status category
  mutate(
    risk_status_category = case_when(
      risk_status <= 1 ~ "good",
      risk_status >= 1 ~ "bad"
    )
  ) %>%
  # Group by the category and summarize counts
  group_by(risk_status_category) %>%
  summarise(
    total = n(),
    .groups = "drop"  # Avoid grouping in the output
  )


ggplot(counts_long, aes(x = "", y = total, fill = factor(risk_status_category))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("good" = "#2c7a1b", "bad" = "#8f1f1f")) +
  labs(fill = "Risk Status", title = "Good vs Bad alleles") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )+geom_text(aes(label = scales::percent(total/sum(total))), position = position_stack(vjust = 0.5))


View(merged_data_pruned)

merged_data_pruned <- merged_data_pruned %>%
  mutate(final_risk = (log(odds_ratio)*risk_status)-population_score) %>%
  filter(!is.na(final_risk))



aggregated_true <- merged_data_pruned %>%
  group_by(trait) %>%
  summarise(
    rsid = first(rsid),  # Keep the first rsid
    final_risk = sum(final_risk, na.rm = TRUE),  # Product of nonzero `odds_ratio`
    n_ris_snps = n(),  # Count the number of SNPs per trait
    risk_allele = first(risk_allele),  # Keep the first risk_allele
    pmid = first(pmid),  # Keep the first pmid
    chromosome = first(chromosome),  # Keep the first chromosome
    personal_genotype = first(personal_genotype), 
    or = first(or),# Keep the first personal_genotype
    .groups = "drop"  # Ungroup after summarising
  ) %>%
  mutate(percentile = pnorm(final_risk)*100) %>%
  filter(or ==T)

View(aggregated_true)


hist(aggregated_true$percentile)
#final_data <- bind_rows(aggregated_false, aggregated_true)
final_data <- aggregated_true

final_data_or = final_data %>%
  arrange(percentile)



plot_sd_diff <- function(sd_diff, trait) {
  # Create data frame for a standard normal distribution
  df <- data.frame(x = seq(-4, 4, length.out = 1000))
  df$y <- dnorm(df$x, mean = 0, sd = 1)
  
  # Calculate the percentile corresponding to sd_diff
  percentile_value <- pnorm(sd_diff) * 100
  
  ggplot() +
    # Fill area to the left of sd_diff in green
    geom_area(data = subset(df, x <= sd_diff), aes(x = x, y = y), fill = "#F7A399") +
    # Fill area to the right of sd_diff in red
    geom_area(data = subset(df, x >= sd_diff), aes(x = x, y = y), fill = "#B4E197") +
    # Draw the line on top
    geom_line(data = df, aes(x = x, y = y), color = "black", size = 1) +
    # Add vertical line at sd_diff
    geom_vline(xintercept = sd_diff, color = "black", linetype = "dashed") +
    # Add percentile label near the vertical line
    geom_text(
      x = sd_diff,
      y = max(df$y) * 0.9,
      label = paste0(round(percentile_value, 2), "%"),
      color = "red",
      hjust = ifelse(sd_diff > 0, -0.1, 1.1),
      vjust = 0
    ) +
    labs(
      title = paste("Your Assesement for", trait),
      subtitle = paste("Population Risk Percentile:", round(percentile_value, 2), "%"),
      x = "Z-score (SD Units)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plot_sd_diff(final_data_or$final_risk[190], final_data_or$trait[190])


final_data_or <- final_data_or %>%
  select(-confidence_interval)
View(final_data_or)
#write to csv
write.csv(final_data_or, "final_data_or.csv", row.names = FALSE)