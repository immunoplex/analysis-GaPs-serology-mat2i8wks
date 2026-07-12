# ============================================================
# Read Figure 3 (panels A, B, C) digitized data into R
# de Graaf et al., Lancet Microbe 2026
# ============================================================

# ── Option 1: from saved file ─────────────────────────────────
# Save the CSV above as "figure3_ABC_data.csv" then run:
df <- read.csv(here::here("./data/figure3_abc_data.csv"))

# ── Option 2: inline (paste CSV text as string) ───────────────
# df <- read.csv(text = "...full csv text here...",
#                stringsAsFactors = FALSE)

# ── Enforce correct column types ─────────────────────────────
df$subject_id          <- as.integer(df$subject_id)

df$antibody            <- factor(df$antibody,
                                 levels = c("IgG", "IgA",
                                            "WT_IgG", "WT_IgA"))

df$antigen             <- factor(df$antigen,
                                 levels = c("PT", "PRN",
                                            "FHA", "FIM", "WT_Bp"))

df$compartment         <- factor(df$compartment,
                                 levels = c("serum", "MLF"))

df$colonization_status <- factor(df$colonization_status,
                                 levels = c("resistant", "colonized"))

df$concentration       <- as.numeric(df$concentration)

# ── Quick sanity checks ───────────────────────────────────────
str(df)
summary(df)

# How many rows per panel?
# Panel A: IgG  x {PT,PRN,FHA,FIM} x serum        = 4*50 = 200
# Panel B: IgA  x {PT,PRN,FHA,FIM} x serum        = 4*50 = 200
# Panel C: WT_IgG/WT_IgA x WT_Bp x {serum, MLF}  = 2*2*50 = 200
cat("Total rows:", nrow(df), "\n")          # expect 600
cat("Panel A rows:", nrow(subset(df, antibody %in% c("IgG")
                                 & antigen != "WT_Bp")), "\n")
cat("Panel B rows:", nrow(subset(df, antibody %in% c("IgA")
                                 & antigen != "WT_Bp")), "\n")
cat("Panel C rows:", nrow(subset(df, antigen == "WT_Bp")), "\n")

# ── Verify geometric means match Table 1 anchors ─────────────
library(dplyr)

df %>%
  filter(antibody == "IgG", antigen != "WT_Bp") %>%
  group_by(antigen, colonization_status) %>%
  summarise(GM = exp(mean(log(concentration))), .groups = "drop")
# Expected from Table 1:
#   PT   colonized ~1.9,  resistant ~6.2
#   PRN  colonized ~7.3,  resistant ~17.6
#   FHA  colonized ~8.5,  resistant ~30.1
#   FIM  colonized ~28.2, resistant ~25.4  (p=0.97, no diff)

# ── Optional reproduction plot (mirrors Figure 3A) ───────────
# install.packages("ggplot2")   # if needed
library(ggplot2)

panel_A <- df %>%
  filter(antibody == "IgG", antigen != "WT_Bp")

ggplot(panel_A,
       aes(x = colonization_status, y = concentration,
           colour = colonization_status)) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.7) +
  stat_summary(fun = "median", geom = "crossbar",
               width = 0.4, colour = "black", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1000")
  ) +
  scale_colour_manual(
    values = c("resistant" = "#4575b4", "colonized" = "#d73027")
  ) +
  facet_wrap(~ antigen, nrow = 1) +
  labs(
    x = "Colonisation status",
    y = "Serum IgG concentration (IU/mL or AU/mL)",
    colour = NULL,
    title = "Figure 3A – Pre-inoculation serum IgG (digitized estimates)",
    caption = "⚠ Values are digitized approximations anchored to Table 1 GMs\nReplace with exact values from WebPlotDigitizer"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
