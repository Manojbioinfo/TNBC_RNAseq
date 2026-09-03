library(dplyr)
data=readRDS("Result_Combined/Combined_eQTM_Results_ALL_annotated.rds")
head(data)


library(dplyr)
library(ggplot2)

data <- readRDS("Result_Combined/Combined_eQTM_Results_ALL_annotated.rds")

## ============================================================
## 1. Filter to significant CpG-Gene pairs (adjust threshold as needed)
## ============================================================

sig <- data %>% filter(FDR < 0.05)

## ============================================================
## 2. Classify direction: Hypomethylation vs Hypermethylation
##    (based on sign of t-statistic; negative t = hypomethylation
##     driving expression, i.e., negative correlation between
##     methylation and gene expression is the typical convention.
##     Adjust sign logic based on your own definition.)
## ============================================================

sig <- sig %>%
  mutate(Direction = ifelse(t < 0, "Hypomethylation", "Hypermethylation"))

## ============================================================
## 3. Clean up Relation_to_Island factor levels/order
## ============================================================

sig$Relation_to_Island <- factor(
  sig$Relation_to_Island,
  levels = c("Island", "N_Shore", "N_Shelf", "S_Shore", "S_Shelf", "OpenSea")
)

## ============================================================
## 4. Summarise counts and percentages per Relation_to_Island
## ============================================================

summary_df <- sig %>%
  filter(!is.na(Relation_to_Island)) %>%
  group_by(Relation_to_Island, Direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Relation_to_Island) %>%
  mutate(total = sum(n),
         pct = 100 * n / total) %>%
  ungroup()

## Total n per group (for labels above bars)
totals_df <- summary_df %>%
  distinct(Relation_to_Island, total)

## Ensure consistent stacking order (Hypomethylation at bottom)
summary_df$Direction <- factor(summary_df$Direction,
                               levels = c("Hypomethylation", "Hypermethylation"))

## ============================================================
## 5. Plot
## ============================================================

p <- ggplot(summary_df, aes(x = Relation_to_Island, y = pct, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.7, color = "white") +
  geom_text(aes(label = paste0(round(pct), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 5) +
  geom_text(data = totals_df,
            aes(x = Relation_to_Island, y = 103, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 4) +
  scale_fill_manual(values = c("Hypomethylation" = "#1F4E8C",
                               "Hypermethylation" = "#C0392B")) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  labs(x = "Relation to CpG Island",
       y = "Percentage of CpG-Gene Pairs (%)",
       fill = "") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

print(p)

ggsave("Result_Combined/eQTM_Direction_by_IslandRelation.pdf",
       p, width = 6, height = 5)
ggsave("Result_Combined/eQTM_Direction_by_IslandRelation.png",
       p, width = 6, height = 5, dpi = 300)
