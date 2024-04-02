# Packages and Functiosn --------------------------------------------------
library(tidyverse)

theme_set(theme_bw())
pal <- RColorBrewer::brewer.pal(8, "Dark2")


# Data --------------------------------------------------------------------
# Importing and removing fixed variables; Saving as RDS:
if (FALSE) {
  data <- read.csv("Data/ipumsi_data.csv") %>%
    as_tibble()
  
  data <- data %>%
    select(-all_of(map_dbl(data, ~length(unique(.x))) %>% {names(.)[. == 1]}))
  
  saveRDS(data, "Data/data.RDS")
}

# Loading data:
data <- readRDS("Data/data.RDS")

# Transforming missing to NA; Creating new variables:
data <- data %>%
  filter(INCTOT > 0) %>%
  mutate(across(where(is.numeric) & !SERIAL,
    ~ifelse(grepl("99+8?", as.character(.x)), NA, .x))
  ) %>%
  mutate(
    INCREST = INCTOT - INCWEL - INCRET - INCEARN,
    .after = INCTOT
  )

data_house <- data %>%
  group_by(as.factor(SERIAL)) %>%
  summarise(
    #across(everything(), mean),
    SIBLINGS = sum(RELATE == 3),
    TWINS = any(duplicated(AGE[RELATE == 3])),
    INCTOT = sum(INCTOT, na.rm = TRUE),
    INCTOTLOG = log(sum(INCTOT) + 1, base = 10)
  ) %>%
  mutate(
    INCTOTCUT = cut(INCTOT, c(0,25000,50000,100000,Inf)),
  ) %>%
  filter(SIBLINGS > 0)

# Other possible definitions of bins:
#INCTOTCUT = cut(INCTOT, quantile(INCTOT, seq(0, 1, length = 6), na.rm = TRUE))


# Question 2 --------------------------------------------------------------
forms_q1 <- list(
  "Linear" = SIBLINGS ~ INCTOT,
  "Quadratic" = SIBLINGS ~ INCTOT + I(INCTOT^2),
  "Linear + bins" = SIBLINGS ~ INCTOT + INCTOTCUT,
  "Linear * bins" = SIBLINGS ~ INCTOT*INCTOTCUT
)

mods_q1 <- map(forms_q1, ~lm(.x, data_house))

mods_res_q1 <- imap_dfr(mods_q1, function(m, name) {
  tibble(
    model = name,
    fit = fitted(m),
    res = residuals(m),
    INCTOT = m$model[,"INCTOT"]
  )
})

ggplot(data_house, aes(INCTOT, SIBLINGS)) +
  geom_point(alpha = 0.1, position = position_jitter(height = 0.1)) +
  geom_line(aes(y = fit, color = model), mods_res_q1,
    linewidth = 1, alpha = 0.7
  ) +
  ylim(0.9, 4.1) + xlim(0,500000) +
  scale_color_manual(values = pal) +
  labs(
    title = "Siblings versus Total Income",
    x = "Total Income",
    y = "Siblings",
    color = "Model"
  )

