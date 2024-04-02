# Packages and Functiosn --------------------------------------------------
library(tidyverse)
library(stargazer)

theme_set(theme_bw())
pal <- RColorBrewer::brewer.pal(8, "Dark2")

prettify_star <- function(star, low = 0.01, high = 1000, digits = 3, scipen = -7, ...) {
  star <- capture.output(star)
  mark  <- '::::'
  
  replace_numbers <- function(x) {
    x <- gsub(mark, '.', x)
    x.num <- as.numeric(x)
    ifelse(
      (x.num >= low) & (x.num < high), 
      round(x.num, digits = digits), 
      prettyNum(x.num, digits = digits, scientific = scipen, ...)
    )
  }    
  
  reg <- paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)")
  cat(gsubfn::gsubfn(reg, ~replace_numbers(x), star), sep = '\n')
}


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
data_raw <- readRDS("Data/data.RDS")

# Transforming missing to NA; Creating new variables:
data_raw <- data_raw %>%
  filter(INCTOT > 0) %>%
  mutate(across(where(is.numeric) & !SERIAL,
    ~ifelse(grepl("99+8?", as.character(.x)), NA, .x))
  )

data <- data_raw %>%
  group_by(as.factor(SERIAL)) %>%
  summarise(
    #across(everything(), mean),
    Siblings = sum(RELATE == 3),
    HasTwins = any(duplicated(AGE[RELATE == 3])),
    IncTot = sum(INCTOT, na.rm = TRUE)
  ) %>%
  filter(Siblings > 0)

ecdf_inc <- ecdf(data$IncTot)

data <- data %>%
  mutate(
    IncTotLog = log(IncTot + 1, base = 10),
    IncTotCut = cut(IncTot, c(0,25000,50000,100000,Inf)) %>%
      fct_relabel(., \(x) set_names(paste(1:4), levels(.))[x]),
    IncTotQuant = edf_inc(IncTot)
  )

# Other possible definitions of cut and quant:
#IncTotCut = cut(IncTot, quantile(IncTot, seq(0, 1, length = 6), na.rm = TRUE))
#IncTotQuant = IncTot/max(IncTot); m$model[,"IncTotQuant"] * max(data$IncTot)


# Question 2 --------------------------------------------------------------
forms_q1 <- list(
  "Linear" = Siblings ~ IncTot,
  "Linear + bins" = Siblings ~ IncTot + IncTotCut,
  "Linear * bins" = Siblings ~ IncTot*IncTotCut,
  "Quadratic" = Siblings ~ IncTot + I(IncTot^2),
  "Log" = Siblings ~ IncTot + IncTotLog,
  "Quantile" = Siblings ~ IncTotQuant,
  "Quantile Q." = Siblings ~ IncTotQuant + I(IncTotQuant^2)
)

mods_q1 <- map(forms_q1, ~lm(.x, data))

mods_res_q1 <- imap_dfr(mods_q1, function(m, name) {
  tibble(model = name, fit = fitted(m), res = residuals(m)) %>%
    mutate(
      IncTot = if (name %in% c("Quantile", "Quantile Q.")) {
        quantile(data$IncTot, m$model[,"IncTotQuant"])
      } else if (name %in% c("Log")) {
        10^m$model[,"IncTotLog"] - 1
      } else {
        m$model[,"IncTot"]
      },
      model_cat = if (grepl("Linear", name)) "Linear" else "Non-Linear",
      model = factor(model, levels = c(names(forms_q1)[1:3], " ", names(forms_q1)[4:8]))
    )
})

ggplot(data, aes(IncTot, Siblings)) +
  geom_point(alpha = 0.1, position = position_jitter(height = 0.1)) +
  geom_line(aes(y = fit, color = model), mods_res_q1,
    linewidth = 1, alpha = 0.7
  ) +
  facet_wrap(vars(model_cat), nrow = 1) +
  ylim(0.9, 4.1) + xlim(0,450000) +
  scale_color_manual(values = c(pal[1:3], "white", pal[4:7]), drop = FALSE) +
  labs(
    title = "Siblings versus Total Income",
    x = "Total Income",
    y = "Siblings",
    color = "Model"
  )

stargazer(mods_q1,
  #covariate.labels = custom_labels(x, FALSE),
  #dep.var.labels = custom_labels(x, TRUE),
  type = "latex",
  #single.row = TRUE,
  no.space = TRUE,
  omit.stat = c("f"), #"rsq", "adj.rsq",
  df = FALSE,
  omit.table.layout = "n",
  decimal.mark = "::::"
) %>%
  prettify_star()
