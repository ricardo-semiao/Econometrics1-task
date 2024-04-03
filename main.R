# Packages and Functiosn --------------------------------------------------
library(tidyverse)
library(rlang)
library(patchwork)
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

gaze <- function(x) {
  stargazer(x,
            type = "latex",
            no.space = TRUE,
            omit.stat = c("f"),
            df = FALSE,
            omit.table.layout = "n",
            decimal.mark = "::::"
  ) %>%
    prettify_star()
}

plot_smooths <- function(data) {
  create_sm <- function(formula, name, data, method = "lm"){
    list(geom_smooth(aes(Siblings, color = !!name), data,
                     formula = formula, method = method, se = FALSE, na.rm = TRUE, alpha = 0.7 
    ))
  }
  
  list(
    create_sm(y ~ x, "Linear", data = data),
    create_sm(y ~ x + I(x^2), "Quadratic", data = data),
    create_sm(y ~ I(as.factor(trunc(x))), "Dummies", data = data)
  )
}

data_gdensity <- function(...) {
  scalings <- table(data$Siblings) %>% {./max(.)}
  result <- list()
  
  result$data_density <- data %>%
    select(c(Siblings, !!!ensyms(...))) %>%
    group_split(Siblings) %>%
    map_dfr(function(group) {
      scalings_group <- scalings[group$Siblings[1]]
      imap_dfc(select(group, -Siblings), function(var, name) {
        d <- density(var, na.rm = TRUE)
        tibble(
          "{name}==value" := d$x,
          "{name}==dens" := d$y/max(d$y)*0.6*scalings_group
        ) %>%
          arrange("{name}==value")
      }) %>%
        mutate(Siblings = group$Siblings[1], .before = 1)
    }) %>%
    pivot_longer(-Siblings, names_sep = "==", names_to = c("var", ".value"))
  
  result$data_true <- data %>%
    select(c(Siblings, !!!ensyms(...))) %>%
    mutate(Siblings_jitter = Siblings + runif(n(), 0, 0.6)) %>%
    pivot_longer(-starts_with("Siblings"), names_to = c("var"))
  
  result
}

plot_gdensity <- function(graph_data) {
  list(
    geom_point(aes(y = Siblings_jitter), graph_data$data_true, alpha = 0.05),
    geom_path(aes(y = Siblings + dens, group = Siblings), linewidth = 0.75, color = pal[7])
  )
}


# Data --------------------------------------------------------------------
# Importing and removing fixed variables; Saving as RDS:
if (FALSE) {
  data_raw <- read.csv("Data/ipumsi_data.csv") %>%
    as_tibble()
  
  data_raw <- data_raw %>%
    select(-all_of(map_dbl(data_raw, ~length(unique(.x))) %>% {names(.)[. == 1]}))
  
  saveRDS(data_raw, "Data/data_raw.RDS")
}

# Loading data:
data_raw <- readRDS("Data/data_raw.RDS")

# Filtering out negative incomes; transforming missing codes to NA:
data_raw <- data_raw %>%
  filter(INCTOT >= 0) %>%
  mutate(across(where(is.numeric) & !SERIAL,
    ~ifelse(grepl("99+8?", as.character(.x)), NA, .x))
  )

# Setup:
age_max <- 21 #max age to be considered a "child of the household"
age_start <- 5 #age of possible school start
is_child <- expr(RELATE == 3 & AGE <= age_max)

# Aggregating data to the household level:
data <- data_raw %>%
  group_by(HouseID = as.factor(SERIAL)) %>%
  summarise(
    Families = mean(NFAMS),
    Siblings = sum(!!is_child), 
    ChildAges = list(AGE[!!is_child]),
    IncTot = sum(INCTOT, na.rm = TRUE),
    RoomsCapta = ROOMS[1]/PERSONS[1],
    YearsEduc = mean(YRSCHOOL[!!is_child] - age_start, na.rm = TRUE),
    StudentsCapta = sum(SCHOOL[!!is_child] == 1)/Siblings,
    HrsUsual = list(HRSUSUAL1[!!is_child] %>% ifelse(is.na(.), 0, .)),
    HasTwins = any(duplicated(AGE[!!is_child])),
  ) %>%
  ungroup() %>%
  filter(Siblings > 0 & Siblings < 8)

# Creating extra variables:
data <- data %>%
  filter(Families == 1) %>%
  mutate(
    IncTotLog = log(IncTot + 1, base = 10),
    IncTotCut = cut(IncTot, c(0,25000,50000,100000,Inf)) %>%
      fct_relabel(., \(x) set_names(paste(1:4), levels(.))[x]),
    IncTotQuant = ecdf(IncTot)(IncTot)
  ) %>%
  rowwise() %>%
  mutate(
    HasWorkers = any(HrsUsual[[1]] != 0),
    WorkedCapta = sum(HrsUsual[[1]])/Siblings
  ) %>%
  ungroup()

# Other possible variables for parenting quality:
#Own = OWNERSHIP[1]
#PresentParents = (MOMLOC[RELATE == 3][1] != 0) + (POPLOC[RELATE == 3][1] != 0),
#BedroomsCapta = BEDROOMS[1]/PERSONS[1],
#NonBedroomsCapta = (ROOMS[1] - BEDROOMS[1])/PERSONS[1],

# Other possible definitions of cut and quant:
#IncTotCut = cut(IncTot, quantile(IncTot, seq(0, 1, length = 6), na.rm = TRUE))
#IncTotQuant = IncTot/max(IncTot); m$model[,"IncTotQuant"] * max(data$IncTot)


# Exploratory Analysis ----------------------------------------------------
# Chekcing frequency of families per household
table(data$Families) / nrow(data)

# Analyzing ages of "childs":
data_raw %>%
  filter(RELATE == 3) %>%
  pull(AGE) %>%
  hist()

# Analyzing marital status of "childs"
data_raw %>%
  filter(RELATE == 3 & AGE < 21) %>%
  pull(MARST) %>%
  {table(.)/length(.)}

# Analyzing education status of adolescents
map(list("16 to 20" = 16:20, "18 to 20" = 18:20), function(ages) {
  data_raw %>%
    filter(RELATE == 3 & AGE %in% ages) %>%
    pull(EDUCPR) %>%
    {c(sum(. < 410), sum(. == 410), sum(. > 410))/length(.)}
})

# Analyzing years of education in children:
data_raw %>%
  filter(AGE < 21) %>%
  select(AGE, YRSCHOOL, EDUCPR) %>%
  pivot_longer(-AGE) %>%
  ggplot(aes(AGE, value)) +
  geom_count() +
  facet_wrap(vars(name), scales = "free_y")

# Analyzing child labor:
data %>%
  rowwise() %>%
  mutate(
    HasWorkersSub18 = map2_lgl(HrsUsual[[1]] != 0, ChildAges[[1]] < 18, `&`),
    HasWorkersSub16 = map2_lgl(HrsUsual[[1]] != 0, ChildAges[[1]] < 16, `&`)
  ) %>%
  select(starts_with("HasWorkers")) %>%
  map(table)



# Task 1 ------------------------------------------------------------------
# -------------------- Models --------------------
model1 <- list()

model1$formulas <- list(
  "Linear" = Siblings ~ IncTot,
  "Linear + bins" = Siblings ~ IncTot + IncTotCut,
  "Linear * bins" = Siblings ~ IncTot*IncTotCut,
  "Quadratic" = Siblings ~ IncTot + I(IncTot^2),
  "Log" = Siblings ~ IncTot + IncTotLog,
  "Quantile" = Siblings ~ IncTotQuant,
  "Quantile Q." = Siblings ~ IncTotQuant + I(IncTotQuant^2)
)

model1$models <- map(model1$formulas, ~lm(.x, data))

model1$results <- imap_dfr(model1$models, function(m, name) {
  tibble(model = name, fit = fitted(m), res = residuals(m)) %>%
    mutate(
      value = if (name %in% c("Quantile", "Quantile Q.")) {
        quantile(data$IncTot, m$model[,"IncTotQuant"])
      } else if (name %in% c("Log")) {
        10^m$model[,"IncTotLog"] - 1
      } else {
        m$model[,"IncTot"]
      },
      model_cat = if (grepl("Linear", name)) "Linear" else "Non-Linear",
      model = factor(model,
        levels = c(names(model1$formulas)[1:3], " ", names(model1$formulas)[4:8])
      )
    )
})


# -------------------- Graphs --------------------
graph1 <- data_gdensity(IncTot)

ggplot(graph1$data_density, aes(value, Siblings)) +
  plot_gdensity(graph1) +
  geom_line(aes(y = fit, color = model), model1$results,
    linewidth = 1, alpha = 0.7
  ) +
  facet_wrap(vars(model_cat), nrow = 1) +
  scale_color_manual(values = c(pal[1:3], "white", pal[4:7]), drop = FALSE) +
  labs(
    title = "Siblings versus Total Income",
    x = "Total Income",
    y = "Siblings",
    color = "Model"
  )

gaze(model1$models)


# Task 2 ------------------------------------------------------------------
# -------------------- Graphs --------------------
graph2 <- data_gdensity(RoomsCapta, YearsEduc, StudentsCapta, WorkedCapta)

graph2$data_hasworkers <- data %>%
  select(c(Siblings, HasWorkers)) %>%
  mutate(HasWorkers = as.integer(HasWorkers), var = "HasWorkers")

graph2$graph_cont <- ggplot(graph2$data_density, aes(Siblings, value)) +
  plot_gdensity(graph2) +
  plot_smooths(graph2$data_true) +
  facet_wrap(vars(var), scales = "free_y") +
  scale_color_manual(values = pal) +
  labs(y = "Value", color = "Models")

graph2$graph_disc <- ggplot(graph2$data_hasworkers, aes(Siblings, HasWorkers)) +
  geom_count(color = pal[7]) +
  plot_smooths(graph2$data_hasworkers) +
  facet_wrap(vars(var)) +
  scale_size_continuous(breaks = c(100, 1000)) +
  scale_color_manual(values = pal) +
  labs(x = "Siblings", y = "", size = "Count/Density", color = "Models")

graph2$layout <- "AAAA#\nAAAAB\nAAAAB\nAAAAB\nAAAA#"

graph2$graph_cont + graph2$graph_disc +
  plot_layout(widths = c(4,1), guides = "collect", design = graph2$layout) +
  plot_annotation(title = "Quality Measures versus Siblings")


# -------------------- Models --------------------
model2 <- list()

model2$formulas <- expand.grid(
  c("RoomsCapta", "YearsEduc", "StudentsCapta", "HasWorkers", "WorkedCapta"),
  "~",
  c("Siblings", "as.factor(Siblings)")
) %>%
  apply(1, \(x) as.formula(paste(x, collapse = " ")))

model2$models <- map(model2$formulas, ~lm(.x, data))

gaze(model2$models)
