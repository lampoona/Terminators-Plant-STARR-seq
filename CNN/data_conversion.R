library(tidyverse)


### prepare data for CNN modelling
load('data_seq.Rdata')
load('data_seq_proto.Rdata')

set.seed(928)

tobacco <- data_seq %>%
  as_tibble() %>%
  select(id, 'species' = Species, 'sequence' = seq, 'enrichment_tobacco' = mean_enrichment)

maize <- data_seq_proto %>%
  ungroup() %>%
  select(id, 'species' = Species, 'sequence' = seq, 'enrichment_maize' = mean_enrichment)

CNN_data <- inner_join(tobacco, maize) %>%
  filter(nchar(sequence) == 170) %>%
  mutate(
    set = sample(c('train', 'test'), n(), replace = TRUE, prob = c(0.9, 0.1))
  ) %>%
  write_tsv('terminator_data.tsv')


### results
predictions <- read_tsv('terminator_data_predictions.tsv') %>%
  rename_with(
    .cols = ends_with('predictions'),
    .fn = ~ str_replace(.x, 'enrichment_(.*)_predictions', 'prediction_\\1')
  ) %>%
  pivot_longer(
    starts_with(c('enrichment', 'prediction')),
    names_to = c('.value', 'system'),
    names_pattern = '(.*)_(.*)'
  )

correlation <- predictions %>%
  group_by(system) %>%
  summarise(
    n = paste0('~~n == "', prettyNum(n(), big.mark = ','), '"'),
    rsquare = sprintf('~~R^2 == "%.2f"', cor(enrichment, prediction)^2),
    spearman = sprintf('~~rho == "%.2f"', cor(enrichment, prediction, method = 'spearman'))
  ) %>%
  ungroup() %>%
  pivot_longer(
    c(n, spearman, rsquare),
    names_to = 'method',
    values_to = 'label'
  ) %>%
  group_by(system) %>%
  mutate(
    vjust = 1.33 * seq_len(n())
  ) %>%
  ungroup()

ggplot(predictions, aes(x = enrichment, y = prediction)) +
  facet_grid(
    cols = vars(system)
  ) +
  geom_hex(
    bins = 50
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'red',
    linetype = 'dashed'
  ) +
  geom_text(
    data = correlation,
    aes(label = label, vjust = vjust),
    x = -Inf,
    y = Inf,
    hjust = 0,
    size = 14/.pt,
    parse = TRUE
  ) +
  scale_x_continuous(
    name = expression(log[2]('enrichment, measured'))
  ) +
  scale_y_continuous(
    name = expression(log[2]('enrichment, predicted'))
  ) +
  scale_fill_viridis_c(
    trans = 'log2'
  ) +
  theme_bw(
    base_size = 18
  )


### check evolution predictions
evo_data <- read_tsv('evolution_data.tsv')

data_long <- evo_data %>%
  pivot_longer(
    cols = starts_with('prediction'),
    values_to = 'enrichment',
    names_to = 'model',
    names_prefix = 'prediction_'
  ) %>%
  mutate(
    opt_for = if_else(opt_for == 'start', 'tobacco/maize/both', opt_for)
  ) %>%
  separate_rows(
    opt_for,
    sep = '/'
  )

ggplot(data_long, aes(x = round, y = enrichment, color = origin)) +
  facet_grid(
    cols = vars(opt_for),
    rows = vars(model),
    labeller = label_both
  ) +
  geom_line() +
  scale_x_continuous(
    name = 'rounds of in silico evolution',
    breaks = seq(0, 10, 2)
  ) +
  theme_bw(
    base_size = 18
  ) +
  theme(
    legend.position = 'none'
  )

ggplot(data_long %>% filter(round %in% c(0, 3, 10)), aes(x = as.factor(round), y = enrichment, fill = as.factor(round))) +
  facet_grid(
    cols = vars(opt_for),
    rows = vars(model),
    labeller = label_both
  ) +
  geom_violin(
    alpha = 0.5
  ) +
  geom_boxplot(
    fill = 'white',
    width = 0.2,
    outlier.shape = NA
  ) +
  scale_x_discrete(
    name = 'rounds of in silico evolution'
  ) +
  theme_bw(
    base_size = 18
  ) +
  theme(
    legend.position = 'none'
  )
