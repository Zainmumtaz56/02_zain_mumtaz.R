## Assignment 2: ggplot2 precipitation data
## Author: Jan Ahmad Jamal
## Put this script in the same folder as prec_data.rds and run from top to bottom.

## load packages
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

## ------------------------------------------------------------
## Part 1: Load data
## ------------------------------------------------------------

## load data
dta <- readRDS(file = "prec_data.rds")

## convert to data.table
dta_t <- as.data.table(x = dta)

## ------------------------------------------------------------
## Part 2: Basic inspection for plotting
## ------------------------------------------------------------

## basic inspection
class(x = dta)
class(x = dta_t)
head(x = dta_t)
names(x = dta_t)
dim(x = dta_t)

## required columns check
req_cols <- c("DT", "STATION", "ELEMENT", "VALUE")
miss_cols <- setdiff(x = req_cols, y = names(x = dta_t))
if (length(x = miss_cols) > 0) {
  stop(paste("Missing required columns:", paste(miss_cols, collapse = ", ")))
}

## create optional columns if the file does not contain them
if (!("QUALITY" %in% names(x = dta_t))) {
  dta_t[, QUALITY := NA_character_]
}
if (!("FLAG" %in% names(x = dta_t))) {
  dta_t[, FLAG := NA_character_]
}

## make sure DT is usable as time
if (!inherits(x = dta_t$DT, what = c("Date", "POSIXct", "POSIXt", "IDate", "IDateTime"))) {
  dta_t[, DT := as.POSIXct(x = DT)]
}

## inspection values
range(x = dta_t$DT, na.rm = TRUE)
uniqueN(x = dta_t$STATION)
uniqueN(x = dta_t$ELEMENT)
uniqueN(x = dta_t$QUALITY)
uniqueN(x = dta_t$FLAG)
summary(object = dta_t$VALUE)

## Interpretation for plotting:
## DT is the natural x variable because the main question is temporal change.
## VALUE is the natural y variable because it stores the measured precipitation value.
## STATION and ELEMENT are useful for colour, group, and facet because they are categories.
## QUALITY and FLAG are better used as quality indicators, not as main axes.
## STATION should not be mapped directly to colour for all stations because too many colours
## would make the plot unreadable. Faceting or selection is more appropriate.

## ------------------------------------------------------------
## Part 3: Data preparation for visualisation
## ------------------------------------------------------------

## helper variables created with :=
dta_t[, year := as.integer(x = format(x = DT, format = "%Y"))]
dta_t[, month := as.integer(x = format(x = DT, format = "%m"))]
dta_t[, date := as.IDate(x = DT)]
dta_t[, month_start := as.IDate(x = paste0(year, "-", sprintf("%02d", month), "-01"))]
dta_t[, pos := VALUE > 0]
dta_t[, zero := VALUE == 0]
dta_t[, quality_missing := is.na(x = QUALITY) | QUALITY == ""]
dta_t[, flag_present := !(is.na(x = FLAG) | FLAG == "")]

## station-element threshold for unusually high positive values
thr_dt <- dta_t[!is.na(x = VALUE), .(
  q99 = as.numeric(x = quantile(x = VALUE, probs = 0.99, na.rm = TRUE)),
  iqr_v = IQR(x = VALUE, na.rm = TRUE)
), by = .(STATION, ELEMENT)]

thr_dt[, high_thr := q99 + 3 * iqr_v]
dta_t <- merge(
  x = dta_t,
  y = thr_dt[, .(STATION, ELEMENT, high_thr)],
  by = c("STATION", "ELEMENT"),
  all.x = TRUE,
  sort = FALSE
)

## suspicious rule: negative values, extremely high values, or quality/flag problems
## This rule is explicit, reproducible, and conservative.
dta_t[, suspicious_value := VALUE < 0 | VALUE > high_thr | quality_missing | flag_present]
dta_t[is.na(x = suspicious_value), suspicious_value := FALSE]

## event rank within station and element
dta_t[, event_rank := frank(x = -VALUE, ties.method = "dense", na.last = "keep"), by = .(STATION, ELEMENT)]

## long zero runs
setorder(x = dta_t, STATION, ELEMENT, DT)
dta_t[, zero_run_id := rleid(zero), by = .(STATION, ELEMENT)]
dta_t[, zero_run_len := .N, by = .(STATION, ELEMENT, zero_run_id)]
dta_t[, long_zero_run := zero & zero_run_len >= 10]

## Why these variables help:
## year, month, date and month_start make temporal aggregation and faceting easier.
## pos and zero help with zero inflation in precipitation data.
## quality_missing and flag_present turn text quality columns into clear plotting indicators.
## suspicious_value is used to reveal unlikely or lower-quality observations.
## event_rank helps isolate extremes without manually choosing rows.
## long_zero_run helps find long repeated dry periods or possible recording gaps.

## ------------------------------------------------------------
## Part 4: Station selection strategy
## ------------------------------------------------------------

## station summary for reproducible selection
station_stat <- dta_t[, .(
  n_obs = .N,
  n_nonzero = sum(x = VALUE > 0, na.rm = TRUE),
  total_value = sum(x = VALUE, na.rm = TRUE),
  n_flag = sum(x = flag_present, na.rm = TRUE),
  n_missing_quality = sum(x = quality_missing, na.rm = TRUE),
  n_suspicious = sum(x = suspicious_value, na.rm = TRUE),
  n_long_zero = sum(x = long_zero_run, na.rm = TRUE),
  max_value = max(x = VALUE, na.rm = TRUE)
), by = STATION]

station_stat[, activity_score := frank(x = -n_nonzero, ties.method = "first")]
station_stat[, problem_score := frank(x = -(n_flag + n_missing_quality + n_suspicious + n_long_zero), ties.method = "first")]

## choose 3 active stations and 3 problematic stations, then ensure exactly 6 stations
st_active <- station_stat[order(activity_score)][1:3, STATION]
st_problem <- station_stat[order(problem_score)][1:3, STATION]
st_sel <- unique(x = c(st_active, st_problem))

if (length(x = st_sel) < 6) {
  st_extra <- station_stat[!(STATION %in% st_sel)][order(activity_score)][1:(6 - length(x = st_sel)), STATION]
  st_sel <- c(st_sel, st_extra)
}

st_sel <- st_sel[1:6]
st_sel

## Selection explanation:
## I select three stations with many non-zero observations and three stations with many
## quality or suspicious-data problems. This gives a visual comparison between active
## precipitation stations and stations that may have data-quality problems.
## What is lost: stations outside st_sel are not shown in detail, so local patterns from
## unselected stations may be missed. The selection improves readability but reduces completeness.

## ------------------------------------------------------------
## Part 8 helper: Element selection object
## ------------------------------------------------------------

el_stat <- dta_t[, .(
  n_obs = .N,
  n_nonzero = sum(x = VALUE > 0, na.rm = TRUE),
  total_value = sum(x = VALUE, na.rm = TRUE),
  n_suspicious = sum(x = suspicious_value, na.rm = TRUE),
  max_value = max(x = VALUE, na.rm = TRUE)
), by = ELEMENT]

## choose elements with most observations because these are most comparable visually
el_sel <- el_stat[order(-n_obs, -n_nonzero)][1:min(4, .N), ELEMENT]
el_sel

## ------------------------------------------------------------
## Part 5: Temporal plots
## ------------------------------------------------------------

## Plot A: raw or near-raw time series
p1 <- ggplot(
  data = dta_t[STATION %in% st_sel & ELEMENT %in% el_sel[1]],
  mapping = aes(x = DT, y = VALUE, group = STATION)
) +
  geom_line(alpha = 0.45, linewidth = 0.25) +
  facet_wrap(facets = vars(STATION), ncol = 2) +
  scale_x_datetime(date_labels = "%Y", date_breaks = "1 year") +
  labs(
    title = "Raw precipitation time series for selected stations",
    subtitle = paste("Element:", el_sel[1], "| zeros kept visible because dry periods are part of the pattern"),
    x = "Time",
    y = "VALUE"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

## Design note:
## A standard line plot with all stations in one panel would be misleading because
## overplotting would hide station-specific behaviour. Faceting keeps each station readable.
## Zeros stay visible because removing them would exaggerate rainfall frequency.

## Plot B: aggregated monthly time series
monthly_dt <- dta_t[STATION %in% st_sel & ELEMENT %in% el_sel[1], .(
  value_month = sum(x = VALUE, na.rm = TRUE),
  n_obs = .N,
  n_suspicious = sum(x = suspicious_value, na.rm = TRUE)
), by = .(STATION, ELEMENT, month_start)]

p2 <- ggplot(
  data = monthly_dt,
  mapping = aes(x = month_start, y = value_month, group = STATION)
) +
  geom_line(linewidth = 0.35) +
  geom_point(size = 0.8, alpha = 0.65) +
  facet_wrap(facets = vars(STATION), ncol = 2) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(
    title = "Monthly aggregated precipitation",
    subtitle = "Monthly sums reduce timestamp noise and make seasonal or large-event patterns clearer",
    x = "Month",
    y = "Monthly VALUE"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2

## Aggregation explanation:
## Month is appropriate because precipitation can be very noisy at the original timestamp level.
## Monthly aggregation makes seasonal and station-level differences easier to see.
## Lost information: exact timing and short extreme bursts are partly hidden.

## Plot C: extreme-event view with highlighting and labels
extreme_dt <- dta_t[
  STATION %in% st_sel & ELEMENT %in% el_sel[1] & !is.na(x = VALUE)
][order(STATION, -VALUE), .SD[1:min(5, .N)], by = STATION]

p3 <- ggplot(
  data = dta_t[STATION %in% st_sel & ELEMENT %in% el_sel[1]],
  mapping = aes(x = DT, y = VALUE)
) +
  geom_line(alpha = 0.25, linewidth = 0.25) +
  geom_point(
    data = extreme_dt,
    mapping = aes(x = DT, y = VALUE),
    colour = "red",
    size = 2
  ) +
  geom_text(
    data = extreme_dt[event_rank <= 2],
    mapping = aes(x = DT, y = VALUE, label = round(x = VALUE, digits = 1)),
    vjust = -0.6,
    size = 2.8,
    colour = "red"
  ) +
  facet_wrap(facets = vars(STATION), ncol = 2) +
  scale_x_datetime(date_labels = "%Y", date_breaks = "1 year") +
  labs(
    title = "Extreme precipitation events highlighted",
    subtitle = "Red points show the five largest events in each selected station; labels mark the two largest",
    x = "Time",
    y = "VALUE"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3

## Extreme-event explanation:
## The deliberate highlighting step makes extremes visible without losing the background time series.
## A default plot would make it hard to identify the largest observations precisely.

## ------------------------------------------------------------
## Part 6: Distribution plots
## ------------------------------------------------------------

## Distribution across stations, positive values only, transformed scale
p4 <- ggplot(
  data = dta_t[STATION %in% st_sel & ELEMENT %in% el_sel[1] & VALUE > 0],
  mapping = aes(x = reorder(x = STATION, X = VALUE, FUN = median), y = VALUE)
) +
  geom_boxplot(outlier.alpha = 0.35) +
  scale_y_log10(labels = label_number()) +
  coord_flip() +
  labs(
    title = "Positive precipitation distribution by station",
    subtitle = "Zero values removed here because log scale cannot show zeros; this focuses on wet observations",
    x = "Station ordered by median positive value",
    y = "VALUE, log10 scale"
  ) +
  theme_bw()

p4

## Boxplot warning:
## A boxplot can be misleading here because many precipitation observations are zero.
## If zeros were included, the median and lower quartiles would often be zero and the wet-day
## distribution would be hidden. Therefore this plot is only about positive observations.

## Distribution across elements using ECDF on log1p values
p5 <- ggplot(
  data = dta_t[STATION %in% st_sel & ELEMENT %in% el_sel],
  mapping = aes(x = log1p(x = pmax(VALUE, 0)), colour = ELEMENT)
) +
  stat_ecdf(linewidth = 0.7) +
  labs(
    title = "Distribution comparison across elements",
    subtitle = "log1p keeps zero values while reducing extreme skewness",
    x = "log1p(VALUE)",
    y = "Cumulative probability",
    colour = "Element"
  ) +
  theme_bw()

p5

## Scale explanation:
## log1p is useful because precipitation data are zero-inflated and right-skewed.
## A normal histogram on the original scale would mostly show a high bar near zero
## and hide the upper tail.

## ------------------------------------------------------------
## Part 7: Quality and suspicious-data plots
## ------------------------------------------------------------

## monthly quality indicators
quality_month <- dta_t[STATION %in% st_sel & ELEMENT %in% el_sel, .(
  missing_quality_rate = mean(x = quality_missing, na.rm = TRUE),
  flag_rate = mean(x = flag_present, na.rm = TRUE),
  suspicious_rate = mean(x = suspicious_value, na.rm = TRUE)
), by = .(STATION, ELEMENT, month_start)]

quality_long <- melt(
  data = quality_month,
  id.vars = c("STATION", "ELEMENT", "month_start"),
  measure.vars = c("missing_quality_rate", "flag_rate", "suspicious_rate"),
  variable.name = "quality_type",
  value.name = "rate"
)

p6 <- ggplot(
  data = quality_long[ELEMENT %in% el_sel[1]],
  mapping = aes(x = month_start, y = rate, fill = quality_type)
) +
  geom_col(position = "identity", alpha = 0.65) +
  facet_grid(rows = vars(quality_type), cols = vars(STATION), scales = "free_y") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Temporal pattern of quality problems",
    subtitle = "This combines time, quality indicators, and station comparison",
    x = "Month",
    y = "Share of observations",
    fill = "Quality indicator"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

p6

## This plot investigates missing QUALITY, values in FLAG, and suspicious rows over time.
## A single total count table would not reveal whether problems are clustered in specific months.

## long zero run and station-element anomaly plot
zero_anom <- dta_t[STATION %in% st_sel & ELEMENT %in% el_sel, .(
  long_zero_rows = sum(x = long_zero_run, na.rm = TRUE),
  suspicious_rows = sum(x = suspicious_value, na.rm = TRUE),
  repeated_zero_share = mean(x = long_zero_run, na.rm = TRUE)
), by = .(STATION, ELEMENT)]

p7 <- ggplot(
  data = zero_anom,
  mapping = aes(x = ELEMENT, y = STATION, fill = repeated_zero_share)
) +
  geom_tile(colour = "white") +
  geom_text(mapping = aes(label = long_zero_rows), size = 3) +
  scale_fill_gradient(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Long zero runs by station and element",
    subtitle = "Tile colour shows share of observations in long zero runs; text shows row counts",
    x = "Element",
    y = "Station",
    fill = "Long-zero share"
  ) +
  theme_bw()

p7

## This plot investigates suspiciously repeated zero values, station-specific anomalies,
## and element-specific differences. Long zero runs can be real dry periods, but they can
## also indicate missing recording or constant sensor output.

## ------------------------------------------------------------
## Part 8: Element comparison plots
## ------------------------------------------------------------

## Element comparison 1: same station across elements, with pseudo-log scale
monthly_el <- dta_t[STATION %in% st_sel & ELEMENT %in% el_sel, .(
  value_month = sum(x = VALUE, na.rm = TRUE),
  suspicious_rate = mean(x = suspicious_value, na.rm = TRUE)
), by = .(STATION, ELEMENT, month_start)]

p8 <- ggplot(
  data = monthly_el,
  mapping = aes(x = month_start, y = value_month, group = ELEMENT, colour = ELEMENT)
) +
  geom_line(linewidth = 0.35) +
  facet_grid(rows = vars(ELEMENT), cols = vars(STATION), scales = "free_y") +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), labels = label_number()) +
  labs(
    title = "Monthly behaviour across elements and stations",
    subtitle = "Pseudo-log scale keeps zeros but compresses extreme values",
    x = "Month",
    y = "Monthly VALUE, pseudo-log scale",
    colour = "Element"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

p8

## Element comparison 2: reliability by element
el_quality <- dta_t[ELEMENT %in% el_sel, .(
  missing_quality_rate = mean(x = quality_missing, na.rm = TRUE),
  flag_rate = mean(x = flag_present, na.rm = TRUE),
  suspicious_rate = mean(x = suspicious_value, na.rm = TRUE),
  positive_rate = mean(x = VALUE > 0, na.rm = TRUE)
), by = ELEMENT]

el_quality_long <- melt(
  data = el_quality,
  id.vars = "ELEMENT",
  measure.vars = c("missing_quality_rate", "flag_rate", "suspicious_rate", "positive_rate"),
  variable.name = "metric",
  value.name = "rate"
)

p9 <- ggplot(
  data = el_quality_long,
  mapping = aes(x = reorder(x = ELEMENT, X = rate, FUN = mean), y = rate)
) +
  geom_col() +
  facet_wrap(facets = vars(metric), scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  coord_flip() +
  labs(
    title = "Element reliability and measurement comparison",
    subtitle = "Lower missing, flag, and suspicious rates suggest more reliable elements",
    x = "Element",
    y = "Rate"
  ) +
  theme_bw()

p9

## Element interpretation:
## Elements may be on different scales, so free y-scales and transformed scales are useful.
## The most trustworthy element is the one with a low missing-quality rate, low flag rate,
## low suspicious rate, and stable behaviour across stations. The least trustworthy element
## is the one with concentrated suspicious values or many long zero runs.

## ------------------------------------------------------------
## Part 9: Grammar of Graphics explanation for four plots
## ------------------------------------------------------------

## p1 explanation:
## data: dta_t filtered to st_sel and one selected element.
## mapping: DT is mapped to x and VALUE to y; station is separated by facets.
## geom: geom_line shows temporal development.
## stat: default identity statistic is used.
## facet: facet_wrap by STATION reduces overplotting.
## scale: x-axis uses readable yearly date breaks.
## theme: theme_bw and rotated x labels improve readability.

## p2 explanation:
## data: monthly_dt, an aggregated data.table.
## mapping: month_start to x and monthly summed VALUE to y.
## geom: geom_line and geom_point show trend and individual months.
## stat: default identity statistic is used because aggregation was computed manually.
## facet: facet_wrap by STATION makes station comparisons readable.
## scale: x-axis is formatted as years.
## theme: simple theme and rotated labels reduce visual noise.

## p4 explanation:
## data: positive observations only from dta_t.
## mapping: station to x and VALUE to y.
## geom: geom_boxplot compares distributions.
## stat: boxplot statistic is computed by ggplot2.
## facet: no facet because the station comparison is already the main axis.
## scale: log10 y-scale handles skewness but requires removing zeros.
## theme: coord_flip and theme_bw make station names readable.

## p6 explanation:
## data: quality_long, a melted monthly quality table.
## mapping: month_start to x, quality rate to y, quality_type to fill.
## geom: geom_col shows monthly rates.
## stat: default identity statistic is used.
## facet: facet_grid separates quality type and station.
## scale: percent y-axis makes rates easier to interpret.
## theme: bottom legend and rotated x-axis improve readability.

## ------------------------------------------------------------
## Part 10: Faceting and scale decisions
## ------------------------------------------------------------

## Decision 1: facet_wrap versus one-colour-per-station.
## I use facet_wrap for station time series because six stations already create a busy plot.
## Colour alone would cause overlapping lines and a difficult legend.

## Decision 2: fixed versus free scales.
## In p1 and p2 I keep the same basic station structure for comparison, but in p6 and p8
## I use free scales because quality rates and element values can be very different.
## Fixed scales were considered for p8, but rejected because small elements became invisible.

## Decision 3: transformed versus untransformed axis.
## p4 uses log10 for positive observations. p8 uses pseudo-log because zeros must stay visible.
## A normal untransformed axis would be dominated by a few large events.

## Decision 4: alpha transparency for overplotting.
## In p1 and p3 I use alpha so the background time series does not dominate highlighted extremes.

## Decision 5: ordering factor levels.
## In p4 I reorder stations by median positive value so the visual order has meaning.

## ------------------------------------------------------------
## Part 11: Reshape task
## ------------------------------------------------------------

## reshape station summaries into long format for faceted plotting
station_sum <- dta_t[STATION %in% st_sel, .(
  total_value = sum(x = VALUE, na.rm = TRUE),
  positive_rate = mean(x = VALUE > 0, na.rm = TRUE),
  missing_quality_rate = mean(x = quality_missing, na.rm = TRUE),
  flag_rate = mean(x = flag_present, na.rm = TRUE),
  suspicious_rate = mean(x = suspicious_value, na.rm = TRUE),
  long_zero_rate = mean(x = long_zero_run, na.rm = TRUE)
), by = STATION]

dta_m <- melt(
  data = station_sum,
  id.vars = "STATION",
  variable.name = "metric",
  value.name = "value"
)

dta_m

p10 <- ggplot(
  data = dta_m,
  mapping = aes(x = reorder(x = STATION, X = value, FUN = mean), y = value)
) +
  geom_col() +
  facet_wrap(facets = vars(metric), scales = "free_y") +
  coord_flip() +
  labs(
    title = "Reshaped station summary metrics",
    subtitle = "melt() creates a long table so several station metrics can be shown with one grammar",
    x = "Station",
    y = "Metric value"
  ) +
  theme_bw()

p10

## Reshape explanation:
## Reshaping was needed because the station summary originally had several metric columns.
## A long format gives one metric column and one value column, which works naturally with
## facet_wrap. Without melt(), each metric would need a separate plot or repeated code.

## ------------------------------------------------------------
## Part 12: Multi-panel final figure
## ------------------------------------------------------------

p_final <- (p2 + p4) / (p6 + p9) +
  plot_annotation(
    title = "Visual report: precipitation patterns, distributions, quality problems, and element reliability",
    caption = "Notice first whether monthly peaks align with suspicious-quality periods. Remaining uncertainty: some long zero runs may be real weather, not measurement error."
  )

p_final

## Final figure story:
## The figure compares monthly temporal behaviour, positive-value distributions,
## quality/anomaly patterns, and element reliability. The reader should first notice
## whether the strongest precipitation peaks are connected to suspicious quality periods.
## The main uncertainty is that visual evidence can identify suspicious patterns, but it
## cannot prove whether each suspicious point is a real weather event or a recording problem.

## ------------------------------------------------------------
## Part 13: Saving output
## ------------------------------------------------------------

ggsave(
  filename = "fig_station_timeseries.png",
  plot = p1,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = "fig_quality_flags.png",
  plot = p6,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  filename = "fig_final_visual_report.png",
  plot = p_final,
  width = 14,
  height = 10,
  dpi = 300
)

## ------------------------------------------------------------
## Part 14: Final discussion
## ------------------------------------------------------------

## Final discussion:
## The most informative visualisation is p6 because it connects time, station, element,
## and quality indicators. It shows whether quality problems are isolated or clustered.
## The hardest visualisation to design well is p3 because extremes must be highlighted
## without hiding the normal time series.
## The most important readability choices are selecting only six stations, using facets,
## using transformed scales for skewed values, and keeping zeros visible where they matter.
## Suspicious problems that become visible after plotting include quality clusters,
## long zero runs, repeated dry-period records, and extreme values far above the usual range.
## The most unusual station should be identified from p6, p7, and p10 as the station with
## the strongest combination of suspicious rate, flag rate, missing quality, and long-zero rate.
## The most trustworthy element should have low missing-quality, flag, and suspicious rates
## in p9, while the least trustworthy element has the opposite pattern.
## Default ggplot2 choices would be poor in three cases: plotting all stations in one panel,
## using an untransformed y-axis for highly skewed precipitation values, and using a boxplot
## with many zeros without explaining zero inflation.
