library(ggplot2)
library(dplyr)
library(lubridate)

data = read_csv("~/Rec_Class/data/crude_daily_10Y.csv")

# Calculate monthly maximum per year
monthly_max <- data %>%
  mutate(
    date = mdy(date),
    year  = year(date),
    month = month(date, label = TRUE, abbr = TRUE) # Returns "Jan", "Feb", etc.
  ) %>%
  group_by(year, month) %>%
  summarise(max_data = max(price, na.rm = TRUE), .groups = "drop") %>%
  # Replace Inf/-Inf with NA if a month has only missing data
  mutate(max_data = ifelse(is.infinite(max_data), NA, max_data))


x <- monthly_max %>%  pull(max_data)
#x=x/max(x)

# Prepare record dataset matching the new Date axis format
# Assumes L contains matching Date objects (e.g., "2016-05-01")
# If L contains strings, convert it using: L <- ymd(L) or mdy(L)
L <- rec_times(x)
R <- rec_values(x)

data_rec <- data.frame(
  date_axis = L,
  record_value = R,
  time = length(x),
  year = monthly_max$year[L],
  month =monthly_max$month[L]
) %>%mutate(
  date_axis = ymd(paste(year, as.character(month), "01"))
)
data_rec

# 1. Prepare historical time series data by creating a continuous Date column
plot_data <- monthly_max %>%
  mutate(
    # Assumes 'month' is an ordered factor or string like "Apr", "May"
    # Adds Day 1 to create a valid R date object (e.g., "2016-04-01")
    date_axis = ymd(paste(year, as.character(month), "01"))
  )

# 3. Generate the Extremes Template Plot
ggplot() +
  # Baseline historical monthly time series line
  geom_line(
    data = plot_data,
    aes(x = date_axis, y = max_data),
    color = "#7f8c8d",
    linewidth = 0.75
  ) +
  # Baseline historical data points
  geom_point(
    data = plot_data,
    aes(x = date_axis, y = max_data),
    color = "#bdc3c7",
    size = 1.5
  ) +
  # Drop-down vertical reference lines for records down to the timeline
  geom_segment(
    data = data_rec,
    aes(x = date_axis, xend = date_axis, y = 0, yend = record_value),
    color = "#0072B2",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  # Highlighted Record points (Solid center)
  geom_point(
    data = data_rec,
    aes(x = date_axis, y = record_value),
    color = "#0072B2",
    size = 3.5,
    shape = 19
  ) +
  # Visual target ring around the records
  geom_point(
    data = data_rec,
    aes(x = date_axis, y = record_value),
    color = "#007",
    size = 5.0,
    shape = 1
  ) +
  # Labels, Titles, and Scale Formatters
  labs(
    title ="",# "Historical Monthly Crude Oil Price Peaks & Record Highs",
    subtitle = "", # "Grey series tracking monthly maximums; Red targets highlight verified record breaks.",
    x = "Timeline (Year-Month)",
    y = "Maximum Crude Oil Price ($)",
    caption = "" #"Source: Extremes Analysis Template Execution"
  ) +
  # Dynamically formats x-axis to show years while maintaining monthly spacing underneath
  scale_x_date(
    date_labels = "%Y",
    date_breaks = "2 years"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text( color = "black", size = 12),
    plot.subtitle = element_text(color = "#7f8c8d", margin = margin(b = 15)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f0f0f0"),
    axis.title = element_text( color = "black")
  )
