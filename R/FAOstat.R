# Install packages if you haven't already:
# install.packages(c("data.table", "dplyr"))

library(data.table)
library(dplyr)

# 1. Read the 500MB dataset efficiently
# data.table::fread is highly optimized for large files
file_path <- "data/Production_Crops_Livestock_E_All_Data_(Normalized).csv"
fao_data <- fread(file_path)

# 2. Define the list of pulse crops based on FAOSTAT classifications
fao_pulses <- c(
  "Beans, dry", "Broad beans and horse beans, dry", "Chick peas, dry", 
  "Cow peas, dry", "Lentils, dry", "Lupins", "Peas, dry", 
  "Pigeon peas, dry", "Other pulses n.e.c.", "Vetches", "Bambara beans, dry"
)

# 3. Filter the data for the last 10 years and the "Production" element
latest_year <- max(fao_data$Year, na.rm = TRUE)

pulses_data <- fao_data %>%
  filter(
    Element == "Production",
    Item %in% fao_pulses,
    Year >= (latest_year - 5),
    `Area Code` < 500 # Remove aggregate regions (e.g., "World", "Americas") so we don't double count.
  ) %>%
  # Adjust this list if you spot other macro-regions in your specific dataset.
  filter(!grepl("World|Total|Africa|Americas|Asia|Europe|Oceania|Union|mainland", Area, ignore.case = TRUE))

# =========================================================================
# PIPELINE A: Average production of pulses & Position of Lentils
# =========================================================================

lentils_global_rank <- pulses_data %>%
  # Group by crop and year, sum global production for that year, then average across 10 years
  group_by(Item, Year) %>%
  summarise(Global_Yearly_Production = sum(Value, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Item) %>%
  summarise(Avg_10yr_Production = mean(Global_Yearly_Production, na.rm = TRUE), .groups = 'drop') %>%
  # Rank them from highest average production to lowest
  arrange(desc(Avg_10yr_Production)) %>%
  mutate(Pulse_Rank = row_number())

print("--- Rank of Lentils among Pulses (Last 10 Years) ---")
print(lentils_global_rank)

# =========================================================================
# PIPELINE B: Position of Canada as a Lentil Producer
# =========================================================================

canada_lentil_rank <- pulses_data %>%
  filter(Item == "Lentils, dry") %>%
  # Average the production for each country over the last 10 years
  group_by(Area) %>%
  summarise(Avg_10yr_Production = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(Avg_10yr_Production)) %>%
  mutate(Producer_Rank = row_number())

print("--- Top 10 Global Lentil Producers ---")
print(head(canada_lentil_rank, 10))

print("--- Canada's Specific Rank ---")
print(canada_lentil_rank %>% filter(Area == "Canada"))

#barplot####
library(ggplot2)

# Filter the data for production > 1,000,000
high_production_pulses <- lentils_global_rank %>%
  filter(Avg_10yr_Production > 1e6) %>%
  mutate(Item = case_when(
    Item == "Beans, dry" ~ "Beans",
    Item == "Broad beans and horse beans, dry" ~ "Fava beans",
    Item == "Chick peas, dry" ~ "Chickpeas",
    Item == "Peas, dry" ~ "Peas",
    Item == "Cow peas, dry" ~ "Cow peas",
    Item == "Pigeon peas, dry" ~ "Pigeon pea",
    Item == "Lentils, dry" ~ "Lentils",
    Item == "Other pulses n.e.c." ~ "Other pulses",
    # Keep all other names exactly as they are
    TRUE ~ Item
  ))
# Create the bar plot
ggplot(high_production_pulses, aes(x = reorder(Item, Avg_10yr_Production), y = Avg_10yr_Production/1e6,fill = Avg_10yr_Production)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flips the axes so the crop names are easy to read
  scale_y_continuous(labels = scales::comma) + # Formats numbers with commas (e.g., 1,000,000)
  scale_fill_gradient(low = "#a1d99b", high = "#005a32", labels = scales::comma) +
  labs(
    title = "Average Global Production of Pulses",
    subtitle = "2019-2024 period",
    x = NULL,
    y = "Average Production (Million Tonnes)",
    caption = "Showing pulses with an average production > 1 Mt"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )

#genetic gain####
# 1. Get the names of the top pulses (> 1,000,000 tonnes) from our previous table
# We use the original names to match against the raw fao_data
top_pulse_items <- lentils_global_rank %>%
  filter(Avg_10yr_Production > 1e6) %>%
  pull(Item)

# 2. Filter the raw data for Yield, the top crops, and the last 30 years
yield_data <- fao_data %>%
  filter(
    Element == "Yield",
    Item %in% top_pulse_items,
    Year >= (latest_year - 50), # Last 30 years inclusive
    `Area Code` < 500, # Remove aggregate regions (e.g., "World", "Americas") so we don't double count.
    Area == "Canada"
  ) %>%
  # Remove aggregate regions as we did before
  filter(!grepl("World|Total|Africa|Americas|Asia|Europe|Oceania|Union|mainland", Area, ignore.case = TRUE)) %>%
  # Apply the exact same renaming logic for consistency in the plot
  mutate(Item = case_when(
    Item == "Beans, dry" ~ "Beans",
    Item == "Broad beans and horse beans, dry" ~ "Fava beans",
    Item == "Chick peas, dry" ~ "Chickpeas",
    Item == "Peas, dry" ~ "Peas",
    Item == "Cow peas, dry" ~ "Cow peas",
    Item == "Pigeon peas, dry" ~ "Pigeon pea",
    Item == "Lentils, dry" ~ "Lentils",
    Item == "Other pulses n.e.c." ~ "Other pulses",
    TRUE ~ Item
  ))

# 3. Calculate the global average yield per year for each crop
global_yield_trend <- yield_data %>%
  group_by(Item, Year) %>%
  # We take the mean across all producing countries for a given year
  summarise(Avg_Yield = mean(Value, na.rm = TRUE), .groups = 'drop')

# 4. Plot the yield over time with trendlines representing genetic/agronomic gain
ggplot(global_yield_trend |> filter(Item %in% c("Lentils", "Peas")), aes(x = Year, y = Avg_Yield, color = Item)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  # Add a linear trendline to clearly highlight the "gain" (slope) over 30 years
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.7) + 
  scale_color_manual(values = c("#85C086","#5F9E6A"))+
  labs(
    title = "Genetic Gain in Peas and Lentils",
    subtitle = "In Canada across 50 years period",
    x = NULL,
    y = "Average Yield (kg/ha)",
    color = "Pulse Type"
  ) +
  scale_x_continuous(breaks = seq(min(global_yield_trend$Year), max(global_yield_trend$Year), by = 5)) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11),
    legend.position = "none",
  )
