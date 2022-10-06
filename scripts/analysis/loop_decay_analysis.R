library(ggplot2)


# load processed data and convert to dataframes ---------------------------

normalized1_gained <- readRDS("data/processed/hic/gained_normalized_mh_index.rds")
## convert to dataframe
Mean_gained <- as.vector(colMeans(normalized1_gained))
df_gained <- data.frame(Group_gained,Mean_gained)
df_gained

normalized1_nullSet <- readRDS("data/processed/hic/nullSet_cont_normalized_mh_index.rds")
## convert to dataframe
Mean_nullSet <- as.vector(colMeans(normalized1_nullSet))
df_nullSet <- data.frame(Group_nullSet,Mean_nullSet)
df_nullSet

normalized1_control <- readRDS("data/processed/hic/cont_normalized_unfiltered_mh_index.rds")
## convert to dataframe
Mean_control <- as.vector(colMeans(normalized1_control))
df_control <- data.frame(Group_control,Mean_control)
df_control

normalized1_lost <- readRDS("data/processed/hic/lost_normalized_mh_index.rds")
## convert to dataframe
Mean_lost<- as.vector(colMeans(normalized1_lost))
df_lost <- data.frame(Group_lost, Mean_lost)
df_lost


# final dataframe for plotting --------------------------------------------

df <- cbind(df_gained, df_nullSet) |> 
  dplyr::select(Mean_nullSet, Mean_gained)

df$Group <- factor(0:10)


# visualization -----------------------------------------------------------

## might want to add 'n' of loops per point

loop_decay_plot <- ggplot(df,
       aes(x = Group, group = 1)) +
  geom_point(aes(y = Mean_gained),
             color="black") +
  geom_line(aes(y = Mean_gained),
            color="black") +
  # geom_point(aes(y = Mean_lost),
  #           color= "red") +
  # geom_line(aes(y = Mean_lost),
  #             color= "red") +
  geom_point(aes(y = Mean_nullSet),
             color= "blue") +
  geom_line(aes(y = Mean_nullSet),
            color= "blue") +
  xlab("Distance from the center (kb)") + 
  ylab("Signal relative to the center") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(
    axis.text.y=element_text(size=8),
    axis.title.x = element_text(color="steelblue", size=12, face="bold"),
    axis.title.y = element_text(color="steelblue", size=12, face="bold"))

ggsave("plots/loop_decay_plot.pdf", loop_decay_plot, device = "pdf")
