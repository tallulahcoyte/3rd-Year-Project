I found that some columns which were meant to be numeric were character, the columns affected in the NWA73_80 data frame were: 'MEAN_PDL', 'SE_PDL', 'AVE_PDS', 'SE_PDS'. In the NWA81_90 data frame the columns that are character instead of numeric are: 'SE_PDL' and 'SE_PDS'. I will ammend this in the following block of code.

```{changing columns from character to numeric}
NWA73_80$MEAN_PDL <- as.numeric(NWA73_80$MEAN_PDL)
NWA73_80$SE_PDL <- as.numeric(NWA73_80$SE_PDL)
NWA73_80$AVE_PDS <- as.numeric(NWA73_80$AVE_PDS)
NWA73_80$SE_PDS <- as.numeric(NWA73_80$SE_PDS)
NWA81_90$SE_PDL <- as.numeric(NWA81_90$SE_PDL)
NWA81_90$SE_PDS <- as.numeric(NWA81_90$SE_PDS)
```

Now I'm going to delete any rows that have missing values, just to keep my data has tidy as possible.

```{deleting missing value rows}
NWA73_80 <- na.omit(NWA73_80)
NWA81_90 <- na.omit(NWA81_90)
```

In the theme of cleaning up and filtering down the data I am going to delete the 'PDSCINAMO' and 'PSCINAMB' columns, as the 'PDSCINAM' column is the string sum of these two. 

```{deleting 'PDSCINAMO' and 'PDSCINAMB' columns}
NWA73_80$PDSCINAMO <- NULL
NWA73_80$PDSCINAMB <- NULL
NWA81_90$PDSCINAMO <- NULL
NWA81_90$PDSCINAMB <- NULL
```

Now I'm going to check for any typos in the character column PDSCINAM.

```{checking for typos}
unique(NWA73_80$PDSCINAM)
unique(NWA81_90$PDSCINAM)
```


Can't see any typos so I'll assume that I can go ahead and look for any numerical errors. I'll start off with going from left to right in the numeric columns, alternating between NWA73_80 and NWA81_90 on each column comparison.

```{plot}
plot(NWA73_80$N_PDL)
```

Going to look at the summaries for the minimum, mean and max predator lengths so I can remove or alter any observations that don't fit the data trend.

```{sum of min max pred length}
summary(NWA73_80$MIN_PDL)
summary(NWA81_90$MIN_PDL)
summary(NWA73_80$MAX_PDL)
summary(NWA81_90$MAX_PDL)
summary(NWA73_80$MEAN_PDL)
summary(NWA81_90$MEAN_PDL)
```

It looks like some values in the minimum column have been multiplied by 10 or 100 looking at the 3rd Quartile and Max value for the 1981-90 data set in comparison to the MAX_PDL, as common sense tells us that each value in the output of the summary function should be greater for the MAX_PDL than the MIN_PDL. In the following block of code I will divide the values in MIN_PDL for the 1981-90 period by a suitable value so it fits in with the dataset.

```{fixing min values in correspondence with relevant mean}
NWA81_90[ 5, 4] <- NWA81_90[ 5, 4]/100
NWA81_90[ 11, 4] <- NWA81_90[ 11, 4]/100
NWA81_90[ 12, 4] <- NWA81_90[ 12, 4]/10
NWA81_90[ 14, 4] <- NWA81_90[ 14, 4]/100
NWA81_90[ 15, 4] <- NWA81_90[ 15, 4]/10
NWA81_90[ 16, 4] <- NWA81_90[ 16, 4]/100
NWA81_90[ 18, 4] <- NWA81_90[ 18, 4]/100
NWA81_90[ 19, 4] <- NWA81_90[ 19, 4]/10
NWA81_90[ 20, 4] <- NWA81_90[ 20, 4]/10
NWA81_90[ 24, 4] <- NWA81_90[ 24, 4]/10
NWA81_90[ 25, 4] <- NWA81_90[ 25, 4]/1000
NWA81_90[ 26, 4] <- NWA81_90[ 26, 4]/100
NWA81_90[ 26, 4] <- NWA81_90[ 26, 4]/10
NWA81_90[ 27, 4] <- NWA81_90[ 27, 4]/100
NWA81_90[ 28, 4] <- NWA81_90[ 28, 4]/100
NWA81_90[ 29, 4] <- NWA81_90[ 29, 4]/1000
NWA81_90[ 30, 4] <- NWA81_90[ 30, 4]/1000
NWA81_90[ 31, 4] <- NWA81_90[ 31, 4]/100
NWA81_90[ 35, 4] <- NWA81_90[ 35, 4]/10
NWA81_90[ 36, 4] <- NWA81_90[ 36, 4]/10
NWA81_90[ 37, 4] <- NWA81_90[ 37, 4]/100
NWA81_90[ 38, 4] <- NWA81_90[ 38, 4]/100
NWA81_90[ 39, 4] <- NWA81_90[ 39, 4]/100
NWA81_90[ 40, 4] <- NWA81_90[ 40, 4]/100
NWA81_90[ 41, 4] <- NWA81_90[ 41, 4]/10
NWA81_90[ 42, 4] <- NWA81_90[ 42, 4]/100
NWA81_90[ 43, 4] <- NWA81_90[ 43, 4]/10
NWA81_90[ 45, 4] <- NWA81_90[ 45, 4]/10
NWA81_90[ 49, 4] <- NWA81_90[ 49, 4]/10
NWA81_90[ 51, 4] <- NWA81_90[ 51, 4]/100
NWA81_90[ 52, 4] <- NWA81_90[ 52, 4]/10
NWA81_90[ 54, 4] <- NWA81_90[ 54, 4]/10
NWA81_90[ 55, 4] <- NWA81_90[ 55, 4]/100
NWA81_90[ 56, 4] <- NWA81_90[ 56, 4]/1000
NWA81_90[ 58, 4] <- NWA81_90[ 58, 4]/10
NWA81_90[ 59, 4] <- NWA81_90[ 59, 4]/100
NWA81_90[ 62, 4] <- NWA81_90[ 62, 4]/100
NWA81_90[ 63, 4] <- NWA81_90[ 63, 4]/100
NWA81_90[ 64, 4] <- NWA81_90[ 64, 4]/10
NWA81_90[ 73, 4] <- NWA81_90[ 73, 4]/10
NWA81_90[ 74, 4] <- NWA81_90[ 74, 4]/100
NWA81_90[ 83, 4] <- NWA81_90[ 83, 4]/10
NWA81_90[ 82, 4] <- NWA81_90[ 82, 4]*10
NWA81_90[ 85, 4] <- NWA81_90[ 85, 4]*10
NWA81_90[ 87, 4] <- NWA81_90[ 87, 4]/100
NWA81_90[ 88, 4] <- NWA81_90[ 88, 4]/100
NWA81_90[ 89, 4] <- NWA81_90[ 89, 4]*10
NWA81_90[ 90, 4] <- NWA81_90[ 90, 4]*10
NWA81_90[ 91, 4] <- NWA81_90[ 91, 4]*10
NWA81_90[ 92, 4] <- NWA81_90[ 92, 4]*10
NWA81_90[ 100, 4] <- NWA81_90[ 100, 4]*10
NWA81_90[ 101, 4] <- NWA81_90[ 101, 4]*10
```

## Visualising the data

Now that the data is sufficiently cleaned up, I want to see the comparison between predator length and stomach samples, I'll start by plotting minimum predator length against average stomach samples from 1973-80.

```{minpred v av stom samples 73-80}
min_pdl_v_ave_pds_73_80 <- ggplot(NWA73_80, aes(x = MIN_PDL , y = AVE_PDS))
min_pdl_v_ave_pds_73_80 + geom_point()+
     geom_smooth(method = "lm")
min_pdl_v_ave_pds_73_80 + geom_point()+
     geom_smooth(method = "loess")
```

From the two plots above we can see that the relationship between minimum predator length and average stomach samples is positively exponentially distributed (especially in the second plot), so as the minimum predator length increases the average number of stomach samples becomes greater. Now I will look at the 1981-90 time period.

```{av stom samples v minpred 81-90}
min_pdl_v_ave_pds_81_90 <- ggplot(NWA81_90, aes(x = MIN_PDL , y = AVE_PDS))
min_pdl_v_ave_pds_81_90 + geom_point()+
     geom_smooth(method = "lm")
min_pdl_v_ave_pds_81_90 + geom_point()+
     geom_smooth(method = "loess")
```

It doesn't seem like there is any correlation between the minimum predator length and the average number of stomach samples in 1981-90, this may be because there were errors in the data input or because this time frame was subject to a certain amount of overfishing.

Next I am going to compare the maximum predator length with the average number of stomach samples in the 1973-80.

```{maxpred v av stom samples 73-80}
max_pdl_v_ave_pds_73_80 <- ggplot(NWA73_80, aes(x = MAX_PDL , y = AVE_PDS))
max_pdl_v_ave_pds_73_80 + geom_point()+
  geom_smooth(method = "lm")
max_pdl_v_ave_pds_73_80 + geom_point()+
  geom_smooth(method = "loess")
```

The correlation between maximum predator length and average stomach samples from 1973-80 isn't very strong however we can see from the plot on the left that it is somewhat positively correlated due to the positive gradient on the line of best fit.

As I have looked at the 1973-80 dataset I will now observe the 1981-90 set.

```{av stom samples v maxpred 81-90}

max_pdl_v_ave_pds_81_90 <- ggplot(NWA81_90, aes(x = MAX_PDL , y = AVE_PDS))
max_pdl_v_ave_pds_81_90 + geom_point()+
     geom_smooth(method = "lm")
max_pdl_v_ave_pds_81_90 + geom_point()+
     geom_smooth(method = "loess")
```

In comparison to the former time frame we can see that the maximum predator length in relation to the average number of stomach samples is positively exponentially distributed, which is evident in the plot on the right hand side. I would assume that the reason for the lack of correlation in the prior time frame in comparison to the 1981-90 frame is that the bigger predators consumed more of the smaller predators, which would also explain why there is a lack of correlation between the minimum predator length and average stomach samples in 1981-90.

In theme with the previous two observations I will now compare the mean predator length with average stomach samples, starting with the 1973-80 period.

```{meanpred v av stom samples 73-80}
mean_pdl_v_ave_pds_73_80 <- ggplot(NWA73_80, aes(x = MEAN_PDL , y = AVE_PDS))
mean_pdl_v_ave_pds_73_80 + geom_point()+
     geom_smooth(method = "lm")
mean_pdl_v_ave_pds_73_80 + geom_point()+
     geom_smooth(method = "loess")
```

The two variables are somewhat positively correlated, as seen more prominently in the left plot, so we can gather from this that the general trend of the data is that as the predator length increases, so do the average number of stomach samples. This is expected as the maximum predator length has little correlation to stomach samples from 1973-90 however, there is a definite positive correlation between the minimum predator length and stomach samples in the same period.

We will now look at the relation between mean predator length and the average number of stomach samples in the 1981-90 time frame.

```{mean pred v av stom samples 81-90}

mean_pdl_v_ave_pds_81_90 <- ggplot(NWA81_90, aes(x = MEAN_PDL , y = AVE_PDS))
mean_pdl_v_ave_pds_81_90 + geom_point()+
     geom_smooth(method = "lm")
mean_pdl_v_ave_pds_81_90 + geom_point()+
     geom_smooth(method = "loess")
```

The relation between mean predator length and average stomach samples from 1981-90 are strongly positively exponentially distributed, as seen in the right plot most clearly. This makes sense due to there being a strong correlation between minimum and maximum predator length when compared to average stomach samples in the same time frame.

I would assume that there is such a difference in correlation between the two time frames due to non-natural intervention ie. more fishing in 1973-80 than 1981-90, meaning that in the latter time frame there were more prey for the predators to feed on as they weren't being fished as frequently.
