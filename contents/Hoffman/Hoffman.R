df <- read.csv("SAS_CHAPTER2.csv")

df$demgroup <- factor(df$demgroup)
df$sex <- factor(df$sexMW)
df$age.c <- df$age - 85
df$grip.c <- df$grip - 9


summary(lm(cognition ~ age.c + grip.c + sexMW + demgroup, df))

# age = 80, grip = 12