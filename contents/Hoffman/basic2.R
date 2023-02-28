df <- read_xlsx("SPSS_Chapter3b.xlsx")
df %>% group_by(session) %>% 
        summarise(var = var(rt))

summary(mod <- lm(rt ~ session, df))
summary(mod2 <- aov(rt ~ session, df))
