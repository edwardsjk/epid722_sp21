### Code for Lecture 9 ####
library(tidyverse)
library(survival)
library(splines)

setwd("/Users/jkedwar/Box Sync/teaching/epid722/2021/lectures/L9")

set.seed(234)
n <- 1000
tau <- 12
z <- rbinom(n, 1, 0.5)
a <- ifelse(z==1, rbinom(n, 1, 0.75), rbinom(n,1,0.25))
py <- plogis(log(0.2/0.8) + log(10)*.5 - log(5)*.5 - log(10)*a + log(5)*z)
delta <- rbinom(n, 1, py)
t <- ifelse(delta == 0, tau+.01, runif(n, 0, 12))

id = c(1:n)
mydat <- data.frame(id, z, a, delta, t)

tab <- table(a,delta)
prop.table(tab, margin = 1)

tabz1 <- table(a[z==1],delta[z==1])
prop.table(tabz1, margin = 1)

tabz0 <- table(a[z==0],delta[z==0])
prop.table(tabz0, margin = 1)

### KM ----

library(survival)

# natural course ----

nckm <- survfit(Surv(t, delta) 
                ~ 1, 
                data = mydat)


nc <- data.frame(t = nckm$time, s = nckm$surv, r = 1 - nckm$surv)

ncplot <- ggplot() +
  geom_step(data = nc, aes(x = t, y = r))+
  xlab("Time")+
  ylab("Risk") +
  scale_y_continuous(limits = c(0, 0.55), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8,10, 12))+
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.28, 0.9), legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"), 
        panel.grid.major = element_line(colour = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"))

ncplot
ggsave(ncplot, file = "ncplot.png", width = 5, heigh = 5, units = "cm")




crudekm <- survfit(Surv(t, delta) 
                   ~ a, 
                   data = mydat)
crude <- data.frame(t = crudekm$time, s = crudekm$surv, r = 1 - crudekm$surv, 
                    a = c(rep(0, crudekm$strata["a=0"]), rep(1, crudekm$strata["a=1"])))



crude_plot <- ggplot() +
  geom_step(data = crude, aes(x = t, y = r,  group = factor(a), color = factor(a)))+
  xlab("Time")+
  ylab("Risk") +
  scale_color_discrete(labels = c("Unvaccinated", "Vaccinated"))+
  scale_y_continuous(limits = c(0, 0.55), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8,10, 12))+
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.28, 0.9), legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"), 
        panel.grid.major = element_line(colour = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"))

crude_plot

ggsave(crude_plot, file = "crude.png", width = 5, heigh = 5, units = "cm")



# g comp, pooled logistic

mydat_long <- mydat %>% 
  mutate(weeks = ceiling((t/12)*52)) %>% 
  uncount(weeks, .id = "week") %>% 
  group_by(id) %>% 
  mutate(y = ifelse(delta == 1 & week == last(week), 1, 0))

mod <- glm(y ~ a + z + a*z + ns(week, knots = 3), 
           data = mydat_long, 
           family = "binomial"(link = "logit"))

mydat_long$py1 <- predict(mod, 
                          newdata = mydat_long %>% mutate(a = 1), 
                          type = "response")
mydat_long$py0 <- predict(mod, 
                          newdata = mydat_long %>% mutate(a = 0), 
                          type = "response")


risks <- mydat_long %>% 
  group_by(id) %>% 
  mutate(mu1 = 1 - cumprod(1 - py1), 
         mu0 = 1 - cumprod(1 - py0)) %>% 
  group_by(week) %>% 
  summarize(r1 = mean(mu1), 
            r0 = mean(mu0)) %>% 
  select(week, r1, r0)


plrisks <- bind_rows(risks %>% 
                       select(week, r1) %>% 
                       rename(r = r1) %>% 
                       mutate(a = 1), 
                     risks %>% select(week, r0) %>% rename(r = r0) %>% mutate(a = 0))


pl_plot <- ggplot() +
  geom_step(data = plrisks, aes(x = week, y = r,  group = factor(a), color = factor(a)))+
  xlab("Time")+
  ylab("Risk") +
  scale_color_discrete(labels = c("Unvaccinated", "Vaccinated"))+
  scale_y_continuous(limits = c(0, 0.55), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  scale_x_continuous(limits = c(0, 52))+
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.28, 0.9), legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"), 
        panel.grid.major = element_line(colour = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"))

pl_plot


# g comp breslow

# A = 1
# Step 1: fit cox model among those with noart == 1, conditional on covariates, save coefficients "betahat"
coxmod <- coxph(Surv(t, delta) ~ z, data = mydat %>% filter(a == 1), method = "breslow")
betahat <- unlist(ifelse(is.na(coxmod$coef), 0, coxmod$coef))

# Step 2: create dataframe with unique event times and number of events for those in the noart == 1 group
events <- mydat %>% filter(a == 1) %>% group_by(t) %>% summarize(events = sum(delta)) %>% filter(events > 0)
event_times <- events[,1]


# Step 3: compute baseline hazard function
h0 <- rep(NA, length(event_times)) # Create empty vector to hold baseline hazard function
# Use the Breslow estimator to estimate h0 (see Lin DY. On the Breslow estimator. Lifetime data analysis. 2007 Dec 1;13(4):4 (page 474))
for(l in 1:nrow(event_times)){
  h0[l] <- as.numeric(events[l,2]) / # number of events at time l, divided by...
    sum(
      exp(
        as.matrix(
          mydat[mydat$t >= as.numeric(event_times[l,]) & mydat$a == 1, c("z")] # make a vector z for those exposed and in risk set at event time l
        ) 
        %*% betahat #multiply by betahat
      ), 
      na.rm = TRUE)
}

# Step 4: Multiply estimated baseline hazard function by the linear predictor from cox model for each person in the FULL dataset (not just noart == 1)
hw <-h0 %*% exp(betahat %*% t(as.matrix(mydat[, c("z")])))
cumhaz <- apply(hw, 2, cumsum) # estimate cumulative hazard function for each person

#Step 5, estimate risk
rt <- 1 - exp(-cumhaz) # probability of outcome at each event time for each person
r <- apply(rt, 1, FUN = mean) # average probability of outcome at each event time
chp1 <- cbind(event_times, r)
chp1$a <- 1

#a = 0


# Step 1: fit cox model among those with noart == 1, conditional on covariates, save coefficients "betahat"
coxmod <- coxph(Surv(t, delta) ~ z, data = mydat %>% filter(a == 0), method = "breslow")
betahat <- unlist(ifelse(is.na(coxmod$coef), 0, coxmod$coef))

# Step 2: create dataframe with unique event times and number of events for those in the noart == 1 group
events <- mydat %>% filter(a == 0) %>% group_by(t) %>% summarize(events = sum(delta)) %>% filter(events > 0)
event_times <- events[,1]


# Step 3: compute baseline hazard function
h0 <- rep(NA, length(event_times)) # Create empty vector to hold baseline hazard function
# Use the Breslow estimator to estimate h0 (see Lin DY. On the Breslow estimator. Lifetime data analysis. 2007 Dec 1;13(4):4 (page 474))
for(l in 1:nrow(event_times)){
  h0[l] <- as.numeric(events[l,2]) / # number of events at time l, divided by...
    sum(
      exp(
        as.matrix(
          mydat[mydat$t >= as.numeric(event_times[l,]) & mydat$a == 0, c("z")] # make a vector z for those exposed and in risk set at event time l
        ) 
        %*% betahat #multiply by betahat
      ), 
      na.rm = TRUE)
}

# Step 4: Multiply estimated baseline hazard function by the linear predictor from cox model for each person in the FULL dataset (not just noart == 0)
hw <-h0 %*% exp(betahat %*% t(as.matrix(mydat[, c("z")])))
cumhaz <- apply(hw, 2, cumsum) # estimate cumulative hazard function for each person

#Step 5, estimate risk
rt <- 1 - exp(-cumhaz) # probability of outcome at each event time for each person
r <- apply(rt, 1, FUN = mean) # average probability of outcome at each event time
chp0 <- cbind(event_times, r)
chp0$a <-  0

chp <- bind_rows(chp1, chp0)


b_plot <- ggplot() +
  geom_step(data = chp, aes(x = t, y = r,  group = factor(a), color = factor(a)))+
  xlab("Time")+
  ylab("Risk") +
  scale_color_discrete(labels = c("Unvaccinated", "Vaccinated"))+
  scale_y_continuous(limits = c(0, 0.55), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8,10, 12))+
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.28, 0.9), legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"), 
        panel.grid.major = element_line(colour = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"))

b_plot
