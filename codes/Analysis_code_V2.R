####################Logistic Mixed effects model####################
library(tidyverse)
library(data.table)
library(janitor)
library(lme4)
library(sjPlot)
library(webshot2)
library(afex)
library(ggeffects)
library(gtsummary)
library(MuMIn)
library(performance)

##########################Load combined data#######################
data<- read_delim(file = './data/data.csv', delim = ',')


Table2_input<- data %>%
  mutate(synch_use = factor(synch_use,
                            levels = c('Absent', 'Present'))) %>%
  mutate(risk.status = factor(risk.status, levels = c('High Risk',  
                                                      'Low Risk')))

Table2_reduced<- glmer(formula = synch_use~ (1|Participant),
                         data = Table2_input,
                         family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

Table2_three_way_interactions<- glmer(formula = synch_use~ risk.status+
                                          age.in.months+
                                          risk.status*age.in.months+
                                          (1|Participant),
                                        data = Table2_input,
                                        family = 'binomial',
                                        control = glmerControl(optimizer = 'bobyqa'))
options(na.action = 'na.fail')

dredge_Table2_three_way_interactions<- MuMIn::dredge(Table2_three_way_interactions, trace = TRUE)

options(na.action = "na.omit")

#############selected model fit ##################
selected_vars_Table2_three_way_interactions<- dredge_Table2_three_way_interactions %>%
  as.data.frame(.) %>%
  filter(delta < 2) %>%
  select(-c(1, (ncol(.)-4): ncol(.))) %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1))) %>%
  select_if(colSums(.) > 0)

colnames(selected_vars_Table2_three_way_interactions)
select_Table2_three_way_interactions<- glmer(formula = as.formula(
  paste('synch_use~', 
        gsub(pattern = ':', 
             replacement = '*', 
             x = paste(paste0(colnames(selected_vars_Table2_three_way_interactions), 
                              collapse = '+'), 
                       '(1|Participant)', sep = '+')))),
  data = Table2_input,
  family = 'binomial',
  control = glmerControl(optimizer = 'bobyqa'))
summary(select_Table2_three_way_interactions)

aov_res<- anova(Table2_reduced, select_Table2_three_way_interactions, method = 'LRT')

chi_sq<- aov_res$Chisq[2]
perf<- performance::performance(select_Table2_three_way_interactions)
r2_cond<- perf$R2_conditional
r2_marginal<- perf$R2_marginal

select_Table2_temp<- tbl_regression(select_Table2_three_way_interactions, tidy_fun = broom.mixed::tidy,
                                    exponentiate = TRUE, intercept = TRUE,
                                    label = list(risk.status ~ 'Risk Status',
                                                 age.in.months ~ 'Scaled age')) %>%
  bold_p() %>%
  add_n() %>%
  bold_labels() %>%
  tbl_split_by_rows(., variables = 'risk.status')

select_Table2_three_way_interactions<- tbl_stack(list(select_Table2_temp[[1]],select_Table2_temp[[2]]),
                                                 group_header = c('Fixed effects',
                                                                  'Random effects')) %>%
  as_gt() %>%
  gt::tab_style(style = gt::cell_text(weight = 'bold'),
                locations =gt::cells_row_groups(groups = everything())) %>%
  
  gt::tab_footnote(footnote = gt::md(glue::glue("full model vs null model: &Chi;
                                                <sup>2</sup> = {format(chi_sq, nsmall = 2)}, p < 0.0001"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>conditional</sub> = 
                                                {format(r2_cond, nsmall = 2)}"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>marginal</sub> = 
                                                {format(r2_marginal, nsmall = 2)}")))

gt::gtsave(select_Table2_three_way_interactions, 
           filename = './results/select_Table2_three_way_interactions.docx', expand = 10)

##############################


Table3_input<- data %>%
  mutate(synch_use = factor(synch_use,
                            levels = c('Absent', 'Present'))) %>%
  mutate(risk.status = factor(risk.status, levels = c('High Risk',  
                                                      'Low Risk'))) %>%
  mutate(ia_loc_type = ia_loc) %>%
  mutate(ia_loc= ifelse(ia_loc_type == 'Away', 0 ,1),
         ia_loc = factor(ia_loc, levels = c(0,1)))


Table3_reduced<- glmer(formula = ia_loc~ (1|Participant),
                         data = Table3_input,
                         family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

Table3_three_way_interactions<- glmer(formula = ia_loc~ risk.status+
                                          synch_use+
                                          age.in.months+
                                          risk.status*age.in.months+
                                          risk.status*synch_use+
                                          synch_use*age.in.months+
                                          risk.status*synch_use*age.in.months+
                                          (1|Participant),
                                        data = Table3_input,
                                        family = 'binomial',
                                        control = glmerControl(optimizer = 'bobyqa'))

options(na.action = 'na.fail')

dredge_Table3_three_way_interactions<- dredge(Table3_three_way_interactions, 
                                                trace = TRUE)

options(na.action = "na.omit")

###################Selected model fit###############################
selected_vars_Table3_three_way_interactions<- dredge_Table3_three_way_interactions %>%
  as.data.frame(.) %>%
  filter(delta < 2) %>%
  select(-c(1, (ncol(.)-4): ncol(.))) %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1))) %>%
  select_if(colSums(.) > 0)

colnames(selected_vars_Table3_three_way_interactions)
select_Table3_three_way_interactions<- glmer(formula = as.formula(
  paste('ia_loc~', 
        gsub(pattern = ':', 
             replacement = '*', 
             x = paste(paste0(colnames(selected_vars_Table3_three_way_interactions), 
                              collapse = '+'), 
                       '(1|Participant)', sep = '+')))),
  data = Table3_input,
  family = 'binomial',
  control = glmerControl(optimizer = 'bobyqa'))
summary(select_Table3_three_way_interactions)

aov_res<- anova(Table3_reduced, select_Table3_three_way_interactions, method = 'LRT')
chi_sq<- aov_res$Chisq[2]
perf<- performance::performance(select_Table3_three_way_interactions)
r2_cond<- perf$R2_conditional
r2_marginal<- perf$R2_marginal


select_Table3_temp<- tbl_regression(select_Table3_three_way_interactions, tidy_fun = broom.mixed::tidy,
                                    exponentiate = TRUE,intercept = TRUE,
                                    label = list(risk.status ~ 'Risk Status',
                                                 synch_use ~ 'Synchrony',
                                                 age.in.months ~ 'Scaled age')) %>%
  bold_p() %>%
  add_n() %>%
  bold_labels() %>%
  tbl_split_by_rows(., variables = 'age.in.months:risk.status:synch_use')

select_Table3_three_way_interactions<- tbl_stack(list(select_Table3_temp[[1]],select_Table3_temp[[2]]),
                                                 group_header = c('Fixed effects',
                                                                  'Random effects')) %>%
  as_gt() %>%
  gt::tab_style(style = gt::cell_text(weight = 'bold'),
                locations =gt::cells_row_groups(groups = everything())) %>%
  
  gt::tab_footnote(footnote = gt::md(glue::glue("full model vs null model: &Chi;
                                                <sup>2</sup> = {format(chi_sq, nsmall = 2)}, p < 0.0001"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>conditional</sub> = 
                                                {format(r2_cond, nsmall = 2)}"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>marginal</sub> = 
                                                {format(r2_marginal, nsmall = 2)}")))

gt::gtsave(select_Table3_three_way_interactions, 
           filename = './results/select_Table3_three_way_interactions.docx', expand = 10)

###############################



Table4_input <- data %>%
  mutate(synch_use = factor(synch_use,
                            levels = c('Absent', 'Present'))) %>%
  mutate(risk.status = factor(risk.status, levels = c('High Risk',  
                                                      'Low Risk'))) %>%
  dplyr::filter(ia_loc != 'Away') %>%
  droplevels() %>%
  mutate(ia_loc = factor(ia_loc, levels = c('Mother', 'Object')))

Table4_reduced<- glmer(formula = ia_loc~ (1|Participant),
                         data = Table4_input,
                         family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

Table4_three_way_interactions<- glmer(formula = ia_loc~ risk.status+
                                          age.in.months+
                                          synch_use+
                                          risk.status*age.in.months+
                                          risk.status*synch_use+
                                          synch_use*age.in.months+
                                          risk.status*synch_use*age.in.months+
                                          (1|Participant),
                                        data = Table4_input,
                                        family = 'binomial',
                                        control = glmerControl(optimizer = 'bobyqa'))


options(na.action = 'na.fail')

dredge_Table4_three_way_interactions<- dredge(Table4_three_way_interactions, trace = TRUE)

options(na.action = "na.omit")

######################Selected model fit##################################
selected_vars_Table4_three_way_interactions<- dredge_Table4_three_way_interactions %>%
  as.data.frame(.) %>%
  filter(delta < 2) %>%
  select(-c(1, (ncol(.)-4): ncol(.))) %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1))) %>%
  select_if(colSums(.) > 0)

colnames(selected_vars_Table4_three_way_interactions)
select_Table4_three_way_interactions<- glmer(formula = as.formula(
  paste('ia_loc~', 
        gsub(pattern = ':', 
             replacement = '*', 
             x = paste(paste0(colnames(selected_vars_Table4_three_way_interactions), 
                              collapse = '+'), 
                       '(1|Participant)', sep = '+')))),
  data = Table4_input,
  family = 'binomial',
  control = glmerControl(optimizer = 'bobyqa'))
summary(select_Table4_three_way_interactions)

aov_res<- anova(Table4_reduced, select_Table4_three_way_interactions, method = 'LRT')
chi_sq<- aov_res$Chisq[2]
perf<- performance::performance(select_Table4_three_way_interactions)
r2_cond<- perf$R2_conditional
r2_marginal<- perf$R2_marginal

select_Table10_temp<- tbl_regression(select_Table4_three_way_interactions, tidy_fun = broom.mixed::tidy,
                                     exponentiate = TRUE,intercept = TRUE,
                                     label = list(risk.status ~ 'Risk Status',
                                                  age.in.months ~ 'Scaled age',
                                                  synch_use ~ 'Synchrony')) %>%
  bold_p() %>%
  add_n() %>%
  bold_labels() %>%
  tbl_split_by_rows(., variables = 'risk.status:synch_use')

select_Table10<- tbl_stack(list(select_Table10_temp[[1]],select_Table10_temp[[2]]),
                           group_header = c('Fixed effects',
                                            'Random effects')) %>%
  as_gt() %>%
  gt::tab_style(style = gt::cell_text(weight = 'bold'),
                locations =gt::cells_row_groups(groups = everything())) %>%
  
  gt::tab_footnote(footnote = gt::md(glue::glue("full model vs null model: &Chi;
                                                <sup>2</sup> = {format(chi_sq, nsmall = 2)}, p < 0.0001"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>conditional</sub> = 
                                                {format(r2_cond, nsmall = 2)}"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>marginal</sub> = 
                                                {format(r2_marginal, nsmall = 2)}")))

gt::gtsave(select_Table10, filename = './results/select_Table4_respecified_three_way_interactions.docx', expand = 10)

#####################################################


Table4_input<- data %>%
  mutate(synch_use = factor(synch_use,
                            levels = c('Absent', 'Present'))) %>%
  mutate(risk.status = factor(risk.status, levels = c('High Risk',  
                                                      'Low Risk'))) %>%
  dplyr::filter(synch_use != 'Absent') %>%
  droplevels() %>%
  mutate(ia_loc_fac = ifelse(ia_loc != 'Object', 'No Object', 'Object'),
         ia_loc_fac = factor(ia_loc_fac, levels = c('No Object', 'Object')))



Table4_reduced<- glmer(formula = ia_loc_fac~ (1|Participant),
                         data = Table4_input,
                         family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

Table4_three_way_interactions<- glmer(formula = ia_loc_fac~ risk.status+
                                          age.in.months+
                                          risk.status*age.in.months+
                                          (1|Participant),
                                        data = Table4_input,
                                        family = 'binomial',
                                        control = glmerControl(optimizer = 'bobyqa'))

options(na.action = 'na.fail')

dredge_Table4_three_way_interactions<- dredge(Table4_three_way_interactions, trace = TRUE)

options(na.action = "na.omit")




###########################Selected model fit#################################

selected_vars_Table4_three_way_interactions<- dredge_Table4_three_way_interactions %>%
  as.data.frame(.) %>%
  filter(delta < 2) %>%
  select(-c(1, (ncol(.)-4): ncol(.))) %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1))) %>%
  select_if(colSums(.) > 0)

colnames(selected_vars_Table4_three_way_interactions)


select_Table4_three_way_interactions<- glmer(formula = as.formula(
  paste('ia_loc_fac~', 
        gsub(pattern = ':', 
             replacement = '*', 
             x = paste(paste0(colnames(selected_vars_Table4_three_way_interactions), 
                              collapse = '+'), 
                       '(1|Participant)', sep = '+')))),
  data = Table4_input,
  family = 'binomial',
  control = glmerControl(optimizer = 'bobyqa'))
summary(select_Table4_three_way_interactions)

aov_res<- anova(Table4_reduced, select_Table4_three_way_interactions, method = 'LRT')
chi_sq<- aov_res$Chisq[2]
perf<- performance::performance(select_Table4_three_way_interactions)
r2_cond<- perf$R2_conditional
r2_marginal<- perf$R2_marginal


select_Table11_temp<- tbl_regression(select_Table4_three_way_interactions, tidy_fun = broom.mixed::tidy,
                                     exponentiate = TRUE,intercept = TRUE,
                                     label = list(risk.status ~ 'Risk Status',
                                                  age.in.months ~ 'Scaled age')) %>%
  bold_p() %>%
  add_n() %>%
  bold_labels() %>%
  tbl_split(., variables = 'age.in.months:risk.status')

select_Table11<- tbl_stack(list(select_Table11_temp[[1]],select_Table11_temp[[2]]),
                           group_header = c('Fixed effects',
                                            'Random effects')) %>%
  as_gt() %>%
  gt::tab_style(style = gt::cell_text(weight = 'bold'),
                locations =gt::cells_row_groups(groups = everything())) %>%
  
  gt::tab_footnote(footnote = gt::md(glue::glue("full model vs null model: &Chi;
                                                <sup>2</sup> = {format(chi_sq, nsmall = 2)}, p < 0.0001"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>conditional</sub> = 
                                                {format(r2_cond, nsmall = 2)}"))) %>% 
  gt::tab_footnote(footnote = gt::md(glue::glue("R<sup>2</sup><sub>marginal</sub> = 
                                                {format(r2_marginal, nsmall = 2)}")))

gt::gtsave(select_Table11, filename = './results/select_Table4_three_way_interactions.docx', expand = 10)

#################################################################

#############################Parametric bootstrap run##############################

source(file = './codes/Parametric_bootstrap_functions.R')

Suppl_table1<- boot_fixef_ci(select_Table2_three_way_interactions, 
                                                            parallel = 'multicore',
                                                            ncpus = 24, nsim=1000, seed=123)

Suppl_table1_parametric_simulation<- gt::gt(Suppl_table1)

gt::gtsave(Suppl_table1_parametric_simulation,
           filename = './results/Suppl_table1.docx', expand = 10)


Suppl_table3 <- boot_fixef_ci(select_Table3_three_way_interactions, 
                                             parallel = 'multicore',
                                             ncpus = 24, nsim=1000, seed=123)

Suppl_table3_parametric_simulation<- gt::gt(Suppl_table3)

gt::gtsave(Suppl_table3_parametric_simulation,
           filename = './results/Suppl_table3.docx', expand = 10)


Suppl_table4 <- boot_fixef_ci(select_Table4_three_way_interactions, 
                                             parallel = 'multicore',
                                             ncpus = 24, nsim=1000, seed=123)

Suppl_table4_parametric_simulation<- gt::gt(Suppl_table4)

gt::gtsave(Suppl_table4_parametric_simulation,
           filename = './results/Suppl_table4.docx', expand = 10)


Suppl_table5 <- boot_fixef_ci(select_Table4_three_way_interactions, 
                                             parallel = 'multicore',
                                             ncpus = 24, nsim=1000, seed=123)

Suppl_table5_parametric_simulation<- gt::gt(Suppl_table5)

gt::gtsave(Suppl_table5_parametric_simulation,
           filename = './results/Suppl_table5.docx', expand = 10)

####################################################################################

###################Bias corrected parametric simulation plots##############################

####################################Implementation###########
synch_cols<- c('Absent' = "#000000",
               'Present' = "#808080")

risk_cols<- c('High Risk' = "#000000",
              'Low Risk' = "#808080")

###############################################################
figure3 <- boot_glmm_marginal_plot(
  mod = select_Table2_three_way_interactions,
  data = Table2_input,
  age_var = "age.in.months",
  group_vars = 'risk.status',
  color_var = "risk.status",
  condition = list(risk.status = 'High Risk'),
  trim = c(0, 1),               
  n_re = 1500,
  nsim = 1000,
  parallel = TRUE,
  ncpus = 24
)

figure3_plot<- figure3$plot+
  xlab('Scaled age')+ ylab('Predicted probabilities of \n Synchrony use')+
  scale_color_manual(values = color_shades)+
  scale_fill_manual(values = color_shades)+
  scale_y_continuous(breaks = seq(0.25, 1, 0.25), 
                     labels = c("25%", "50%", "75%", "100%"), 
                     limits = c(0, 1))+
  theme(legend.position = 'none')


ggsave(filename = './results/figure3_plot.png',
       plot = figure3_plot, width = 6, 
       height = 5, 
       units = 'in',
       device = 'png')


###############################################################

#################Pairwise comparison##########################
pred_Table3_age_synch_risk <- ggpredict(
  select_Table3_three_way_interactions,
  terms = c("age.in.months [all]", "risk.status", "synch_use")
)

pred_Table3_age_synch_risk_df <- as.data.frame(pred_Table3_age_synch_risk) 


Table3_pw_age_risk_synch<- emmeans::emtrends(select_Table3_three_way_interactions,
                                            specs = ~ synch_use | risk.status,
                                            var = 'age.in.months')
slope_Table3_age_risk_synch_pairs<- pairs(Table3_pw_age_risk_synch)

slope_Table3_age_risk_synch_tbl<- slope_Table3_age_risk_synch_pairs %>% 
  broom::tidy() %>% 
  mutate(stars = case_when(
    p.value < .001 ~ "***",
    p.value < .01  ~ "**",
    p.value < .05  ~ "*",
    TRUE           ~ " "))


ann_slope_Table3_age_risk_synch <- slope_Table3_age_risk_synch_tbl %>%
  
  mutate(risk.status = factor(risk.status, levels = c('High Risk', 'Low Risk'))) %>% 
  mutate(
    x = max(pred_Table3_age_synch_risk_df$x) * 0.9,                            
    y = 0.9,
    label = stars
  )
###################################################################

color_shades<- synch_cols

figure4 <- boot_glmm_marginal_plot(
  mod = select_Table3_three_way_interactions,
  data = Table3_input,
  age_var = "age.in.months",
  group_vars = c("synch_use", 'risk.status'),
  facet_var = "risk.status",          
  color_var = "synch_use",        
  trim = c(0, 1),               
  n_re = 1500,
  nsim = 1000,
  parallel = TRUE,
  ncpus = 24
)

figure4_plot<- figure4$plot+
  xlab('Scaled age')+ ylab('Predicted probabilities of \n infant attention response')+
  scale_color_manual(values = color_shades)+
  scale_fill_manual(values = color_shades)+
  scale_y_continuous(breaks = seq(0.25, 1, 0.25), 
                     labels = c("25%", "50%", "75%", "100%"), 
                     limits = c(0, 1))+
  labs(title = '',
       color = 'Synchrony',
       fill = 'Synchrony')+
  geom_text(data = ann_slope_Table3_age_risk_synch,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 10)+
  theme(strip.text = element_text(face = 'bold', size = 19))


ggsave(filename = './results/figure4_plot.png',
       plot = figure4_plot, width = 8, 
       height = 5, 
       units = 'in',
       device = 'png')



##########################################################################

color_shades<- risk_cols

figure5a <- boot_glmm_marginal_plot(
  mod = select_Table4_three_way_interactions,
  data = Table4_input,
  age_var = "age.in.months",
  group_vars = 'risk.status',
  color_var = "risk.status",       
  trim = c(0, 1),               
  n_re = 1500,
  nsim = 1000,
  parallel = TRUE,
  ncpus = 24
)

figure5a_plot<- figure5a$plot+
  xlab('Scaled age')+ ylab('Predicted probabilities of \n infant attention (to Object)')+
  scale_color_manual(values = color_shades)+
  scale_fill_manual(values = color_shades)+
  scale_y_continuous(breaks = seq(0.25, 1, 0.25), 
                     labels = c("25%", "50%", "75%", "100%"), 
                     limits = c(0, 1))+
  labs(title = '',
       color = 'Risk Status',
       fill = 'Risk Status')


ggsave(filename = './results/figure5a_plot.png',
       plot = figure5a_plot, width = 6, 
       height = 5, 
       units = 'in',
       device = 'png')



##########################################################################



color_shades<- risk_cols

figure5b <- boot_glmm_marginal_plot(
  mod = select_Table4_three_way_interactions,
  data = Table4_input,
  age_var = "age.in.months",
  group_vars = 'risk.status',
  color_var = "risk.status",
  condition = list(risk.status = 'High Risk'),
  trim = c(0, 1),               
  n_re = 1500,
  nsim = 1000,
  parallel = TRUE,
  ncpus = 24
)

figure5b_plot<- figure5b$plot+
  xlab('Scaled age')+ ylab('Predicted probabilities of \n infant attention (to Object)')+
  scale_color_manual(values = color_shades)+
  scale_fill_manual(values = color_shades)+
  scale_y_continuous(breaks = seq(0.25, 1, 0.25), 
                     labels = c("25%", "50%", "75%", "100%"), 
                     limits = c(0, 1))+
  theme(legend.position = 'none')


ggsave(filename = './results/figure5b_plot.png',
       plot = figure5b_plot, width = 6, 
       height = 5, 
       units = 'in',
       device = 'png')
############################################################################


figure6b_boot <- boot_glmm_pointrange(
  mod        = select_Table4_three_way_interactions,
  data       = Table4_input,
  focal_var  = "synch_use",
  condition  = list(age.in.months = 0,          
                    risk.status   = "High Risk"),
  n_re       = 1500,
  nsim       = 1000,
  parallel   = TRUE,
  ncpus      = 24
)

figure6b_df <- figure6b_boot$plot_df

synch_use_cols <- c("Absent" = "#000000", "Present" = "#808080")

figure6b_plt <- ggplot() +
  geom_pointrange(
    data    = figure6b_df,
    mapping = aes(x = x, y = predicted,
                  ymin = conf.low, ymax = conf.high,
                  group = x, color = x, fill = x),
    lwd = 0.9, size = 2, shape = 16
  ) +
  scale_y_continuous(expand = c(0, 0.2)) +
  scale_color_manual(values = synch_use_cols) +
  ggpubr::geom_bracket(
    xmin        = c("Absent"),
    xmax        = c("Present"),
    label       = c("***"),
    y.position  = c(0.9),
    label.size  = 9,
    size        = 0.75
  ) +
  xlab("Synchrony") +
  ylab("Predicted probabilities of \n infant attention (to Object)") +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x     = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.y     = element_text(face = "bold", size = 15),
    axis.title      = element_text(face = "bold", size = 17),
    panel.grid      = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = "./results/figure6b_plot.png",
  plot     = figure6b_plt,
  width = 6, height = 5, units = "in", device = "png"
)


##################### Risk plot######################################
figure6a_boot <- boot_glmm_pointrange(
  mod        = select_Table4_three_way_interactions,
  data       = Table4_input,
  focal_var  = "risk.status",
  condition  = list(age.in.months = 0,          
                    synch_use   = "Absent"),
  n_re       = 1500,
  nsim       = 1000,
  parallel   = TRUE,
  ncpus      = 24
)

figure6a_df <- figure6a_boot$plot_df

risk_use_cols <- c("High Risk" = "#000000", "Low Risk" = "#808080")

figure6a_plt <- ggplot() +
  geom_pointrange(
    data    = figure6a_df,
    mapping = aes(x = x, y = predicted,
                  ymin = conf.low, ymax = conf.high,
                  group = x, color = x, fill = x),
    lwd = 0.9, size = 2, shape = 16
  ) +
  scale_y_continuous(expand = c(0, 0.2)) +
  scale_color_manual(values = risk_use_cols) +
  ggpubr::geom_bracket(
    xmin        = c("High Risk"),
    xmax        = c("Low Risk"),
    label       = c("*"),
    y.position  = c(0.9),
    label.size  = 9,
    size        = 0.75
  ) +
  xlab("Risk Status") +
  ylab("Predicted probabilities of \n infant attention (to Object)") +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x     = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.y     = element_text(face = "bold", size = 15),
    axis.title      = element_text(face = "bold", size = 17),
    panel.grid      = element_blank(),
    legend.position = "none"
  )


figure6a_plt

ggsave(
  filename = "./results/figure6a_plot.png",
  plot     = figure6a_plt,
  width = 6, height = 5, units = "in", device = "png"
)


####################Model-2e risk plot##########################################

##################### Risk plot######################################
figure5c <- boot_glmm_pointrange(
  mod        = select_Table4_three_way_interactions,
  data       = Table4_input,
  focal_var  = "risk.status",
  condition  = list(age.in.months = 0,          
                    synch_use   = "Absent"),
  n_re       = 1500,
  nsim       = 1000,
  parallel   = TRUE,
  ncpus      = 24
)

figure5c_df <- figure5c$plot_df

risk_use_cols <- c("High Risk" = "#000000", "Low Risk" = "#808080")

figure5c_plot <- ggplot() +
  geom_pointrange(
    data    = figure5c_df,
    mapping = aes(x = x, y = predicted,
                  ymin = conf.low, ymax = conf.high,
                  group = x, color = x, fill = x),
    lwd = 0.9, size = 2, shape = 16
  ) +
  scale_y_continuous(expand = c(0, 0.2)) +
  scale_color_manual(values = risk_use_cols) +
  ggpubr::geom_bracket(
    xmin        = c("High Risk"),
    xmax        = c("Low Risk"),
    label       = c("*"),
    y.position  = c(0.9),
    label.size  = 9,
    size        = 0.75
  ) +
  xlab("Risk Status") +
  ylab("Predicted probabilities of \n infant attention (to Object)") +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x     = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.y     = element_text(face = "bold", size = 15),
    axis.title      = element_text(face = "bold", size = 17),
    panel.grid      = element_blank(),
    legend.position = "none"
  )




ggsave(
  filename = "./results/figure5c_plot.png",
  plot     = figure5c_plot,
  width = 6, height = 5, units = "in", device = "png"
)

