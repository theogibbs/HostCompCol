source("0_Functions.R")
library(doFuture)
library(progressr)
library(lhs)
library(shiny)
library(shinydashboard)
# Set up progress bar
handlers(global = TRUE)


ch <- 1
chm <- 7.5
dh <- 1
dhp <- 0

cp <- 4
dp <- 1

cm <- 6.5
dm <- 1


in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

end_time <- 10000
time_step <- 1

cur_params <- iterated_params
cur_pars <- as.list(cur_params)

root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
h_soln <- root_soln$root
mp_soln <- GetPM(h_soln, cur_pars)
p_soln <- mp_soln[1]
m_soln <- mp_soln[2]

cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)

J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
cur_eig <- GetEig(J)

# Create trajectory from random initial condition

ini_cond_sds <- c(0.001, 0.05, 0.5)
ini_cond_sd = 0.01
ini_resolution = 100
LHS_mat = randomLHS(ini_resolution, 3)
ini_states <- tibble(h = LHS_mat[,1], p = LHS_mat[,2], m = LHS_mat[,3])

plan(multisession)



all_trajs <-  data.frame()
all_trajs <- foreach(
  i = 1:nrow(ini_states),
  .combine = 'rbind',
  .options.future = list(seed = TRUE)) %dofuture% {
    
    # Cover the entire state space of initial conditions to see if any admit 2-cycles
    ini_state = ini_states[i,]
    
    out_dyn <- IntegrateDynamics(as.double(ini_state), cur_pars,
                                 end_time, time_step,
                                 fn = CoralMutPathDynamics)
    plot_df <- out_dyn %>%
      filter(time > 0) %>% 
      filter(time %% 10 == 0) %>% 
      # filter(time > (end_time - 101)) %>%
      melt(id.vars = c("time")) %>% 
      mutate(iterate = i)
    plot_df
  }


trajectory_plot <- all_trajs %>% 
  ggplot(aes(x = time, y = value, group = iterate)) +
  facet_wrap(~ variable, nrow = 1) +
  geom_line(alpha = 0.5) +
  # geom_point(size = 2, color = "black", alpha = 0.25) +
  scale_x_log10() +
  theme_classic()
trajectory_plot

hp_plot <- all_trajs %>%
  pivot_wider(names_from = variable, values_from = value) %>% 
  ggplot(aes(x = V2, y = p)) +
  geom_point(size = 2, color = "black", alpha = 0.5) + 
  theme_classic()
hp_plot

hm_plot <- all_trajs %>%
  pivot_wider(names_from = variable, values_from = value) %>% 
  ggplot(aes(x = V2, y = m)) +
  geom_point(size = 2, color = "black", alpha = 0.5) + 
  theme_classic()
hm_plot

pm_plot <- all_trajs %>%
  pivot_wider(names_from = variable, values_from = value) %>% 
  ggplot(aes(x = p, y = m)) +
  geom_point(size = 2, color = "black", alpha = 0.5) + 
  theme_classic()
pm_plot



# Plot eigenvalues
ch <- 0.05
chm <- seq(0.05, 10, length.out = 1000)
dh <- 1
dhp <- 0

cp <- 4
dp <- 1

cm <- 7.5
dm <- 1


in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

eig_df = data.frame()
eig_df <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind') %dofuture% {
    cur_params <- iterated_params[i,]
    cur_pars <- as.list(cur_params)
    
    root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
    h_soln <- root_soln$root
    mp_soln <- GetPM(h_soln, cur_pars)
    p_soln <- mp_soln[1]
    m_soln <- mp_soln[2]
    
    cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
    
    J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
    eigenvalues = eigen(J)$values
    
    temp_df <- data.frame(
      real = Re(eigenvalues),
      imag = Im(eigenvalues),
      label = c(1,2,3),
      # real_1 = Re(eigenvalues[1]),
      # imag_1 = Im(eigenvalues[1]),
      # real_2 = Re(eigenvalues[2]),
      # imag_2 = Im(eigenvalues[2]),
      # real_3 = Re(eigenvalues[3]),
      # imag_3 = Im(eigenvalues[3]),
      chm = cur_pars$chm
    )
    temp_df
  }

eig_plot <- eig_df %>% 
  # filter(label != 1) %>% 
  ggplot(aes(x = real, y = imag, color = chm)) +
  geom_vline(xintercept = 0, lwd = 1) +
  # geom_hline(yintercept = 0, lwd = 1) + 
  geom_point(shape = 16, size = 3) +
  # geom_line() +
  facet_wrap(~label,
             # scales = "free",
             ncol = 1) +
  theme_minimal()
eig_plot

chm_Re_plot <- eig_df %>% 
  # filter(label != 1) %>% 
  ggplot(aes(x = chm, y = real)) +
  geom_vline(xintercept = 0, lwd = 1) +
  # geom_hline(yintercept = 0, lwd = 1) + 
  geom_point(shape = 16, size = 1) +
  # geom_line() +
  theme_minimal()
chm_Re_plot

chm_eig_plot <- eig_df %>% 
  pivot_longer(cols = c(real, imag)) %>% 
  # filter(label != 1) %>% 
  ggplot(aes(x = chm, y = value, color = chm)) +
  geom_hline(yintercept = 0, lwd = 1) +
  geom_point(shape = 16, size = 1) +
  scale_color_viridis_c() +
  # geom_line() +
  facet_wrap(~name) +
  theme_minimal()
chm_eig_plot


handlers("progress")
eig_func <- function(chm_vec, ch_vec, cm_vec) {
  dh <- 1
  dhp <- 0
  cp <- 5
  dp <- 1
  dm <- 1
  
  full_factorial = expand.grid(chm = chm_vec, cm = cm_vec, ch = ch_vec) %>% 
    filter(chm > ch)
  
  plan(multisession, workers = 20)
  p = progressor(along = 1:dim(full_factorial)[1])
  
  real_out = foreach(
    index = 1:nrow(full_factorial),
    # j = 1:length(ch),
    # k = 1:length(cm),
    .inorder = FALSE,
    .combine = 'rbind') %dofuture% {
      # for (index in 1:nrow(full_factorial)) {
      
      # Iterate progress bar
      p()
      
      cur_pars = list(ch = full_factorial$ch[index], 
                      chm = full_factorial$chm[index], 
                      dh = dh, 
                      dhp = dhp,
                      cp = cp, 
                      dp = dp,
                      cm = full_factorial$cm[index], 
                      dm = dm)
      
      h_soln = uniroot(f = PredEq, interval = c(0, 1), cur_pars)$root
      mp_soln = GetPM(h_soln, cur_pars)
      p_soln = mp_soln[1]
      m_soln = mp_soln[2]
      
      J = BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
      eigenvalues = eigen(J)$values
      
      temp_df = data.frame(
        real = Re(eigenvalues),
        imag = Im(eigenvalues),
        label = c(1,2,3),
        chm = cur_pars$chm,
        cm = cur_pars$cm,
        ch = cur_pars$ch
      )
      temp_df
    }
}

step_size = 0.1

chm_min = 2; chm_max = 10
chm_vec = seq(chm_min, chm_max, step_size)

cm_min = 2; cm_max = 16
cm_vec = seq(cm_min, cm_max, step_size)

ch_min = 0; ch_max = 2
ch_vec = seq(ch_min, ch_max, step_size)

full_df = eig_func(chm_vec, ch_vec, cm_vec)


write.csv(full_df, "eigenvalues.csv")

full_df <- read.csv("eigenvalues.csv")
full_df2 <- full_df %>% 
  mutate(chm = signif(chm, 2),
         ch = signif(ch, 3),
         cm = signif(cm, 3))

library(shinyWidgets)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(sliderInput("slider_cm","cm", min = cm_min, max = cm_max, step = step_size, value = cm_min),
                   sliderInput("slider_ch","ch", min = ch_min, max = 2, step = step_size, value = ch_min)
                   # sliderTextInput(
                   #   inputId = "myslider",
                   #   label = "Choose a value:", 
                   #   choices = unique(full_df2$ch),
                   #   grid = TRUE
                   # )
                   ),
  dashboardBody(
    fluidRow(column(6, plotOutput('chm_eig_plot')))
  ))

server <- function(input, output, session) { 
  output$chm_eig_plot <- renderPlot({
    chm <- seq(chm_min, chm_max, step_size)
    yfxn <- function(chm) {
      filter(full_df2,
             chm == chm,
             cm == input$slider_cm,
             ch == input$slider_ch) %>% 
        select(chm, real, label)
    }
    eig_real_part <- yfxn(chm)
    
    ggplot(eig_real_part, aes(x = chm, y = real, color = as.factor(label))) +
      geom_hline(yintercept = 0, lwd = 1) +
      geom_point(shape = 16, size = 1) +
      scale_color_viridis_d() +
      theme_minimal(12)
  })
}

shinyApp(ui, server)