# call functions file
source("0_Functions.R")

# Bifurcation diagrams as the benefit to host colonization from the mutualist is increased ----

# Choose parameter values. Choose a range of values for the added colonization rate of host due to mutualist settlement. Choose three values for colonizer (mutualist) colonization rate
ch <- 0.25
chm <- seq(0.25+1e-10, 25, length.out = 200) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- c(10, 11, 15)
dm <- 1

# get all possible combinations of parameter values
in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# choose initial condition standard deviation
ini_cond_sds <- c(0.001, 0.5)

# indiciate how long to run simulation and time step size
end_time <- 4000
time_step <- 0.5

# initialize output data frame
out_data <- data.frame()

out_data <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind') %:%
  foreach(
    ini_cond_sd = ini_cond_sds,
    .combine = 'rbind',
    .inorder = FALSE) %dofuture% {
      cur_params <- iterated_params[i,]
      cur_pars <- as.list(cur_params)
      
      # solve for equilibrium for each population
      root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
      h_soln <- root_soln$root
      mp_soln <- GetPM(h_soln, cur_pars)
      p_soln <- mp_soln[1]
      m_soln <- mp_soln[2]
      
      # determine feasibility
      cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
      
      # choose positive initial conditions to run simulation
      ini_state <- c(h_soln, p_soln, m_soln) + rnorm(n = 3, mean = 0, sd = ini_cond_sd)
      ini_state[ini_state <= 0] <- 0.0001
      names(ini_state) <- c("Host", "Pathogen", "Mutualist")
      
      # simulate the dynamics
      out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                                   end_time, time_step,
                                   fn = CoralMutPathDynamics)
      
      # extract last 100 time steps of the output
      out_dyn <- out_dyn %>%
        filter(time > (end_time - 100)) %>%
        melt(id.vars = c("time"))
      
      # extract leading eigenvalue for the Jacobian matrix to determine stability
      J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
      cur_eig <- GetEig(J)
      
      cur_params$IniCondSd <- ini_cond_sd
      
      # label outcomes accordingly with respective population, feasibility, and stability.
      cur_dyn <- cbind(cur_params, out_dyn) %>%
        mutate(pred_soln = case_when(variable == "Host" ~ h_soln,
                                     variable == "Mutualist" ~ m_soln,
                                     variable == "Pathogen" ~ p_soln)) %>%
        mutate(Eigenvalue = cur_eig,
               Stable = ifelse(cur_feas, ifelse(cur_eig < 0,
                                                "Feasible and stable",
                                                "Feasible but unstable"),
                               "Not feasible"))
      cur_dyn
    }

proc_data <- out_data %>%
  mutate(pred_soln = ifelse(pred_soln > 0, pred_soln, NA))

# Reorder data to have labels in descending order
proc_data <- mutate(proc_data, cm_label = paste0("c[m] == ", cm))
proc_data$cm_label <- factor(proc_data$cm_label,
                             levels = c("c[m] == 15", "c[m] == 11", "c[m] == 10", "c[m] == 5"))

# plot output dataframe. X-axis represents additional host colonization rate from mutualist settlement while y-axis represents frequency of the respective population
# Blue indicates feasible and stable, yellow indicates feasible and unstable, and red indicates infeasible equilbria
# Each panel is a different population: host, mutualist (colonizer), or pathogen (competitor)
plMutBenefit <- ggplot(proc_data,
                       aes(x = chm, y = value)) +
  facet_grid(
    rows = vars(cm_label),
    cols = vars(variable),
    labeller = label_parsed
  ) +
  geom_point(size = 2, alpha = 0.005) + theme_classic() +
  geom_line(linewidth = 1, aes(x = chm, y = pred_soln, color = Stable)) +
  scale_color_manual(breaks = c("Feasible and stable", "Feasible but unstable", "Not feasible"),
                     values = c("#0072B2", "#F0E442", "#D55E00")) +
  labs(x = expression("Added Host Colonization from the Mutualist" ~ (c[hm])),
       y = "Frequency",
       color = "",
       linetype = "") +
  ggtitle("A") +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        legend.position = "none",
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plMutBenefit # Figure 3 A in the main text

# Feasibility and stability of coexistence as mutualist colonization rate and benefit to host colonization are varied ----

# initialize parameter values. Choose range of values for additional host colonization rate due to mutualist settlement and mutualist colonization rate. Choose 2 values for intrinsic host colonization rate. 
ch <- c(0.25, 1.25)
chm <- seq(1.75, 25, length.out = 200) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- seq(4, 16, length.out = 200)
dm <- 1

# Find all possible combinations of parameters
in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

# initialize output dataframe
out_data <- data.frame()

out_data <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind') %dofuture% {
    
    # initialize current iteration's parameter values
    cur_params <- iterated_params[i,]
    cur_pars <- as.list(cur_params)
    
    # solve for equilibria of each population
    root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
    h_soln <- root_soln$root
    mp_soln <- GetPM(h_soln, cur_pars)
    p_soln <- mp_soln[1]
    m_soln <- mp_soln[2]
    
    # extract leading eigenvalue from the Jacobian matrix to determine stability
    J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
    cur_eig <- GetEig(J)
    
    # Determine feasibility and stability of equilbria
    cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
    cur_stable <- cur_eig < 0
    
    # Label the feasibility and stability accordingly
    cur_outcome <- ifelse(!cur_feas, "\nNot feasible\n",
                          ifelse(cur_stable, "\nFeasible\nand stable\n",
                                 "\nFeasible\nbut unstable\n"))
    
    cur_dyn <- cbind(cur_params, data.frame(Outcome = cur_outcome))
    cur_dyn
  }

# Reorder data to have labels in descending order
out_data <- mutate(out_data, ch_label = paste0("c[h] == ", ch))
out_data$ch_label <- factor(out_data$ch_label,
                            levels = c("c[h] == 1.75", "c[h] == 1.25", "c[h] == 0.75", "c[h] == 0.25"))

# create heatmap to visualize output data. X-axis represents additional host colonization due to mutualist settlement while the y-axis represents mutualist (colonizer) colonization rate
# Each panel indicates a different intrinsic host colonization rate
# Blue indices feasible and stable, yellow indicates feasible and unstable, and red indicates infeasible. 
plHeatMap <- ggplot(out_data,
                    aes(x = chm, y = cm, fill = Outcome)) +
  geom_tile() + theme_classic() +
  scale_fill_manual("Coexistence\nstatus:", values = c("#0072B2", "#F0E442", "#D55E00")) +
  labs(x = expression(atop("Added Host Colonization from the Mutualist" ~ (c[hm]))),
       y = expression("Mutualist Colonization" ~ (c[m])),
       fill = "") +
  facet_grid(
    rows = vars(ch_label),
    labeller = label_parsed
  ) +
  ggtitle("B") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text( size = 10 ),
        panel.background = element_rect(fill = NA),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plHeatMap # Figure 3 B in the main text

# Example trajectories where a limit cycle is present ----

# initialize parameters
ch <- 0.5
chm <- 10 # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- 10
dm <- 1

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

# set initial conditions
ini_state <- c(1, 0.5, 0.5) * 0.125
names(ini_state) <- c("Host", "Pathogen", "Mutualist")

# set length of simulation and time step size
end_time <- 50
time_step <- 0.001

# simulate dynamics
out_dyn <- IntegrateDynamics(ini_state, pars,
                             end_time, time_step,
                             fn = CoralMutPathDynamics)

melt_dyn <- melt(out_dyn, id.vars = c("time"))

# plot output. X-axis is time and y-axis is proportion of space occupied by the respective population. 
plDyn <- ggplot(melt_dyn, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1) + 
  theme_classic() + 
  labs(x = "Time", y = "Frequency", color = "", linetype = "") +
  ggtitle("C") +
  scale_color_viridis_d() +
  theme(text = element_text(size=15),
        legend.position = "top",
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plDyn # Figure 3 C in the main text

# Phase-space trajectories as the limit cycle is approached ----

# initialize parameters values. Choose 4 different values for additional host colonization rate from mutualist settlement
ch <- 0.5
chm <- c(3, 3.56, 3.57, 10) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- 10
dm <- 1

# find all combinations of parameter values
in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# Set length of simulation and time step size
end_time <- 2e3
time_step <- 0.1

# initialize output dataframe
out_data <- data.frame()

for(i in 1:nrow(iterated_params)) {
  
  # intiailize current iteration's parameters
  cur_params <- iterated_params[i,]
  cur_pars <- as.list(cur_params)
  
  # solve for equilibria of all populations
  root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
  h_soln <- root_soln$root
  mp_soln <- GetPM(h_soln, cur_pars)
  p_soln <- mp_soln[1]
  m_soln <- mp_soln[2]
  
  # Choose initial condition
  ini_state <- c(h_soln, p_soln, m_soln) + rnorm(n = 3, mean = 0, sd = 0.001)
  ini_state[ini_state <= 0] <- 0.001
  names(ini_state) <- c("host", "pathogen", "mutualist")
  
  # simulate the dynamics
  out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                               end_time, time_step,
                               fn = CoralMutPathDynamics)
  
  cur_dyn <- cbind(cur_params, out_dyn)
  
  out_data <- rbind(out_data, cur_dyn)
  
}

# indicate final frequencies at the end of the simulation
end_points <- out_data %>%
  filter(time == max(time)) %>%
  mutate(PatEnd = pathogen) %>% 
  mutate(MutEnd = mutualist) %>% 
  select(chm, PatEnd, MutEnd) %>% 
  distinct()

plot_data <- out_data %>%
  right_join(end_points, by = "chm")

# plot output dataframe. X-axis is pathogen frequency while y-axis is mutualist frequency. Color indicates time at which the respective frequencies occurred. 
plBifurcation <- ggplot(plot_data, aes(x = pathogen, y = mutualist, color = time)) +
  geom_path(linewidth = 0.5) +
  geom_point(aes(x = PatEnd, y = MutEnd), color = "red", size = 2) +
  facet_wrap(~chm, 
             scales = "free",
             labeller = label_bquote(cols = c[hm] == .(chm)), nrow = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0.005), ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0.005)) +
  scale_color_viridis_c(option = "G", direction = -1) +
  labs(x = "Pathogen",
       y = "Mutualist",
       color = "Time") +
  ggtitle("D") +
  theme(text = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.key.width = unit(2, 'cm'),
        legend.position = "top",
        legend.justification = "center", 
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plBifurcation # Figure 3 D in the main text

# Regions of stability for initial conditions when there are multiple stable states ----

# Demonstrates the allee effect in the model

# initialize parameter values
ch <- 0.5
chm <- 10 # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- 10
dm <- 1

in_pars <- data.frame(ch = ch, chm = chm, dh = dh, dhp = dhp,
                      cp = cp, dp = dp,
                      cm = cm, dm = dm)

# set duration of simulationa, time step size, and minimum population frequency
end_time <- 5e2
time_step <- 1
freq_cutoff <- 1e-7

# initialize a range of intiial conditions for host and mutualist populations
ini_hosts <- seq(0.001, 0.9, length.out = 200)
ini_muts <- seq(0.0001, 0.15, length.out = 200)

# choose one initial condition for pathogen population
ini_path <- 0.00005 # should be less than the smaller value above

ini_vals <- expand_grid(host = ini_hosts, mut = ini_muts) %>% 
  filter(mut < host - ini_path)

# initialize output dataframe
out_data <- data.frame()
i=1
plan(multisession)

out_data <- foreach(
  j = 1:nrow(ini_vals),
  .combine = 'rbind'
) %dofuture% {
  ini_host = ini_vals$host[j]
  ini_mut = ini_vals$mut[j]  
  
  # initialize parameters
  cur_params <- in_pars[i,]
  cur_pars <- as.list(cur_params)
  
  # set initial conditions for the simulation
  ini_state <- c(ini_host, ini_path, ini_mut)
  names(ini_state) <- c("Host", "Pathogen", "Mutualist")
  
  # simulate the dynamics for the iteration's respective initial conditions
  out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                               end_time, time_step,
                               fn = CoralMutPathDynamics)
  
  # extract final densities of each population at the end of the simulation
  end_state <- as.numeric(out_dyn[nrow(out_dyn),2:4])
  
  cur_params$IniMut <- ini_mut
  cur_params$IniHost <- ini_host
  cur_params$IniPath <- ini_path
  
  # Determine coexistence comparative to the frequency cutoff set
  out_coexist <- data.frame(Coexist = ifelse(prod(end_state > freq_cutoff),
                                             "Coexistence",
                                             "No coexistence"))
  
  cur_dyn <- cbind(cur_params, out_coexist)
}

out_data <- out_data %>%
  mutate(Coexist = ifelse(IniMut + IniPath > IniHost, NA, Coexist))

# plot output dataframe. X-axis represents initial condition of the host population while the y-axis represents the initial condition of the mutualist (colonizer) population
# Blue represents coexistence while grey represents no coexistence
plIniCond <- ggplot(out_data,
                    aes(x = IniHost, y = IniMut, fill = Coexist)) +
  geom_tile() + theme_classic() +
  scale_fill_manual(values = c("blue", "gray"), na.value="white", na.translate = F) +
  labs(x = "Initial Frequency of the Host",
       y = "Initial Frequency of the Mutualist",
       fill = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("E") +
  theme(axis.text = element_text( size = 10),
        legend.position = "top",
        panel.background = element_rect(fill = NA),
        strip.background = element_blank(),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plIniCond # Figure 3 E in the main text

# Set up and save figure ----
layout_mat <- matrix(0, nrow = 5, ncol = 10)
layout_mat[1,] <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
layout_mat[2,] <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
layout_mat[3,] <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
layout_mat[4,] <- c(3, 3, 3, 4, 4, 4, 4, 5, 5, 5)
layout_mat[5,] <- c(3, 3, 3, 4, 4, 4, 4, 5, 5, 5)

jpeg("./figs/Fig3MutBenefit.jpeg",
     width = 4100, height = 3250, res = 300)
grid.arrange(plMutBenefit, plHeatMap,
             plDyn, plBifurcation, plIniCond,
             layout_matrix = layout_mat)
dev.off()





