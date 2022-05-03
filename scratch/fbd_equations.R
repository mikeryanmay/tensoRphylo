library(ape)
library(dispRity)

format = function(probs) {
  cat(paste0("(",sprintf("%.10e", probs),")"), sep="")
}

factor = function(tree) {

  # compute the ages
  ages = tree.age(tree, digits=10)

  # compute the number of extinct and extant species
  taxa        = tree$tip.label
  num_taxa    = length(taxa)
  num_extinct = sum(ages[1:num_taxa,]$ages > 0)
  num_extant  = num_taxa - num_extinct

  # compute the factor
  f = (num_taxa - 1) * log(2) - lfactorial(num_extant)

  return(f)

}

p_condition = function( condition, age, lambda, mu, phi, delta, rho ) {

  if ( delta != 0 ) {
    stop("Analytical conditions not implemented for delta > 0.")
  }

  if ( condition == "time" ) {
    return(0)
  } else if ( condition == "survivalRoot" ) {
    p_survival = 1 - p0t(age, lambda, mu, 0, rho)
    return(2 * log(p_survival))
  } else if ( condition == "MRCA" ) {
    p_survival = 1 - p0t(age, lambda, mu, 0, rho)
    return(log(lambda) + 2 * log(p_survival))
  } else if ( condition == "sampledRoot" ) {
    p_survival = 1 - p0t(age, lambda, mu, phi, rho)
    return(2 * log(p_survival))
  } else if ( condition == "survivalStem" ) {
    p_survival = 1 - p0t(age, lambda, mu, 0, rho)
    return(log(p_survival))
  } else if ( condition == "sampledStem" ) {
    p_survival = 1 - p0t(age, lambda, mu, phi, rho)
    return(log(p_survival))
  } else if ( condition == "twoExtant" ) {
    p_survival  = 1 - p0t(age, lambda, mu, 0, rho)
    p_singleton = p1t(age, lambda, mu, 0, rho)
    return(log(p_survival - p_singleton))
  } else if ( condition == "twoSampled" ) {
    p_survival  = 1 - p0t(age, lambda, mu, phi, rho)
    p_singleton = p1t(age, lambda, mu, phi, rho)
    return(log(p_survival - p_singleton))
  } else if ( condition == "sampledMRCA" ) {
    p_survival = 1 - p0t(age, lambda, mu, phi, rho)
    return(log(lambda) + 2 * log(p_survival))
  } else {
    stop("Invalid condition selected.")
  }

}

c1 = function(lambda, mu, phi, rho) {

  abs(sqrt((lambda - mu - phi)^2 + 4 * lambda * phi))

}

c2 = function(lambda, mu, phi, rho) {

  -(lambda - mu - 2 * lambda * rho - phi) / c1(lambda, mu, phi, rho)

}

p0t = function(t, lambda, mu, phi, rho) {

  c_one = c1(lambda, mu, phi, rho)
  c_two = c2(lambda, mu, phi, rho)

  # res = lambda + mu + phi
  # res = res + c_one * ( exp(-c_one * t) * (1 - c_two) - (1 + c_two) ) / ( exp(-c_one * t) * (1 - c_two) + (1 + c_two) )
  # res = res / (2 * lambda)

  res = (lambda + mu + phi + c_one * ( ( exp(-c_one * t) * (1 - c_two) - (1 + c_two) ) / ( exp(-c_one * t) * (1 - c_two) + (1 + c_two) ) )) / (2 * lambda)
  return(res)

}

p1t = function(t, lambda, mu, phi, rho) {

  c_one = c1(lambda, mu, phi, rho)
  c_two = c2(lambda, mu, phi, rho)

  res = 4 * rho / ( 2 * (1 - c_two^2) + exp(-c_one * t) * (1 - c_two)^2 + exp(c_one * t) * (1 + c_two)^2 )
  return(res)

}

qt = function(t, lambda, mu, phi, rho) {

  c_one = c1(lambda, mu, phi, rho)
  c_two = c2(lambda, mu, phi, rho)

  res = 2 * (1 - c_two^2) + exp(-c_one * t) * (1 - c_two)^2 + exp(c_one * t) * (1 + c_two)^2
  return(res)

}

# eq. 3 from Stadler 2010
fbd_likelihood_stem = function(tree, lambda, mu, phi, rho) {

  # compute the ages
  ages = tree.age(tree, digits=5)

  # compute the number of extinct and extant species
  taxa        = tree$tip.label
  num_taxa    = length(taxa)
  num_extinct = sum(ages[1:num_taxa,]$ages > 0)
  num_extant  = num_taxa - num_extinct

  # compute the internal node times
  node_ages = ages[-c(1:num_taxa),]$ages

  # compute the origin time
  origin_time = max(node_ages) + tree$root.edge

  # compute end times for fossils
  if ( num_extinct > 0 ) {
    end_times = ages$ages[which(ages[1:num_taxa,]$ages > 0)]
  } else {
    end_times = c()
  }

  # compute the probability density
  if ( phi > 0 ) {
    ll = (num_extant + num_extinct - 1) * log(lambda) + (num_extinct) * log(phi) + num_extant * log(4 * rho) - log(qt(origin_time, lambda, mu, phi, rho))
  } else {
    ll = (num_extant + num_extinct - 1) * log(lambda) + num_extant * log(4 * rho) - log(qt(origin_time, lambda, mu, phi, rho))
  }

  # the probability of nodes
  ll = ll - sum(log(qt(node_ages, lambda, mu, phi, rho)))

  # the probability of end times
  ll = ll + sum(log(p0t(end_times, lambda, mu, phi, rho) * qt(end_times, lambda, mu, phi, rho)))

  return(ll)

}

# eq. 5 from Stadler 2010
fbd_likelihood_crown = function(tree, lambda, mu, phi, rho) {

  # compute the ages
  ages = tree.age(tree, digits=5)

  # compute the number of extinct and extant species
  taxa        = tree$tip.label
  num_taxa    = length(taxa)
  num_extinct = sum(ages[1:num_taxa,]$ages > 0)
  num_extant  = num_taxa - num_extinct

  # compute the internal node times
  node_ages = ages[-c(1:num_taxa),]$ages

  # compute end times for fossils
  if ( num_extinct > 0 ) {
    end_times = ages$ages[1:num_extinct + num_extant]
  } else {
    end_times = c()
  }

  # compute the probability density
  if ( phi > 0 ) {
    ll = (num_extant + num_extinct - 1) * log(lambda) + (num_extinct) * log(phi) + log(p1t(max(node_ages), lambda, mu, phi, rho))
  } else {
    ll = (num_extant + num_extinct - 1) * log(lambda) + log(p1t(max(node_ages), lambda, mu, phi, rho))
  }

  # the probability of nodes
  ll = ll + sum(log(p1t(node_ages, lambda, mu, phi, rho)))

  # the probability of end times
  ll = ll + sum(log(p0t(end_times, lambda, mu, phi, rho) / p1t(end_times, lambda, mu, phi, rho)))

  return(ll)

}

#########################################
# functions for singleton probabilities #
#########################################

ds = function(s, lambda, mu, phi, delta, rho, t, dt) {
  u = p0t(t, lambda, mu, phi, rho)
  d = (-(lambda + mu + phi + delta) * s + 2 * lambda * s * u + phi * u + delta) * dt
  return(d)
}

st = function(lambda, mu, phi, delta = 0, rho, t, nt = 1001) {

  # time
  dt    = t / nt
  times = seq(0, t, dt)

  # initial values
  s = rho

  # integrate
  for(i in 2:nt) {
    d = ds(s, lambda, mu, phi, delta, rho, times[i], dt)
    s = s + d
  }

  return(s)

}

simulate_bdsp = function(lambda, mu, phi, delta, rho, t) {

  # initial number of lineages
  num_lineages = 1
  num_samples  = 0

  # rates
  total_rate = lambda + mu + phi + delta
  event_prob = c(lambda, mu, phi, delta) / total_rate

  # time
  current_time = t

  # simulate
  while ( num_lineages > 0 ) {

    # draw the next event time
    waiting_time = rexp(1, num_lineages * total_rate)

    # increment time
    current_time = current_time - waiting_time

    # if past the present, break
    if ( current_time < 0 ) {
      break
    }

    # choose the type of event
    event_type = sample.int(4, size=1, prob=event_prob)

    if ( event_type == 1 ) {
      # speciation
      num_lineages = num_lineages + 1
    } else if ( event_type == 2 ) {
      # extinction
      num_lineages = num_lineages - 1
    } else if ( event_type == 3 ) {
      # sampling
      num_samples = num_samples + 1
    } else if ( event_type == 4 ) {
      # destructive sampling
      num_lineages = num_lineages - 1
      num_samples  = num_samples + 1
    }

    # cat(sprintf("%.3f  %i  %i  %i ", current_time, num_lineages, num_samples, event_type), "\n")

  }

  # simulate sampling at the present
  if ( current_time < 0 ) {
    num_samples = num_samples + rbinom(1, size=num_lineages, prob=rho)
  }

  return(num_samples)

}





