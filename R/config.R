# Study parameters. Concept IDs are OMOP standard concepts; ingredient-level
# IDs pick up every product containing the drug via CONCEPT_ANCESTOR.

study <- list(
  target_concept_id     = 1118084,  # celecoxib (ingredient)
  comparator_concept_id = 1124300,  # diclofenac (ingredient)
  outcome_concept_id    = 192671,   # gastrointestinal hemorrhage (and descendants)

  washout_days = 365,               # prior observation required before index
  min_age      = 18,

  # Covariates are assessed from day -365 through day 0 relative to index.
  covariate_window         = c(-365, 0),
  min_covariate_prevalence = 0.01,

  # Time-at-risk windows in days after index. Exposure end dates in this
  # source are not usable (end = start for every record), so follow-up uses
  # fixed windows instead of an on-treatment definition.
  # The intention-to-treat window runs to the end of observation; it is also
  # the window for negative controls, which are too rare in 90 days.
  tar_primary     = c(1, 90),
  tar_sensitivity = list(c(1, 30)),
  tar_itt         = c(1, 36500),

  caliper  = 0.2,                   # SD of the logit propensity score
  min_cell = 11,                    # small-cell suppression threshold for shared tables

  # Negative control outcomes: conditions with no plausible causal link to
  # either NSAID, used to estimate residual systematic error.
  negative_controls = c(
    81151,    # sprain of ankle
    78272,    # sprain of wrist
    378001,   # concussion with no loss of consciousness
    4218389,  # whiplash injury to neck
    4113008,  # laceration of hand
    4156265,  # facial laceration
    4155034,  # laceration of forearm
    4152936,  # laceration of thigh
    4109685,  # laceration of foot
    4278672,  # fracture of forearm
    4237458,  # fracture of clavicle
    4059173,  # fracture of ankle
    4142905,  # fracture of rib
    4296204,  # first degree burn
    4296205,  # second degree burn
    372328,   # otitis media
    28060,    # streptococcal sore throat
    4116491,  # Escherichia coli urinary tract infection
    195588    # cystitis
  )
)

dataset_url <- "https://raw.githubusercontent.com/OHDSI/EunomiaDatasets/main/datasets/GiBleed/GiBleed_5.3.zip"
