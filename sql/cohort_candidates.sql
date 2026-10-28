-- New-user cohorts: one row per person per study drug at their first-ever
-- exposure, with the flags used to apply eligibility criteria in R so that
-- each exclusion step can be counted for the attrition table.
CREATE OR REPLACE TABLE cohort_candidates AS
WITH first_exposure AS (
  SELECT de.person_id,
         ca.ancestor_concept_id          AS ingredient_concept_id,
         MIN(de.drug_exposure_start_date) AS index_date
  FROM drug_exposure de
  JOIN concept_ancestor ca
    ON ca.descendant_concept_id = de.drug_concept_id
  WHERE ca.ancestor_concept_id IN (@target_concept_id, @comparator_concept_id)
  GROUP BY 1, 2
),
outcome_history AS (
  SELECT co.person_id, MIN(co.condition_start_date) AS first_outcome_date
  FROM condition_occurrence co
  JOIN concept_ancestor ca
    ON ca.descendant_concept_id = co.condition_concept_id
  WHERE ca.ancestor_concept_id = @outcome_concept_id
  GROUP BY 1
)
SELECT fe.person_id,
       CASE WHEN fe.ingredient_concept_id = @target_concept_id THEN 1 ELSE 0 END AS treatment,
       fe.index_date,
       op.observation_period_end_date                                    AS obs_end_date,
       TRY_CAST(d.death_date AS DATE)                                     AS death_date,
       date_diff('day', op.observation_period_start_date, fe.index_date)  AS prior_obs_days,
       year(fe.index_date) - p.year_of_birth                              AS age,
       COUNT(*) OVER (PARTITION BY fe.person_id) > 1                      AS used_both_drugs,
       COALESCE(oh.first_outcome_date <= fe.index_date, FALSE)            AS prior_outcome
FROM first_exposure fe
JOIN person p
  ON p.person_id = fe.person_id
LEFT JOIN observation_period op
  ON op.person_id = fe.person_id
 AND fe.index_date BETWEEN op.observation_period_start_date AND op.observation_period_end_date
LEFT JOIN death d
  ON d.person_id = fe.person_id
LEFT JOIN outcome_history oh
  ON oh.person_id = fe.person_id;
