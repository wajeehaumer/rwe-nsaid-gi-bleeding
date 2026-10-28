-- Follow-up for one outcome and one time-at-risk window, per person in the
-- study population. Follow-up is censored at the end of the window, the end
-- of observation, or death, whichever comes first. Reused for the primary
-- outcome and every negative control outcome.
WITH outcome_events AS (
  SELECT co.person_id, co.condition_start_date AS event_date
  FROM condition_occurrence co
  JOIN concept_ancestor ca
    ON ca.descendant_concept_id = co.condition_concept_id
  WHERE ca.ancestor_concept_id = @outcome_concept_id
),
windows AS (
  SELECT person_id,
         index_date,
         index_date + @tar_start                                   AS tar_start,
         LEAST(index_date + @tar_end, obs_end_date,
               COALESCE(death_date, obs_end_date))                 AS tar_end
  FROM study_population
)
SELECT w.person_id,
       w.index_date,
       w.tar_start,
       w.tar_end,
       COALESCE(BOOL_OR(o.event_date <= w.index_date), FALSE)        AS prior_outcome,
       MIN(o.event_date) FILTER (WHERE o.event_date BETWEEN w.tar_start AND w.tar_end) AS outcome_date
FROM windows w
LEFT JOIN outcome_events o
  ON o.person_id = w.person_id
GROUP BY w.person_id, w.index_date, w.tar_start, w.tar_end;
