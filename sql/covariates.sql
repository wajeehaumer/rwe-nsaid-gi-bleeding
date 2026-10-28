-- Baseline covariates in long format (person_id, covariate_id, covariate_name, value).
-- Conditions and drug ingredients are binary indicators for any record in the
-- covariate window; the study drugs themselves are excluded.
CREATE OR REPLACE TABLE covariates_long AS
SELECT DISTINCT
       sp.person_id,
       'condition_' || co.condition_concept_id      AS covariate_id,
       'Condition: ' || c.concept_name              AS covariate_name,
       1.0                                          AS value
FROM study_population sp
JOIN condition_occurrence co
  ON co.person_id = sp.person_id
 AND co.condition_start_date BETWEEN sp.index_date + (@window_start) AND sp.index_date + (@window_end)
JOIN concept c
  ON c.concept_id = co.condition_concept_id

UNION ALL

SELECT DISTINCT
       sp.person_id,
       'drug_' || ca.ancestor_concept_id,
       'Drug: ' || c.concept_name,
       1.0
FROM study_population sp
JOIN drug_exposure de
  ON de.person_id = sp.person_id
 AND de.drug_exposure_start_date BETWEEN sp.index_date + (@window_start) AND sp.index_date + (@window_end)
JOIN concept_ancestor ca
  ON ca.descendant_concept_id = de.drug_concept_id
JOIN concept c
  ON c.concept_id = ca.ancestor_concept_id
 AND c.concept_class_id = 'Ingredient'
WHERE ca.ancestor_concept_id NOT IN (@target_concept_id, @comparator_concept_id)

UNION ALL

SELECT sp.person_id,
       'visit_count',
       'Healthcare visits in prior year',
       COUNT(vo.visit_occurrence_id)
FROM study_population sp
LEFT JOIN visit_occurrence vo
  ON vo.person_id = sp.person_id
 AND vo.visit_start_date BETWEEN sp.index_date + (@window_start) AND sp.index_date - 1
GROUP BY sp.person_id;
