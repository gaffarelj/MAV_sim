-- SQLite


-- SELECT COUNT(*) FROM solutions_multi_fin WHERE h_p_score IS NULL AND dv_used = "all" ORDER BY id DESC;

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE h_p_score IS NULL AND dv_used = "init_angle_only" ORDER BY id DESC;

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE h_p_score IS NULL AND dv_used = "TVC_only" ORDER BY id DESC;

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE h_p_score IS NULL AND dv_used = "SRM_only" ORDER BY id DESC;

-- SELECT id, final_time, h_a/1e3, h_p/1e3, t_b_1, t_b_2, mass FROM solutions_multi_fin WHERE h_a > 0 AND h_p > 0 AND h_a < 800e3 AND h_p < 800e3 ORDER BY h_p DESC;

-- SELECT * FROM solutions_multi_fin WHERE h_a_score > 1e2 OR h_p_score > 1e2;

-- SELECT COUNT(*) FROM solutions_multi_fin;

-- SELECT id, h_a_score, h_p_score, mass_score FROM solutions_multi_fin ORDER BY h_a_score ASC LIMIT 10;

-- SELECT id, h_a_score, h_p_score, mass_score FROM solutions_multi_fin ORDER BY h_p_score ASC LIMIT 10;

-- SELECT id, h_a_score, h_p_score, mass_score FROM solutions_multi_fin ORDER BY mass_score ASC LIMIT 10;

-- SELECT id, h_a_score+h_p_score+mass_score, dv_used, h_a_score, h_p_score, mass_score FROM solutions_multi_fin ORDER BY h_a_score+h_p_score+mass_score ASC LIMIT 500;

-- SELECT * FROM solutions_multi_fin WHERE id=7559 OR id=13157;

-- SELECT h_a_score+h_p_score, * FROM solutions_multi_fin WHERE dv_used = 'SRM_only' ORDER BY h_a_score+h_p_score+0*mass_score ASC LIMIT 10;

-- SELECT id, dv_used, h_p_score, h_a_score, mass_score FROM solutions_multi_fin WHERE dv_used LIKE "%_tuning_%" ORDER BY id DESC;

-- DELETE FROM solutions_multi_fin  WHERE dv_used LIKE "%_tuning_%";

-- SELECT id, dv_used, * FROM solutions_multi_fin WHERE dv_used LIKE "%_tuning_%";

-- SELECT * FROM solutions_multi_fin WHERE dv_used LIKE "NSPSO_%_tuning_%";

-- SELECT id, dv_used, h_p_score, h_a_score, mass_score FROM solutions_multi_fin WHERE id = 56827;--dv_used = "NSGA2_42_tuning_1";

-- SELECT id, dv_used, h_p_score, h_a_score, mass_score FROM solutions_multi_fin WHERE h_a_score is NULL OR h_p_score is NULL OR mass_score is NULL;

-- SELECT * FROM solutions_multi_fin LIMIT 1;

-- SELECT id, h_p FROM solutions_multi_fin LIMIT 10--ORDER BY h_p_score+h_a_score+10*mass_score DESC LIMIT 1000 -- ORDER BY RANDOM() LIMIT 1000

-- DELETE FROM solutions_multi_fin  WHERE dv_used LIKE "%_pt_%";

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE dv_used LIKE "%_pt_%" AND h_p_score IS NOT NULL;

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE h_p_score IS NOT NULL AND angle_1 IS NOT NULL;

-- SELECT DISTINCT(dv_used) FROM solutions_multi_fin WHERE dv_used LIKE "%_pt_%" AND h_p_score IS NOT NULL ORDER BY ID DESC LIMIT 5;

-- SELECT DISTINCT(dv_used) FROM solutions_multi_fin WHERE dv_used LIKE "%_pt_1" AND h_p_score IS NOT NULL ORDER BY ID DESC;

-- SELECT * FROM solutions_multi_fin WHERE h_p_score = 86.76264100955248 AND h_a_score = 4.50567081948868 AND mass_score = 0.8031542570936888 AND angle_1 IS NOT NULL LIMIT 1;

-- SELECT COUNT(h_p_score), * FROM solutions_multi_fin WHERE angle_1 IS NOT NULL GROUP BY h_p_score HAVING COUNT(h_p_score) > 1

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE h_p_score IS NOT NULL AND angle_1 IS NOT NULL;

-- SELECT dv_used FROM solutions_multi_fin WHERE dv_used LIKE "%_pt_%";

-- SELECT dv_used, h_p_score, angle_1 FROM solutions_multi_fin WHERE dv_used LIKE "opti_%";

-- SELECT id, h_a_score, h_p_score, mass_score FROM solutions_multi_fin WHERE h_a_score + h_p_score < 1.5 AND mass_score < 1.5;

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE dv_used LIKE "refinement3_%%" AND h_a_score IS NOT NULL;

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE dv_used LIKE "refinement3_%%" AND h_a_score+h_p_score = 0;

-- SELECT id, dv_used, mass_score, h_a_score+h_p_score FROM solutions_multi_fin WHERE dv_used LIKE "refinement3_%%" AND h_a_score+h_p_score = 0 ORDER BY mass_score ASC LIMIT 3;

-- SELECT id, dv_used, mass_score, h_a_score+h_p_score FROM solutions_multi_fin WHERE dv_used NOT LIKE "refinement3_%%" AND h_a_score+h_p_score = 0 ORDER BY mass_score ASC LIMIT 3;

-- SELECT id, dv_used, mass_score, h_a_score+h_p_score FROM solutions_multi_fin WHERE dv_used LIKE "refinement3_%%" AND h_a_score IS NOT NULL ORDER BY id DESC LIMIT 3;

-- SELECT DISTINCT dv_used FROM solutions_multi_fin WHERE dv_used LIKE "refinement%" LIMIT 3;

-- SELECT DISTINCT dv_used FROM solutions_multi_fin WHERE dv_used LIKE "refinement2_%" LIMIT 3;

-- SELECT DISTINCT dv_used FROM solutions_multi_fin WHERE dv_used LIKE "refinement3_%" LIMIT 3;

-- SELECT * FROM solutions_multi_fin WHERE angle_1 IS NOT NULL and h_a_score IS NOT NULL AND h_a_score + h_p_score = 0 AND mass_score < 1.0 ORDER BY mass_score ASC LIMIT 5;

-- SELECT id, h_a, h_p FROM solutions_multi_fin WHERE dv_used = "extra_accurate";

-- DELETE FROM solutions_multi_fin WHERE dv_used = "extra_accurate";

-- SELECT h_a, h_p FROM solutions_multi_fin WHERE angle_1 IS NOT NULL AND h_a_score + h_p_score = 0 ORDER BY mass ASC LIMIT 20;

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE h_a > 200e3 AND h_p > 200e3 AND h_a < 500e3 AND h_p < 500e3 AND dv_used = "TVC_only";

-- SELECT COUNT(*) FROM solutions_multi_fin WHERE dv_used = "NSGA2_42_tuning_14";

-- SELECT * FROM solutions_multi_fin WHERE h_p_score + h_a_score = 0 ORDER BY mass_score ASC LIMIT 1;
-- SELECT * FROM solutions_multi_fin WHERE h_p_score + h_a_score = 0 ORDER BY mass_score DESC LIMIT 1;

-- SELECT * FROM solutions_multi_fin WHERE h_a IS NOT NULL AND h_a_score + h_p_score = 0 ORDER BY t_b_1 + t_b_2 ASC LIMIT 1;
-- SELECT * FROM solutions_multi_fin WHERE h_a IS NOT NULL AND h_a_score + h_p_score = 0 ORDER BY t_b_1 ASC LIMIT 1;
-- SELECT * FROM solutions_multi_fin WHERE h_a IS NOT NULL AND h_a_score + h_p_score = 0 ORDER BY t_b_2 ASC LIMIT 1;

-- SELECT COUNT(*) FROM solutions_anchor WHERE angle_1 IS NOT NULL AND h_a_score IS NULL;
-- SELECT COUNT(*) FROM solutions_anchor WHERE angle_1 IS NOT NULL AND h_a_score IS NOT NULL;
-- SELECT DISTINCT(dv_used) FROM solutions_anchor;
-- SELECT id, h_a_score, h_p_score, h_a, h_p, mass, dv_used FROM solutions_anchor WHERE h_a_score IS NOT NULL ORDER BY h_a_score + h_p_score ASC LIMIT 5;

-- SELECT "multi_fin", COUNT(*) FROM solutions_multi_fin;

-- SELECT "tubular", COUNT(*) FROM solutions_tubular;
-- SELECT DISTINCT(dv_used) FROM solutions_tubular ORDER BY id DESC LIMIT 1;
-- SELECT id, h_a_score, h_p_score, mass, dv_used, h_a, h_p FROM solutions_tubular WHERE h_a_score IS NOT NULL ORDER BY h_a_score + h_p_score + mass_score ASC LIMIT 5;

-- SELECT "rod_and_tube", COUNT(*) FROM solutions_rod_and_tube;

-- SELECT "anchor", COUNT(*) FROM solutions_anchor;

-- DELETE FROM solutions_multi_fin WHERE dv_used LIKE "extra%";
-- DELETE FROM solutions_tubular WHERE dv_used LIKE "extra%";
-- DELETE FROM solutions_rod_and_tube WHERE dv_used LIKE "extra%";
-- DELETE FROM solutions_anchor WHERE dv_used LIKE "extra%";

SELECT COUNT(*) from solutions_multi_fin WHERE dv_used LIKE "extra%" AND h_a_score IS NOT NULL;
SELECT COUNT(*) from solutions_tubular WHERE dv_used LIKE "extra%" AND h_a_score IS NOT NULL;
SELECT COUNT(*) from solutions_rod_and_tube WHERE dv_used LIKE "extra%" AND h_a_score IS NOT NULL;
SELECT COUNT(*) from solutions_anchor WHERE dv_used LIKE "extra%" AND h_a_score IS NOT NULL;

-- SELECT * FROM solutions_multi_fin WHERE h_a_score = 0 AND h_p_score = 0 ORDER BY mass_score ASC LIMIT 1;
-- SELECT * FROM solutions_tubular WHERE h_a_score = 0 AND h_p_score = 0 ORDER BY mass_score ASC LIMIT 1;
-- SELECT * FROM solutions_rod_and_tube WHERE h_a_score = 0 AND h_p_score = 0 ORDER BY mass_score ASC LIMIT 1;
-- SELECT * FROM solutions_anchor WHERE h_a_score = 0 AND h_p_score = 0 ORDER BY mass_score ASC LIMIT 1;
-- SELECT * FROM solutions_anchor WHERE h_a_score = 0 AND h_p > 305e3 ORDER BY mass_score ASC LIMIT 1;

-- SELECT * FROM solutions_anchor WHERE h_a_score + h_p_score < 0.7 AND mass_score < 1.0 AND h_a_score IS NOT NULL and angle_1 IS NOT NULL AND dv_used LIKE '%_tuning_%'

-- SELECT * FROM solutions_multi_fin WHERE h_a_score = 0 AND h_p > 305e3 AND mass_score <= 1 ORDER by t_b_2 ASC LIMIT 1;
-- SELECT t_b_2 FROM solutions_tubular WHERE h_a_score = 0 AND h_p > 305e3 AND mass_score <= 1 ORDER by t_b_2 ASC LIMIT 1;
-- SELECT t_b_2 FROM solutions_rod_and_tube WHERE h_a_score = 0 AND h_p > 305e3 AND mass_score <= 1 ORDER by t_b_2 ASC LIMIT 1;
-- SELECT t_b_2 FROM solutions_anchor WHERE h_a_score = 0 AND h_p > 305e3 AND mass_score <= 1 ORDER by t_b_2 ASC LIMIT 1;

-- SELECT * FROM solutions_rod_and_tube WHERE h_a_score + h_p_score = 0 ORDER BY mass_score ASC LIMIT 1;