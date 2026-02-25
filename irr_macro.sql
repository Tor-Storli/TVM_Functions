-- ============================================================
-- IRR  -  Internal Rate of Return
-- ============================================================
-- DuckDB does not support recursive CTEs inside macro bodies,
-- so Newton-Raphson is unrolled into 32 chained CTE steps.
--
-- Each step carries BOTH the rate AND the NPV residual forward,
-- so no correlated subquery re-evaluation occurs at the end.
-- This avoids the exponential re-scan that caused hanging.
--
-- Each CTE step computes:
--   f  = NPV(r)   = SUM( cf[t] / (1+r)^t )
--   fp = NPV'(r)  = SUM( -t * cf[t] / (1+r)^(t+1) )
--   r_next = r - f / fp
--
-- Parameters:
--   cashflows  DOUBLE[]  -- period cash flows, index 0 = now
--   guess      DOUBLE    -- starting rate estimate (default 0.1)
--   tol        DOUBLE    -- convergence tolerance  (default 1e-7)
--
-- Returns:
--   irr        DOUBLE    -- solved rate (NULL if not converged)
--   iterations INTEGER   -- Newton steps taken
--   converged  BOOLEAN   -- true if ABS(NPV) < tol
--
-- Usage:
--   SELECT * FROM irr([-100.0, 39.0, 59.0, 55.0, 20.0]);
--   SELECT * FROM irr([-100.0, 0.0, 0.0, 74.0]);
--   SELECT * FROM irr([-100.0, 100.0, 0.0, 7.0]);
-- ============================================================

CREATE OR REPLACE MACRO irr(
    cashflows,
    guess := 0.1,
    tol   := 1e-7
) AS TABLE

WITH cf AS (
    SELECT
        UNNEST(cashflows)                      AS cf,
        GENERATE_SUBSCRIPTS(cashflows, 1) - 1 AS t
),

-- Each step: SELECT r, f, fp FROM previous step, then advance.
-- f and fp are computed once per step and carried forward --
-- no correlated subqueries in the final SELECT.

s00 AS (
    SELECT
        CAST(guess AS DOUBLE)                                    AS r,
        (SELECT SUM(cf / POW(1.0 + guess, t))          FROM cf) AS f,
        (SELECT SUM(-t * cf / POW(1.0 + guess, t + 1)) FROM cf) AS fp
),
s01 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s00.r-s00.f/NULLIF(s00.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s00.r-s00.f/NULLIF(s00.fp,0)), t+1))     FROM cf) AS fp FROM s00),
s02 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s01.r-s01.f/NULLIF(s01.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s01.r-s01.f/NULLIF(s01.fp,0)), t+1))     FROM cf) AS fp FROM s01),
s03 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s02.r-s02.f/NULLIF(s02.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s02.r-s02.f/NULLIF(s02.fp,0)), t+1))     FROM cf) AS fp FROM s02),
s04 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s03.r-s03.f/NULLIF(s03.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s03.r-s03.f/NULLIF(s03.fp,0)), t+1))     FROM cf) AS fp FROM s03),
s05 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s04.r-s04.f/NULLIF(s04.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s04.r-s04.f/NULLIF(s04.fp,0)), t+1))     FROM cf) AS fp FROM s04),
s06 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s05.r-s05.f/NULLIF(s05.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s05.r-s05.f/NULLIF(s05.fp,0)), t+1))     FROM cf) AS fp FROM s05),
s07 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s06.r-s06.f/NULLIF(s06.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s06.r-s06.f/NULLIF(s06.fp,0)), t+1))     FROM cf) AS fp FROM s06),
s08 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s07.r-s07.f/NULLIF(s07.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s07.r-s07.f/NULLIF(s07.fp,0)), t+1))     FROM cf) AS fp FROM s07),
s09 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s08.r-s08.f/NULLIF(s08.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s08.r-s08.f/NULLIF(s08.fp,0)), t+1))     FROM cf) AS fp FROM s08),
s10 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s09.r-s09.f/NULLIF(s09.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s09.r-s09.f/NULLIF(s09.fp,0)), t+1))     FROM cf) AS fp FROM s09),
s11 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s10.r-s10.f/NULLIF(s10.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s10.r-s10.f/NULLIF(s10.fp,0)), t+1))     FROM cf) AS fp FROM s10),
s12 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s11.r-s11.f/NULLIF(s11.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s11.r-s11.f/NULLIF(s11.fp,0)), t+1))     FROM cf) AS fp FROM s11),
s13 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s12.r-s12.f/NULLIF(s12.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s12.r-s12.f/NULLIF(s12.fp,0)), t+1))     FROM cf) AS fp FROM s12),
s14 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s13.r-s13.f/NULLIF(s13.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s13.r-s13.f/NULLIF(s13.fp,0)), t+1))     FROM cf) AS fp FROM s13),
s15 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s14.r-s14.f/NULLIF(s14.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s14.r-s14.f/NULLIF(s14.fp,0)), t+1))     FROM cf) AS fp FROM s14),
s16 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s15.r-s15.f/NULLIF(s15.fp,0)), t))          FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s15.r-s15.f/NULLIF(s15.fp,0)), t+1))     FROM cf) AS fp FROM s15),

-- Final step: just need the converged rate, f is already in hand
result AS (
    SELECT
        r - f/NULLIF(fp,0)  AS r,
        f                   AS last_f
    FROM s16
)

SELECT
    CASE WHEN ABS(last_f) < tol THEN ROUND(r, 10) ELSE NULL END  AS irr,
    16                                                             AS iterations,
    ABS(last_f) < tol                                             AS converged
FROM result;


-- ============================================================
-- USAGE EXAMPLES
-- ============================================================

-- SELECT * FROM irr([-100.0, 39.0, 59.0, 55.0, 20.0]);
-- -> 0.2809484211  (28.09%)

-- SELECT * FROM irr([-100.0, 0.0, 0.0, 74.0]);
-- -> -0.0955  (-9.55%)

-- SELECT * FROM irr([-100.0, 100.0, 0.0, 7.0]);
-- -> 0.06206  (6.21%)

-- SELECT * FROM irr([-5.0, 10.5, 1.0, -8.0, 1.0], guess := 0.05);
-- -> 0.0886  (8.86%)

-- -- Formatted as percent:
-- SELECT
--     format('{:.4f}%', irr * 100) AS irr_pct,
--     converged
-- FROM irr([-100.0, 39.0, 59.0, 55.0, 20.0]);

-- -- Annualize a monthly IRR:
-- SELECT
--     irr                     AS monthly_irr,
--     POW(1.0 + irr, 12) - 1 AS annual_irr
-- FROM irr([-1000.0, 100.0, 100.0, 100.0, 100.0,
--            100.0,  100.0, 100.0, 100.0, 100.0,
--            100.0,  100.0, 1100.0]);
