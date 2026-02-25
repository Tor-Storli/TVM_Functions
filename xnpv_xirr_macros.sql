-- ============================================================
-- XNPV and XIRR Macros for DuckDB  (fixed)
-- Companion to tvm_macros.sql
-- ============================================================
-- XNPV and XIRR operate on irregular (date-based) cash flows.
-- Cash flow dates do NOT need to be evenly spaced.
--
-- Day-count convention: Actual/365 (matches Excel, numpy-financial)
--   t[i] = (date[i] - date[0]) / 365.0
--
-- Sign convention: negative = cash out, positive = cash in.
-- The first cash flow is typically the initial investment (negative).
-- ============================================================


-- ============================================================
-- 1. XNPV  -  Net Present Value for irregular cash flows
-- ============================================================
-- Formula:
--   XNPV = SUM( cf[i] / (1 + rate)^((date[i] - date[0]) / 365) )
--
-- Reference date (t0) is isolated in its own CTE to avoid
-- mixing window functions inside aggregates.
--
-- Parameters:
--   rate      DOUBLE   -- discount rate (e.g. 0.10 for 10%)
--   cashflows DOUBLE[] -- cash flow amounts
--   dates     DATE[]   -- corresponding dates (same length, ascending)
--
-- Usage:
--   SELECT xnpv(
--       0.10,
--       [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--       ['2025-01-01','2025-07-01','2026-01-01',
--        '2026-07-01','2027-01-01','2027-07-01']::DATE[]
--   ) AS xnpv;
-- ============================================================
CREATE OR REPLACE MACRO xnpv(rate, cashflows, dates) AS (
    WITH unpacked AS (
        SELECT
            UNNEST(cashflows)  AS cf,
            UNNEST(dates)      AS d
    ),
    ref AS (
        SELECT MIN(d) AS t0 FROM unpacked
    )
    SELECT SUM(
        cf / POW(1.0 + rate, DATEDIFF('day', ref.t0, d) / 365.0)
    )
    FROM unpacked, ref
);


-- ============================================================
-- 2. XIRR  -  Internal Rate of Return for irregular cash flows
-- ============================================================
-- XIRR is the rate `r` such that XNPV(r, cashflows, dates) = 0.
--
-- Uses the same chained CTE Newton-Raphson pattern as irr()
-- to avoid DuckDB's recursive CTE limitation inside macros.
--
-- Each step carries r, f=XNPV(r), and fp=XNPV'(r) forward
-- so no correlated subquery re-scanning occurs at the end.
--
-- Formula:
--   f(r)  = SUM( cf[i] / (1+r)^t[i] )
--   f'(r) = SUM( -t[i] * cf[i] / (1+r)^(t[i]+1) )
--   r_next = r - f(r) / f'(r)
--
-- Parameters:
--   cashflows  DOUBLE[]  -- cash flow amounts
--   dates      DATE[]    -- corresponding dates (ascending)
--   guess      DOUBLE    -- starting rate estimate (default 0.1)
--   tol        DOUBLE    -- convergence tolerance  (default 1e-7)
--
-- Returns:
--   xirr       DOUBLE    -- solved rate (NULL if not converged)
--   iterations INTEGER   -- steps taken
--   converged  BOOLEAN   -- true if ABS(XNPV) < tol
--
-- Usage:
--   SELECT * FROM xirr(
--       [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--       ['2025-01-01','2025-07-01','2026-01-01',
--        '2026-07-01','2027-01-01','2027-07-01']::DATE[]
--   );
-- ============================================================
CREATE OR REPLACE MACRO xirr(
    cashflows,
    dates,
    guess := 0.1,
    tol   := 1e-7
) AS TABLE

WITH

-- Unpack and compute time fractions once
unpacked AS (
    SELECT
        UNNEST(cashflows)  AS cf,
        UNNEST(dates)      AS d
),
ref AS (
    SELECT MIN(d) AS t0 FROM unpacked
),
cf AS (
    SELECT
        cf,
        DATEDIFF('day', ref.t0, d) / 365.0  AS t
    FROM unpacked, ref
),

-- Seed: evaluate f and fp at initial guess
s00 AS (
    SELECT
        CAST(guess AS DOUBLE)                                    AS r,
        (SELECT SUM(cf / POW(1.0 + guess, t))          FROM cf) AS f,
        (SELECT SUM(-t * cf / POW(1.0 + guess, t + 1)) FROM cf) AS fp
),

-- 16 chained Newton-Raphson steps.
-- Each step carries r, f, fp forward -- no re-scan of cf at the end.
s01 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s00.r-s00.f/NULLIF(s00.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s00.r-s00.f/NULLIF(s00.fp,0)), t+1)) FROM cf) AS fp FROM s00),
s02 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s01.r-s01.f/NULLIF(s01.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s01.r-s01.f/NULLIF(s01.fp,0)), t+1)) FROM cf) AS fp FROM s01),
s03 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s02.r-s02.f/NULLIF(s02.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s02.r-s02.f/NULLIF(s02.fp,0)), t+1)) FROM cf) AS fp FROM s02),
s04 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s03.r-s03.f/NULLIF(s03.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s03.r-s03.f/NULLIF(s03.fp,0)), t+1)) FROM cf) AS fp FROM s03),
s05 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s04.r-s04.f/NULLIF(s04.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s04.r-s04.f/NULLIF(s04.fp,0)), t+1)) FROM cf) AS fp FROM s04),
s06 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s05.r-s05.f/NULLIF(s05.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s05.r-s05.f/NULLIF(s05.fp,0)), t+1)) FROM cf) AS fp FROM s05),
s07 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s06.r-s06.f/NULLIF(s06.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s06.r-s06.f/NULLIF(s06.fp,0)), t+1)) FROM cf) AS fp FROM s06),
s08 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s07.r-s07.f/NULLIF(s07.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s07.r-s07.f/NULLIF(s07.fp,0)), t+1)) FROM cf) AS fp FROM s07),
s09 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s08.r-s08.f/NULLIF(s08.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s08.r-s08.f/NULLIF(s08.fp,0)), t+1)) FROM cf) AS fp FROM s08),
s10 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s09.r-s09.f/NULLIF(s09.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s09.r-s09.f/NULLIF(s09.fp,0)), t+1)) FROM cf) AS fp FROM s09),
s11 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s10.r-s10.f/NULLIF(s10.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s10.r-s10.f/NULLIF(s10.fp,0)), t+1)) FROM cf) AS fp FROM s10),
s12 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s11.r-s11.f/NULLIF(s11.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s11.r-s11.f/NULLIF(s11.fp,0)), t+1)) FROM cf) AS fp FROM s11),
s13 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s12.r-s12.f/NULLIF(s12.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s12.r-s12.f/NULLIF(s12.fp,0)), t+1)) FROM cf) AS fp FROM s12),
s14 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s13.r-s13.f/NULLIF(s13.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s13.r-s13.f/NULLIF(s13.fp,0)), t+1)) FROM cf) AS fp FROM s13),
s15 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s14.r-s14.f/NULLIF(s14.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s14.r-s14.f/NULLIF(s14.fp,0)), t+1)) FROM cf) AS fp FROM s14),
s16 AS (SELECT r - f/NULLIF(fp,0) AS r,
    (SELECT SUM(cf / POW(1.0+(s15.r-s15.f/NULLIF(s15.fp,0)), t))      FROM cf) AS f,
    (SELECT SUM(-t*cf / POW(1.0+(s15.r-s15.f/NULLIF(s15.fp,0)), t+1)) FROM cf) AS fp FROM s15),

-- Final step: advance once more, use last_f (already computed) for convergence check
result AS (
    SELECT
        r - f/NULLIF(fp,0)  AS r,
        f                   AS last_f
    FROM s16
)

SELECT
    CASE WHEN ABS(last_f) < tol THEN ROUND(r, 10) ELSE NULL END  AS xirr,
    16                                                             AS iterations,
    ABS(last_f) < tol                                             AS converged
FROM result;


-- ============================================================
-- USAGE EXAMPLES
-- ============================================================

-- -- XNPV: $1000 investment, five $250 returns every ~6 months at 10%
-- SELECT xnpv(
--     0.10,
--     [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--     ['2025-01-01','2025-07-01','2026-01-01',
--      '2026-07-01','2027-01-01','2027-07-01']::DATE[]
-- ) AS xnpv;

-- -- XIRR: find the rate that makes XNPV = 0
-- SELECT * FROM xirr(
--     [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--     ['2025-01-01','2025-07-01','2026-01-01',
--      '2026-07-01','2027-01-01','2027-07-01']::DATE[]
-- );

-- -- Verify: XNPV at the XIRR rate should be ~0
-- SELECT xnpv(
--     (SELECT xirr FROM xirr(
--         [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--         ['2025-01-01','2025-07-01','2026-01-01',
--          '2026-07-01','2027-01-01','2027-07-01']::DATE[]
--     )),
--     [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--     ['2025-01-01','2025-07-01','2026-01-01',
--      '2026-07-01','2027-01-01','2027-07-01']::DATE[]
-- ) AS xnpv_at_xirr;   -- should be ~0

-- -- Classic Excel XIRR example (~37.34%):
-- SELECT * FROM xirr(
--     [-10000.0, 2750.0, 4250.0, 3250.0, 2750.0],
--     ['2008-01-01','2008-03-01','2008-10-30',
--      '2009-02-15','2009-04-01']::DATE[]
-- );

-- -- Formatted output:
-- SELECT
--     format('${:,.2f}', xnpv(
--         0.10,
--         [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--         ['2025-01-01','2025-07-01','2026-01-01',
--          '2026-07-01','2027-01-01','2027-07-01']::DATE[]
--     )) AS xnpv,
--     format('{:.4f}%', (SELECT xirr * 100 FROM xirr(
--         [-1000.0, 250.0, 250.0, 250.0, 250.0, 250.0],
--         ['2025-01-01','2025-07-01','2026-01-01',
--          '2026-07-01','2027-01-01','2027-07-01']::DATE[]
--     ))) AS xirr_pct;
