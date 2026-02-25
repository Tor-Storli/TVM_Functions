-- ============================================================
-- Time Value of Money (TVM) Macros for DuckDB
-- Based on numpy-financial (https://github.com/numpy/numpy-financial)
-- ============================================================
-- Conventions (follow numpy-financial / spreadsheet sign convention):
--   rate  = interest rate per period (decimal, e.g. 0.05 for 5%)
--   nper  = total number of periods
--   pmt   = payment per period (negative = cash out)
--   pv    = present value     (positive = cash received, e.g. loan principal)
--   fv    = future value      (default 0)
--   when  = 0 end-of-period (default), 1 beginning-of-period
-- ============================================================


-- ============================================================
-- 1. FV  -  Future Value
-- ============================================================
-- Equation: fv + pv*(1+r)^n + pmt*(1+r*w)/r*((1+r)^n - 1) = 0
-- When r=0: fv = -(pv + pmt*n)
--
-- Example: fv(0.05/12, 10*12, -100, -100) -> 15692.929
-- ============================================================
CREATE OR REPLACE MACRO fv(rate, nper, pmt, pv, "when" := 0) AS
    CASE
        WHEN rate = 0
            THEN -(pv + pmt * nper)
        ELSE
              -pv  * POW(1 + rate, nper)
            - pmt  * (1 + rate * "when") / rate * (POW(1 + rate, nper) - 1)
    END;


-- ============================================================
-- 2. PV  -  Present Value
-- ============================================================
-- Equation solved for pv.
-- When r=0: pv = -(fv + pmt*n)
--
-- Example: pv(0.05/12, 10*12, -100, 15692.93) -> -100.001
-- ============================================================
CREATE OR REPLACE MACRO pv(rate, nper, pmt, fv := 0.0, "when" := 0) AS
    CASE
        WHEN rate = 0
            THEN -(fv + pmt * nper)
        ELSE
            -(fv + pmt * (1 + rate * "when") / rate * (POW(1 + rate, nper) - 1))
            / POW(1 + rate, nper)
    END;


-- ============================================================
-- 3. PMT  -  Periodic Payment
-- ============================================================
-- Equation solved for pmt.
-- When r=0: pmt = -(fv + pv) / n
--
-- Example: pmt(0.075/12, 12*15, 200000) -> -1854.025
-- ============================================================
CREATE OR REPLACE MACRO pmt(rate, nper, pv, fv := 0.0, "when" := 0) AS
    CASE
        WHEN rate = 0
            THEN -(fv + pv) / nper
        ELSE
            -(fv + pv * POW(1 + rate, nper))
            / ((1 + rate * "when") * (POW(1 + rate, nper) - 1) / rate)
    END;


-- ============================================================
-- 4. NPER  -  Number of Periods
-- ============================================================
-- Solved analytically via logarithms.
-- When r=0: nper = -(fv + pv) / pmt
--
-- Example: nper(0.07/12, -150, 8000) -> 64.073
-- ============================================================
CREATE OR REPLACE MACRO nper(rate, pmt, pv, fv := 0.0, "when" := 0) AS
    CASE
        WHEN rate = 0
            THEN -(fv + pv) / pmt
        ELSE
            LN(
                (pmt * (1 + rate * "when") / rate - fv)
                / (pmt * (1 + rate * "when") / rate + pv)
            )
            / LN(1 + rate)
    END;


-- ============================================================
-- 5. IPMT  -  Interest Portion of a Payment
-- ============================================================
-- ipmt(per) = remaining_balance(per-1) * rate
-- remaining_balance after k periods = -fv(rate, k, pmt, pv)
--
-- Special cases:
--   per < 1          -> NULL  (invalid period)
--   when=1 & per=1   -> 0.0   (no interest on first begin-of-period payment)
--   when=1 & per > 1 -> discount result by one period
--
-- Example: ipmt(0.0824/12, 1, 12, 2500) -> -17.167
-- ============================================================
CREATE OR REPLACE MACRO ipmt(rate, per, nper, pv, fv := 0.0, "when" := 0) AS
    CASE
        WHEN per < 1
            THEN NULL
        WHEN "when" = 1 AND per = 1
            THEN 0.0
        WHEN "when" = 1 AND per > 1
            THEN
                -fv(rate, per - 1, pmt(rate, nper, pv, fv, "when"), pv, "when")
                * rate / (1 + rate)
        ELSE
            -fv(rate, per - 1, pmt(rate, nper, pv, fv, "when"), pv, "when")
            * rate
    END;


-- ============================================================
-- 6. PPMT  -  Principal Portion of a Payment
-- ============================================================
-- ppmt = total_pmt - ipmt
--
-- Example: ppmt(0.0824/12, 1, 12, 2500) -> -200.582
-- ============================================================
CREATE OR REPLACE MACRO ppmt(rate, per, nper, pv, fv := 0.0, "when" := 0) AS
    pmt(rate, nper, pv, fv, "when")
    - ipmt(rate, per, nper, pv, fv, "when");


-- ============================================================
-- 7. NPV  -  Net Present Value  (scalar macro, takes a DOUBLE[])
-- ============================================================
-- NPV = SUM( cf[t] / (1+rate)^t )  for t = 0 .. N-1
--
-- Pass cashflows as a typed array literal, e.g.:
--   SELECT npv(0.08, [-40000.0, 5000.0, 8000.0, 12000.0, 30000.0]);
--
-- Example: -> 3065.223
-- ============================================================
CREATE OR REPLACE MACRO npv(rate, cashflows) AS (
    SELECT SUM(v / POW(1 + rate, idx))
    FROM (
        SELECT UNNEST(cashflows)                        AS v,
               GENERATE_SUBSCRIPTS(cashflows, 1) - 1   AS idx
    ) t
);


-- ============================================================
-- 8. MIRR  -  Modified Internal Rate of Return  (table macro)
-- ============================================================
-- Matches numpy-financial formula:
--   numer = NPV of positive cashflows at reinvest_rate
--   denom = NPV of negative cashflows at finance_rate (absolute value)
--   MIRR  = (numer / denom)^(1/(n-1)) * (1+reinvest_rate) - 1
--
-- Usage:
--   SELECT mirr FROM mirr([-100.0, 50.0, -60.0, 70.0], 0.10, 0.12);
--   -> -0.03909
-- ============================================================
CREATE OR REPLACE MACRO mirr(cashflows, fin_rate, reinv_rate) AS TABLE
WITH cf AS (
    SELECT
        UNNEST(cashflows)                       AS v,
        GENERATE_SUBSCRIPTS(cashflows, 1) - 1   AS t,
        ARRAY_LENGTH(cashflows, 1)              AS n
),
parts AS (
    SELECT
        MAX(n)                                                              AS n,
        SUM(CASE WHEN v > 0 THEN  v / POW(1 + reinv_rate, t) ELSE 0 END)  AS npv_pos,
        SUM(CASE WHEN v < 0 THEN  v / POW(1 + fin_rate,   t) ELSE 0 END)  AS npv_neg
    FROM cf
)
SELECT
    POW(npv_pos / ABS(npv_neg), 1.0 / (n - 1)) * (1 + reinv_rate) - 1  AS mirr
FROM parts;


-- ============================================================
-- 9. AMORTIZATION_SCHEDULE  -  Table Macro
-- ============================================================
-- Recursive CTE: each row advances balance by one period.
-- Balance recurrence (end-of-period):
--   balance[k] = balance[k-1] * (1 + rate) + pmt
--   interest[k] = balance[k-1] * rate
--   principal[k] = pmt - interest[k]
--
-- Columns returned:
--   period | payment | interest | principal | remaining_balance
--
-- Usage:
--   SELECT * FROM amortization_schedule(0.075/12, 180, 200000);
-- ============================================================
CREATE OR REPLACE MACRO amortization_schedule(rate, nper, pv, fv := 0.0, "when" := 0)
AS TABLE
WITH RECURSIVE sched AS (

    -- Seed: period 0 carries opening balance.
    -- Interest and principal are computed in the recursive step
    -- directly from prev_balance, so LAG is never needed.
    SELECT
        0                                            AS period,
        CAST(pv AS DOUBLE)                           AS prev_balance,
        CAST(pv AS DOUBLE)                           AS balance,
        pmt(rate, nper, pv, fv, "when")              AS _pmt,
        CAST(rate AS DOUBLE)                         AS _rate,
        CAST(nper AS INTEGER)                        AS _nper,
        CAST(0.0 AS DOUBLE)                          AS interest,
        CAST(0.0 AS DOUBLE)                          AS principal

    UNION ALL

    SELECT
        period + 1,
        balance                                      AS prev_balance,
        balance * (1.0 + _rate) + _pmt               AS balance,
        _pmt,
        _rate,
        _nper,
        -- Interest = opening balance * rate  (positive)
        balance * _rate                              AS interest,
        -- Principal = |payment| - interest   (positive)
        ABS(_pmt) - balance * _rate                  AS principal

    FROM sched
    WHERE period < _nper
)

SELECT
    period,
    CAST(ROUND(ABS(_pmt),                                2) AS DOUBLE)  AS payment,
    CAST(ROUND(interest,                                 2) AS DOUBLE)  AS interest,
    CAST(ROUND(principal,                                2) AS DOUBLE)  AS principal,
    CAST(ROUND(SUM(interest)  OVER (ORDER BY period),   2) AS DOUBLE)  AS cum_interest,
    CAST(ROUND(SUM(principal) OVER (ORDER BY period),   2) AS DOUBLE)  AS cum_principal,
    CAST(ROUND(balance,                                  2) AS DOUBLE)  AS remaining_balance

FROM sched
WHERE period > 0
ORDER BY period;


-- ============================================================
-- QUICK SMOKE-TEST QUERIES
-- ============================================================

-- Scalar macros:
-- SELECT fv(0.05/12, 10*12, -100, -100)                              AS fv;    -- 15692.929
-- SELECT pv(0.05/12, 10*12, -100, 15692.93)                          AS pv;    -- -100.001
-- SELECT pmt(0.075/12, 12*15, 200000)                                AS pmt;   -- -1854.025
-- SELECT nper(0.07/12, -150, 8000)                                   AS nper;  -- 64.073
-- SELECT ipmt(0.0824/12, 1, 12, 2500)                                AS ipmt;  -- -17.167
-- SELECT ppmt(0.0824/12, 1, 12, 2500)                                AS ppmt;  -- -200.582
-- SELECT npv(0.08, [-40000.0, 5000.0, 8000.0, 12000.0, 30000.0])    AS npv;   -- 3065.223

-- Table macros:
-- SELECT * FROM mirr([-100.0, 50.0, -60.0, 70.0], 0.10, 0.12);                -- -0.039
-- SELECT * FROM amortization_schedule(0.0325/12, 180, 350000) LIMIT 16;
