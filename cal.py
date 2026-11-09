#!/usr/bin/env python3
"""Bikram Samvat calendar from rules (no month tables, no hand-tuned cutoffs), like Unix `cal`.

Month-start convention (which civil day is day 1 of a solar month):
  sunrise    (default) day 1 = first day whose sunrise finds the Sun in the new sign.
                       Same principle as taking the tithi at sunrise.
  sunset     same, but tested at sunset
  noon       sankranti before local noon -> same day, else next day
  civil_date month starts on the civil date of the sankranti (midnight-based)
Eras (ayanamsa, tz, convention per BS-year range) can be set via --rules rules.json.
AD years are astronomical (1 BC = 0). Dates before 1582-10-15 are Julian.
"""
import argparse
import datetime
import json
import math
import sys
from dataclasses import dataclass, replace
from functools import lru_cache

import swisseph as swe

BS_MONTH_NAMES = ["Baishakh", "Jestha", "Ashadh", "Shrawan", "Bhadra", "Ashwin",
                  "Kartik", "Mangsir", "Poush", "Magh", "Falgun", "Chaitra"]
DAYS = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"]
AYANAMSAS = {"lahiri": "SIDM_LAHIRI", "raman": "SIDM_RAMAN",
             "krishnamurti": "SIDM_KRISHNAMURTI", "fagan": "SIDM_FAGAN_BRADLEY",
             "suryasiddhanta": "SIDM_SURYASIDDHANTA",
             "suryasiddhanta_msun": "SIDM_SURYASIDDHANTA_MSUN"}
CONVENTIONS = ["sunrise", "sunset", "noon", "civil_date"]
GREG_START_JDN = 2299161           # 1582-10-15
EPH = swe.FLG_MOSEPH               # built-in ephemeris, no data files needed
KTM = (85.324, 27.7172, 1300.0)    # lon, lat, alt (m)


# ============================ RULES ======================================
@dataclass(frozen=True)
class Rule:
    from_bs: int = -10**9
    ayanamsa: str = "lahiri"
    tz: float = 5.75               # hours from UT for the civil day (Nepal: +5:45)
    convention: str = "sunrise"


RULES = [Rule()]


def load_rules(path):
    """JSON list of era dicts; each inherits unspecified fields from the previous one."""
    global RULES
    with open(path) as f:
        eras = sorted(json.load(f), key=lambda d: d.get("from_bs", -10**9))
    base, out = Rule(), []
    for d in eras:
        base = replace(base, **d)
        if base.convention not in CONVENTIONS:
            sys.exit(f"unknown convention {base.convention!r}; choose from {CONVENTIONS}")
        out.append(base)
    RULES = out or [Rule()]
    month_starts.cache_clear()


def rule_for(bs_year):
    r = RULES[0]
    for x in RULES:
        if x.from_bs <= bs_year:
            r = x
    return r


# ============================ DATES ======================================
def jdn_from_ymd(y, m, d):
    cal = swe.GREG_CAL if (y, m, d) >= (1582, 10, 15) else swe.JUL_CAL
    return int(round(swe.julday(y, m, d, 12.0, cal)))


def ymd_from_jdn(jdn):
    cal = swe.GREG_CAL if jdn >= GREG_START_JDN else swe.JUL_CAL
    y, m, d, _ = swe.revjul(float(jdn), cal)
    return y, m, d


def fmt_ad(jdn):
    y, m, d = ymd_from_jdn(jdn)
    return f"{y if y > 0 else f'{1 - y} BC'}-{m:02d}-{d:02d}"


def dow(jdn):
    return (jdn + 1) % 7       # 0 = Sunday


# ============================ ASTRONOMY ==================================
def sun_sid(jd, ayan):
    swe.set_sid_mode(getattr(swe, AYANAMSAS[ayan]))
    pos, _ = swe.calc_ut(jd, swe.SUN, EPH | swe.FLG_SIDEREAL | swe.FLG_SPEED)
    return pos[0], pos[3]


def solve_crossing(target, guess_jd, ayan):
    jd = guess_jd
    for _ in range(40):
        lon, speed = sun_sid(jd, ayan)
        delta = ((target - lon + 180.0) % 360.0) - 180.0
        jd += delta / speed
        if abs(delta) < 1e-9:
            break
    return jd


@lru_cache(maxsize=None)
def sankrantis(bs_year, ayan):
    """13 UT Julian days: Sun's entry into each sidereal sign; [12] = next Mesha."""
    ad = bs_year - 57
    cal = swe.GREG_CAL if ad >= 1583 else swe.JUL_CAL
    jd = solve_crossing(0.0, swe.julday(ad, 4, 14, 0.0, cal), ayan)
    out = [jd]
    for k in range(1, 13):
        jd = solve_crossing((30.0 * k) % 360.0, jd + 29.5, ayan)
        out.append(jd)
    return tuple(out)


_warned = [False]


@lru_cache(maxsize=None)
def sun_event(jdn, tz, kind):
    """UT Julian day of sunrise/sunset at Kathmandu on civil day jdn. kind: 'rise' or 'set'."""
    midnight_ut = jdn - 0.5 - tz / 24.0
    flag = swe.CALC_RISE if kind == "rise" else swe.CALC_SET
    try:
        res = swe.rise_trans(midnight_ut, swe.SUN, flag, KTM, 1013.25, 15.0, EPH)
        return res[1][0]
    except Exception as e:
        if not _warned[0]:
            print(f"warning: sunrise/sunset failed ({e}); using 06:00 / 18:00 local",
                  file=sys.stderr)
            _warned[0] = True
        return midnight_ut + (0.25 if kind == "rise" else 0.75)


def sunrise_jd(jdn, tz):
    return sun_event(jdn, tz, "rise")


def civil_start(jd, r):
    """Civil day (JDN) on which the month beginning at sankranti jd starts."""
    local = jd + r.tz / 24.0 + 0.5
    n = math.floor(local)
    c = r.convention
    if c == "sunrise":
        return n if jd <= sun_event(n, r.tz, "rise") else n + 1
    if c == "sunset":
        return n if jd <= sun_event(n, r.tz, "set") else n + 1
    if c == "noon":
        return n if (local - n) * 24.0 < 12.0 else n + 1
    return n   # civil_date


def year_start(bs_year):
    r = rule_for(bs_year)
    return civil_start(sankrantis(bs_year, r.ayanamsa)[0], r)


@lru_cache(maxsize=None)
def month_starts(bs_year):
    """13 JDNs: 12 month starts + next Baishakh 1 (years never overlap or leave gaps)."""
    r = rule_for(bs_year)
    s = sankrantis(bs_year, r.ayanamsa)
    t = [civil_start(s[i], r) for i in range(12)]
    t.append(year_start(bs_year + 1))
    return tuple(t)


def month_info(y, idx):
    t = month_starts(y)
    return t[idx], t[idx + 1] - t[idx]


def ad_to_bs(jdn):
    y = ymd_from_jdn(jdn)[0]
    for by in (y + 56, y + 57, y + 58):
        t = month_starts(by)
        if t[0] <= jdn < t[12]:
            for i in range(12):
                if t[i] <= jdn < t[i + 1]:
                    return {"year": by, "month_idx": i, "day": jdn - t[i] + 1}
    return None


def bs_to_jdn(y, idx, day):
    start, total = month_info(y, idx)
    if not 1 <= day <= total:
        raise ValueError(f"{BS_MONTH_NAMES[idx]} {y} has {total} days, not {day}")
    return start + day - 1


def shift_month(y, idx, delta):
    n = y * 12 + idx + delta
    return n // 12, n % 12


# ============================ TITHI ======================================
TITHI_NAMES = ["Pratipada", "Dwitiya", "Tritiya", "Chaturthi", "Panchami", "Shashthi",
               "Saptami", "Ashtami", "Navami", "Dashami", "Ekadashi", "Dwadashi",
               "Trayodashi", "Chaturdashi"]
HIGHLIGHT = {8, 11, 15, 23, 26, 30}


def tithi_at(jd):
    s = swe.calc_ut(jd, swe.SUN, EPH)[0][0]
    m = swe.calc_ut(jd, swe.MOON, EPH)[0][0]
    return int(((m - s) % 360.0) // 12.0) + 1


def tithi_label(n):
    if n == 15:
        return "Purnima"
    if n == 30:
        return "Amavasya"
    return f"{'Shukla' if n < 15 else 'Krishna'} {TITHI_NAMES[(n - 1) % 15]}"


def skipped(a, b):
    """Tithis that start and end between two consecutive sunrises (kshaya)."""
    d = (b - a) % 30
    return [(a + i - 1) % 30 + 1 for i in range(1, d)] if d >= 2 else []


def elong(jd):
    s = swe.calc_ut(jd, swe.SUN, EPH | swe.FLG_SPEED)[0]
    m = swe.calc_ut(jd, swe.MOON, EPH | swe.FLG_SPEED)[0]
    return (m[0] - s[0]) % 360.0, m[3] - s[3]


def tithi_end_jd(n, from_jd):
    target = (12.0 * n) % 360.0
    e, rate = elong(from_jd)
    jd = from_jd + ((target - e) % 360.0) / rate
    for _ in range(30):
        e, rate = elong(jd)
        delta = ((target - e + 180.0) % 360.0) - 180.0
        jd += delta / rate
        if abs(delta) < 1e-7:
            break
    return jd


def local_hm(jd, tz, civil_jdn=None):
    local = jd + tz / 24.0 + 0.5
    n = math.floor(local)
    h = (local - n) * 24.0
    tag = f" {n - civil_jdn:+d}d" if civil_jdn is not None and n != civil_jdn else ""
    return f"{int(h):02d}:{int(h % 1 * 60):02d}{tag}"


# ============================ RENDERING ==================================
RESET, INV = "\033[0m", "\033[1;7m"


def format_month_lines(y, idx, today):
    start, total = month_info(y, idx)
    lines = [f"{BS_MONTH_NAMES[idx]} {y}".center(20), "Su Mo Tu We Th Fr Sa"]
    col, row = dow(start), ["   "] * 7
    for day in range(1, total + 1):
        is_today = today and (today["year"], today["month_idx"], today["day"]) == (y, idx, day)
        row[col] = f"{INV}{day:2d}{RESET} " if is_today else f"{day:2d} "
        col += 1
        if col > 6:
            lines.append("".join(row).rstrip())
            row, col = ["   "] * 7, 0
    if any(c.strip() for c in row):
        lines.append("".join(row).rstrip())
    while len(lines) < 8:
        lines.append("")
    return lines


def pad(s, w):
    return s + " " * max(0, w - len(s.replace(RESET, "").replace(INV, "")))


def render_three(y, idx, today):
    cols = [format_month_lines(*shift_month(y, idx, d), today) for d in (-1, 0, 1)]
    print(f"\n====================== {y} ======================\n")
    for a, b, c in zip(*cols):
        print(f"{pad(a, 22)}   {pad(b, 22)}   {pad(c, 22)}")
    print()


def render_month(y, idx, today):
    for line in format_month_lines(y, idx, today):
        if line:
            print(line)
    s, n = month_info(y, idx)
    print(f"({fmt_ad(s)} .. {fmt_ad(s + n - 1)})")


def render_year(y, today):
    print(f"\n================ {y} ================\n")
    for m in range(12):
        for line in format_month_lines(y, m, today):
            if line:
                print(line)
        print()


def render_detail(y, idx, today):
    start, total = month_info(y, idx)
    tz = rule_for(y).tz
    print(f"\n{BS_MONTH_NAMES[idx]} {y}   ({fmt_ad(start)} .. {fmt_ad(start + total - 1)})\n")
    print(f"{'BS':<4} {'Day':<4} {'AD date':<14} Tithi at sunrise (Kathmandu)")
    th = [tithi_at(sunrise_jd(start + i, tz)) for i in range(total + 1)]
    for i in range(total):
        sk = skipped(th[i], th[i + 1])
        extra = f"  (skips {', '.join(tithi_label(x) for x in sk)})" if sk else ""
        star = " *" if th[i] in HIGHLIGHT else ""
        here = today and (today["year"], today["month_idx"], today["day"]) == (y, idx, i + 1)
        print(f"{i + 1:<4} {DAYS[dow(start + i)]:<4} {fmt_ad(start + i):<14} "
              f"{tithi_label(th[i])}{star}{extra}{'  <-- today' if here else ''}")


# ============================ REFERENCE / FIT (optional) =================
def load_ref(path):
    """CSV lines: bs_year,month(1-12),ad_year,ad_month,ad_day. '#' starts a comment."""
    ref = {}
    with open(path) as f:
        for line in f:
            parts = [p.strip() for p in line.split("#")[0].split(",")]
            try:
                by, m, y, mo, d = map(int, parts[:5])
            except ValueError:
                continue
            ref[(by, m - 1)] = jdn_from_ymd(y, mo, d)
    return ref


def cmd_fit(ref):
    """Score every (ayanamsa, tz, convention) combination against the reference."""
    res = []
    for a in AYANAMSAS:
        if not hasattr(swe, AYANAMSAS[a]):
            continue
        for tz in (5.5, 5.75):
            for c in CONVENTIONS:
                r = Rule(ayanamsa=a, tz=tz, convention=c)
                ok = sum(civil_start(sankrantis(by, a)[m], r) == j for (by, m), j in ref.items())
                res.append((ok, a, tz, c))
    res.sort(key=lambda x: -x[0])
    print(f"{len(ref)} reference month starts:")
    for ok, a, tz, c in res[:10]:
        print(f"  {ok:>4}/{len(ref)}  ayanamsa={a:<20} tz={tz:<5} convention={c}")
    _, a, tz, c = res[0]
    r = Rule(ayanamsa=a, tz=tz, convention=c)
    bad = {}
    for (by, m), j in sorted(ref.items()):
        d = j - civil_start(sankrantis(by, a)[m], r)
        if d:
            bad.setdefault(by, []).append(f"{BS_MONTH_NAMES[m][:3]}{d:+d}")
    print(f"\nMismatches for the best rule ({a}, tz={tz}, {c}):")
    for by, v in bad.items():
        print(f"  {by}: {' '.join(v)}")
    if not bad:
        print("  none")


# ============================ DEBUG ======================================
def compare_year(y):
    """Show month starts under every convention. No data needed: rows that differ are the
    only months where the convention matters."""
    r = rule_for(y)
    s = sankrantis(y, r.ayanamsa)
    print(f"BS {y}: ayanamsa={r.ayanamsa} tz={r.tz}")
    print(f"{'Month':<9} {'Sankranti (local)':<19} " +
          " ".join(f"{c:<12}" for c in CONVENTIONS) + " Differ?")
    for i in range(12):
        local = s[i] + r.tz / 24.0 + 0.5
        n = math.floor(local)
        h = (local - n) * 24.0
        starts = [civil_start(s[i], replace(r, convention=c)) for c in CONVENTIONS]
        print(f"{BS_MONTH_NAMES[i]:<9} {fmt_ad(n)} {int(h):02d}:{int(h % 1 * 60):02d}  " +
              " ".join(f"{fmt_ad(x):<12}" for x in starts) +
              ("  <-- convention-sensitive" if len(set(starts)) > 1 else ""))


def debug_year(y, ref):
    r = rule_for(y)
    s = sankrantis(y, r.ayanamsa)
    t = month_starts(y)
    print(f"BS {y}: ayanamsa={r.ayanamsa} tz={r.tz} convention={r.convention}")
    print(f"{'Month':<9} {'Sankranti (local)':<19} {'Start':<13} {'Dow':<4} {'Len':<4} "
          f"{'vs sunrise':<14} Ref")
    for i in range(12):
        local = s[i] + r.tz / 24.0 + 0.5
        n = math.floor(local)
        h = (local - n) * 24.0
        margin = (s[i] - sunrise_jd(n, r.tz)) * 1440.0   # + = after that day's sunrise
        refc = ""
        if ref and (y, i) in ref:
            d = ref[(y, i)] - t[i]
            refc = "ok" if d == 0 else f"ref is {d:+d} day"
        print(f"{BS_MONTH_NAMES[i]:<9} {fmt_ad(n)} {int(h):02d}:{int(h % 1 * 60):02d}   "
              f"{fmt_ad(t[i]):<13} {DAYS[dow(t[i])]:<4} {t[i + 1] - t[i]:<4} "
              f"{margin:>+8.0f} min   {refc}")
    print(f"Next Baishakh 1: {fmt_ad(t[12])}")


def debug_month(y, idx, ref=None):
    r = rule_for(y)
    tz = r.tz
    start, total = month_info(y, idx)
    sj = sankrantis(y, r.ayanamsa)[idx]
    n = math.floor(sj + tz / 24.0 + 0.5)
    print(f"BS {y} {BS_MONTH_NAMES[idx]}: ayanamsa={r.ayanamsa} tz={tz} convention={r.convention}")
    print(f"Sankranti (local): {fmt_ad(n)} {local_hm(sj, tz)}   "
          f"sunrise that day {local_hm(sunrise_jd(n, tz), tz)}   "
          f"->  month starts {fmt_ad(start)} {DAYS[dow(start)]}, {total} days")
    if ref and (y, idx) in ref:
        d = ref[(y, idx)] - start
        print(f"Reference start: {fmt_ad(ref[(y, idx)])}  ->  {'ok' if d == 0 else f'ref is {d:+d} day'}")

    rises = [sunrise_jd(start + i, tz) for i in range(total + 1)]
    th = [tithi_at(j) for j in rises]
    t = datetime.date.today()
    today_jdn = jdn_from_ymd(t.year, t.month, t.day)
    print(f"\n{'BS':<3} {'Dow':<4} {'AD date':<13} {'Sunrise':<8} {'Tithi at sunrise':<22} "
          f"{'Ends (local)':<14} Notes")
    for i in range(total):
        k = th[i]
        end = local_hm(tithi_end_jd(k, rises[i]), tz, start + i)
        notes = []
        if k in HIGHLIGHT:
            notes.append("*")
        if i > 0 and th[i - 1] == k:
            notes.append("repeat (vriddhi)")
        sk = skipped(k, th[i + 1])
        if sk:
            notes.append("skips " + ", ".join(tithi_label(x) for x in sk))
        if start + i == today_jdn:
            notes.append("<-- today")
        print(f"{i + 1:<3} {DAYS[dow(start + i)]:<4} {fmt_ad(start + i):<13} "
              f"{local_hm(rises[i], tz, start + i):<8} {tithi_label(k):<22} {end:<14} "
              f"{' '.join(notes)}")
    print("\n* = Ashtami/Ekadashi/Purnima/Amavasya at sunrise.  "
          "'Ends' is when that tithi finishes (+1d = next civil day).")


# ============================ CLI ========================================
def main():
    p = argparse.ArgumentParser(description="Bikram Samvat calendar (rule-based)")
    p.add_argument("-m", type=int, metavar="MONTH")
    p.add_argument("-y", action="store_true", help="whole year")
    p.add_argument("--three", action="store_true", help="previous, current, next month")
    p.add_argument("-d", "--detail", action="store_true", help="detailed month with tithi")
    p.add_argument("--ad", nargs=3, type=int, metavar=("Y", "M", "D"))
    p.add_argument("--bc", action="store_true", help="with --ad: Y is a BC year")
    p.add_argument("--bs", nargs=3, type=int, metavar=("Y", "M", "D"))
    p.add_argument("--convention", choices=CONVENTIONS, help="override month-start convention")
    p.add_argument("--compare", type=int, metavar="BS_YEAR",
                   help="month starts under every convention")
    p.add_argument("--debug", type=int, metavar="BS_YEAR")
    p.add_argument("--debug-month", nargs=2, type=int, metavar=("BS_YEAR", "MONTH"))
    p.add_argument("--rules", metavar="FILE", help="JSON era rules")
    p.add_argument("--ref", metavar="FILE", help="optional reference CSV")
    p.add_argument("--fit", action="store_true", help="score rules against --ref")
    p.add_argument("args", nargs="*", type=int, help="[month] [year] or [year]")
    a = p.parse_args()

    global RULES
    if a.rules:
        load_rules(a.rules)
    if a.convention:
        RULES = [replace(r, convention=a.convention) for r in RULES]
        month_starts.cache_clear()
    ref = load_ref(a.ref) if a.ref else None

    if a.fit:
        if not ref:
            sys.exit("--fit needs --ref FILE")
        return cmd_fit(ref)
    if a.compare is not None:
        return compare_year(a.compare)
    if a.debug is not None:
        return debug_year(a.debug, ref)
    if a.debug_month:
        return debug_month(a.debug_month[0], a.debug_month[1] - 1, ref)
    if a.ad:
        y, m, d = a.ad
        r = ad_to_bs(jdn_from_ymd(1 - y if a.bc else y, m, d))
        return print(f"{r['year']} {BS_MONTH_NAMES[r['month_idx']]} {r['day']}")
    if a.bs:
        y, m, d = a.bs
        try:
            return print(fmt_ad(bs_to_jdn(y, m - 1, d)))
        except ValueError as e:
            sys.exit(str(e))

    t = datetime.date.today()
    today = ad_to_bs(jdn_from_ymd(t.year, t.month, t.day))
    year, month = today["year"], None
    if len(a.args) == 1:
        year, a.y = a.args[0], True
    elif len(a.args) >= 2:
        month, year = a.args[0], a.args[1]
    if a.m:
        month = a.m

    if a.detail:
        render_detail(year, (month or today["month_idx"] + 1) - 1, today)
    elif a.three:
        render_three(year, (month or today["month_idx"] + 1) - 1, today)
    elif a.y:
        render_year(year, today)
    elif month:
        render_month(year, month - 1, today)
    else:
        render_three(year, today["month_idx"], today)


if __name__ == "__main__":
    main()
